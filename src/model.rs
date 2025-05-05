use crate::math::numerical_methods::{
    forward_time_centered_space_linear, forward_time_centered_space_radial, ftcs_stable,
};
use crate::math::utils::arcsinh;
use crate::ocv;
use crate::Simulate;

use std::f64::consts::PI;
use std::fs::File;
use std::io::BufWriter;
use std::io::Write;
use std::fs::OpenOptions;

pub const PARTICLE_DISCRETISATION: usize = 20;
const ELECTROLYTE_DISCRETISATION: usize = 20;
const FARADAY: f64 = 96485.3321233100184; // C/mol (=As/mol), 2019 SI revision definition
const GAS_CONSTANT: f64 = 8.31446261815324; // J/(mol*K), 2019 SI revision definition
const STANDARD_TEMPERATURE: f64 = 298.15; // Kelvin
const CATION_TRANSFERENCE_NUMBER: f64 = 0.2594; // dimensionless

#[derive(Debug, Clone, Copy)]
pub struct Particle {
    pub radius: f64,
    pub dr: f64,
    pub diffusion_coeff: f64,
    pub concentration: [f64; PARTICLE_DISCRETISATION],
    pub concentration_max: f64,
    pub concentration_init: f64,
}

impl Particle {
    pub fn new(
        radius: f64,
        diffusion_coeff: f64,
        concentration_max: f64,
        concentration_init: f64,
    ) -> Self {
        let dr = radius / (PARTICLE_DISCRETISATION as f64);
        let concentration = [concentration_init; PARTICLE_DISCRETISATION];

        Particle {
            radius,
            dr,
            diffusion_coeff,
            concentration,
            concentration_max,
            concentration_init,
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct Electrolyte {
    // This is used in place of separator for the time being
    pub concentration: [f64; ELECTROLYTE_DISCRETISATION],
    pub conductivity: f64,
    pub diffusion_coeff: f64,
    pub thickness: f64,
}

#[derive(Debug, Clone, Copy)]
pub struct Electrode {
    pub height: f64,
    pub width: f64,
    pub thickness: f64,
    pub particle: Particle,
    pub active_material_volume_fraction: f64,
}

#[derive(Debug, Clone)]
pub struct SPMeModel {
    pub negative_electrode: Electrode,
    pub positive_electrode: Electrode,
    pub electrolyte: Electrolyte,
    pub concentration: Vec<[f64; ELECTROLYTE_DISCRETISATION]>,
}

impl Default for SPMeModel {
    // Default parameters for an LG MJ1 18650 cylindrical cell
    fn default() -> Self {
        let model = SPMeModel {
            negative_electrode: Electrode {
                height: 0.059,      // meters
                width: 1.22,        // meters
                thickness: 86.7e-6, // meters
                particle: Particle::new(
                    6.1e-6,  // meters
                    5e-14,   // m^2/s
                    34684.0, // mol/m^3
                    1000.0,
                ),
                active_material_volume_fraction: 0.694,
            },
            positive_electrode: Electrode {
                height: 0.059,      // meters
                width: 1.22,        // meters
                thickness: 66.2e-6, // meters
                particle: Particle::new(
                    3.8e-6,  // meters
                    5e-14,   // m^2/s
                    50060.0, // mol/m^3
                    49000.0, // mol/m^3
                ),
                active_material_volume_fraction: 0.754,
            },
            electrolyte: Electrolyte {
                concentration: [1000.0; ELECTROLYTE_DISCRETISATION],
                conductivity: 0.8,      // S/m
                diffusion_coeff: 1.7e-10, // m^2/s avg of Nyman et al. (2008) (fluctuates between 2.2e-10-1.3e-10 between 800-1200mol/m^3)
                thickness: 12e-6,      // meters
            },
            concentration: vec![[1000.0; ELECTROLYTE_DISCRETISATION]; 1],
        };
        model
    }
}

impl SPMeModel {
    fn assert_ftcs_stability(&self, dt: f64) {
        // Check stability of numerical method in particles and electrolyte
        assert!(
            ftcs_stable(
                dt,
                self.negative_electrode.particle.radius / PARTICLE_DISCRETISATION as f64,
                self.negative_electrode.particle.diffusion_coeff
            ),
            "FTCS method not stable for negative particle"
        );
        assert!(
            ftcs_stable(
                dt,
                self.positive_electrode.particle.radius / PARTICLE_DISCRETISATION as f64,
                self.positive_electrode.particle.diffusion_coeff
            ),
            "FTCS method not stable for positive particle"
        );
        assert!(
            ftcs_stable(
                dt,
                self.electrolyte.thickness / ELECTROLYTE_DISCRETISATION as f64,
                self.electrolyte.diffusion_coeff
            ),
            "FTCS method not stable for electrolyte"
        );
    }


    fn butler_volmer_overpotential(&self, current_density: f64, electrode: &Electrode) -> f64 {
        // Butler volmer overpotential, eta
        // TODO: make B-V function parameters tunable
        let alpha: f64 = 0.5; // charge transfer coefficient
        let reaction_rate_constant: f64 = 1e-3; // reaction rate constant
        let c_e: f64 = 1000.0; // electrolyte concentration mol/m^3, TODO: change this to real value
        let exchange_current_density: f64 = 
            reaction_rate_constant 
            * c_e.powf(alpha) 
            * electrode.particle.concentration[PARTICLE_DISCRETISATION - 1].powf(alpha)
            * (1.0 - ( electrode.particle.concentration[PARTICLE_DISCRETISATION - 1] / electrode.particle.concentration_max ) ).powf(alpha);

        // Symmetric butler volmer, only valid for alpha = 0.5
        let butler_volmer: f64 = ( 2.0 * GAS_CONSTANT * STANDARD_TEMPERATURE / FARADAY )
            * arcsinh(
                current_density / ( 2.0* exchange_current_density * self.specific_interfacial_surface_area(electrode) * electrode.thickness )
            );
        - butler_volmer // negative sign since we define positive current as charge.
    }

    fn electrolyte_concentration_overpotential(&self) -> f64 {
        // The electrolyte concentration overpotential is the voltage induced by the concentration gradient in the electrolyte.
        let electrolyte_concentration_n: f64 = self.electrolyte.concentration[0];
        let electrolyte_concentration_p: f64 = self.electrolyte.concentration[ELECTROLYTE_DISCRETISATION - 1];
        let eta_c: f64 = 
            2.0 // Accounts for potential drop at both sides
            * ( 1.0 - CATION_TRANSFERENCE_NUMBER ) // Describes how much of the current is carried by cations (Li+)
            * (GAS_CONSTANT * STANDARD_TEMPERATURE / FARADAY) // Nernst potential part 1
            * (electrolyte_concentration_p - electrolyte_concentration_n); // Nernst potential part 2
        eta_c
    }

    fn cell_potential(&self, current: f64) -> f64 {
        let cell_area: f64 = self.negative_electrode.height * self.negative_electrode.width;
        let current_density: f64 = current / cell_area; // A/m^2

        // Open circuit voltages, U(c)
        ocv::open_circuit_voltage_nmc811(&self.positive_electrode.particle)
            - ocv::open_circuit_voltage_graphite_si(&self.negative_electrode.particle)
            // Reaction/charge transfer overpotential, eta_r
            - self.butler_volmer_overpotential(current_density, &self.negative_electrode)
            - self.butler_volmer_overpotential(current_density, &self.positive_electrode)
            // Electrolyte concentration overpotential, eta_c
            + self.electrolyte_concentration_overpotential()
            // TODO: Ohmic overpotential in solid
            // TODO: Ohmic overpotential in electrolyte
    }

    fn specific_interfacial_surface_area(&self, electrode: &Electrode) -> f64 {
        // The specific interfacial surface area is the surface area per unit volume, and it's use
        // assumes a uniform distribution of monodisperse spherical particles.
        // $a = 3 \cdot \frac{\epsilon}{r}$
        let a: f64 = 3.0 * electrode.active_material_volume_fraction / electrode.particle.radius;
        a
    }

    fn particle_surface_flux(&self, current: f64, electrode: &Electrode) -> f64 {
        // The calculation of flux at the particle surface is the current density (cell current divided by electrode area) 
        // divided by Faraday's constant, $F$, (conversion of current to moles), the 
        // specific interfacial surface area, $a$, and the thickness of the electrode, $L$.
        let current_density: f64 = current / (electrode.height * electrode.width); // A/m^2
        let a: f64 = self.specific_interfacial_surface_area(electrode);
        let flux: f64 = current_density / ( FARADAY * a * electrode.thickness ); // mol/(s*m^2)
        flux
    }

    fn electrolyte_boundary_flux(&self, current: f64) -> f64 {
        // The flux at the electrolyte boundary is the current density (cell current divided by electrode area) 
        // divided by Faraday's constant, $F$, (conversion of current to moles), the 
        // specific interfacial surface area, $a$, and the thickness of the electrolyte, $L$.
        let electrode_area: f64 = self.negative_electrode.height * self.negative_electrode.width;
        let current_density: f64 = current / electrode_area; // A/m^2
        let flux: f64 = current_density / FARADAY; // mol/(s*m^2)
        flux
    }

    fn save_to_file(&self, vec: &Vec<[f64;20]>, filename: &str) -> std::io::Result<()> {
        let file = OpenOptions::new()
            .create(true) // Create the file if it doesn't exist
            .append(true) // Open the file in append mode
            .open(filename)?;
        let mut writer = BufWriter::new(file);
    
        // Iterate over each line (inner Vec<f64>)
        for line_vec in vec {
            let line = line_vec
                .iter()
                .map(|value| value.to_string()) // Convert each value to a string
                .collect::<Vec<String>>()
                .join(","); // Join values with commas
    
            // Write the entire line to the file
            writeln!(writer, "{}", line)?;
        }
    
        Ok(())
    }

    fn save_model_state(&self) -> std::io::Result<()> {
        self.save_to_file(
            &self.concentration.to_vec(),
            "electrolyte_concentration.csv",
        )?;
        Ok(())
    }

}

impl Simulate for SPMeModel {
    fn simulate(&mut self, time: &[f64], current: &[f64]) -> Vec<f64> {
        // Check that the time and current vectors are the same length
        assert_eq!(
            time.len(),
            current.len(),
            "Time and current vectors must be the same length"
        );
        // Check that the time vector is sorted
        assert!(
            time.windows(2).all(|w| w[0] < w[1]),
            "Time vector must be sorted"
        );
        // Check that the current vector is not empty
        assert!(!current.is_empty(), "Current vector must not be empty");
        // Check that the time vector has a constant timestep
        let dt: f64 = time[1] - time[0];
        assert!(
            time.windows(2).all(|w| (w[1] - w[0] - dt).abs() < 1e-9),
            "Time vector must have a constant timestep"
        );

        self.assert_ftcs_stability(dt);

        // Set up cell potential over time
        let mut cell_potential: Vec<f64> = vec![0.0; time.len()];
        let electrolyte_dx: f64 = self.electrolyte.thickness / ELECTROLYTE_DISCRETISATION as f64;

        for i in 0..time.len() {
            // Step the electrolyte concentration in time
            let flux_e: f64 = self.electrolyte_boundary_flux(current[i]);
            forward_time_centered_space_linear(
                &mut self.electrolyte.concentration,
                electrolyte_dx,
                dt,
                self.electrolyte.diffusion_coeff,
                flux_e/1000.0,
            );

            // Step the particles' concentration in time
            let flux_n: f64 = -self.particle_surface_flux(current[i], &self.negative_electrode);
            forward_time_centered_space_radial(
                &mut self.negative_electrode.particle.concentration,
                self.negative_electrode.particle.dr,
                dt,
                self.negative_electrode.particle.diffusion_coeff,
                self.negative_electrode.particle.radius,
                flux_n, // flux
            );
            let flux_p: f64 = self.particle_surface_flux(current[i], &self.positive_electrode);
            forward_time_centered_space_radial(
                &mut self.positive_electrode.particle.concentration,
                self.positive_electrode.particle.dr,
                dt,
                self.positive_electrode.particle.diffusion_coeff,
                self.positive_electrode.particle.radius,
                flux_p, // flux
            );

            // Calculate cell potential
            cell_potential[i] = self.cell_potential(current[i]);
            
            // TODO: Find a better way to save timeseries model state.
            if std::env::var("WRITE_MODEL_OUTPUT").is_ok() {
                self.concentration.push(self.electrolyte.concentration);
            }
        }
        if std::env::var("WRITE_MODEL_OUTPUT").is_ok() {
            self.save_model_state().unwrap();
        }
        cell_potential
        
    }
}
