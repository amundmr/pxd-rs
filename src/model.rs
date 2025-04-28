use crate::math::numerical_methods::{
    forward_time_centered_space_linear, forward_time_centered_space_radial, ftcs_stable,
};
use crate::ocv;
use crate::Simulate;

use std::f64::consts::PI;

pub const PARTICLE_DISCRETISATION: usize = 20;
const ELECTROLYTE_DISCRETISATION: usize = 20;
const FARADAY: f64 = 96485.3321233100184; // C/mol (=As/mol), 2019 SI revision definition
const STANDARD_TEMPERATURE: f64 = 298.15; // Kelvin

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
}

#[derive(Debug, Clone)]
pub struct SPMeModel {
    pub negative_electrode: Electrode,
    pub positive_electrode: Electrode,
    pub electrolyte: Electrolyte,
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
            },
            electrolyte: Electrolyte {
                concentration: [1000.0; ELECTROLYTE_DISCRETISATION],
                conductivity: 0.8,      // S/m
                diffusion_coeff: 1e-11, // m^2/s
                thickness: 100e-6,      // meters
            },
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

    fn cell_potential(&self) -> f64 {
        ocv::open_circuit_voltage_nmc811(&self.positive_electrode.particle)
            - ocv::open_circuit_voltage_graphite_si(&self.negative_electrode.particle)
    }

    fn get_number_of_particles(&self, electrode: &Electrode) -> usize {
        // TODO: this scaling concept doesn't work at all. Dividing by 3_000_000.0 is just a guess.
        // Calculate number of particles
        // (electrode.thickness * electrode.height * electrode.width
        // / (4.0/3.0 * PI * electrode.particle.radius.powi(3))) as usize
        ((electrode.thickness / electrode.particle.radius)
            * (electrode.height / electrode.particle.radius)
            * (electrode.width / electrode.particle.radius)
            / 3_000_00.0) as usize
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

        // Calculate number of particles
        // TODO: Scaling the current by the number of particles doesn't work at all.
        let n_negative_particles: usize = self.get_number_of_particles(&self.negative_electrode);
        let n_positive_particles: usize = self.get_number_of_particles(&self.positive_electrode);

        // Set up cell potential over time
        let mut cell_potential: Vec<f64> = vec![0.0; time.len()];
        let electrolyte_dx: f64 = self.electrolyte.thickness / ELECTROLYTE_DISCRETISATION as f64;

        for i in 0..time.len() {
            // Step the electrolyte concentration in time
            forward_time_centered_space_linear(
                &mut self.electrolyte.concentration,
                electrolyte_dx,
                dt,
                self.electrolyte.diffusion_coeff,
            );
            // Step the particles' concentration in time
            forward_time_centered_space_radial(
                &mut self.negative_electrode.particle.concentration,
                self.negative_electrode.particle.dr,
                dt,
                self.negative_electrode.particle.diffusion_coeff,
                self.negative_electrode.particle.radius,
                -current[i] / n_negative_particles as f64, // flux
            );
            forward_time_centered_space_radial(
                &mut self.positive_electrode.particle.concentration,
                self.positive_electrode.particle.dr,
                dt,
                self.positive_electrode.particle.diffusion_coeff,
                self.positive_electrode.particle.radius,
                current[i] / n_positive_particles as f64, // flux
            );

            // Calculate cell potential
            cell_potential[i] = self.cell_potential();
        }
        cell_potential
    }
}
