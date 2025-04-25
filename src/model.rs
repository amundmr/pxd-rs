use crate::input::CurrentFunction;

use std::f64::consts::PI;

const PARTICLE_DISCRETISATION: usize = 20;
const ELECTROLYTE_DISCRETISATION: usize = 20;
const FARADAY: f64 = 96485.3321233100184; // C/mol (=As/mol), 2019 SI revision definition
const STANDARD_TEMPERATURE: f64 = 298.15; // Kelvin

#[derive(Debug, Clone, Copy)]
pub struct Particle {
    pub radius: f64,
    pub diffusion_coeff: f64,
    pub concentration: [f64; PARTICLE_DISCRETISATION],
    pub concentration_max: f64,
    pub open_circuit_voltage: fn(f64, f64)->f64,
}

#[derive(Debug, Clone, Copy)]
pub struct Electrolyte {
    pub concentration: [f64; ELECTROLYTE_DISCRETISATION],
    pub conductivity: f64,
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
    pub current_function: CurrentFunction,
}

impl SPMeModel {

    /// Yields an SPMe model with default parameters for an LG MJ1 18650 cylindrical cell
    pub fn default() -> Self {
        fn open_circuit_voltage_graphite_si(c_max: f64, c: f64) -> f64 {
            let x: f64 = c / c_max;
    
            let p = [
                1.20912055e+00, 5.62297420e+01, -1.11020020e-01, -2.53458213e-01,
                4.92581391e+01, 1.22046522e-02, 4.73538620e-02, 1.79631246e+01,
                1.75283209e-01, 1.88038929e-02, 3.03255334e+01, 4.66328034e-01,
            ];
    
            p[0] * (-p[1] * x).exp()
            + p[2]
            - p[3] * (p[4] * (x - p[5])).tanh()
            - p[6] * (p[7] * (x - p[8])).tanh()
            - p[9] * (p[10] * (x - p[11])).tanh()
        }
        pub fn open_circuit_voltage_nmc811(c_max: f64, c: f64) -> f64 {
            let x = c / c_max;
    
            let p = [
                0.74041974, 4.39107343, 0.03434767, 18.16841489, 0.53463176,
                17.68283504, 14.59709162, 0.28835348, 17.58474971, 14.69911523,
                0.28845641,
            ];
    
            -p[0] * x
            + p[1]
            - p[2] * (p[3] * (x - p[4])).tanh()
            - p[5] * (p[6] * (x - p[7])).tanh()
            + p[8] * (p[9] * (x - p[10])).tanh()
        }

        SPMeModel {
            negative_electrode: Electrode{
                height: 0.059, // meters
                width: 1.22,    // meters
                thickness: 86.7e-6, // meters
                particle: Particle{
                    radius: 6.1e-6, // meters
                    diffusion_coeff: 5e-14, // m^2/s
                    concentration: [1000.0; PARTICLE_DISCRETISATION], // mol/m^3
                    concentration_max: 34684.0, // mol/m^3
                    open_circuit_voltage: open_circuit_voltage_graphite_si,
                },
            },
            positive_electrode: Electrode {
                height: 0.059, // meters
                width: 1.22, // meters
                thickness: 66.2e-6, // meters
                particle: Particle{
                    radius: 3.8e-6, // meters
                    diffusion_coeff: 5e-14, // m^2/s
                    concentration: [49000.0; PARTICLE_DISCRETISATION], // mol/m^3
                    concentration_max: 50060.0, // mol/m^3
                    open_circuit_voltage: open_circuit_voltage_nmc811,
                }
            },
            electrolyte: Electrolyte{
                concentration: [1000.0; ELECTROLYTE_DISCRETISATION],
                conductivity: 0.8, // S/m
            },
            current_function: CurrentFunction::default(),
        }
    }

    // /// Returns cell voltage given the particle surface lithium concentration of the negative and positive electrode
    // fn u_cell(c_n_s: f64, c_p_s: f64) -> f64 {
    //     self.p.u_p(c_p_s) - self.p.u_n(c_n_s)
    // }

    // fn sphere_volume(r: f64) -> f64 {
    //     4.0/3.0 * self.p.pi*r**3
    // }

    // /// Returns the change in concentration for a given time, concentration, electrode domain and particle radius
    // fn dcdt(&self, time: f64, c: f64, domain: u8, r: f64) -> f64 {

    //     let mut polarity = 1;
    //     if domain == 0 { // negative electrode
    //         polarity = 1;
    //     } else if domain == 2 { // positive electrode
    //         polarity = -1;
    //     }

    //     polarity * self.cf.current(time) / (self.p.n * self.p.F * self.sphere_volume(r))
    // }
}