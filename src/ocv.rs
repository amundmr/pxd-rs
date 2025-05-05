use crate::{PARTICLE_DISCRETISATION, model::Particle};
// Holds the open circuit voltage functions for LG MJ1 18650 cell

pub fn open_circuit_voltage_graphite_si(particle: &Particle) -> f64 {
    let c: f64 = particle.concentration[PARTICLE_DISCRETISATION - 1];
    let c_max: f64 = particle.concentration_max;
    let x: f64 = c / c_max;

    let p = [
        1.20912055e+00,
        5.62297420e+01,
        -1.11020020e-01,
        -2.53458213e-01,
        4.92581391e+01,
        1.22046522e-02,
        4.73538620e-02,
        1.79631246e+01,
        1.75283209e-01,
        1.88038929e-02,
        3.03255334e+01,
        4.66328034e-01,
    ];

    p[0] * (-p[1] * x).exp() + p[2]
        - p[3] * (p[4] * (x - p[5])).tanh()
        - p[6] * (p[7] * (x - p[8])).tanh()
        - p[9] * (p[10] * (x - p[11])).tanh()
}
pub fn open_circuit_voltage_nmc811(particle: &Particle) -> f64 {
    let c: f64 = particle.concentration[PARTICLE_DISCRETISATION - 1];
    let c_max: f64 = particle.concentration_max;
    let x = c / c_max;

    let p = [
        0.74041974,
        4.39107343,
        0.03434767,
        18.16841489,
        0.53463176,
        17.68283504,
        14.59709162,
        0.28835348,
        17.58474971,
        14.69911523,
        0.28845641,
    ];

    -p[0] * x + p[1] - p[2] * (p[3] * (x - p[4])).tanh() - p[5] * (p[6] * (x - p[7])).tanh()
        + p[8] * (p[9] * (x - p[10])).tanh()
}
