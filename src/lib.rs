const RADIAL_DISCRETISATION: usize = 20;
const AXIAL_DISCRETISATION: usize = 20;

#[derive(Debug, Clone, Copy)]
pub struct Particle {
    pub radius: f64,
    pub diffusion_coeff: f64,
    pub concentration: [f64; RADIAL_DISCRETISATION],
}

#[derive(Debug, Clone, Copy)]
pub struct Electrolyte {
    pub concentration: [f64; AXIAL_DISCRETISATION],
    pub conductivity: f64,
    pub diffusion_coeff: f64,
}

#[derive(Debug, Clone, Copy)]
pub struct Potential {
    pub phi_s: [f64; AXIAL_DISCRETISATION],
    pub phi_e: [f64; AXIAL_DISCRETISATION],
}

#[derive(Debug, Clone, Copy)]
pub struct SPMeModel {
    pub time: f64,
    pub delta_x: f64,
    pub negative_electrode_particle: Particle,
    pub positive_electrode_particle: Particle,
    pub electrolyte: Electrolyte,
    pub potential: Potential,
    pub current_density: f64,
}
