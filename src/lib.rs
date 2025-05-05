// Todo: Build an actual API
pub mod math;
pub mod model;
pub mod ocv;

pub const PARTICLE_DISCRETISATION: usize = 20;
pub const ELECTROLYTE_DISCRETISATION: usize = 20;

pub struct SimulationStepOutput {
    pub cell_potential: f64,
    pub electrolyte_concentration: [f64; ELECTROLYTE_DISCRETISATION],
}

pub trait Simulate {
    fn simulate(
        &mut self,
        time: &[f64],
        current: &[f64],
    ) -> impl Iterator<Item = SimulationStepOutput>;
}
