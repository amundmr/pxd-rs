// Todo: Build an actual API
pub mod math;
pub mod model;
pub mod ocv;

pub trait Simulate {
    fn simulate(&mut self, time: &[f64], current: &[f64]) -> Vec<f64>;
}
