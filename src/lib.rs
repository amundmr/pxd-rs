// Todo: Build an actual API
pub mod model;
pub mod math;


pub trait Simulate {
    fn simulate(&mut self, time: &[f64], current: &[f64]) -> Vec<f64>;
}