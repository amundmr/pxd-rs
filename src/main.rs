use pxd::Simulate;
use pxd::model::SPMeModel;

use std::error::Error;
use std::fs::OpenOptions;
use std::io::{BufWriter, Write};

fn main() -> Result<(), Box<dyn Error>> {
    // Generate input
    // 1C charge and 1C discharge for LG MJ1 18650 cell
    let time: f64 = 60.0 * 60.0; // seconds
    let dt: f64 = 0.001; // seconds
    let current: f64 = 3.2; // Amperes
    let n_steps: usize = (2.0 * time / dt).ceil() as usize;
    let mut t: Vec<f64> = vec![0.0; n_steps];
    let mut i: Vec<f64> = vec![0.0; n_steps];

    for step in 0..n_steps {
        let time_val: f64 = step as f64 * dt;
        t[step] = time_val;

        let current_val: f64 = if time_val < time {
            current // charge
        } else {
            -current
        };
        i[step] = current_val;
    }
    let mut model = SPMeModel::default();

    let mut w = BufWriter::new(
        OpenOptions::new()
            .create(true) // Create the file if it doesn't exist
            .append(true) // Open the file in append mode
            .open("cell_voltage.csv")?,
    );

    Ok(model
        .simulate(&t, &i)
        .try_for_each(|out| writeln!(w, "{}", out.cell_potential))?)
}
