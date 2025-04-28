use pxd::model::SPMeModel;
use pxd::Simulate;

use std::fs::File;
use std::io::{Write, BufWriter};

fn save_vec_to_file(vec: &Vec<f64>, filename: &str) -> std::io::Result<()> {
    let file = File::create(filename)?;
    let mut writer = BufWriter::new(file);
    for value in vec {
        writeln!(writer, "{}", value)?;
    }
    Ok(())
}

#[test]
fn simulate_default_cycle() {
    // Generate input
    // 1C charge and 1C discharge for LG MJ1 18650 cell
    let time: f64 = 60.0 * 60.0; // seconds
    let dt: f64 = 0.1; // seconds
    let current: f64 = 3.2; // Amperes
    let n_steps: usize = (2.0 * time / dt).ceil() as usize;
    let mut t: Vec<f64> = vec![0.0; n_steps];
    let mut i: Vec<f64> = vec![0.0; n_steps];

    for step in 0..n_steps {
        let time_val = step as f64 * dt;
        t[step] = time_val;

        let current_val = if time_val < time {
                current // charge
            } else {
                -current
            };
        i[step] = current_val;
    }
    let mut model = SPMeModel::default();

    let cell_potential = model.simulate(&t, &i);
    println!("{:?}", cell_potential);
    save_vec_to_file(&cell_potential, "cell_voltage.csv").unwrap();
}