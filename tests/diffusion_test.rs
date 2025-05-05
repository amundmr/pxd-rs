use pxd::math::numerical_methods::{ftcs_stable, forward_time_centered_space_linear, forward_time_centered_space_radial};
use std::fs::File;
use std::io::{Write, BufWriter};

#[test]
fn test_fcts() {
    // Set up space system
    let length: f64 = 1e-3; // meters
    let dx: f64 = 1e-6; // meters
    let nx: usize = (length / dx) as usize;
    let initial_concentration: f64 = 1000.0; // mol/m^3
    let mut concentration: Vec<f64> = vec![initial_concentration; nx];

    // Set a gradient so we can see the diffusion
    concentration[0] = 1300.0; // mol/m^3

    // Set up time
    let time_span: f64 = 60.0*20.0; // seconds
    let dt: f64 = 0.1; // seconds
    let nt: usize = (time_span / dt) as usize;
    let diffusion_coeff: f64 = 1e-12; // ??

    // Check stability with given conditions
    assert!(ftcs_stable(dt, dx, diffusion_coeff), "The parameters used makes the FTCS method unstable.");

    // Optionally create a file writer if environment variable is set
    let mut maybe_writer = if std::env::var("WRITE_TEST_OUTPUT").is_ok() {
        let file = File::create("concentration_over_time.csv").expect("Could not create file");
        Some(BufWriter::new(file))
    } else {
        None
    };

    println!("The first 60 timesteps of the leftmost and second leftmost points:");
    println!("{:?}, {:?}", concentration[0], concentration[1]);
    let mut x0: f64 = concentration[0];
    let x_rest_init: Vec<f64> = concentration[1..].to_vec();

    for i in 0..nt {
        // Step the method one step
        forward_time_centered_space_linear(&mut concentration, dx, dt, diffusion_coeff, 0.0);

        // Print if test fails
        if i < 30 {
            println!("{:?}, {:?}", concentration[0], concentration[1]);
        }

        // Maybe save to file
        if let Some(writer) = maybe_writer.as_mut() {
            let line = concentration
                .iter()
                .map(|c| c.to_string())
                .collect::<Vec<_>>()
                .join(",");
            writeln!(writer, "{line}").expect("Could not write to file");
        }

        // Assert monotonic decrease of peak
        assert!(x0 >= concentration[0], "The concentration peak is not monotonicaly decreasing.");
        x0 = concentration[0];

        // Assert that all the other points have increased compared to t = 0
        let all_less_or_equal = x_rest_init
            .iter()
            .zip(&concentration[1..])
            .all(|(a, b)| a <= b);
        assert!(all_less_or_equal, "The concentration of some x-point is below the initial value.");

    }
}