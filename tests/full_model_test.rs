use pxd::model::SPMeModel;
use pxd::Simulate;

#[test]
fn simulate_default_cycle() {
    // Generate input
    let time: f64 = 60.0 * 20.0; // seconds
    let dt: f64 = 0.1; // seconds
    let current: f64 = 1.0; // Amperes
    let n_steps: usize = (2.0 * time / dt).ceil() as usize;
    let mut t: Vec<f64> = vec![0.0; n_steps];
    let mut i: Vec<f64> = vec![0.0; n_steps];

    for step in 0..n_steps {
        let time_val = step as f64 * dt;
        t[step] = time_val;

        let current_val = if time_val < time {
                current
            } else {
                -current
            };
        i[step] = current_val;
    }
    let mut model = SPMeModel::default();

    model.simulate(&t, &i);
}