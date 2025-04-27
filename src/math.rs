pub mod numerical_methods {

    pub fn forward_time_centered_space(y: &mut Vec<f64>, dx: f64, dt: f64, a: f64) {
        // Forward Time Centered Space (FTCS) is an integration method that utilises
        // the finite difference method (FDM). Refer to https://en.wikipedia.org/wiki/FTCS_scheme.
        // It ises forward euler in time, and central difference in space.
        //
        // borrows mutable y-vector, x-step, timestep, the diffusion coefficient, and steps the y-vector forward by one timestep.
        let n: usize = y.len();
        let mut y_new: Vec<f64> = vec![0.0; n];

        let dx2 = dx * dx;
        let adt = a * dt;

        // Interior points
        for i in 1..n-1 {
            let d2y_dx2 = (y[i-1] - 2.0 * y[i] + y[i+1]) / dx2; // Central difference
            y_new[i] = y[i] + adt * d2y_dx2; // Forward Euler
        }

        // Left boundary by forward difference and flux
        let dy_dx_left = (y[1] - y[0]) / dx;
        let flux_left = 0.0;
        y_new[0] = y[0] + adt * (2.0 * (dy_dx_left - flux_left)/dx);

        // right boundary by backward difference and flux
        let dy_dx_right = (y[n-1] - y[n-2]) / dx;
        let flux_right = 0.0;
        y_new[n-1] = y[n-1] + adt * (2.0 * (dy_dx_right - flux_right)/dx);

        y.copy_from_slice(&y_new);
    }

    pub fn ftcs_stable(dt: f64, dx: f64, alpha: f64) -> bool {
        // Returns stability (bool) of forward time centered space method
        // for the heat equation (equaling fickian diffusion)
        // given the timestep, x-step and heat transfer coefficient (diffusion coeff).
        dt <= dx * dx / (2.0 * alpha)
    }
}