pub mod numerical_methods {

    pub fn forward_time_centered_space_linear(y: &mut [f64], dx: f64, dt: f64, a: f64) {
        // Forward Time Centered Space (FTCS) is an integration method that utilises
        // the finite difference method (FDM). Refer to https://en.wikipedia.org/wiki/FTCS_scheme.
        // It ises forward euler in time, and central difference in space.
        // This is for the linear case (e.g. lithium in electrolyte)
        //
        // borrows mutable y-vector, x-step, timestep, the diffusion coefficient, and steps the y-vector forward by one timestep.
        let n: usize = y.len();

        let dx2 = dx * dx;
        let adt = a * dt;

        let mut y_prev: f64 = y[0];

        // Left boundary by forward difference (before y[1] is iterated)
        let dy_dx_left = (y[1] - y[0]) / dx;
        // Right boundary by backward difference (before y[n-2] is iterated)
        let dy_dx_right = (y[n-1] - y[n-2]) / dx;

        // Interior points
        for i in 1..n-1 {
            let d2y_dx2 = (y_prev - 2.0 * y[i] + y[i+1]) / dx2; // Central difference
            y_prev = y[i]; // Make sure we store the value before overwriting
            y[i] = y[i] + adt * d2y_dx2; // Forward Euler
        }

        // Left boundary by forward difference
        let flux_left = 0.0;
        y[0] = y[0] + adt * (2.0 * (dy_dx_left - flux_left)/dx);

        // Right boundary by backward difference
        let flux_right = 0.0;
        y[n-1] = y[n-1] + adt * (2.0 * (dy_dx_right - flux_right)/dx);
    }

    pub fn ftcs_stable(dt: f64, dx: f64, alpha: f64) -> bool {
        // Returns stability (bool) of forward time centered space method
        // for the heat equation (equaling fickian diffusion)
        // given the timestep, x-step and heat transfer coefficient (diffusion coeff).
        dt <= dx * dx / (2.0 * alpha)
    }

    pub fn forward_time_centered_space_radial(y: &mut [f64], dx: f64, dt: f64, a: f64) {
        // Forward Time Centered Space (FTCS) is an integration method that utilises
        // the finite difference method (FDM). Refer to https://en.wikipedia.org/wiki/FTCS_scheme.
        // It ises forward euler in time, and central difference in space.
        // This is for the linear case (e.g. lithium in electrolyte)
        //
        // borrows mutable y-vector, x-step, timestep, the diffusion coefficient, and steps the y-vector forward by one timestep.
        let n: usize = y.len();

        let dx2 = dx * dx;
        let adt = a * dt;

        let mut y_prev: f64 = y[0];

        // Left boundary by forward difference (before y[1] is iterated)
        let dy_dx_left = (y[1] - y[0]) / dx;
        // Right boundary by backward difference (before y[n-2] is iterated)
        let dy_dx_right = (y[n-1] - y[n-2]) / dx;

        // Interior points
        for i in 1..n-1 {
            let d2y_dx2 = (y_prev - 2.0 * y[i] + y[i+1]) / dx2; // Central difference
            y_prev = y[i]; // Make sure we store the value before overwriting
            y[i] = y[i] + adt * d2y_dx2; // Forward Euler
        }

        // Left boundary by forward difference
        let flux_left = 0.0;
        y[0] = y[0] + adt * (2.0 * (dy_dx_left - flux_left)/dx);

        // Right boundary by backward difference
        let flux_right = 0.0;
        y[n-1] = y[n-1] + adt * (2.0 * (dy_dx_right - flux_right)/dx);
    }
}