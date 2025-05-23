pub mod numerical_methods {


    pub fn forward_time_centered_space_linear(y: &mut [f64], dx: f64, dt: f64, a: f64, flux: f64) {
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

        // Left boundary by central difference (before y[1] is iterated)
        // let dy_dx_left = (y[1] - y[0]) / dx;
        let dy_dx_left: f64 = -flux / a;
        let left_ghost: f64 = y[0] + dx * dy_dx_left;
        let d2y_dx2_left: f64 = (y[1] - 2.0 * y[0] + left_ghost) / dx2;

        // Right boundary by central difference (before y[n-2] is iterated)
        // let dy_dx_right = (y[n - 1] - y[n - 2]) / dx;
        let dy_dx_right: f64 = flux / a;
        let right_ghost: f64 = y[n-1] + dx * dy_dx_right;
        let d2y_dx2_right: f64 = (right_ghost - 2.0 * y[n-1] + y[n-2]) / dx2;

        // Interior points
        for i in 1..n - 1 {
            let d2y_dx2 = (y_prev - 2.0 * y[i] + y[i + 1]) / dx2; // Central difference
            y_prev = y[i]; // Make sure we store the value before overwriting
            y[i] = y[i] + adt * d2y_dx2; // Forward Euler
        }

        // Left boundary by forward euler
        y[0] = y[0] + adt * d2y_dx2_left;

        // Right boundary by forward euler
        y[n - 1] = y[n - 1] + adt * d2y_dx2_right;
    }

    pub fn ftcs_stable(dt: f64, dx: f64, alpha: f64) -> bool {
        // Returns stability (bool) of forward time centered space method
        // for the heat equation (equaling fickian diffusion)
        // given the timestep, x-step and heat transfer coefficient (diffusion coeff).
        println!("dt: {:?}, dx: {:?}, alpha: {:?}", dt, dx, alpha);
        dt <= dx * dx / (2.0 * alpha)
    }

    pub fn forward_time_centered_space_radial(
        y: &mut [f64],
        dr: f64,
        dt: f64,
        a: f64,
        r: f64,
        flux: f64,
    ) {
        // Forward Time Centered Space (FTCS) is an integration method that utilises
        // the finite difference method (FDM). Refer to https://en.wikipedia.org/wiki/FTCS_scheme.
        // It uses forward euler in time, and central difference in space.
        // This is for spherical symmetry (e.g. lithium in particle)
        //
        // borrows mutable y-vector, r-step, timestep, the diffusion coefficient, and steps the y-vector forward by one timestep.
        let n: usize = y.len();

        let dr2 = dr * dr;
        let adt = a * dt;

        let mut y_prev: f64 = y[0];

        // Boundaries found (before y is iterated) using ghost points, central difference and Neumann boundary conditions

        // Spherical symmetry yields dy/dr=0 at r=0
        let dy_dr_center: f64 = 0.0;
        let center_ghost: f64 = y[0] + dr * dy_dr_center;
        let d2y_dr2_center = (y[1] - 2.0 * y[0] + center_ghost) / dr2;

        // The surface takes all the flux, yo
        let dy_dr_surface: f64 = -flux / a;
        let surface_ghost: f64 = y[n - 1] + dr * dy_dr_surface;
        let d2y_dr2_surface = (surface_ghost - 2.0 * y[n - 1] + y[n - 2]) / dr2;

        // Interior points
        for i in 1..n - 1 {
            let rn: f64 = i as f64 * dr; // r is the distance from the center of the particle
            let d2y_dr2 = (y_prev - 2.0 * y[i] + y[i + 1]) / dr2; // Central difference
            let dy_dr = (y[i + 1] - y[i - 1]) / (2.0 * dr); // Central difference
            y_prev = y[i]; // Make sure we store the value before overwriting
            y[i] = y[i] + adt * (d2y_dr2 + 2.0 / rn * dy_dr); // Forward Euler
        }

        // Left boundary by fwd euler
        y[0] = y[0] + adt * (d2y_dr2_center); //  + 2.0/r * dy_dr_center is ommitted since divisio by zero is naughty

        // Right boundary by fwd euler
        y[n - 1] = y[n - 1] + adt * (d2y_dr2_surface + 2.0 / r * dy_dr_surface);
        // r is the surface so need to do n*dr here
    }
}

pub mod utils {
    pub fn arcsinh(x: f64) -> f64 {
        // Returns the inverse hyperbolic sine of x
        // arcsinh(x) = ln(x + sqrt(x^2 + 1))
        (x + (x * x + 1.0).sqrt()).ln()
    }
}