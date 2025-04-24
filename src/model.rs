mod parameters;
mod input;

struct SPMe {
    p: parameters::ParameterSet;
    cf: input::CurrentFunction;
}

impl SPMe {
    /// Returns cell voltage given the particle surface lithium concentration of the negative and positive electrode
    fn u_cell(c_n_s: f64, c_p_s: f64) -> f64 {
        self.p.u_p(c_p_s) - self.p.u_n(c_n_s)
    }

    fn sphere_volume(r: f64) -> f64 {
        4.0/3.0 * self.p.pi*r**3
    }

    /// Returns the change in concentration for a given time, concentration, electrode domain and particle radius
    fn dcdt(&self, time: f64, c: f64, domain: u8, r: f64) -> f64 {

        let mut polarity = 1;
        if domain == 0 { // negative electrode
            polarity = 1;
        } else if domain == 2 { // positive electrode
            polarity = -1;
        }

        polarity * self.cf.current(time) / (self.p.n * self.p.F * self.sphere_volume(r))
    }
}