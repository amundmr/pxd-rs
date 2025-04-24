struct CurrentFunction {
    pub t: Vec<f64>;  // Time in seconds
    pub i: Vec<f64>;  // Current in ampere
}

impl CurrentFunction {
    /// Generates a CurrentFunction that is a full cycle given:
    /// * current in ampere
    /// * half cycle time in seconds
    /// * timestep in seconds
    pub fn cycle_from_current(current: f64, time: f64, dt: f64) -> CurrentFunction {
        let n_steps = (2.0 * time / dt).ceil() as usize;
        let mut t = Vec::with_capacity(n_steps);
        let mut i = Vec::with_capacity(n_steps);

        for step in 0..n_steps {
            let time_val = step as f64 * dt;
            t.push(time_val);

            let current_val = if time_val < time {
                current
            } else {
                -current
            };
            i.push(current_val);
        }

        CurrentFunction { t, i }
    }

    /// Returns the current given a specific time
    pub fn current(&self, time: f64) -> f64 {
        // Find the index where the time is less than or equal to the provided time
        for (i, &t_val) in self.t.iter().enumerate() {
            if t_val >= time {
                return self.i[i]; // Return the corresponding current
            }
        }
        // If no time is found (i.e., time is greater than the last element), return the last current
        *self.i.last().unwrap()
    }
}