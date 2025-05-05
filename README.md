# pxd-rs
Rust implementation of the SPMe battery model.

Numerical method is central difference and forward Euler.

## Status

- [x] Empirical open circuit voltage functions
- [x] Fickian diffusion in electrolyte NB: Requires 1ms time-step for 1.7e-10 diffusion coeff. A better timestepping algo is required to allow increased timestep and increased performance.
- [ ] Migration in electrolyte
- [ ] Bruggeman correction
- [x] Fickian diffusion in particles
- [x] Active material volume fraction
- [x] Butler-Volmer kinetics
- [ ] Double layer capacitance


![Current Status](current_status.png)
