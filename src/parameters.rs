

#[derive(Debug, Clone, PartialEq)]
pub struct ParameterSet {
    // Max Lithium concentrations in electrodes
    let c_n_max: f64 = 34684.; // mol/m^3
    let c_p_max: f64 = 50060.; // mol/m^3

    // Initial parameters
    let c_n_init: f64 = 1000.; // mol/m^3
    let c_p_init: f64 = 49000.; // mol/m^3

    // Radiuses, calculated from c_max and theoretical capacity of 2Ah
    let r_n: f64 = 0.2009/2.; // m
    let r_p: f64 = 0.01778/2.; // m

    // Constants
    let n: f64 = 1.; // unitless, number of elementary charges of the ion. float type for ease of operating with other floats
    const F: f64 = 96485.3321233100184; // C/mol (=As/mol), 2019 SI revision definition
    let pi: f64 = std::f64::consts::PI;
}

impl ParameterSet {

    ///Returns the negative electrode open circuit voltage given a lithium concentration in mol/m^3
    pub fn u_n(&self, c: f64) -> f64 {
        let x: f64 = c / self.c_max_n;

        let p = [
            1.20912055e+00, 5.62297420e+01, -1.11020020e-01, -2.53458213e-01,
            4.92581391e+01, 1.22046522e-02, 4.73538620e-02, 1.79631246e+01,
            1.75283209e-01, 1.88038929e-02, 3.03255334e+01, 4.66328034e-01,
        ];

        p[0] * (-p[1] * x).exp()
        + p[2]
        - p[3] * (p[4] * (x - p[5])).tanh()
        - p[6] * (p[7] * (x - p[8])).tanh()
        - p[9] * (p[10] * (x - p[11])).tanh()
    }

    /// Returns the positive electrode open circuit voltage given a lithium concentration in mol/m^3
    pub fn u_p(&self, c: f64) -> f64 {
        let x = c / self.c_max_p;

        let p = [
            0.74041974, 4.39107343, 0.03434767, 18.16841489, 0.53463176,
            17.68283504, 14.59709162, 0.28835348, 17.58474971, 14.69911523,
            0.28845641,
        ];

        -p[0] * x
        + p[1]
        - p[2] * (p[3] * (x - p[4])).tanh()
        - p[5] * (p[6] * (x - p[7])).tanh()
        + p[8] * (p[9] * (x - p[10])).tanh()
    }
}