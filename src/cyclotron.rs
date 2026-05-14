use std::f64::consts::{FRAC_PI_2, FRAC_PI_4};


#[derive(Debug, Clone, Copy)]
pub struct Cyclotron {
    half_apex_angle: f64,
    angular_width: f64,
}

impl Cyclotron {
    // Default.
    pub fn new() -> Self {
        Self {
            half_apex_angle: 0.0,
            angular_width: 0.0,
        }
    }

    // Standard constructor
    pub fn with_params(
        half_apex_angle: f64,
        angular_width: f64,
    ) -> Self {
        Self {
            half_apex_angle,
            angular_width,
        }
    }


    /// Computes I(mu)
    pub fn i(&self, theta: f64) -> f64 {
        let phi: f64 = (FRAC_PI_2 / self.angular_width) * (theta - self.half_apex_angle + FRAC_PI_4) - ((FRAC_PI_2 / self.angular_width) - 1.0) * FRAC_PI_4;
        if phi < 0.0 || phi > FRAC_PI_2 {
            return 0.0
        }
        let (sinp, cosp) = phi.sin_cos();
        let i = sinp*cosp;
        if i < 0.0 {
            return 0.0
        }
        i
    }
}

impl Default for Cyclotron {
    fn default() -> Self {
        Self::new()
    }
}
