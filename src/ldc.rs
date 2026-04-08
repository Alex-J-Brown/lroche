use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum LDCType {
    Poly,
    Claret,
}

#[derive(Debug, Clone, Copy)]
pub struct LDC {
    ldc1: f64,
    ldc2: f64,
    ldc3: f64,
    ldc4: f64,
    mucrit: f64,
    ltype: LDCType,
}

impl LDC {
    // Default. Sets all to zero and type to POLY.
    pub fn new() -> Self {
        Self {
            ldc1: 0.0,
            ldc2: 0.0,
            ldc3: 0.0,
            ldc4: 0.0,
            mucrit: 0.0,
            ltype: LDCType::Poly,
        }
    }

    // Standard constructor
    pub fn with_params(
        ldc1: f64,
        ldc2: f64,
        ldc3: f64,
        ldc4: f64,
        mucrit: f64,
        ltype: LDCType,
    ) -> Self {
        Self {
            ldc1,
            ldc2,
            ldc3,
            ldc4,
            mucrit,
            ltype,
        }
    }

    /// Computes I(mu)
    pub fn imu(&self, mu: f64) -> f64 {
        if mu <= 0.0 {
            0.0
        } else {
            let mu = mu.min(1.0);
            let ommu = 1.0 - mu;
            let mut im = 1.0;

            match self.ltype {
                LDCType::Poly => {
                    im -= ommu
                        * (self.ldc1 + ommu * (self.ldc2 + ommu * (self.ldc3 + ommu * self.ldc4)));
                }
                LDCType::Claret => {
                    im -= self.ldc1 + self.ldc2 + self.ldc3 + self.ldc4;
                    let msq = mu.sqrt();
                    im +=
                        msq * (self.ldc1 + msq * (self.ldc2 + msq * (self.ldc3 + msq * self.ldc4)));
                }
            }

            im
        }
    }

    /// To help applying mucrit
    pub fn see(&self, mu: f64) -> bool {
        mu > self.mucrit
    }
}

impl Default for LDC {
    fn default() -> Self {
        Self::new()
    }
}
