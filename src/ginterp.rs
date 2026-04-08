pub struct Ginterp {
    // Start phase of coarse grid 0 -- 0.5
    pub phase1: f64,
    // End phase of coarse grid 0 -- 0.5
    pub phase2: f64,
    // Scale factor star 1 at phase1
    pub scale11: f64,
    // Scale factor star 1 at 1-phase1
    pub scale12: f64,
    // Scale factor star 2 at -phase2
    pub scale21: f64,
    // Scale factor star 2 at phase2
    pub scale22: f64,
}

impl Ginterp {
    // returns scale factor for star1 at given phase
    pub fn scale1(&self, phase: f64) -> f64 {
        // assume coarse grid outside -phase1 to +phase1
        let pnorm: f64 = phase - phase.floor();
        if pnorm <= self.phase1 || pnorm >= 1.0 - self.phase1 {
            1.0
        } else {
            (self.scale11 * (1.0 - self.phase1 - pnorm) + self.scale12 * (pnorm - self.phase1))
                / (1.0 - 2.0 * self.phase1)
        }
    }

    // Returns scale factor for star 2 at a given phase
    pub fn scale2(&self, phase: f64) -> f64 {
        let pnorm: f64 = phase - phase.floor();
        if pnorm >= self.phase2 || pnorm <= 1.0 - self.phase2 {
            1.0
        } else if pnorm < 0.5 {
            (self.scale22 * (self.phase2 - pnorm) + self.scale12 * (pnorm + self.phase2))
                / (2.0 * self.phase2)
        } else {
            (self.scale21 * (1.0 + self.phase2 - pnorm)
                + self.scale22 * (pnorm - 1.0 + self.phase2))
                / (2.0 * self.phase2)
        }
    }

    pub fn interp_type(&self, phase: f64) -> i32 {
        let pnorm: f64 = phase - phase.floor();
        if pnorm <= self.phase1 || pnorm >= 1.0 - self.phase1 {
            // coarse grid for star 2, fine for star 1
            1
        } else if (pnorm > self.phase1 && pnorm < self.phase2)
            || (pnorm > 1.0 - self.phase2 && pnorm < 1.0 - self.phase1)
        {
            // coarse grid for both stars
            2
        } else {
            // coarse grid for star 1, fine for star 2
            3
        }
    }
}
