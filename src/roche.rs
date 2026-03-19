use std::panic;
use std::f64::consts::{PI, TAU};
use std::cmp;
use bulirsch::{self, Integrator};
// use numpy::ndarray;
use ndarray;

use crate::constants::{C, H, K};
use crate::vec3::Vec3;


pub struct Xy {
    pub x: f64,
    pub y: f64,
}

#[derive(Debug, PartialEq, Eq, Clone, Copy)]
pub enum Star {
    Primary,
    Secondary,
}

pub struct RocheContext {
    pub q: f64,
    pub star: Star,
    pub spin: f64,
    pub xl1: f64,
}

impl RocheContext {

    pub fn new(q: f64, star: Star, spin: f64) -> Self {
        if q <= 0. {
            panic!("q = {} <= 0", q);
        }
        let xl1: f64 = match star {
            Star::Primary => xl11(q, spin),
            Star::Secondary => xl12(q, spin),
        };
        Self { q, star, spin, xl1 }
    }

    pub fn potential(&self, earth: &Vec3, p:&Vec3, lam: f64) -> f64 {
        rpot_val(self.q, self.star, self.spin, earth, p, lam)
    }


    pub fn gradient(&self, earth: &Vec3, p: &Vec3, lam: f64) -> (f64, f64) {
        let (dp, dl) = rpot_grad(self.q, self.star, self.spin, earth, p, lam);
        (dp, dl)
    }


    pub fn potential_grad(&self, earth: &Vec3, p: &Vec3, lam: f64) -> (f64, f64, f64) {
        let f = self.potential(earth, p, lam);
        let (dp, dl) = rpot_grad(self.q, self.star, self.spin, earth, p, lam);
        (f, dp, dl)
    }


    pub fn ref_sphere(&self, ffac: f64) -> (f64, f64) {
        let tref: f64;
        let rref: f64;
        let pref: f64;
        if self.star == Star::Primary {
            tref = self.xl1;
            rref = tref * 1.0_f64.min(1.001*ffac);
            pref = rpot1(self.q, self.spin, &Vec3 { x:ffac*tref, y: 0.0, z: 0.0 });
            (rref, pref)
        } else if self.star == Star::Secondary {
            tref = 1. - self.xl1;
            rref = tref * 1.0_f64.min(1.001*ffac);
            pref = rpot2(self.q, self.spin, &Vec3 { x:1. - ffac*tref, y: 0.0, z: 0.0 });
            (rref, pref)
        } else {
            panic!("star is not an instance of Star")
        }
    }


    pub fn fblink(&self, ffac: f64, acc: f64, earth: &Vec3, p: &Vec3) -> Result<bool, &'static str> {

        let (rref, pref) = self.ref_sphere(ffac);

        let cofm: Vec3 = match self.star {
            Star::Primary => Vec3::cofm1(),
            Star::Secondary => Vec3::cofm2(),
        };
        
        // First compute the multipliers cutting the reference sphere (if any)
        let mut lam1 = 0.0;
        let mut lam2 = 0.0;
        if !sphere_eclipse_vector(earth, p, &cofm, rref, &mut lam1, &mut lam2) {
            return Ok(false);
        }
        if lam1 == 0.0 {
            return Ok(true);
        }

        // Create function objects for 1D minimisation in lambda direction
        let func = |lam: f64| {
            self.potential(earth, p, lam)
        };

        // Now try to bracket a minimum. We just crudely compute function at regularly spaced intervals filling in the
        // gaps until the step size between the points drops below the threshold. Take every opportunity to jump out early
        // either if the potential is below the threshold or if we have bracketed a minimum.
        let mut nstep: i32 = 1;
        let mut step: f64 = lam2 - lam1;

        let mut f1: f64 = 0.0;
        let mut f2: f64 = 0.0;
        let mut flam: f64 = 1.0;
        let mut lam: f64 = lam1;

        while step > acc {

            lam = lam1 + step/2.0;

            for _ in 0..nstep {

                flam = func(lam);
                if flam <= pref {
                    return Ok(true);
                }

                // Calculate these as late as possible because they may often not be needed
                if nstep == 1 {
                    f1 = func(lam1);
                    f2 = func(lam2);
                }

                if flam < f1 && flam < f2 {
                    break;
                }

                lam += step;
            }
            if flam < f1 && flam < f2 {
                break;
            }
            step /= 2.0;
            nstep *= 2;
        }

        if flam < f1 && flam < f2 {

            // OK, minimum bracketted, so finally pin it down accurately
            // Possible that multiple minima could cause problems but I have
            // never seen this in practice.
            let dfunc = |lam: f64| {
                let (_dp, dl) = self.gradient(earth, p, lam);
                dl
            };

            let (_xmin, flam) = dbrent(lam1, lam, lam2, |x| func(x), |x| dfunc(x), acc, true, pref)?;

            Ok(flam < pref)
        } else {
            // Not bracketted even after a detailed search, and we have not jumped 
	        // out either, so assume no eclipse
            Ok(false)
        }


    }


    pub fn linmin(&self, cosi: f64, sini: f64, p: &Vec3, phi: &mut f64, lam: &mut f64, mut dphi: f64, mut dlam: f64, phi1: f64, phi2: f64, lam1: f64, lam2: f64, pref: f64, acc: f64) -> (f64, bool) {

        let mut jammed = false;

        let make_func = |phi0: f64, lam0: f64, dphi0: f64, dlam0: f64| {
            
            let func = move |x: f64| {
                let earth: Vec3 = set_earth(cosi, sini, phi0 + dphi0 * x);
                self.potential(&earth, p, lam0 + dlam0 * x)
            };

            let dfunc = move |x: f64| {
                let earth: Vec3 = set_earth(cosi, sini, phi0 + dphi0 * x);
                let (dp, dl) = self.gradient(&earth, p, lam0 + dlam0 * x);
                dp * dphi0 + dl * dlam0
            };

            (func, dfunc)
        };

        let (mut func, mut dfunc) = make_func(*phi, *lam, dphi, dlam);

        // --- determine xmax from boundaries ---
        let mut xmax: f64 = 1.0e30;
        let mut nbound: i32 = 0;

        let mut check = |bound: f64, val: f64, d: f64, id: i32| {
            if d != 0.0 {
                let x: f64 = (bound - val) / d;
                if x > 0.0 && x < xmax {
                    xmax = x;
                    nbound = id;
                }
            }
        };

        check(phi1, *phi, dphi, 1);
        check(phi2, *phi, dphi, 2);
        check(lam1, *lam, dlam, 3);
        check(lam2, *lam, dlam, 4);

        // --- initial bracketing ---
        // Now the aim is to bracket the minimum, while accounting for the maximum
        // possible step so that we can then apply dbrent.
        let xa: f64 = 0.0;
        let fa: f64 = func(xa);

        let mut xb: f64 = 1e-8 * xmax;
        let mut fb: f64 = func(xb);

        // for _ in 0..7 {
        //     if fb < fa && xa != xb {
        //         break;
        //     }
        //     xb *= 10.0;
        //     fb = func(xb);
        // }
        let mut nten: i32 = 0;
        const NTEN: i32 = 7;
        while (fb >= fa || xa == xb) && nten < NTEN {
            nten += 1;
            xb *= 10.0;
            fb = func(xb);
        }

        if fb <= pref {
            *phi += dphi * xb;
            *lam += dlam * xb;
            return (fb, jammed);
        }

        if fb >= fa {
            // Let's hope that we have not stepped past the minimum without
            // knowing it
            return (fa, jammed);
        }

        // --- bracket other side ---
        // OK, so fb < fa so we are heading downhill at least. Now try
        // to find other side starting from xb, looking for a point when
        // we go up or dip below the critical potential
        let mut bracketted = false;
        
        let mut xc: f64 = 0.0;
        let mut fc: f64 = 0.0;

        let xbold: f64 = xb;
        let fbold: f64 = fb;
        xmax -= xb;

        const NTRY: i32 = 5;

        let xmax_rem: f64 = xmax - xb;

        for n in 1..=NTRY {
            xc = xbold + xmax_rem * (n as f64) / (NTRY as f64);
            fc = func(xc);

            if fc <= pref {
                *phi += dphi * xc;
                *lam += dlam * xc;
                return (fc, jammed);
            }

            if fc < fb {
                xb = xc;
                fb = fc;
            } else {
                bracketted = true;
                break;
            }
        }

        jammed = false;

        // --- boundary logic ---
        if !bracketted {

            let dc: f64 = dfunc(xc);

            if dc > 0.0 {

                // We have crashed into the end stop without crossing the minimum
                // but the derivative says that we are going up. Damn!  Go back to
                // old xb and check that derivative was going down there it really
                // ought to have been ...
                xb = xbold;
                let db: f64 = dfunc(xb);
                if db < 0.0 {

                    // OK, let's try to zero in on the point at which the derivative
                    // switches sign.
                    let mut xm: f64 = (xb+xc)/2.0;
                    while (xc-xb) > 1.0e-8*xc {
                        xm = (xb+xc)/2.0;
                        if dfunc(xm) > 0.0 {
                            xc = xm;
                        } else {
                            xb = xm;
                        }
                    }
                    let fm: f64 = func(xm);
                    if fm <= pref {
                        *phi += dphi*xm;
                        *lam += dlam*xm;
                        return (fm, jammed);
                    }
                    if fm < fc && fm < fbold {
                        xb = xm;
                        bracketted = true;
                    } else {
                        panic!("Roche::linmin: failed to bracket minimum, error 3");
                    }
                } else {
                    panic!("Roche::linmin: failed to bracket minimum, error 1");
                }
            } else {
                // We are trapped on a boundary; re-define line minimisation functions.
                jammed = true;
                *phi += dphi*xmax;
                *lam += dlam*xmax;
                xmax = 1.0;
                (dphi, dlam) = rpot_grad(self.q, self.star, self.spin, &set_earth(cosi, sini, *phi), p, *lam);

                match nbound {
                    1 => {
                        *phi = phi1;
                        dphi = 0.0;
                        if dlam > 0.0 {
                            dlam = lam1 - *lam;
                        } else {
                            dlam = lam2 - *lam;
                        }
                    }
                    2 => {
                        *phi = phi2;
                        dphi = 0.0;
                        if dlam > 0.0 {
                            dlam = lam1 - *lam;
                        } else {
                            dlam = lam2 - *lam;
                        }
                    }
                    3 => {
                        *lam = lam1;
                        dlam = 0.0;
                        if dphi > 0.0 {
                            dphi = phi1 - *phi;
                        } else {
                            dphi = phi2 - *phi;
                        }
                    }
                    4 => {
                        *lam = lam1;
                        dlam = 0.0;
                        if dphi > 0.0 {
                            dphi = phi1 - *phi;
                        } else {
                            dphi = phi2 - *phi;
                        }
                    }
                    _ => {}
                }
                (func, dfunc) = make_func(*phi, *lam, dphi, dlam);

                // Again try to bracket minimum
                nten = 0;
                xb = 1.0e-6;
                fb = func(xb);

                while fb >= fa && nten < NTEN {
                    xb *= 10.0;
                    nten += 1;
                    fb = func(xb);
                }
                if fb <= pref {
                    *phi += dphi*xb;
                    *lam += dlam*xb;
                    return (fb, jammed);
                }
                if fb >= fa {
                    return (fa, jammed);
                }

                // Now the bracketting steps.
                bracketted = false;
                for n in 1..=NTRY {
                    xc = xmax*n as f64/NTRY as f64;
                    fc = func(xc);
                    if fc <= pref {
                        *phi += dphi*xc;
                        *lam += dlam*xc;
                        return (fc, jammed)
                    } else {
                        bracketted = true;
                        break;
                    }
                }

                if !bracketted {

                    if dfunc(xc) > 0.0 {
                        panic!("Roche::linmin: failed to bracket minimum, error 2");
                    } else {
                        *phi += dphi*xmax;
                        *lam += dlam*xmax;
                        return (fc, jammed);
                    }
                }

            }

            // if dc <= 0.0 {

            //     jammed = true;

            //     *phi += dphi * xmax;
            //     *lam += dlam * xmax;
                
            //     let earth = set_earth(cosi, sini, *phi);
            //     let (ndphi, ndlam) = self.gradient(&earth, p, *lam);

            //     dphi = ndphi;
            //     dlam = ndlam;

            //     match nbound {
            //         1 => {
            //             *phi = phi1;
            //             dphi = 0.0;
            //         }
            //         2 => {
            //             *phi = phi2;
            //             dphi = 0.0;
            //         }
            //         3 => {
            //             *lam = lam1;
            //             dlam = 0.0;
            //         }
            //         4 => {
            //             *lam = lam2;
            //             dlam = 0.0;
            //         }
            //         _ => {}
            //     }

            //     let (f, df) = make_func(*phi, *lam, dphi, dlam);
            //     func = f;
            //     dfunc = df;

            //     xb = 1e-6;
            //     fb = func(xb);

            //     for _ in 0..7 {
            //         if fb < fa {
            //             break;
            //         }
            //         xb *= 10.0;
            //         fb = func(xb);
            //     }

            //     if fb <= pref {
            //         *phi += dphi * xb;
            //         *lam += dlam * xb;
            //         return (fb, jammed);
            //     }

            //     if fb >= fa {
            //         return (fa, jammed);
            //     }
                
            // } else {
            //     panic!("linmin failed to bracket minimum");
            // }
        }

        // --- Brent refinement ---

        let xacc: f64 = acc / ((TAU*dphi).powi(2) + dlam*dlam).sqrt();

        let (xmin, pmin) = dbrent(xa, xb, xc, &func, &dfunc, xacc, true, pref).unwrap();

        *phi += dphi * xmin;
        *lam += dlam * xmin;

        (pmin, jammed)
    }


    pub fn pot_min(&self, cosi: f64, sini: f64, p: &Vec3, phi1: f64, phi2: f64, lam1: f64, lam2: f64, rref: f64, pref: f64, acc: f64, phi: &mut f64, lam: &mut f64) -> bool {

        *phi = (phi1 + phi2)/2.;
        *lam = (lam1 + lam2)/2.;

        let mut rp: f64 = TAU * *phi;
        let (sinp, cosp) = rp.sin_cos();

        let mut earth: Vec3 = Vec3::new(sini*cosp, -sini*sinp, cosi);
        let mut pot: f64;
        let mut dphi: f64;
        let mut dlam: f64;
        (pot, dphi, dlam) = self.potential_grad(&earth, p, *lam);
        if pot <= pref {
            return true
        };

        let mut gdphi: f64 = -dphi;
        let mut gdlam: f64 = -dlam;
        dphi = gdphi;
        let mut hdphi: f64 = gdphi;
        dlam = gdlam;
        let mut hdlam: f64 = gdlam;

        const ITMAX: i32 = 200;

        let delphi: f64 = self.q/(1.0 + self.q)*((acc*acc)/(rref*rref))/2.0;
        let mut pmin: f64;
        let mut gam: f64;
        let mut dgg: f64;
        let mut gg: f64;
        let mut jammed: bool;
        for _ in 0..ITMAX {

            (pmin, jammed) = self.linmin(cosi, sini, p, phi, lam, dphi, dlam, phi1, phi2, lam1, lam2, pref, acc);

            if pmin <= pref {
                return true
            }
            if jammed || (pmin - pot).abs() < delphi {
                return false
            }

            pot = pmin;
            rp = TAU * *phi;
            let (sinp, cosp) = rp.sin_cos();
            earth.set(sini*cosp, -sini*sinp, cosi);
            (dphi, dlam) = self.gradient(&earth, p, *lam);

            gg = gdphi*gdphi + gdlam*gdlam;
            if gg == 0. {
                return false
            }

            dgg = (dphi+gdphi)*dphi + (dlam+gdlam)*dlam;
            gam = dgg/gg;

            gdphi = -dphi;
            gdlam = -dlam;
            dphi = gdphi + gam*hdphi;
            hdphi = gdphi + gam*hdphi;
            dlam = gdlam + gam*hdlam;
            hdlam = gdlam + gam*hdlam;

        }
        panic!("Too many iterations.")
    }


    pub fn face(&self, direction: Vec3, rref: f64, pref: f64, acc: f64) -> (Vec3, Vec3, f64, f64) {

        let mut pvec: Vec3;
        let mut r: f64;

        let cofm: Vec3 = match self.star {
            Star::Primary => Vec3::cofm1(),
            Star::Secondary => Vec3::cofm2(),
        };

        let rp: fn(f64, f64, &Vec3) -> f64 = match self.star {
            Star::Primary => rpot1,
            Star::Secondary => rpot2,
        };

        let drp: fn(f64, f64, &Vec3) -> Vec3 = match self.star {
            Star::Primary => drpot1,
            Star::Secondary => drpot2,
        };

        let mut tref: f64 = rp(self.q, self.spin, &(cofm + rref*direction));
        if tref < pref {
            panic!("stuff")
        }

        let mut r1: f64 = rref/2.;
        let mut r2: f64 = rref;
        tref = pref + 1.;

        const MAXSEARCH: i32 = 30;
        let mut i: i32 = 0;
        while i < MAXSEARCH && tref > pref {
            r1 = r2/2.;
            tref = rp(self.q, self.spin, &(cofm + r1*direction));
            if tref > pref {
                r2 = r1;
            }
            i+=1;
        }
        if tref > pref {
            panic!("other stuff");
        }

        const MAXCHOP: i32 = 100;
        let mut nchop: i32 = 0;
        while r2 - r1 > acc && nchop < MAXCHOP {
            r = (r1 + r2)/2.;
            pvec = cofm + r*direction;
            if rp(self.q, self.spin, &pvec) < pref {
                r1 = r;
            }else {
                r2 = r;
            }
            nchop += 1;
        }
        if nchop == MAXCHOP {
            panic!("even more stuff");
        }
        r = (r1 + r2)/2.;
        pvec = cofm + r*direction;
        let mut dvec: Vec3 = drp(self.q, self.spin, &pvec);
        let g = dvec.length();
        dvec /= g;
        return (pvec, dvec, r, g)
    }


    pub fn ingress_egress(&self, ffac: f64, iangle: f64, delta: f64, r: &Vec3, ingress: &mut f64, egress: &mut f64) -> bool {
        let rref: f64;
        let pref: f64;
        (rref, pref) = self.ref_sphere(ffac);
        let ri: f64 = iangle.to_radians();
        let (sini, cosi) = ri.sin_cos();

        let cofm: Vec3 = match self.star {
            Star::Primary => Vec3::cofm1(),
            Star::Secondary => Vec3::cofm2(),
        };

        let mut phi1: f64 = 0.0;
        let mut phi2: f64 = 0.0;
        let mut lam1: f64 = 0.0;
        let mut lam2: f64 = 0.0;
        let mut phi: f64 = 0.0;
        let mut lam: f64 = 0.0;

        if sphere_eclipse(cosi, sini, r, &cofm, rref, &mut phi1, &mut phi2, &mut lam1, &mut lam2) {
            
            let acc: f64 = 2.*(2.*TAU*(lam2 - lam1)*delta).sqrt();

            if self.pot_min(cosi, sini, r, phi1, phi2, lam1, lam2, rref, pref, acc, &mut phi, &mut lam) {

                let mut pin: f64 = phi;
                let mut pout: f64 = phi1;
                let mut pmid: f64;

                while (pin - pout).abs() > delta {
                    pmid = (pin + pout)/2.0;
                    if self.fblink(ffac, acc, &set_earth(cosi, sini, pmid), r).unwrap() {
                        pin = pmid;
                    } else {
                        pout = pmid;
                    }
                }
                *ingress = (pin+pout)/2.0;
                *ingress = *ingress - ingress.floor();

                pin = phi;
                pout = phi2;
                while (pin-pout).abs() > delta {
                    pmid = (pin+pout)/2.;
                    if self.fblink(ffac, acc, &set_earth(cosi, sini, pmid), r).unwrap() {
                        pin = pmid;
                    } else {
                        pout = pmid;
                    }
                }
                *egress = (pin+pout)/2.0;
                *egress = *egress - egress.floor();
                if *egress < *ingress {
                    *egress += 1.0;
                }
                return true;
            } else {
                return false;
            }
        } else {
            return false;
        }

    }


    pub fn x_l1(&self) -> f64 {
        xl1(self.q)
    }


    pub fn x_l1_asyncronous(&self) -> f64 {
        match self.star {
            Star::Primary => xl11(self.q, self.spin),
            Star::Secondary => xl12(self.q, self.spin),
        }
    }


    pub fn x_l2(&self) -> f64 {
        xl2(self.q)
    }


    pub fn x_l3(&self) -> f64 {
        xl3(self.q)
    }

}


pub fn sphere_eclipse(cosi: f64, sini: f64, r: &Vec3, c: &Vec3, rsphere: f64, phi1: &mut f64, phi2: &mut f64, lam1: &mut f64, lam2: &mut f64) -> bool {

    let d: Vec3 = *r - *c;

    let pdist: f64 = (d.x*d.x + d.y*d.y).sqrt();
    let bquad: f64 = d.z*cosi - pdist*sini;
    if bquad >= 0. {
        return false;
    }
    let cquad: f64 = d.sqr() - rsphere*rsphere;

    let mut fac: f64 = bquad*bquad - cquad;
    if fac <= 0.0 {
        return false
    }
    fac = fac.sqrt();

    *lam2 = -bquad + fac;
    *lam1 = 0_f64.max(cquad/(*lam2));

    if cquad < 0. {
        *phi1 = 0.;
        *phi2 = 1.
    }else {
        let delta: f64 = ((cosi*d.z + cquad.sqrt())/(sini*pdist)).acos();
        let phi: f64 = d.y.atan2(-d.x);
        *phi1 = (phi - delta)/TAU;
        *phi1 -= phi1.floor();
        *phi2 = *phi1 + 2.*delta/TAU;
    }
    return true;

}


pub fn sphere_eclipse_vector(earth: &Vec3, r: &Vec3, c: &Vec3, rsphere: f64, lam1: &mut f64, lam2: &mut f64) -> bool {

    let d: Vec3 = *r - *c;

    let bquad: f64 = earth.dot(&d);
    if bquad >= 0. {
        return false;
    }
    let cquad: f64 = d.sqr() - rsphere*rsphere;

    let mut fac: f64 = bquad*bquad - cquad;
    if fac <= 0. {
        return false
    }
    fac = fac.sqrt();

    *lam2 = -bquad + fac;
    *lam1 = 0_f64.max(cquad/(*lam2));
    return true;
}


pub fn rpot_val(q: f64, star: Star, spin: f64, earth: &Vec3, p: &Vec3, lam: f64) -> f64 {

    let r: Vec3 = *p + lam* *earth;
    match star {
        Star::Primary => rpot1(q, spin, &r),
        Star::Secondary => rpot2(q, spin, &r),
        }
}


pub fn rpot_val_grad(q: f64, star: Star, spin: f64, earth: &Vec3, p: &Vec3, lam: f64) -> (f64, f64, f64) {

        let r: Vec3 = *p + lam* *earth;
        let d: Vec3 = match star {
            Star::Primary => drpot1(q, spin, &r),
            Star::Secondary => drpot2(q, spin, &r),
        };
        let rpot: f64 = match star {
            Star::Primary => rpot1(q, spin, &r),
            Star::Secondary => rpot2(q, spin, &r),
        };
        let ed: Vec3 = Vec3::new(earth.y, -earth.x, 0.);
        let dphi: f64 = 2.*PI*lam*d.dot(&ed);
        let dlam: f64 = d.dot(earth);

        (rpot, dphi, dlam)
    }


pub fn rpot_grad(
    q: f64,
    star: Star,
    spin: f64,
    earth: &Vec3,
    p: &Vec3,
    lam: f64
) -> (f64, f64) {

        let r: Vec3 = *p + lam * *earth;

        let d: Vec3 = match star {
            Star::Primary => drpot1(q, spin, &r),
            Star::Secondary => drpot2(q, spin, &r),
        };

        let ed: Vec3 = Vec3::new(earth.y, -earth.x, 0.);

        let dphi: f64 = 2.*PI * lam*d.dot(&ed);
        let dlam: f64 = d.dot(earth);

        (dphi, dlam)
}


pub fn dbrent<F, G>(
    ax: f64,
    bx: f64,
    cx: f64,
    func: F,
    dfunc: G,
    acc: f64,
    stopfast: bool,
    fref: f64,
) -> Result<(f64, f64), &'static str>
where
    F: Fn(f64) -> f64,
    G: Fn(f64) -> f64,
{
    const ITMAX: usize = 100;

    let mut a: f64 = ax.min(cx);
    let mut b: f64 = ax.max(cx);

    let mut x: f64 = bx;
    let mut w: f64 = bx;
    let mut v: f64 = bx;

    let mut fx: f64 = func(x);
    let mut fw: f64 = fx;
    let mut fv: f64 = fx;

    if stopfast && fx < fref {
        return Ok((x, fx));
    }

    let mut dx: f64 = dfunc(x);
    let mut dw: f64 = dx;
    let mut dv: f64 = dx;

    let mut e: f64 = 0.0;
    let mut d: f64 = 0.0;

    for _ in 0..ITMAX {

        let xm: f64 = 0.5 * (a + b);
        let tol1: f64 = acc;
        let tol2: f64 = 2.0 * tol1;

        if (x - xm).abs() <= (tol2 - 0.5 * (b - a)) {
            return Ok((x, fx));
        }

        let mut d1: f64;
        let mut d2: f64;

        if e.abs() > tol1 {

            d1 = 2.0 * (b - a);
            d2 = d1;

            if dw != dx {
                d1 = (w - x) * dx / (dx - dw);
            }

            if dv != dx {
                d2 = (v - x) * dx / (dx - dv);
            }

            let u1: f64 = x + d1;
            let u2: f64 = x + d2;

            let ok1: bool = (a - u1) * (u1 - b) > 0.0 && dx * d1 <= 0.0;
            let ok2: bool = (a - u2) * (u2 - b) > 0.0 && dx * d2 <= 0.0;

            let olde: f64 = e;
            e = d;

            if ok1 || ok2 {

                d = if ok1 && ok2 {
                    if d1.abs() < d2.abs() { d1 } else { d2 }
                } else if ok1 {
                    d1
                } else {
                    d2
                };

                if d.abs() <= 0.5 * olde.abs() {

                    let u: f64 = x + d;

                    if (u - a) < tol2 || (b - u) < tol2 {
                        d = tol1.copysign(xm - x);
                    }

                } else {
                    e = if dx >= 0.0 { a - x } else { b - x };
                    d = 0.5 * e;
                }

            } else {

                e = if dx >= 0.0 { a - x } else { b - x };
                d = 0.5 * e;

            }

        } else {

            e = if dx >= 0.0 { a - x } else { b - x };
            d = 0.5 * e;

        }

        let u: f64;
        let fu: f64;

        if d.abs() >= tol1 {

            u = x + d;
            fu = func(u);

            if stopfast && fu < fref {
                return Ok((u, fu));
            }

        } else {

            u = x + tol1.copysign(d);
            fu = func(u);

            if stopfast && fu < fref {
                return Ok((u, fu));
            }

            if fu > fx {
                return Ok((x, fx));
            }
        }

        let du: f64 = dfunc(u);

        if fu <= fx {

            if u >= x { a = x } else { b = x };

            v = w; fv = fw; dv = dw;
            w = x; fw = fx; dw = dx;
            x = u; fx = fu; dx = du;

        } else {

            if u < x { a = u } else { b = u };

            if fu <= fw || w == x {

                v = w; fv = fw; dv = dw;
                w = u; fw = fu; dw = du;

            } else if fu < fv || v == x || v == w {

                v = u; fv = fu; dv = du;

            }
        }
    }

    Err("dbrent: too many iterations")
}


pub fn rpot(q: f64, p: &Vec3) -> f64 {
    if q <= 0. {
        panic!("q = {} <= 0", q);
    }
    let mu: f64 = q/(1. + q);
    let comp: f64 = 1. - mu;
    let x2y2: f64 = p.x*p.x + p.y*p.y;
    let z2: f64 = p.z*p.z;
    let r1sq: f64 = x2y2+z2;
    let r1: f64 = r1sq.sqrt();
    let r2: f64 = (r1sq + 1. - 2.*p.x).sqrt();
    -comp/r1 - mu/r2 - (x2y2 + mu*(mu - 2.*p.x))/2.
}


pub fn rpot1(q: f64, spin: f64, p: &Vec3) -> f64 {
    if q <= 0. {
        panic!("q = {} <= 0", q);
    }
    let mu: f64 = q/(1. + q);
    let comp: f64 = 1. - mu;
    let x2y2: f64 = p.x*p.x + p.y*p.y;
    let z2: f64 = p.z*p.z;
    let r1sq: f64 = x2y2+z2;
    let r1: f64 = r1sq.sqrt();
    let r2: f64 = (r1sq + 1. - 2.*p.x).sqrt();
    -comp/r1 - mu/r2 - spin*spin*x2y2/2. + mu*p.x
}


pub fn rpot2(q: f64, spin: f64, p: &Vec3) -> f64 {
    if q <= 0. {
        panic!("q = {} <= 0", q);
    }
    let mu: f64 = q/(1. + q);
    let comp: f64 = 1. - mu;
    let x2y2: f64 = p.x*p.x + p.y*p.y;
    let z2: f64 = p.z*p.z;
    let r1sq: f64 = x2y2+z2;
    let r1: f64 = r1sq.sqrt();
    let r2: f64 = (r1sq + 1. - 2.*p.x).sqrt();
    -comp/r1 - mu/r2 - spin*spin*(0.5 + 0.5*x2y2 - p.x) - comp*p.x
}


pub fn drpot(q: f64, p: &Vec3) -> Vec3 {
    if q <= 0. {
        panic!("q = {} <= 0", q);
    }
    let r1sq: f64 = p.sqr();
    let r1: f64 = r1sq.sqrt();
    let r2sq: f64 = r1sq + 1. - 2.*p.x;
    let r2: f64 = r2sq.sqrt();
    let mu: f64 = q/(1. + q);
    let mu1: f64 = mu/r2/r2sq;
    let comp: f64 = (1. - mu)/r1/r1sq;
    Vec3::new(comp*p.x + mu1*(p.x - 1.) - p.x + mu,
              comp*p.y + mu1*p.y - p.y,
              comp*p.z + mu1*p.z)
}


pub fn drpot1(q: f64, spin: f64, p: &Vec3) -> Vec3 {
    if q <= 0. {
        panic!("q = {} <= 0", q);
    }
    let r1sq: f64 = p.sqr();
    let r1: f64 = r1sq.sqrt();
    let r2sq: f64 = r1sq + 1. - 2.*p.x;
    let r2: f64 = r2sq.sqrt();
    let mu: f64 = q/(1. + q);
    let mu1: f64 = mu/r2/r2sq;
    let comp: f64 = (1. - mu)/r1/r1sq;
    let ssq: f64 = spin*spin;
    Vec3::new(comp*p.x + mu1*(p.x - 1.) - ssq*p.x + mu,
              comp*p.y + mu1*p.y - ssq*p.y,
              comp*p.z + mu1*p.z)
}


pub fn drpot2(q: f64, spin: f64, p: &Vec3) -> Vec3 {
    if q <= 0. {
        panic!("q = {} <= 0", q);
    }
    let r1sq: f64 = p.sqr();
    let r1: f64 = r1sq.sqrt();
    let r2sq: f64 = r1sq + 1. - 2.*p.x;
    let r2: f64 = r2sq.sqrt();
    let mu: f64 = q/(1. + q);
    let mu1: f64 = mu/r2/r2sq;
    let comp: f64 = (1. - mu)/r1/r1sq;
    let ssq: f64 = spin*spin;
    Vec3::new(comp*p.x + mu1*(p.x - 1.) - ssq*(p.x - 1.) + mu - 1.,
              comp*p.y + mu1*p.y - ssq*p.y,
              comp*p.z + mu1*p.z)
}


pub fn numface(nlat: u32, infill: bool, thelo: f64, thehi: f64, nlatfill: u32, nlngfill: u32) -> u32 {
    let mut nface: u32 = 0;
    let mut theta: f64;
    let mut dphi: f64;
    let dtheta: f64 = PI/nlat as f64;

    if infill {
        let nl1: u32 = (thelo/dtheta).ceil() as u32;
        for i in 0..nl1 {
            theta = thelo*(i as f64 + 0.5)/nl1 as f64;
            dphi = dtheta/theta.sin();
            nface += cmp::max(16, (2.*PI/dphi) as u32);
        }
        let nl2: u32 = ((1. + nlatfill as f64)*(thehi - thelo)/dtheta).ceil() as u32;
        for i in 0..nl2 {
            theta = thelo + (thehi - thelo)*(i as f64 + 0.5)/nl2 as f64;
            dphi = dtheta/theta.sin()/(1. + nlngfill as f64);
            nface += cmp::max(8, (PI/dphi) as u32);
        }
        let nl3: u32 = ((thehi - thelo)/dtheta).ceil() as u32;
        for i in 0..nl3 {
            theta = thelo + (thehi - thelo)*(i as f64 + 0.5)/nl3 as f64;
            dphi = dtheta/theta.sin();
            nface += cmp::max(8, (PI/dphi) as u32);
        }
        let nl4: u32 = ((PI - thehi)/dtheta).ceil() as u32;
        for i in 0..nl4 {
            theta = thehi + (PI - thehi)*(i as f64 + 0.5)/nl4 as f64;
            dphi = dtheta/theta.sin();
            nface += cmp::max(16, (2.*PI/dphi) as u32)
        }
    } else {
        let nlat: u32 = (PI/dtheta).ceil() as u32;
        for i in 0..nlat {
            theta = PI*(i as f64 + 0.5)/nlat as f64;
            dphi = dtheta/theta.sin();
            nface += cmp::max(16, (2.*PI/dphi) as u32)
        }
    }
    nface
}


pub fn set_earth_iangle(iangle: f64, phase: f64) -> Vec3 {
    let iangle_rad: f64 = iangle.to_radians();
    let phase_rad: f64 = 2.*PI*phase;
    let (sini, cosi) = iangle_rad.sin_cos();
    let (sinp, cosp) = phase_rad.sin_cos();
    Vec3{
        x: sini * cosp,
        y: -sini * sinp,
        z: cosi,
    }
}


pub  fn set_earth(cosi: f64, sini: f64, phase: f64) -> Vec3 {
    let phase_rad: f64 = 2.*PI*phase;
    let (sinp, cosp) = phase_rad.sin_cos();
    Vec3{
        x: sini*cosp,
        y: -sini*sinp,
        z: cosi
    }
}


pub fn envelope(rangle: f64, lambda: f64, r1: f64) -> Xy {
    let (sini, cosi) = rangle.sin_cos();
    let (sinl, cosl) = lambda.sin_cos();
    let norm: f64 = (cosi*cosi + sini*sini*cosl*cosl).sqrt();
    Xy { x: (sinl + r1*cosi*sinl/norm), y: (-cosi*cosl - r1*cosl/norm) }
}


pub fn strinit(q: f64) -> (Vec3, Vec3) {
    // strinit sets a particle just inside the L1 point with the 
    // correct velocity as given in Lubow and Shu.
    //
    // \param q mass ratio = M2/M1
    // \param r start position returned
    // \param v start velocity returned
    //
    
    const SMALL: f64 = 1.0e-5;
    let rl1: f64 = xl1(q);
    let mu: f64 = q/(1.0+q);
    let a: f64 = (1.0-mu)/rl1.powi(3)+mu/(1.0-rl1).powi(3);
    let lambda1: f64 = (((a-2.0) + (a*(9.0*a-8.0)).sqrt())/2.0).sqrt();
    let m1: f64 = (lambda1*lambda1-2.0*a-1.0)/2.0/lambda1;

    let r: Vec3 = Vec3::new(rl1-SMALL, -m1*SMALL, 0.0);
    let v: Vec3 = Vec3::new(-lambda1*SMALL, -lambda1*m1*SMALL, 0.0);

    (r, v)

}


pub fn stradv(q: f64, r: &mut Vec3, v: &mut Vec3, rad: f64, acc: f64, smax: f64) -> f64 {

    // stradv advances a particle of given position and velocity until
    // it reaches a specified radius. It then returns with updated position and
    // velocity. It is up to the user not to request a value that cannot be reached.
    //
    // \param q    mass ratio = M2/M1
    // \param r    Initial and final position
    // \param v    Initial and final velocity
    // \param rad  Radius to aim for
    // \param acc  Accuracy with which to place output point at rad.
    // \param smax Largest time step allowed. It is possible that the
    // routine could take such a large step that it misses
    // the point when the stream is inside the requested
    // radius. This allows one to control this. Typical
    // value = 1.e-3.
    //
    // Returns time step taken
    //

    const TMAX: f64 = 10.0;
    let t_next: f64 = 1.0e-2;

    let mut time: f64 = 0.0;

    // let to: f64;
    let mut ro = *r;
    let mut vo = *v;
    
    // Store initial radius
    let rinit: f64 = r.length();
    let mut rnow: f64 = rinit;


    // set up Bulirsch-Stoer integrator
    let system = OrbitalSystem{ q: q };
    let mut integrator = Integrator::default().with_abs_tol(1.0e-8).with_rel_tol(1.0e-8).into_adaptive();
    // Initialise arrays
    let mut y = ndarray::array![r.x, r.y, r.z, v.x, v.y, v.z];
    let mut y_next = ndarray::Array::zeros(y.raw_dim());
    
    let mut yo = y.clone();
    let mut delta_t = t_next.min(smax);
    // Step until radius crossed
    while (rinit > rad && rnow > rad) || (rinit < rad && rnow < rad) {
        ro = *r;
        vo = *v;
        yo = y.clone();
        integrator
            .step(&system, delta_t, y.view(), y_next.view_mut())
            .unwrap();
        y.assign(&y_next);
        r.set(y[0], y[1], y[2]);
        v.set(y[3], y[4], y[5]);
        rnow = r.length();
        time += delta_t;
        
        if time > TMAX {
            panic!("roche::stradv taken too long without crossing given radius.")
        }
    }

    // Now refine by reinitialising and binary chopping until
    // close enough to requested radius.

    let mut lo: f64 = 0.0;
    let mut hi: f64 = delta_t;
    let mut rlo: f64 = ro.length();
    let mut rhi: f64 = rnow;
    let to: f64 = time;

    while (rhi-rlo).abs() > acc {
        delta_t = (lo+hi)/2.0;
        y = yo.clone();
        *r = ro;
        *v = vo;
        time = to;

        integrator
            .step(&system, delta_t, y.view(), y_next.view_mut())
            .unwrap();
        y.assign(&y_next);

        r.set(y[0], y[1], y[2]);
        v.set(y[3], y[4], y[5]);
        rnow = r.length();

        if (rhi > rad && rnow > rad) || (rhi < rad && rnow < rad) {
            rhi = rnow;
            hi = delta_t;
        } else {
            rlo = rnow;
            lo = delta_t;
        }
    }

    time

}



pub fn xl1(q: f64) -> f64 {

    const NMAX: i32 = 1000;
    const EPS: f64 = 1.0e-12;

    if q <= 0. {
        panic!("q = {} <= 0", q);
    }

    let mu: f64 = q/(1. + q);
    let a1: f64 = -1.0 + mu;
    let a2: f64 =  2.0 - 2.*mu;
    let a3: f64 = -1.0 + mu;
    let a4: f64 =  1.0 + 2.*mu;
    let a5: f64 = -2.0 - mu;
    let a6: f64 = 1.0;
    let d1: f64 = 1.0*a2;
    let d2: f64 = 2.0*a3;
    let d3: f64 = 3.0*a4;
    let d4: f64 = 4.0*a5;
    let d5: f64 = 5.0*a6;

    let mut n: i32 = 0;
    let mut xold: f64 = 0.;
    let mut x: f64 = 1./(1. + q);
    let mut f: f64;
    let mut df: f64;
    while n < NMAX && (x - xold).abs() > EPS*x.abs() {
        xold = x;
        f    = x*(x*(x*(x*(x*a6+a5)+a4)+a3)+a2)+a1;
        df   = x*(x*(x*(x*d5+d4)+d3)+d2)+d1;
        x   -= f/df;
        n+=1;
    }
    return x
}


pub fn xl11(q: f64, spin: f64) -> f64 {

    const NMAX: i32 = 1000;
    const EPS: f64 = 1.0e-12;

    if q <= 0. {
        panic!("q = {} <= 0", q);
    }

    let spin_squared: f64 = spin*spin;
    let mu: f64 = q/(1. + q);
    let a1: f64 = -1. + mu;
    let a2: f64 =  2. - 2.*mu;
    let a3: f64 = -1. + mu;
    let a4: f64 =  spin_squared + 2.*mu;
    let a5: f64 = -2.*spin_squared - mu;
    let a6: f64 = spin_squared;
    let d1: f64 = 1.*a2;
    let d2: f64 = 2.*a3;
    let d3: f64 = 3.*a4;
    let d4: f64 = 4.*a5;
    let d5: f64 = 5.*a6;

    let mut n: i32 = 0;
    let mut xold: f64 = 0.;
    let mut x: f64 = 1./(1. + q);
    let mut f: f64;
    let mut df: f64;
    while n < NMAX && (x - xold).abs() > EPS*x.abs() {
        xold = x;
        f    = x*(x*(x*(x*(x*a6+a5)+a4)+a3)+a2)+a1;
        df   = x*(x*(x*(x*d5+d4)+d3)+d2)+d1;
        x   -= f/df;
        n+=1;
    }
    return x
}


#[inline]
pub fn xl12(q: f64, spin: f64) -> f64 {

    const NMAX: i32 = 1000;
    const EPS: f64 = 1.0e-12;

    if q <= 0. {
        panic!("q = {} <= 0", q);
    }

    let spin_squared: f64 = spin*spin;
    let mu: f64 = q/(1. + q);
    let a1: f64 = -1. + mu;
    let a2: f64 =  2. - 2.*mu;
    let a3: f64 = -spin_squared + mu;
    let a4: f64 =  3.*spin_squared + 2.*mu - 2.;
    let a5: f64 = 1. - mu - 3.*spin_squared;
    let a6: f64 = spin_squared;
    let d1: f64 = 1.*a2;
    let d2: f64 = 2.*a3;
    let d3: f64 = 3.*a4;
    let d4: f64 = 4.*a5;
    let d5: f64 = 5.*a6;

    let mut n: i32 = 0;
    let mut xold: f64 = 0.;
    let mut x: f64 = 1./(1. + q);
    let mut f: f64;
    let mut df: f64;
    while n < NMAX && (x - xold).abs() > EPS*x.abs() {
        xold = x;
        f    = x*(x*(x*(x*(x*a6+a5)+a4)+a3)+a2)+a1;
        df   = x*(x*(x*(x*d5+d4)+d3)+d2)+d1;
        x   -= f/df;
        n+=1;
    }
    return x
}


pub fn xl2(q: f64) -> f64 {

    const NMAX: i32 = 1000;
    const EPS: f64 = 1.0e-12;

    if q <= 0. {
        panic!("q = {} <= 0", q);
    }

    let mu = q/(1. + q);
    let a1 = -1. + mu;
    let a2 =  2. - 2.*mu;
    let a3 = -1. - mu;
    let a4 =  1. + 2.*mu;
    let a5 = -2. - mu;
    let a6 = 1.;
    let d1 = 1.*a2;
    let d2 = 2.*a3;
    let d3 = 3.*a4;
    let d4 = 4.*a5;
    let d5 = 5.*a6;

    let mut n: i32 = 0;
    let mut xold: f64 = 0.;
    let mut x: f64 = 1.5;
    let mut f: f64;
    let mut df: f64;
    while n < NMAX && (x - xold).abs() > EPS*x.abs() {
        xold = x;
        f    = x*(x*(x*(x*(x*a6+a5)+a4)+a3)+a2)+a1;
        df   = x*(x*(x*(x*d5+d4)+d3)+d2)+d1;
        x   -= f/df;
        n+=1;
    }
    return x
}


pub fn xl3(q: f64) -> f64 {

    const NMAX: i32 = 1000;
    const EPS: f64 = 1.0e-12;

    if q <= 0. {
        panic!("q = {} <= 0", q);
    }

    let mu = q/(1. + q);
    let a1 = -1. + mu;
    let a2 =  -2. + 2.*mu;
    let a3 = 1. - mu;
    let a4 =  1. + 2.*mu;
    let a5 = -2. - mu;
    let a6 = 1.;
    let d1 = 1.*a2;
    let d2 = 2.*a3;
    let d3 = 3.*a4;
    let d4 = 4.*a5;
    let d5 = 5.*a6;

    let mut n: i32 = 0;
    let mut xold: f64 = 0.;
    let mut x: f64 = -1.;
    let mut f: f64;
    let mut df: f64;
    while n < NMAX && (x - xold).abs() > EPS*x.abs() {
        xold = x;
        f    = x*(x*(x*(x*(x*a6+a5)+a4)+a3)+a2)+a1;
        df   = x*(x*(x*(x*d5+d4)+d3)+d2)+d1;
        x   -= f/df;
        n+=1;
    }
    return x
}


pub fn zeta_rlobe_eggleton(q: f64) -> f64 {
    let q1 = q.powf(1./3.);
    let loneq = (1. + q1).ln();
    (1. + q)/3.*(2.*loneq-q1/(1.+q1))/(0.6*q1*q1*loneq)
}


pub fn dzetadq_rlobe_eggleton(q: f64) -> f64 {
    let q1 = q.powf(1./3.);
    let q2 = q1*q1;
    let opq1 = 1. + q1;
    let loneq = opq1.ln();
    let denom = 0.6*q2 + loneq;
    let numer = 2.*loneq - q1/opq1;
    numer/denom/3. + (1. + q)/3.*((1. + 2.*q1)/3./(q1*opq1).powi(2) - numer*(0.4/q1 + 1./(3.*q2*(1. + q1)))/denom)/denom
}


pub fn planck(wave: f64, temp: f64) -> f64 {

    let fac1: f64 = 2.0e27*H*C;
    let fac2: f64 = 1.0e9*H*C / K;

    let exponent: f64 = fac2/(wave*temp);
    if exponent > 40.0 {
        fac1*-exponent.exp()/(wave*wave*wave)
    } else {
        fac1/exponent.exp_m1()/(wave*wave*wave)
    }
}


pub fn dplanck(wave: f64, temp: f64) -> f64 {

    let fac2: f64 = 1.0e9*H*C / K;

    let exponent: f64 = fac2/(wave*temp);
    exponent/(1.0 - -exponent.exp()) - 3.0
}


pub fn dlpdlt(wave: f64, temp: f64) -> f64 {

    let fac2: f64 = 1.0e9*H*C / K;

    let exponent: f64 = fac2/(wave*temp);
    exponent /(1.0 - (-exponent).exp())
}




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
        if pnorm <= self.phase1 || pnorm >= 1.0-self.phase1 {
            return 1.0;
        } else {
            return (self.scale11*(1.0-self.phase1-pnorm)+self.scale12*(pnorm-self.phase1))/(1.0-2.0*self.phase1);
        }
    }

    // Returns scale factor for star 2 at a given phase
    pub fn scale2(&self, phase: f64) -> f64 {
        let pnorm: f64 = phase - phase.floor();
        if pnorm >= self.phase2 || pnorm <= 1.0-self.phase2 {
            return 1.0;
        } else if pnorm < 0.5 {
            return (self.scale22*(self.phase2-pnorm)+self.scale12*(pnorm+self.phase2))/(2.0*self.phase2);
        } else {
            return (self.scale21*(1.0+self.phase2-pnorm)+self.scale22*(pnorm-1.0+self.phase2))/(2.0*self.phase2);
        }
    }
    
    
    pub fn interp_type(&self, phase: f64) -> i32 {
        let pnorm: f64 = phase - phase.floor();
        if pnorm <= self.phase1 || pnorm >= 1.0-self.phase1 {
            // coarse grid for star 2, fine for star 1
            return 1;
        } else if (pnorm > self.phase1 && pnorm < self.phase2) || (pnorm > 1.0-self.phase2 && pnorm < 1.0-self.phase1){
            // coarse grid for both stars
            return 2;
        } else {
            // coarse grid for star 1, fine for star 2
            return 3;
        }
    }
}


pub fn rocacc(q: f64, r: &Vec3, v: &Vec3) -> (f64, f64, f64) {

    //
    // rocacc calculates and returns the acceleration (in the rotating frame)
    // in a Roche potential of a particle of given position and velocity.
    //
    // \param q mass ratio = M2/M1
    // \param r position, scaled in units of separation.
    // \param v velocity, scaled in units of separation
    //

    let f1: f64 = 1.0 / (1.0+q);
    let f2: f64 = f1*q;

    let yzsq: f64 = r.y*r.y + r.z*r.z;
    let r1sq: f64 = r.x*r.x + yzsq;
    let r2sq: f64 = (r.x-1.0)*(r.x-1.0) + yzsq;
    let fm1: f64 = f1/(r1sq*(r1sq.sqrt()));
    let fm2: f64 = f2/(r2sq*(r2sq.sqrt()));
    let fm3 = fm1+fm2;

    let x: f64 = -fm3*r.x + fm2 + 2.0*v.y + r.x - f2;
    let y: f64 = -fm3*r.y       - 2.0*v.x + r.y;
    let z: f64 = -fm3*r.z;
    (x, y, z)
}


struct OrbitalSystem {
    q: f64,
}

impl bulirsch::System for OrbitalSystem {
    type Float = f64;
    
    fn system(&self, y: bulirsch::ArrayView1<Self::Float>, mut dydt: bulirsch::ArrayViewMut1<Self::Float>) {
        dydt[[0]] = y[[3]];
        dydt[[1]] = y[[4]];
        dydt[[2]] = y[[5]];
        let r = Vec3::new(y[[0]], y[[1]], y[[2]]);
        let v = Vec3::new(y[[3]], y[[4]], y[[5]]);
        (dydt[[3]], dydt[[4]], dydt[[5]]) = rocacc(self.q, &r, &v);
    }
}

