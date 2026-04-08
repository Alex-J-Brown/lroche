use std::cmp;
use std::f64::consts::PI;

pub fn numface(
    nlat: u32,
    infill: bool,
    thelo: f64,
    thehi: f64,
    nlatfill: u32,
    nlngfill: u32,
) -> u32 {
    let mut nface: u32 = 0;
    let mut theta: f64;
    let mut dphi: f64;
    let dtheta: f64 = PI / nlat as f64;

    if infill {
        let nl1: u32 = (thelo / dtheta).ceil() as u32;
        for i in 0..nl1 {
            theta = thelo * (i as f64 + 0.5) / nl1 as f64;
            dphi = dtheta / theta.sin();
            nface += cmp::max(16, (2. * PI / dphi) as u32);
        }
        let nl2: u32 = ((1. + nlatfill as f64) * (thehi - thelo) / dtheta).ceil() as u32;
        for i in 0..nl2 {
            theta = thelo + (thehi - thelo) * (i as f64 + 0.5) / nl2 as f64;
            dphi = dtheta / theta.sin() / (1. + nlngfill as f64);
            nface += cmp::max(8, (PI / dphi) as u32);
        }
        let nl3: u32 = ((thehi - thelo) / dtheta).ceil() as u32;
        for i in 0..nl3 {
            theta = thelo + (thehi - thelo) * (i as f64 + 0.5) / nl3 as f64;
            dphi = dtheta / theta.sin();
            nface += cmp::max(8, (PI / dphi) as u32);
        }
        let nl4: u32 = ((PI - thehi) / dtheta).ceil() as u32;
        for i in 0..nl4 {
            theta = thehi + (PI - thehi) * (i as f64 + 0.5) / nl4 as f64;
            dphi = dtheta / theta.sin();
            nface += cmp::max(16, (2. * PI / dphi) as u32)
        }
    } else {
        let nlat: u32 = (PI / dtheta).ceil() as u32;
        for i in 0..nlat {
            theta = PI * (i as f64 + 0.5) / nlat as f64;
            dphi = dtheta / theta.sin();
            nface += cmp::max(16, (2. * PI / dphi) as u32)
        }
    }
    nface
}
