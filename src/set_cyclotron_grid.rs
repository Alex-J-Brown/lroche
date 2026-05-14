use crate::model::Model;
use crate::numface::numface;
use crate::set_star_grid::{band_geometry};
use rayon::prelude::*;
use roche::errors::RocheError;
use roche::{self, Etype, Point, RocheContext, Star, Vec3};
use std::f64::consts::PI;
use std::panic;


pub fn set_cycotron_grid(model: &Model, star: Star, fine: bool) -> Result<Vec<Point>, RocheError> {
    let (mut r1, mut r2) = model.get_r1r2();
    r1 *= model.cyclotron_radfac.value;

    let eclipse: bool = match star {
        Star::Primary => model.eclipse1,
        Star::Secondary => model.eclipse2,
    };
    let nlat: u32 = match (star, fine) {
        (Star::Primary, true) => model.nlat1f,
        (Star::Primary, false) => model.nlat1c,
        (Star::Secondary, true) => model.nlat2f,
        (Star::Secondary, false) => model.nlat2c,
    };


    let roche_context1 = RocheContext::new(model.q.value, Star::Primary, model.spin1.value)?;
    let roche_context2 = RocheContext::new(model.q.value, Star::Secondary, model.spin2.value)?;

    let rl1: f64 = roche_context1.x_l1;
    if r1 < 0.0 {
        r1 = rl1;
    } else if r1 > rl1 {
        panic!("set_star_grid: the primary star is larger than its Roche lobe!");
    }

    let rl2: f64 = 1.0 - roche_context2.x_l1;
    if r2 < 0.0 {
        r2 = rl2;
    } else if r2 > rl2 {
        panic!("set_star_grid: the secondary star is larger than its Roche lobe!");
    }

    if model.glens1 && star == Star::Secondary && model.roche1 && model.eclipse2 {
        panic!(
            "set_star_grid: cannot have gravitational lensing, eclipse and Roche lobe geometry at the same time"
        );
    }


    // Calculate a reference radius and potential for the two stars
    let ffac1: f64 = r1 / rl1;
    let (rref1, pref1) = roche_context1.ref_sphere(ffac1)?;

    let ffac2: f64 = r2 / rl2;

    // Compute latitude range over which extra points will be added. Only enabled
    // when setting the secondary grid and when the grid North pole is the genuine
    // North pole and when r2 > r1

    // Compute number of faces needed
    let nface: u32 = numface(nlat, false, 0.0, 0.0, 0, 0);

    // Generate arrays over the star's face
    let mut cyclotron_grid: Vec<Point> = Vec::with_capacity(nface as usize);

    let acc: f64 = model.delta_phase / 10.0;

    // let mut dirn: Vec3;
    let _posn: Vec3;
    let _dvec: Vec3;
    let _rad: f64;
    let gref: f64;
    // Compute reference gravity value, from the side of the star opposite from the L1 point
    // to ensure a non-zero value. Set to 1 if Roche distortion being ignored.
    if star == Star::Primary && model.roche1 {
        let dirn = Vec3::new(-1.0, 0.0, 0.0);
        (_posn, _dvec, _rad, gref) = roche_context1.face(dirn, rref1, pref1, acc)?;
    } else {
        gref = 1.0;
    }

    // The grid starts at the North pole and ends at the South, proceeding in
    // a series of equi-latitudinal rings.  The polar axis is parallel to the
    // x-axis which points from one star to the other. The North pole is
    // defined to be the point closest to the other star (or the genuine North
    // Pole if npole is true). The angle theta is measured away from the North
    // pole, the angle phi is measured from the Y axis towards the Z axis.

    let dtheta: f64 = PI / (nlat as f64);

    
    add_cyclotron_faces(
        &mut cyclotron_grid,
        0.0,
        PI,
        dtheta,
        0,
        0,
        model.npole,
        star,
        &roche_context1,
        &roche_context2,
        model.iangle.value,
        r1,
        r2,
        rref1,
        model.roche1,
        model.roche2,
        eclipse,
        gref,
        pref1,
        ffac1,
        ffac2,
        model.delta_phase,
    );
    Ok(cyclotron_grid)
}

pub fn add_cyclotron_faces(
    cyclotron_grid: &mut Vec<Point>,
    tlo: f64,
    thi: f64,
    dtheta: f64,
    nlatfill: u32,
    nlngfill: u32,
    npole: bool,
    star: Star,
    roche_context1: &RocheContext,
    roche_context2: &RocheContext,
    iangle: f64,
    r1: f64,
    r2: f64,
    rref1: f64,
    roche1: bool,
    roche2: bool,
    eclipse: bool,
    gref: f64,
    pref1: f64,
    ffac1: f64,
    ffac2: f64,
    delta: f64,
) {
    let ri: f64 = iangle.to_radians();
    let (sini, cosi) = ri.sin_cos();

    // Centre of masses of the stars
    let cofm1 = Vec3::cofm1();
    let cofm2 = Vec3::cofm2();

    // Can afford to be pretty careful on the location of faces as it is a fast computation
    let acc: f64 = delta / 10.0;

    let infill: bool = (nlatfill > 0) || (nlngfill > 0);

    let nlat: usize = ((thi - tlo) / dtheta).ceil() as usize;
    let nlat1: usize = 0;
    let nlat2: usize = 0;

    let bands: Vec<Vec<Point>> = (0..nlat)
        .into_par_iter()
        .map(|nt| {
            let g = band_geometry(
                nt, tlo, thi, dtheta, nlat, nlat1, nlat2, infill, nlngfill, star,
            );

            let mut band = Vec::with_capacity(g.nphi);

            for np in 0..g.nphi {
                let phi = g.phi1 + (g.phi2 - g.phi1) * (np as f64 + 0.5) / g.nphi as f64;

                let (sinp, cosp) = phi.sin_cos();

                let mut dirn = Vec3::new(0.0, 0.0, 0.0);
                if npole {
                    dirn.set(g.sint * cosp, g.sint * sinp, g.cost);
                } else {
                    dirn.set(g.cost, g.sint * cosp, g.sint * sinp);
                }

                let posn: Vec3;
                let dvec: Vec3;
                let rad: f64;
                let gravity: f64;

                let mut lam1: f64 = 0.0;
                let mut lam2: f64 = 0.0;
                let mut ingress: f64 = 0.0;
                let mut egress: f64 = 0.0;

                // Direction is now defined, so calculate radius and thus the
                // position according to whether we are accounting for Roche
                // geometry or not.

                if roche1 {
                    (posn, dvec, rad, gravity) =
                        roche_context1.face(dirn, rref1, pref1, acc).unwrap();
                } else {
                    // Ignore Roche distortion
                    rad = r1;
                    posn = cofm1 + rad * dirn;
                    dvec = dirn;
                    gravity = 1.0;
                }

                // Area, accounting for angle of face
                let area: f64 = ((g.phi2 - g.phi1) / (g.nphi as f64) * rad * g.sint)
                    * ((thi - tlo) / (g.nl as f64) * rad)
                    / dirn.dot(&dvec);

                // Eclipse computation. We calculate whether a point is
                // eclipsed, and, if it is, its ingress and egress
                // phases. Account for spherical or Roche geometry of other
                // star. Since the cyclotron emission may originate from above
                // the photosphere, eclipse phases should also be calculated
                // for the star itself.

                let mut eclipses = Etype::new();

                if eclipse {
                    if roche2
                             && roche_context2
                                .ingress_egress(
                                    ffac2,
                                    iangle,
                                    delta,
                                    &posn,
                                    &mut ingress,
                                    &mut egress,
                                )
                                .unwrap()
                            || (!roche2
                                && roche::sphere_eclipse(
                                    cosi,
                                    sini,
                                    &posn,
                                    &cofm2,
                                    r2,
                                    &mut ingress,
                                    &mut egress,
                                    &mut lam1,
                                    &mut lam2,
                                )) {
                        eclipses.push((ingress, egress));
                    }
                    if roche1
                             && roche_context1
                                .ingress_egress(
                                    ffac1,
                                    iangle,
                                    delta,
                                    &posn,
                                    &mut ingress,
                                    &mut egress,
                                )
                                .unwrap()
                            || (!roche1
                                && roche::sphere_eclipse(
                                    cosi,
                                    sini,
                                    &posn,
                                    &cofm1,
                                    r1,
                                    &mut ingress,
                                    &mut egress,
                                    &mut lam1,
                                    &mut lam2,
                                )) {
                        eclipses.push((ingress, egress));
                    }
                }
                band.push(Point::new(posn, dvec, area, gravity / gref, eclipses));
            }

            band
        })
        .collect();

    cyclotron_grid.clear();

    for band in bands {
        cyclotron_grid.extend(band);
    }
}

pub fn star_eclipse(
    roche_context: &RocheContext,
    r: f64,
    ffac: f64,
    iangle: f64,
    posn: &Vec3,
    delta: f64,
    roche: bool,
    star: Star,
    eclipses: &mut Etype,
) -> Result<(), RocheError> {
    let ri: f64 = iangle.to_radians();
    let (sini, cosi) = ri.sin_cos();
    let cofm = match star {
        Star::Primary => Vec3::cofm1(),
        Star::Secondary => Vec3::cofm2(),
    };
    let mut lam1: f64 = 0.0;
    let mut lam2: f64 = 0.0;
    let mut ingress: f64 = 0.0;
    let mut egress: f64 = 0.0;
    // let mut eclipses = Etype::new();
    if (roche
        && roche_context.ingress_egress(ffac, iangle, delta, posn, &mut ingress, &mut egress)?)
        || (!roche
            && roche::sphere_eclipse(
                cosi,
                sini,
                posn,
                &cofm,
                r,
                &mut ingress,
                &mut egress,
                &mut lam1,
                &mut lam2,
            ))
    {
        eclipses.push((ingress, egress));
    }
    Ok(())
}
