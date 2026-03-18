// use std::panic;
// use std::f64::consts::{FRAC_PI_2, TAU};
// use crate::model::{Point, Model, Etype};
// use crate::roche::{self, RocheContext, Star};
// use crate::set_star_grid::{self, star_eclipse};
// use crate::vec3::Vec3;

//
// set_bright_spot_grid sets up the elements needed to define the bright spot
// at the edge of the disc. This is modelled as two straight lines of elements
// coincident in position, but with one oriented parallel to the plane of the
// disc and the other tilted by an arbitrary angle. This is to allow the
// elements to make a hump while allowing flexibility in the step heights of
// the bright-spot ingress and egress features. The brightness of the spot
// along the lines is modelled as a function of distance x from the start of
// each line as (x/l)**p1*exp(-(x/l)**p2) where l is a scale length and p1 and
// p2 are power law exponents. The maximum of this occurs at x =
// l*(p1/p2)**(1/p2), and is set to be the intersection of the gas stream with
// the defined spot radius from the white dwarf.
//
// \param mdl       Model object defining parameters
// \param spot      array representing the spot
// \exception Exceptions are thrown if the specified radii over-fill the Roche lobes.
//


// pub fn set_bright_spot_grid(model: &Model) -> Vec<Point> {

//     let (mut r1, mut r2) = model.get_r1r2();

//     let roche_context1 = RocheContext{q: model.q.value, star: Star::Primary, spin: model.spin1.value};
//     let roche_context2 = RocheContext{q: model.q.value, star: Star::Secondary, spin: model.spin2.value};

//     let rl1 = roche_context1.x_l1_asyncronous();
//     if r1 < 0.0 {
//         r1 = rl1;
//     } else if r1 > rl1 {
//         panic!("Primary star is larger than its Roche Lobe!");
//     }

//     let rl2 = 1.0 - roche_context2.x_l1_asyncronous();
//     if r2 < 0.0 {
//         r2 = rl2;
//     } else if r2 > rl2 {
//         panic!("Secondary star is larger than its Roche Lobe!");
//     }

//     let ffac1: f64 = r1/rl1;
//     let ffac2: f64 = r2/rl2;

//     let (bspot, v) = roche::strinit(model.q.value);
//     roche::stradv(model.q.value, &mut bspot, &mut v, model.radius_spot.value, 1.0e-10, 1.0e-3);

//     // Now measure bright-spot angle relative to tangent to disc edge so we need
//     // to add 90 + angle of bright-spot to the input value.
//     let theta: f64 = bspot.y.atan2(bspot.x) + FRAC_PI_2 + model.angle_spot.value.to_radians();

//     let alpha: f64 = model.yaw_spot.value.to_radians();
//     let tilt: f64 = model.tilt_spot.value.to_radians();

//     // The direction of the line of elements is set by angle_spot, but the
//     // beaming direction adds in yaw_spot as well.

// }

