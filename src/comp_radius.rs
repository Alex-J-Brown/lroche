use rust_roche::{Star, Vec3, Point};

//*
// comp_radius1 computes a volume-averaged scaled radius for star 1.
//
// \param star1 -- grid for star1
// \return the value of volume-averaged R1/a
///

pub fn comp_radius(star_grid: &Vec<Point>, star: Star) -> f64 {

    let mut vec: Vec3;
    let cofm = match star {
        Star::Primary => Vec3::cofm1(),
        Star::Secondary => Vec3::cofm2(),
    };

    let mut r: f64;
    let mut rcosa: f64;
    let mut sum_solidangle = 0.0;
    let mut sum_volume = 0.0;

    for point in star_grid {
        vec = point.position - cofm;
        r = vec.length();
        rcosa = vec.dot(&point.direction);

        // sum solid angle and 3x volume of all elements
        sum_solidangle += (point.area as f64) * rcosa / (r*r*r);
        sum_volume += (point.area as f64) * rcosa;
    }

    (sum_volume/sum_solidangle).powf(1.0/3.0)
}
