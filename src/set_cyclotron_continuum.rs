use roche::constants::EFAC;
use roche::errors::RocheError;
use roche::{Point, Vec3};
use crate::model::Model;

pub fn set_cyclotron_continuum(model: &Model, cyclotron_grid: &mut Vec<Point>) -> Result<(), RocheError> {
    
    let (slong, clong) = model.cyclotron_long.value.to_radians().sin_cos();
    let (slat, clat) = model.cyclotron_lat.value.to_radians().sin_cos();
    let cyclotron_spot: Vec3 = Vec3::new(clat * clong, clat * slong, slat);

    let mut t_cyclotron: f64;
    for point in cyclotron_grid {
        let dist: f64 = (cyclotron_spot.dot(&point.position) / point.position.length())
            .acos()
            .to_degrees();
        let exponent: f64 = -(dist / (model.cyclotron_fwhm.value / EFAC)).powi(2) / 2.0;
        t_cyclotron = model.cyclotron_tcen.value * exponent.exp();
        
        let flux: f32 = ((point.area as f64) * roche::planck(model.wavelength, t_cyclotron)) as f32;
        point.set_flux(flux);
    }
    Ok(())
}