use roche::{self, Point, Vec3};
use std::f64::consts::TAU;
use pyo3::prelude::*;

/// 
/// Grid is a struct to hold the vector of `Vec<Points>` defining a
/// component grid along with methods to act on this.
/// 
#[pyclass(skip_from_py_object)]
#[derive(Clone, Debug)]
pub struct Grid {
    #[pyo3(get)]
    pub points: Vec<Point>,
    pub q: f64,
    pub iangle: f64,
}

impl Grid {

    pub fn new(points: Vec<Point>, q: f64, iangle: f64) -> Self {
        Self {
            points,
            q,
            iangle,
        }
    }
}

#[pymethods]
impl Grid {

    #[pyo3(signature = (phase=None))]
    pub fn area(&self, phase: Option<f64>) -> Vec<f32> {
        
        let area: Vec<f32> = match phase {
            Some(phase) => {
                let mut area: Vec<f32> = vec![];
                let earth = roche::set_earth_iangle(self.iangle, phase);
                for point in &self.points {
                    if earth.dot(&point.direction) > 0.0 && point.is_visible(phase) {
                        area.push(point.area);
                    }
                }
                area
            },
            None => {
                let mut area: Vec<f32> = vec![];
                for point in &self.points {
                    area.push(point.area);
                }
                area
            }
        };

        area
    }

    #[pyo3(signature = (phase=None))]
    pub fn flux(&self, phase: Option<f64>) -> Vec<f32> {
        
        let flux: Vec<f32> = match phase {
            Some(phase) => {
                let mut flux: Vec<f32> = vec![];
                let earth = roche::set_earth_iangle(self.iangle, phase);
                for point in &self.points {
                    if earth.dot(&point.direction) > 0.0 && point.is_visible(phase) {
                        flux.push(point.flux);
                    }
                }
                flux
            },
            None => {
                let mut flux: Vec<f32> = vec![];
                for point in &self.points {
                    flux.push(point.flux);
                }
                flux
            }
        };

        flux
    }

    ///
    /// Projects the grid onto a 2D plane as seen at the model inclination
    /// at the supplied phase with the binary centre of mass as the origin.
    /// 
    /// Arguments
    /// * `phase` - Orbital phase at which to project the grid
    /// 
    /// Returns
    /// (x, y) - Arrays of projected grid point positions
    /// 
    #[pyo3(signature = (phase))]
    pub fn project_2d(&self, phase: f64) -> (Vec<f64>, Vec<f64>) {

        let mut x_arr: Vec<f64> = vec![];
        let mut y_arr: Vec<f64> = vec![];
        let cofm = Vec3::new(self.q/(1.0+self.q), 0.0, 0.0);
        let (sinp, cosp) = (TAU*phase).sin_cos();
        let earth = roche::set_earth_iangle(self.iangle, phase);
        let xsky = Vec3::new(sinp, cosp, 0.0);
        let ysky = earth.cross(&xsky);
        
        for point in &self.points {
            if earth.dot(&point.direction) > 0.0 && point.is_visible(phase) {
                let r = point.position - cofm;
                x_arr.push(r.dot(&xsky));
                y_arr.push(r.dot(&ysky));
            }
        }
        (x_arr, y_arr)
    }
    
}
