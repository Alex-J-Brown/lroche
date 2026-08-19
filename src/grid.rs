use roche::{self, Point, Vec3};
use std::f64::consts::TAU;
use pyo3::prelude::*;

#[pyclass(skip_from_py_object)]
#[derive(Clone, Debug)]
pub struct Grid {
    #[pyo3(get)]
    pub points: Vec<Point>,
    #[pyo3(get)]
    pub q: f64,
    #[pyo3(get)]
    pub iangle: f64,
}

#[pymethods]
impl Grid {
    #[pyo3(signature = (phase, iangle=None))]
    pub fn project_2d(&self, phase: f64, iangle: Option<f64>) -> (Vec<f64>, Vec<f64>) {
        let iangle = match iangle {
            Some(iangle) => iangle,
            None => self.iangle,
        };
        let mut x_arr: Vec<f64> = vec![];
        let mut y_arr: Vec<f64> = vec![];
        let cofm = Vec3::new(self.q/(1.0+self.q), 0.0, 0.0);
        let (sinp, cosp) = (TAU*phase).sin_cos();
        let earth = roche::set_earth_iangle(iangle, phase);
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
