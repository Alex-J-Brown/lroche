use std::collections::HashMap;
use std::f64::consts::TAU;
use rayon::prelude::*;
use pyo3::prelude::*;
use numpy::{PyReadonlyArray1, PyReadwriteArray1};
use pyo3::types::{PyDict, PyDictMethods, PyFloat};
use crate::comp_light::{comp_bright_spot, comp_disc, comp_disc_edge, comp_star1, comp_star2};
use crate::constants::{C, DAY};
use crate::roche::{Ginterp, Star, planck, xl12};
use crate::model::{Entry, Etype, LDC, Model, Point};
use crate::set_disc_continuum::{set_disc_continuum, set_edge_continuum};
use crate::set_disc_grid::{set_disc_edge_grid, set_disc_grid};
use crate::set_bright_spot_grid::set_bright_spot_grid;
use crate::set_star_continuum::set_star_continuum;
use crate::set_star_grid::{disc_eclipse, set_star_grid};

#[pyclass]
pub struct BinaryModel {
    star1_coarse_grid: Vec<Point>,
    star2_coarse_grid: Vec<Point>,
    star1_fine_grid: Vec<Point>,
    star2_fine_grid: Vec<Point>,
    disc_grid: Vec<Point>,
    disc_edge_grid: Vec<Point>,
    bright_spot_grid: Vec<Point>,
    gint: Ginterp,
    rlens1: f64,
    model: Model
}

#[pymethods]
impl BinaryModel {
    
    #[new]
    pub fn new(dict: Bound<'_, PyDict>) -> PyResult<Self> {

        let map = map_from_pydict(dict)?;
        let model = Model::from_map(map).map_err(|e| pyo3::exceptions::PyValueError::new_err(e))?;

        let mut star1_fine_grid = set_star_grid(&model, Star::Primary, true);
        let mut star2_fine_grid = set_star_grid(&model, Star::Secondary, true);
        let mut star1_coarse_grid: Vec<Point>;
        let mut star2_coarse_grid: Vec<Point>;

        let (r1, mut r2) = model.get_r1r2();
        let rl2: f64 = 1.0 - xl12(model.q.value, model.spin2.value);
        if r2 < 0.0 {
            r2 = rl2;
        } else if r2 > rl2 {
            panic!("Secondary is larger than Roche Lobe.")
        }

        set_star_continuum(&model, &mut star1_fine_grid, &mut star2_fine_grid);
        

        if model.nlat1f == model.nlat1c {
            star1_coarse_grid = star1_fine_grid.clone();
        } else {
            star1_coarse_grid = set_star_grid(&model, Star::Primary, false);
        }
        let copy2: bool = (model.nlat2f == model.nlat2c) &&
                        (!model.npole || r1 >= r2 || (model.nlatfill == 0 && model.nlngfill == 0));
        if copy2 {
            star2_coarse_grid = star2_fine_grid.clone();
        } else {
            star2_coarse_grid = set_star_grid(&model, Star::Secondary, false)
        }

        if model.nlat1c != model.nlat1f || !copy2 {
            set_star_continuum(&model, &mut star1_coarse_grid, &mut star2_coarse_grid);
        }

        let disc_grid: Vec<Point> = vec![];
        let disc_edge_grid: Vec<Point> = vec![];
        let bright_spot_grid: Vec<Point> = vec![];

        let gint: Ginterp = Ginterp{ phase1: model.phase1, phase2: model.phase2, scale11: 1.0, scale12: 1.0, scale21: 1.0, scale22: 1.0};


        let mut rlens1 = 0.0;
        if model.glens1 {
            let gm: f64 = (1000.0*model.velocity_scale.value).powi(3)*model.tperiod*DAY/TAU;
            let a: f64 = (gm/((TAU/DAY/model.tperiod)*(TAU/DAY/model.tperiod))).powf(1.0/3.0);
            rlens1 = 4.0*gm/(1.0+model.q.value)/a/(C*C);
        }

        // if model.nlat1c != model.nlat1f


        Ok(Self {
            star1_coarse_grid,
            star2_coarse_grid,
            star1_fine_grid,
            star2_fine_grid,
            disc_grid,
            disc_edge_grid,
            bright_spot_grid,
            gint,
            rlens1,
            model
        })
    }

    #[staticmethod]
    pub fn from_file(filename: &str) -> PyResult<Self> {
        let model = Model::from_file(filename)
            .map_err(|e| pyo3::exceptions::PyIOError::new_err(e))?;
        let mut star1_fine_grid = set_star_grid(&model, Star::Primary, true);
        let mut star2_fine_grid = set_star_grid(&model, Star::Secondary, true);
        let mut star1_coarse_grid: Vec<Point>;
        let mut star2_coarse_grid: Vec<Point>;

        let (r1, mut r2) = model.get_r1r2();
        let rl2: f64 = 1.0 - xl12(model.q.value, model.spin2.value);
        if r2 < 0.0 {
            r2 = rl2;
        } else if r2 > rl2 {
            panic!("Secondary is larger than Roche Lobe.")
        }
        
        set_star_continuum(&model, &mut star1_fine_grid, &mut star2_fine_grid);
        
        
        if model.nlat1f == model.nlat1c {
            star1_coarse_grid = star1_fine_grid.clone();
        } else {
            star1_coarse_grid = set_star_grid(&model, Star::Primary, false);
        }
        
        let copy2: bool = (model.nlat2f == model.nlat2c) &&
        (!model.npole || r1 >= r2 || (model.nlatfill == 0 && model.nlngfill == 0));
        
        if copy2 {
            star2_coarse_grid = star2_fine_grid.clone();
        } else {
            star2_coarse_grid = set_star_grid(&model, Star::Secondary, false)
        }
        
        if model.nlat1c != model.nlat1f || !copy2 {
            set_star_continuum(&model, &mut star1_coarse_grid, &mut star2_coarse_grid);
        }
        
        let mut disc_grid: Vec<Point> = vec![];
        let mut disc_edge_grid: Vec<Point> = vec![];
        let mut bright_spot_grid: Vec<Point> = vec![];


        let mut gint: Ginterp = Ginterp{ phase1: model.phase1, phase2: model.phase2, scale11: 1.0, scale12: 1.0, scale21: 1.0, scale22: 1.0};


        let mut rlens1 = 0.0;
        if model.glens1 {
            let gm: f64 = (1000.0*model.velocity_scale.value).powi(3)*model.tperiod*DAY/TAU;
            let a: f64 = (gm/((TAU/DAY/model.tperiod)*(TAU/DAY/model.tperiod))).powf(1.0/3.0);
            rlens1 = 4.0*gm/(1.0+model.q.value)/a/(C*C);
        }

        let ldc1: LDC = model.get_ldc1();
        let ldc2: LDC = model.get_ldc2();

        if model.nlat1c != model.nlat1f {
            let ff: f64 = comp_star1(model.iangle.value, &ldc1, 0.9999999999*model.phase1, 0.0, 1, model.q.value, model.beam_factor1.value, model.velocity_scale.value, &gint, &star1_fine_grid, &star1_coarse_grid);
            let fc: f64 = comp_star1(model.iangle.value, &ldc1, 1.0000000001*model.phase1, 0.0, 1, model.q.value, model.beam_factor1.value, model.velocity_scale.value, &gint, &star1_fine_grid, &star1_coarse_grid);
            gint.scale11 = ff/fc;

            let ff: f64 = comp_star1(model.iangle.value, &ldc1, 1.0-0.9999999999*model.phase1, 0.0, 1, model.q.value, model.beam_factor1.value, model.velocity_scale.value, &gint, &star1_fine_grid, &star1_coarse_grid);
            let fc: f64 = comp_star1(model.iangle.value, &ldc1, 1.0-1.0000000001*model.phase1, 0.0, 1, model.q.value, model.beam_factor1.value, model.velocity_scale.value, &gint, &star1_fine_grid, &star1_coarse_grid);
            gint.scale12 = ff/fc;
        }

        if !copy2 {
            let ff: f64 = comp_star2(model.iangle.value, &ldc2, 1.0-1.0000000001*model.phase2, 0.0, 1, model.q.value, model.beam_factor2.value, model.velocity_scale.value, model.glens1, rlens1, &gint, &star2_fine_grid, &star2_coarse_grid);
            let fc: f64 = comp_star2(model.iangle.value, &ldc2, 1.0-0.9999999999*model.phase2, 0.0, 1, model.q.value, model.beam_factor2.value, model.velocity_scale.value, model.glens1, rlens1, &gint, &star2_fine_grid, &star2_coarse_grid);
            gint.scale21 = ff/fc;

            let ff: f64 = comp_star2(model.iangle.value, &ldc2, 1.0000000001*model.phase2, 0.0, 1, model.q.value, model.beam_factor2.value, model.velocity_scale.value, model.glens1, rlens1, &gint, &star2_fine_grid, &star2_coarse_grid);
            let fc: f64 = comp_star2(model.iangle.value, &ldc2, 0.9999999999*model.phase2, 0.0, 1, model.q.value, model.beam_factor2.value, model.velocity_scale.value, model.glens1, rlens1, &gint, &star2_fine_grid, &star2_coarse_grid);
            gint.scale22 = ff/fc;
        }

        if model.add_disc {
            disc_grid = set_disc_grid(&model);
            disc_edge_grid = set_disc_edge_grid(&model, true, false);

            let rdisc1 = if model.rdisc1.value > 0.0 {
                model.rdisc1.value
            } else {
                r1
            };
            let rdisc2 = if model.rdisc2.value > 0.0 {
                model.rdisc2.value
            } else {
                model.radius_spot.value
            };

            let mut eclipses: Etype;
            if model.opaque {
                for point in &mut star1_fine_grid {
                    eclipses = disc_eclipse(model.iangle.value, rdisc1, rdisc2, model.beta_disc.value, model.height_disc.value, &point.position);
                    for i in 0..eclipses.len() {
                        point.eclipse.push(eclipses[i]);
                    }
                }
                for point in &mut star1_coarse_grid {
                    eclipses = disc_eclipse(model.iangle.value, rdisc1, rdisc2, model.beta_disc.value, model.height_disc.value, &point.position);
                    for i in 0..eclipses.len() {
                        point.eclipse.push(eclipses[i]);
                    }
                }
                for point in &mut star2_fine_grid {
                    eclipses = disc_eclipse(model.iangle.value, rdisc1, rdisc2, model.beta_disc.value, model.height_disc.value, &point.position);
                    for i in 0..eclipses.len() {
                        point.eclipse.push(eclipses[i]);
                    }
                }
                for point in &mut star2_coarse_grid {
                    eclipses = disc_eclipse(model.iangle.value, rdisc1, rdisc2, model.beta_disc.value, model.height_disc.value, &point.position);
                    for i in 0..eclipses.len() {
                        point.eclipse.push(eclipses[i]);
                    }
                }
            }
            
            // Set the surface brightness of the disc
            set_disc_continuum(rdisc2, model.temp_disc.value, model.texp_disc.value, model.wavelength, &mut disc_grid);

            // Set the surface brightness of outer edge, accounting for
            // irradiation by star 2
            set_edge_continuum(model.temp_edge.value, r2, model.t2.value.abs(), model.absorb_edge.value, model.wavelength, &mut disc_edge_grid);

        }

        if model.add_spot {

            bright_spot_grid = set_bright_spot_grid(&model);
        }

        
        Ok(Self {
            star1_coarse_grid,
            star2_coarse_grid,
            star1_fine_grid,
            star2_fine_grid,
            disc_grid,
            disc_edge_grid,
            bright_spot_grid,
            gint,
            rlens1,
            model
        })
    }


    pub fn update_wavelength_dependent(&mut self, new_model: Bound<'_, PyDict>) -> () {
        let map = map_from_pydict(new_model).unwrap();
        let new_model = Model::from_map(map).map_err(|e| pyo3::exceptions::PyValueError::new_err(e)).unwrap();
        self.model.t1.value = new_model.t1.value;
        self.model.t2.value = new_model.t2.value;
        self.model.ldc1_1.value = new_model.ldc1_1.value;
        self.model.ldc1_2.value = new_model.ldc1_2.value;
        self.model.ldc1_3.value = new_model.ldc1_3.value;
        self.model.ldc1_4.value = new_model.ldc1_4.value;
        self.model.ldc2_1.value = new_model.ldc2_1.value;
        self.model.ldc2_2.value = new_model.ldc2_2.value;
        self.model.ldc2_3.value = new_model.ldc2_3.value;
        self.model.ldc2_4.value = new_model.ldc2_4.value;
        self.model.gravity_dark1.value = new_model.gravity_dark1.value;
        self.model.gravity_dark2.value = new_model.gravity_dark2.value;
        self.model.beam_factor1.value = new_model.beam_factor1.value;
        self.model.beam_factor2.value = new_model.beam_factor2.value;
        self.model.absorb.value = new_model.absorb.value;
        self.model.slope.value = new_model.slope.value;
        self.model.quad.value = new_model.quad.value;
        self.model.cube.value = new_model.cube.value;
        self.model.third.value = new_model.third.value;
        self.model.wavelength = new_model.wavelength;
        set_star_continuum(&self.model, &mut self.star1_coarse_grid, &mut self.star2_coarse_grid);
        set_star_continuum(&self.model, &mut self.star1_fine_grid, &mut self.star2_fine_grid);
        
    }


    pub fn planck(
        &self,
        wavelength: Bound<'_, PyFloat>,
        temp: Bound<'_, PyFloat>
    ) -> PyResult<f64> {
        let wavelength = wavelength.extract::<f64>().unwrap();
        let temp = temp.extract::<f64>().unwrap();
        let flux = planck(wavelength, temp);
        Ok(flux)
    }


    pub fn compute_light_curve(
        &self,
        time: PyReadonlyArray1<f64>,
        t_exp: PyReadonlyArray1<f64>,
        flux: PyReadonlyArray1<f64>,
        flux_err: PyReadonlyArray1<f64>,
        weight: PyReadonlyArray1<f64>,
        n_div: PyReadonlyArray1<f64>,
        mut star1: PyReadwriteArray1<f64>,
        mut star2: PyReadwriteArray1<f64>,
        mut disc: PyReadwriteArray1<f64>,
        mut disc_edge: PyReadwriteArray1<f64>,
        mut bright_spot: PyReadwriteArray1<f64>,
    ) -> PyResult<()> {


        let time: &[f64] = time.as_slice()?;
        let t_exp: &[f64] = t_exp.as_slice()?;
        let flux: &[f64] = flux.as_slice()?;
        let flux_err: &[f64] = flux_err.as_slice()?;
        let weight: &[f64] = weight.as_slice()?;
        let n_div: &[f64] = n_div.as_slice()?;

        let star1: &mut [f64] = star1.as_slice_mut()?;
        let star2: &mut[f64] = star2.as_slice_mut()?;
        let disc: &mut[f64] = disc.as_slice_mut()?;
        let disc_edge: &mut[f64] = disc_edge.as_slice_mut()?;
        let bright_spot: &mut[f64] = bright_spot.as_slice_mut()?;

        let ldc1: LDC = self.model.get_ldc1();
        let ldc2: LDC = self.model.get_ldc2();

        star1
        .par_iter_mut()
        .enumerate()
        .for_each(|(i, out)| {
            let phase = (time[i] - self.model.t0.value) / self.model.period.value;
            let expose: f64 = t_exp[i]/self.model.period.value;
            *out = self.compute_star1_flux(phase, &ldc1, expose, n_div[i] as i32);
        });

        // let scale = rescale(flux, flux_err, weight, star1);
        // for val in star1 {
        //     *val *= scale;

        // }

        star2
        .par_iter_mut()
        .enumerate()
        .for_each(|(i, out)| {
            let phase = (time[i] - self.model.t0.value) / self.model.period.value;
            let expose: f64 = t_exp[i]/self.model.period.value;
            *out = self.compute_star2_flux(phase, &ldc2, expose, n_div[i] as i32);
        });


        disc
        .par_iter_mut()
        .enumerate()
        .for_each(|(i, out)| {
            let phase = (time[i] - self.model.t0.value) / self.model.period.value;
            let expose: f64 = t_exp[i]/self.model.period.value;
            *out = self.compute_disc_flux(phase, expose, n_div[i] as i32);
        });

        disc_edge
        .par_iter_mut()
        .enumerate()
        .for_each(|(i, out)| {
            let phase = (time[i] - self.model.t0.value) / self.model.period.value;
            let expose: f64 = t_exp[i]/self.model.period.value;
            *out = self.compute_disc_edge_flux(phase, expose, n_div[i] as i32);
        });

        bright_spot
        .par_iter_mut()
        .enumerate()
        .for_each(|(i, out)| {
            let phase = (time[i] - self.model.t0.value) / self.model.period.value;
            let expose: f64 = t_exp[i]/self.model.period.value;
            *out = self.compute_bright_spot_flux(phase, expose, n_div[i] as i32);
        });

        Ok(())

    }
}


impl BinaryModel {
    fn compute_star1_flux(&self, phase: f64, ldc1: &LDC, expose: f64, n_div: i32) -> f64 {

        let flux: f64;

        // flux = comp_light(self.model.iangle.value, ldc1, ldc2, phase, expose, n_div, self.model.q.value, self.model.beam_factor1.value, self.model.beam_factor2.value, self.model.spin1.value, self.model.spin2.value, self.model.velocity_scale.value, self.model.glens1, self.rlens1, &self.gint, &self.star1_fine_grid, &self.star2_fine_grid, &self.star1_coarse_grid, &self.star2_coarse_grid);
        flux = comp_star1(self.model.iangle.value, ldc1, phase, expose, n_div, self.model.q.value, self.model.beam_factor1.value, self.model.velocity_scale.value, &self.gint, &self.star1_fine_grid, &self.star1_coarse_grid);

        flux
    }

    fn compute_star2_flux(&self, phase: f64, ldc2: &LDC, expose: f64, n_div: i32) -> f64 {

        let star2_flux: f64;

        // (_, star2_flux) = comp_light(self.model.iangle.value, ldc1, ldc2, phase, expose, n_div, self.model.q.value, self.model.beam_factor1.value, self.model.beam_factor2.value, self.model.spin1.value, self.model.spin2.value, self.model.velocity_scale.value, self.model.glens1, self.rlens1, &self.gint, &self.star1_fine_grid, &self.star2_fine_grid, &self.star1_coarse_grid, &self.star2_coarse_grid);
        star2_flux = comp_star2(self.model.iangle.value, ldc2, phase, expose, n_div, self.model.q.value, self.model.beam_factor2.value, self.model.velocity_scale.value, self.model.glens1, self.rlens1, &self.gint, &self.star2_fine_grid, &self.star2_coarse_grid);

        star2_flux
    }

    fn compute_disc_flux(&self, phase: f64, expose: f64, n_div: i32) -> f64 {

        let disc_flux: f64;

        // (_, star2_flux) = comp_light(self.model.iangle.value, ldc1, ldc2, phase, expose, n_div, self.model.q.value, self.model.beam_factor1.value, self.model.beam_factor2.value, self.model.spin1.value, self.model.spin2.value, self.model.velocity_scale.value, self.model.glens1, self.rlens1, &self.gint, &self.star1_fine_grid, &self.star2_fine_grid, &self.star1_coarse_grid, &self.star2_coarse_grid);
        disc_flux = comp_disc(self.model.iangle.value, self.model.lin_limb_disc.value, self.model.quad_limb_disc.value, phase, expose, n_div, &self.disc_grid);

        disc_flux
    }

    fn compute_disc_edge_flux(&self, phase: f64, expose: f64, n_div: i32) -> f64 {

        let disc_edge_flux: f64;

        // (_, star2_flux) = comp_light(self.model.iangle.value, ldc1, ldc2, phase, expose, n_div, self.model.q.value, self.model.beam_factor1.value, self.model.beam_factor2.value, self.model.spin1.value, self.model.spin2.value, self.model.velocity_scale.value, self.model.glens1, self.rlens1, &self.gint, &self.star1_fine_grid, &self.star2_fine_grid, &self.star1_coarse_grid, &self.star2_coarse_grid);
        disc_edge_flux = comp_disc_edge(self.model.iangle.value, self.model.lin_limb_disc.value, self.model.quad_limb_disc.value, phase, expose, n_div, &self.disc_edge_grid);

        disc_edge_flux
    }

    fn compute_bright_spot_flux(&self, phase: f64, expose: f64, n_div: i32) -> f64 {

        let bright_spot_flux: f64;

        // (_, star2_flux) = comp_light(self.model.iangle.value, ldc1, ldc2, phase, expose, n_div, self.model.q.value, self.model.beam_factor1.value, self.model.beam_factor2.value, self.model.spin1.value, self.model.spin2.value, self.model.velocity_scale.value, self.model.glens1, self.rlens1, &self.gint, &self.star1_fine_grid, &self.star2_fine_grid, &self.star1_coarse_grid, &self.star2_coarse_grid);
        bright_spot_flux = comp_bright_spot(self.model.iangle.value, phase, expose, n_div, &self.bright_spot_grid);

        bright_spot_flux
    }

}


pub fn map_from_pydict(dict: Bound<'_, PyDict>) -> PyResult<HashMap<String, Entry>> {

    let mut map: HashMap<String, Entry> = HashMap::new();

    for (key, value) in dict.iter() {

        let name: String = key.extract::<String>()?;

        if let Ok(v) = value.extract::<f64>() {
            map.insert(name, Entry::Scalar(v.to_string()));
        }
        else if let Ok(v) = value.extract::<bool>() {
            map.insert(name, Entry::Scalar(if v { "1".to_string() } else { "0".to_string() }));
        }
        else if let Ok(v) = value.extract::<String>() {
        map.insert(name, Entry::Scalar(v));
        }
        else {
            return Err(PyErr::new::<pyo3::exceptions::PyTypeError, _>(
            format!("Unsupported type for key {}", name)
            ));
        }
    }

    Ok(map)
}


// pub fn rescale(flux: &[f64], flux_err: &[f64], weight: &[f64], model_flux: &mut [f64]) -> f64 {
//     let mut sdy: f64 = 0.0;
//     let mut syy: f64 = 0.0;
//     for i in 0..flux.len() {
//         if weight[i] > 0.0 {
//             let wgt: f64 = weight[i] / (flux_err[i]*flux_err[i]);
//             sdy += wgt*flux[i]*model_flux[i];
//             syy += wgt*model_flux[i]*model_flux[i];
//         }
//     }
//     let scale: f64 = sdy/syy;
//     scale
// }