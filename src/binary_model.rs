use crate::comp_gravity::{comp_gravity1, comp_gravity2};
use crate::comp_light::{comp_bright_spot, comp_disc, comp_disc_edge, comp_star1, comp_star2};
use crate::comp_radius::comp_radius;
use crate::constants::{C, DAY};
use crate::ginterp::Ginterp;
use crate::ldc::LDC;
use crate::model::{Entry, Model, ModelUpdate};
use crate::set_bright_spot_grid::set_bright_spot_grid;
use crate::set_disc_continuum::{set_disc_continuum, set_edge_continuum};
use crate::set_disc_grid::{set_disc_edge_grid, set_disc_grid};
use crate::set_star_continuum::set_star_continuum;
use crate::set_star_grid::set_star_grid;
use numpy::{IntoPyArray, PyArray1, PyReadonlyArray1};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyDictMethods};
use rayon::prelude::*;
use roche::errors::RocheError;
use roche::{self, Etype, Point, Star, disc_eclipse};
use serde_pyobject::from_pyobject;
use std::collections::HashMap;
use std::f64::consts::TAU;

#[pyclass]
pub struct LightCurve {
    #[pyo3(get)]
    pub star1: Py<PyArray1<f64>>,

    #[pyo3(get)]
    pub star2: Py<PyArray1<f64>>,

    #[pyo3(get)]
    pub disc: Py<PyArray1<f64>>,

    #[pyo3(get)]
    pub disc_edge: Py<PyArray1<f64>>,

    #[pyo3(get)]
    pub bright_spot: Py<PyArray1<f64>>,

    #[pyo3(get)]
    pub total: Py<PyArray1<f64>>,

    #[pyo3(get)]
    pub star1_contribution: f64,

    #[pyo3(get)]
    pub logg1: f64,

    #[pyo3(get)]
    pub logg2: f64,

    #[pyo3(get)]
    pub rva1: f64,

    #[pyo3(get)]
    pub rva2: f64,

    #[pyo3(get)]
    pub chi2: Option<f64>,

    #[pyo3(get)]
    pub log_prob: Option<f64>,
}

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

    #[pyo3(get)]
    pub model: Model,
}

#[pymethods]
impl BinaryModel {

    #[staticmethod]
    pub fn from_file(filename: &str) -> PyResult<Self> {
        let model = Model::from_file(filename).map_err(pyo3::exceptions::PyIOError::new_err)?;
        let (
            star1_coarse_grid,
            star2_coarse_grid,
            star1_fine_grid,
            star2_fine_grid,
            disc_grid,
            disc_edge_grid,
            bright_spot_grid,
            gint,
            rlens1,
        ) = build_grids(&model)?;
        
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
            model,
        })
    }

    #[staticmethod]
    pub fn from_model(model: Model) -> PyResult<Self> {
        let (
            star1_coarse_grid,
            star2_coarse_grid,
            star1_fine_grid,
            star2_fine_grid,
            disc_grid,
            disc_edge_grid,
            bright_spot_grid,
            gint,
            rlens1,
        ) = build_grids(&model)?;

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
            model,
        })
    }

    pub fn update(&mut self, _py: Python, dict: &Bound<'_, PyAny>) -> PyResult<()> {
        let upd: ModelUpdate = from_pyobject(dict.clone())?;
        let grid_changed = upd.grid_changed();
        self.model.apply_update(upd);
        if grid_changed {
            (
                self.star1_coarse_grid,
                self.star2_coarse_grid,
                self.star1_fine_grid,
                self.star2_fine_grid,
                self.disc_grid,
                self.disc_edge_grid,
                self.bright_spot_grid,
                self.gint,
                self.rlens1,
            ) = build_grids(&self.model)?;
        } else {
            self.reset_grid_continuum()?;
        }
        Ok(())
    }

    #[pyo3(signature = (
        time,
        t_exp,
        n_div=None,
        flux=None,
        flux_err=None,
        weight=None,
        autoscale=true
    ))]
    pub fn compute_light_curve(
        &self,
        py: Python,
        time: PyReadonlyArray1<f64>,
        t_exp: PyReadonlyArray1<f64>,
        n_div: Option<PyReadonlyArray1<f64>>,
        flux: Option<PyReadonlyArray1<f64>>,
        flux_err: Option<PyReadonlyArray1<f64>>,
        weight: Option<PyReadonlyArray1<f64>>,
        autoscale: bool,
    ) -> PyResult<LightCurve> {


        let time: &[f64] = time.as_slice()?;
        let t_exp: &[f64] = t_exp.as_slice()?;
        // let n_div: &[f64] = n_div.as_slice()?;
        let n: usize = time.len();

        let n_div_default;
        let n_div = if let Some(ref ndiv) = n_div {
            ndiv.as_slice()?
        } else {
            n_div_default = vec![1.0_f64; n];
            &n_div_default
        };
        let flux = match &flux {
            Some(f) => Some(f.as_slice()?),
            None => None,
        };
        let flux_err = match &flux_err {
            Some(f) => Some(f.as_slice()?),
            None => None,
        };
        let weight = match &weight {
            Some(f) => Some(f.as_slice()?),
            None => None,
        };

        let mut star1 = vec![0.0; n];
        let mut star2 = vec![0.0; n];
        let mut disc = vec![0.0; n];
        let mut disc_edge = vec![0.0; n];
        let mut bright_spot = vec![0.0; n];
        let mut total = vec![0.0; n];

        let ldc1: LDC = self.model.get_ldc1();
        let ldc2: LDC = self.model.get_ldc2();

        star1.par_iter_mut().enumerate().for_each(|(i, out)| {
            let mut phase = (time[i] - self.model.t0.value) / self.model.period.value;
            // small Newton-Raphson iteration
            for _ in 0..4 {
                phase -= (self.model.t0.value
                    + phase * (self.model.period.value + self.model.pdot.value * phase)
                    - time[i])
                    / (self.model.period.value + 2.0 * self.model.pdot.value * phase);
            }
            // advance/retard by time offset between primary & secondary eclipse
            phase += self.model.deltat.value / self.model.period.value / 2.0
                * ((TAU * phase).cos() - 1.0);
            let expose: f64 = t_exp[i] / self.model.period.value;
            *out = self.compute_star1_flux(phase, &ldc1, expose, n_div[i] as i32);
        });

        star2.par_iter_mut().enumerate().for_each(|(i, out)| {
            let mut phase = (time[i] - self.model.t0.value) / self.model.period.value;
            // small Newton-Raphson iteration
            for _ in 0..4 {
                phase -= (self.model.t0.value
                    + phase * (self.model.period.value + self.model.pdot.value * phase)
                    - time[i])
                    / (self.model.period.value + 2.0 * self.model.pdot.value * phase);
            }
            // advance/retard by time offset between primary & secondary eclipse
            phase += self.model.deltat.value / self.model.period.value / 2.0
                * ((TAU * phase).cos() - 1.0);
            let expose: f64 = t_exp[i] / self.model.period.value;
            *out = self.compute_star2_flux(phase, &ldc2, expose, n_div[i] as i32);
        });
        
        disc.par_iter_mut().enumerate().for_each(|(i, out)| {
            let mut phase = (time[i] - self.model.t0.value) / self.model.period.value;
            // small Newton-Raphson iteration
            for _ in 0..4 {
                phase -= (self.model.t0.value
                    + phase * (self.model.period.value + self.model.pdot.value * phase)
                    - time[i])
                    / (self.model.period.value + 2.0 * self.model.pdot.value * phase);
            }
            // advance/retard by time offset between primary & secondary eclipse
            phase += self.model.deltat.value / self.model.period.value / 2.0
                * ((TAU * phase).cos() - 1.0);
            let expose: f64 = t_exp[i] / self.model.period.value;
            *out = self.compute_disc_flux(phase, expose, n_div[i] as i32);
        });
        
        disc_edge.par_iter_mut().enumerate().for_each(|(i, out)| {
            let mut phase = (time[i] - self.model.t0.value) / self.model.period.value;
            // small Newton-Raphson iteration
            for _ in 0..4 {
                phase -= (self.model.t0.value
                    + phase * (self.model.period.value + self.model.pdot.value * phase)
                    - time[i])
                    / (self.model.period.value + 2.0 * self.model.pdot.value * phase);
            }
            // advance/retard by time offset between primary & secondary eclipse
            phase += self.model.deltat.value / self.model.period.value / 2.0
                * ((TAU * phase).cos() - 1.0);
            let expose: f64 = t_exp[i] / self.model.period.value;
            *out = self.compute_disc_edge_flux(phase, expose, n_div[i] as i32);
        });

        bright_spot.par_iter_mut().enumerate().for_each(|(i, out)| {
            let mut phase = (time[i] - self.model.t0.value) / self.model.period.value;
            // small Newton-Raphson iteration
            for _ in 0..4 {
                phase -= (self.model.t0.value
                    + phase * (self.model.period.value + self.model.pdot.value * phase)
                    - time[i])
                    / (self.model.period.value + 2.0 * self.model.pdot.value * phase);
            }
            // advance/retard by time offset between primary & secondary eclipse
            phase += self.model.deltat.value / self.model.period.value / 2.0
                * ((TAU * phase).cos() - 1.0);
            let expose: f64 = t_exp[i] / self.model.period.value;
            *out = self.compute_bright_spot_flux(phase, expose, n_div[i] as i32);
        });

        let mut star1_contribution: f64 = self.compute_star1_flux(0.5, &ldc1, 0.0, 1);
        
        for i in 0..time.len() {
            total[i] = star1[i] + star2[i] + disc[i] + disc_edge[i] + bright_spot[i];
        }

        let (chisq, log_prob) = if flux.is_some() && flux_err.is_some() && autoscale {
            let scale_factor = rescale(flux.unwrap(), flux_err.unwrap(), weight, &total);
            for i in 0..time.len() {
                star1[i] *= scale_factor;
                star2[i] *= scale_factor;
                disc[i] *= scale_factor;
                disc_edge[i] *= scale_factor;
                bright_spot[i] *= scale_factor;
                total[i] *= scale_factor;
            }

            star1_contribution *= scale_factor;

            let (chisq, log_prob) =
                chisq_log_prob(flux.unwrap(), flux_err.unwrap(), weight, &total);
            (Some(chisq), Some(log_prob))
        } else {
            (None, None)
        };

        let logg1: f64 = comp_gravity1(&self.model, &self.star1_fine_grid)?;
        let logg2: f64 = comp_gravity2(&self.model, &self.star2_fine_grid)?;
        let rva1: f64 = if self.model.roche1 {
            comp_radius(&self.star1_coarse_grid, Star::Primary)
            } else {
                self.model.r1.value
            };
        let rva2: f64 = comp_radius(&self.star2_coarse_grid, Star::Secondary);
        
        Ok(LightCurve {
            star1: star1.into_pyarray(py).unbind(),
            star2: star2.into_pyarray(py).unbind(),
            disc: disc.into_pyarray(py).unbind(),
            disc_edge: disc_edge.into_pyarray(py).unbind(),
            bright_spot: bright_spot.into_pyarray(py).unbind(),
            total: total.into_pyarray(py).unbind(),
            star1_contribution,
            logg1,
            logg2,
            rva1,
            rva2,
            chi2: chisq,
            log_prob,
        })
    }
}

impl BinaryModel {
    fn compute_star1_flux(&self, phase: f64, ldc1: &LDC, expose: f64, n_div: i32) -> f64 {
        comp_star1(
            self.model.iangle.value,
            ldc1,
            phase,
            expose,
            n_div,
            self.model.q.value,
            self.model.beam_factor1.value,
            self.model.velocity_scale.value,
            &self.gint,
            &self.star1_fine_grid,
            &self.star1_coarse_grid,
        )
    }

    fn compute_star2_flux(&self, phase: f64, ldc2: &LDC, expose: f64, n_div: i32) -> f64 {
        comp_star2(
            self.model.iangle.value,
            ldc2,
            phase,
            expose,
            n_div,
            self.model.q.value,
            self.model.beam_factor2.value,
            self.model.velocity_scale.value,
            self.model.glens1,
            self.rlens1,
            &self.gint,
            &self.star2_fine_grid,
            &self.star2_coarse_grid,
        )
    }

    fn compute_disc_flux(&self, phase: f64, expose: f64, n_div: i32) -> f64 {
        comp_disc(
            self.model.iangle.value,
            self.model.lin_limb_disc.value,
            self.model.quad_limb_disc.value,
            phase,
            expose,
            n_div,
            &self.disc_grid,
        )
    }

    fn compute_disc_edge_flux(&self, phase: f64, expose: f64, n_div: i32) -> f64 {
        comp_disc_edge(
            self.model.iangle.value,
            self.model.lin_limb_disc.value,
            self.model.quad_limb_disc.value,
            phase,
            expose,
            n_div,
            &self.disc_edge_grid,
        )
    }

    fn compute_bright_spot_flux(&self, phase: f64, expose: f64, n_div: i32) -> f64 {
        comp_bright_spot(
            self.model.iangle.value,
            phase,
            expose,
            n_div,
            &self.bright_spot_grid,
        )
    }

    fn reset_grid_continuum(&mut self) -> Result<(), RocheError> {
        let (r1, mut r2) = self.model.get_r1r2();
        let rl2: f64 = 1.0 - roche::x_l1_2(self.model.q.value, self.model.spin2.value)?;
        if r2 < 0.0 {
            r2 = rl2;
        } else if r2 > rl2 {
            panic!("Secondary is larger than Roche Lobe.")
        }

        let ldc1: LDC = self.model.get_ldc1();
        let ldc2: LDC = self.model.get_ldc2();

        set_star_continuum(
            &self.model,
            &mut self.star1_fine_grid,
            &mut self.star2_fine_grid,
        )?;
        set_star_continuum(
            &self.model,
            &mut self.star1_coarse_grid,
            &mut self.star2_coarse_grid,
        )?;

        let copy2: bool = (self.model.nlat2f == self.model.nlat2c)
            && (!self.model.npole
                || r1 >= r2
                || (self.model.nlatfill == 0 && self.model.nlngfill == 0));
        if self.model.nlat1c != self.model.nlat1f {
            let ff: f64 = comp_star1(
                self.model.iangle.value,
                &ldc1,
                0.9999999999 * self.model.phase1,
                0.0,
                1,
                self.model.q.value,
                self.model.beam_factor1.value,
                self.model.velocity_scale.value,
                &self.gint,
                &self.star1_fine_grid,
                &self.star1_coarse_grid,
            );
            let fc: f64 = comp_star1(
                self.model.iangle.value,
                &ldc1,
                1.0000000001 * self.model.phase1,
                0.0,
                1,
                self.model.q.value,
                self.model.beam_factor1.value,
                self.model.velocity_scale.value,
                &self.gint,
                &self.star1_fine_grid,
                &self.star1_coarse_grid,
            );
            self.gint.scale11 = ff / fc;

            let ff: f64 = comp_star1(
                self.model.iangle.value,
                &ldc1,
                1.0 - 0.9999999999 * self.model.phase1,
                0.0,
                1,
                self.model.q.value,
                self.model.beam_factor1.value,
                self.model.velocity_scale.value,
                &self.gint,
                &self.star1_fine_grid,
                &self.star1_coarse_grid,
            );
            let fc: f64 = comp_star1(
                self.model.iangle.value,
                &ldc1,
                1.0 - 1.0000000001 * self.model.phase1,
                0.0,
                1,
                self.model.q.value,
                self.model.beam_factor1.value,
                self.model.velocity_scale.value,
                &self.gint,
                &self.star1_fine_grid,
                &self.star1_coarse_grid,
            );
            self.gint.scale12 = ff / fc;
        }

        if !copy2 {
            let ff: f64 = comp_star2(
                self.model.iangle.value,
                &ldc2,
                1.0 - 1.0000000001 * self.model.phase2,
                0.0,
                1,
                self.model.q.value,
                self.model.beam_factor2.value,
                self.model.velocity_scale.value,
                self.model.glens1,
                self.rlens1,
                &self.gint,
                &self.star2_fine_grid,
                &self.star2_coarse_grid,
            );
            let fc: f64 = comp_star2(
                self.model.iangle.value,
                &ldc2,
                1.0 - 0.9999999999 * self.model.phase2,
                0.0,
                1,
                self.model.q.value,
                self.model.beam_factor2.value,
                self.model.velocity_scale.value,
                self.model.glens1,
                self.rlens1,
                &self.gint,
                &self.star2_fine_grid,
                &self.star2_coarse_grid,
            );
            self.gint.scale21 = ff / fc;

            let ff: f64 = comp_star2(
                self.model.iangle.value,
                &ldc2,
                1.0000000001 * self.model.phase2,
                0.0,
                1,
                self.model.q.value,
                self.model.beam_factor2.value,
                self.model.velocity_scale.value,
                self.model.glens1,
                self.rlens1,
                &self.gint,
                &self.star2_fine_grid,
                &self.star2_coarse_grid,
            );
            let fc: f64 = comp_star2(
                self.model.iangle.value,
                &ldc2,
                0.9999999999 * self.model.phase2,
                0.0,
                1,
                self.model.q.value,
                self.model.beam_factor2.value,
                self.model.velocity_scale.value,
                self.model.glens1,
                self.rlens1,
                &self.gint,
                &self.star2_fine_grid,
                &self.star2_coarse_grid,
            );
            self.gint.scale22 = ff / fc;
        }

        if self.model.add_disc {
            let rdisc2 = if self.model.rdisc2.value > 0.0 {
                self.model.rdisc2.value
            } else {
                self.model.radius_spot.value
            };
            
            // Set the surface brightness of the disc
            set_disc_continuum(
                rdisc2,
                self.model.temp_disc.value,
                self.model.texp_disc.value,
                self.model.wavelength,
                &mut self.disc_grid,
            );

            // Set the surface brightness of outer edge, accounting for
            // irradiation by star 2
            set_edge_continuum(
                self.model.temp_edge.value,
                r2,
                self.model.t2.value.abs(),
                self.model.absorb_edge.value,
                self.model.wavelength,
                &mut self.disc_edge_grid,
            );
        }

        if self.model.add_spot {
            self.bright_spot_grid = set_bright_spot_grid(&self.model)?;
        }
        Ok(())
    }
}

fn build_grids(
    model: &Model,
) -> Result<
    (
        Vec<Point>,
        Vec<Point>,
        Vec<Point>,
        Vec<Point>,
        Vec<Point>,
        Vec<Point>,
        Vec<Point>,
        Ginterp,
        f64,
    ),
    RocheError,
> {
    let mut star1_fine_grid = set_star_grid(model, Star::Primary, true)?;
    let mut star2_fine_grid = set_star_grid(model, Star::Secondary, true)?;
    let mut star1_coarse_grid: Vec<Point>;
    let mut star2_coarse_grid: Vec<Point>;

    let (r1, mut r2) = model.get_r1r2();
    let rl2: f64 = 1.0 - roche::x_l1_2(model.q.value, model.spin2.value)?;
    if r2 < 0.0 {
        r2 = rl2;
    } else if r2 > rl2 {
        panic!("Secondary is larger than Roche Lobe.")
    }
    
    set_star_continuum(model, &mut star1_fine_grid, &mut star2_fine_grid)?;
    
    if model.nlat1f == model.nlat1c {
        star1_coarse_grid = star1_fine_grid.clone();
    } else {
        star1_coarse_grid = set_star_grid(model, Star::Primary, false)?;
    }
    
    let copy2: bool = (model.nlat2f == model.nlat2c)
        && (!model.npole || r1 >= r2 || (model.nlatfill == 0 && model.nlngfill == 0));
    
    if copy2 {
        star2_coarse_grid = star2_fine_grid.clone();
    } else {
        star2_coarse_grid = set_star_grid(model, Star::Secondary, false)?
    }
    
    if model.nlat1c != model.nlat1f || !copy2 {
        set_star_continuum(model, &mut star1_coarse_grid, &mut star2_coarse_grid)?;
    }
    
    let mut disc_grid: Vec<Point> = vec![];
    let mut disc_edge_grid: Vec<Point> = vec![];
    let mut bright_spot_grid: Vec<Point> = vec![];

    let mut gint: Ginterp = Ginterp {
        phase1: model.phase1,
        phase2: model.phase2,
        scale11: 1.0,
        scale12: 1.0,
        scale21: 1.0,
        scale22: 1.0,
    };

    let mut rlens1 = 0.0;
    if model.glens1 {
        let gm: f64 = (1000.0 * model.velocity_scale.value).powi(3) * model.tperiod * DAY / TAU;
        let a: f64 =
            (gm / ((TAU / DAY / model.tperiod) * (TAU / DAY / model.tperiod))).powf(1.0 / 3.0);
        rlens1 = 4.0 * gm / (1.0 + model.q.value) / a / (C * C);
    }

    let ldc1: LDC = model.get_ldc1();
    let ldc2: LDC = model.get_ldc2();

    if model.nlat1c != model.nlat1f {
        let ff: f64 = comp_star1(
            model.iangle.value,
            &ldc1,
            0.9999999999 * model.phase1,
            0.0,
            1,
            model.q.value,
            model.beam_factor1.value,
            model.velocity_scale.value,
            &gint,
            &star1_fine_grid,
            &star1_coarse_grid,
        );
        let fc: f64 = comp_star1(
            model.iangle.value,
            &ldc1,
            1.0000000001 * model.phase1,
            0.0,
            1,
            model.q.value,
            model.beam_factor1.value,
            model.velocity_scale.value,
            &gint,
            &star1_fine_grid,
            &star1_coarse_grid,
        );
        gint.scale11 = ff / fc;

        let ff: f64 = comp_star1(
            model.iangle.value,
            &ldc1,
            1.0 - 0.9999999999 * model.phase1,
            0.0,
            1,
            model.q.value,
            model.beam_factor1.value,
            model.velocity_scale.value,
            &gint,
            &star1_fine_grid,
            &star1_coarse_grid,
        );
        let fc: f64 = comp_star1(
            model.iangle.value,
            &ldc1,
            1.0 - 1.0000000001 * model.phase1,
            0.0,
            1,
            model.q.value,
            model.beam_factor1.value,
            model.velocity_scale.value,
            &gint,
            &star1_fine_grid,
            &star1_coarse_grid,
        );
        gint.scale12 = ff / fc;
    }

    if !copy2 {
        let ff: f64 = comp_star2(
            model.iangle.value,
            &ldc2,
            1.0 - 1.0000000001 * model.phase2,
            0.0,
            1,
            model.q.value,
            model.beam_factor2.value,
            model.velocity_scale.value,
            model.glens1,
            rlens1,
            &gint,
            &star2_fine_grid,
            &star2_coarse_grid,
        );
        let fc: f64 = comp_star2(
            model.iangle.value,
            &ldc2,
            1.0 - 0.9999999999 * model.phase2,
            0.0,
            1,
            model.q.value,
            model.beam_factor2.value,
            model.velocity_scale.value,
            model.glens1,
            rlens1,
            &gint,
            &star2_fine_grid,
            &star2_coarse_grid,
        );
        gint.scale21 = ff / fc;

        let ff: f64 = comp_star2(
            model.iangle.value,
            &ldc2,
            1.0000000001 * model.phase2,
            0.0,
            1,
            model.q.value,
            model.beam_factor2.value,
            model.velocity_scale.value,
            model.glens1,
            rlens1,
            &gint,
            &star2_fine_grid,
            &star2_coarse_grid,
        );
        let fc: f64 = comp_star2(
            model.iangle.value,
            &ldc2,
            0.9999999999 * model.phase2,
            0.0,
            1,
            model.q.value,
            model.beam_factor2.value,
            model.velocity_scale.value,
            model.glens1,
            rlens1,
            &gint,
            &star2_fine_grid,
            &star2_coarse_grid,
        );
        gint.scale22 = ff / fc;
    }

    if model.add_disc {
        disc_grid = set_disc_grid(model)?;
        disc_edge_grid = set_disc_edge_grid(model, true, false)?;

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
                eclipses = disc_eclipse(
                    model.iangle.value,
                    rdisc1,
                    rdisc2,
                    model.beta_disc.value,
                    model.height_disc.value,
                    &point.position,
                )?;
                for eclipse_pair in eclipses {
                    point.eclipse.push(eclipse_pair);
                }
            }
            for point in &mut star1_coarse_grid {
                eclipses = disc_eclipse(
                    model.iangle.value,
                    rdisc1,
                    rdisc2,
                    model.beta_disc.value,
                    model.height_disc.value,
                    &point.position,
                )?;
                for eclipse_pair in eclipses {
                    point.eclipse.push(eclipse_pair);
                }
            }
            for point in &mut star2_fine_grid {
                eclipses = disc_eclipse(
                    model.iangle.value,
                    rdisc1,
                    rdisc2,
                    model.beta_disc.value,
                    model.height_disc.value,
                    &point.position,
                )?;
                for eclipse_pair in eclipses {
                    point.eclipse.push(eclipse_pair);
                }
            }
            for point in &mut star2_coarse_grid {
                eclipses = disc_eclipse(
                    model.iangle.value,
                    rdisc1,
                    rdisc2,
                    model.beta_disc.value,
                    model.height_disc.value,
                    &point.position,
                )?;
                for eclipse_pair in eclipses {
                    point.eclipse.push(eclipse_pair);
                }
            }
        }
        
        // Set the surface brightness of the disc
        set_disc_continuum(
            rdisc2,
            model.temp_disc.value,
            model.texp_disc.value,
            model.wavelength,
            &mut disc_grid,
        );

        // Set the surface brightness of outer edge, accounting for
        // irradiation by star 2
        set_edge_continuum(
            model.temp_edge.value,
            r2,
            model.t2.value.abs(),
            model.absorb_edge.value,
            model.wavelength,
            &mut disc_edge_grid,
        );
    }

    if model.add_spot {
        bright_spot_grid = set_bright_spot_grid(model)?;
    }
    
    Ok((
        star1_coarse_grid,
        star2_coarse_grid,
        star1_fine_grid,
        star2_fine_grid,
        disc_grid,
        disc_edge_grid,
        bright_spot_grid,
        gint,
        rlens1,
    ))
}

pub fn map_from_pydict(dict: Bound<'_, PyDict>) -> PyResult<HashMap<String, Entry>> {
    let mut map: HashMap<String, Entry> = HashMap::new();

    for (key, value) in dict.iter() {
        let name: String = key.extract::<String>()?;

        if let Ok(v) = value.extract::<f64>() {
            map.insert(name, Entry::Scalar(v.to_string()));
        } else if let Ok(v) = value.extract::<bool>() {
            map.insert(
                name,
                Entry::Scalar(if v { "1".to_string() } else { "0".to_string() }),
            );
        } else if let Ok(v) = value.extract::<String>() {
        map.insert(name, Entry::Scalar(v));
        } else {
            return Err(PyErr::new::<pyo3::exceptions::PyTypeError, _>(format!(
                "Unsupported type for key {}",
                name
            )));
        }
    }

    Ok(map)
}

pub fn rescale(flux: &[f64], flux_err: &[f64], weight: Option<&[f64]>, model_flux: &[f64]) -> f64 {
    let mut sdy: f64 = 0.0;
    let mut syy: f64 = 0.0;
    for i in 0..flux.len() {
        if flux_err[i] < 0.0 {
            continue;
        } else if let Some(weight) = weight {
                if weight[i] > 0.0 {
                let wgt = weight[i] / (flux_err[i] * flux_err[i]);
                sdy += wgt * flux[i] * model_flux[i];
                syy += wgt * model_flux[i] * model_flux[i];
                }
            } else {
                let wgt = 1.0 / (flux_err[i] * flux_err[i]);
            sdy += wgt * flux[i] * model_flux[i];
            syy += wgt * model_flux[i] * model_flux[i];
            }
        }
    let scale: f64 = sdy / syy;
    scale
}

pub fn chisq_log_prob(
    flux: &[f64],
    flux_err: &[f64],
    weight: Option<&[f64]>,
    model_flux: &[f64],
) -> (f64, f64) {
    let mut chisq_i: f64 = 0.0;
    let mut chisq_sum: f64 = 0.0;
    let mut log_prob: f64 = 0.0;
    for i in 0..flux.len() {
        if flux_err[i] < 0.0 {
            chisq_i = 0.0
        } else if let Some(weight) = weight {
                if weight[i] > 0.0 {
                chisq_i = weight[i] * ((flux[i] - model_flux[i]) / flux_err[i]).powi(2);
                } else {
                chisq_i = ((flux[i] - model_flux[i]) / flux_err[i]).powi(2);
                }
            }

        chisq_sum += chisq_i;
        log_prob += -0.5 * (chisq_i + (TAU * flux_err[i] * flux_err[i]).ln())
    }
    (chisq_sum, log_prob)
}
