use pyo3::prelude::*;

pub mod model;
pub mod numface;
pub mod set_star_grid;
pub mod set_star_continuum;
pub mod set_disc_grid;
pub mod set_disc_continuum;
pub mod set_bright_spot_grid;
pub mod ginterp;
pub mod comp_light;
pub mod constants;
pub mod comp_gravity;
pub mod comp_radius;
pub mod binary_model;


#[pymodule]
fn lroche(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<binary_model::BinaryModel>()?;
    m.add_function(wrap_pyfunction!(binary_model::x_l1, m)?)?;
    Ok(())
}

// fn main() {
//     let model = Model::from_file("/home/alex/Astro/Analysis/photometry/eclipse_MCMC/model_files/WDdM.mod").unwrap();
//     let mut star1_fine_grid: Vec<model::Point> = set_star_grid::set_star_grid(&model, roche::Star::Primary, true);
//     let mut star2_fine_grid: Vec<model::Point> = set_star_grid::set_star_grid(&model, roche::Star::Secondary, true);
//     set_star_continuum::set_star_continuum(&model, &mut star1_fine_grid, &mut star2_fine_grid);
// }