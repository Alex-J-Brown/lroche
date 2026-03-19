use pyo3::prelude::*;

mod model;
mod vec3;
mod roche;
mod set_star_grid;
mod set_star_continuum;
mod set_disc_grid;
mod set_disc_continuum;
mod set_bright_spot_grid;
mod comp_light;
mod constants;
mod comp_gravity;
mod binary_model;

#[pymodule]
fn lroche(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<binary_model::BinaryModel>()?;
    Ok(())
}

// fn main() {
//     let model = Model::from_file("/home/alex/Astro/Analysis/photometry/eclipse_MCMC/model_files/WDdM.mod").unwrap();
//     let mut star1_fine_grid: Vec<model::Point> = set_star_grid::set_star_grid(&model, roche::Star::Primary, true);
//     let mut star2_fine_grid: Vec<model::Point> = set_star_grid::set_star_grid(&model, roche::Star::Secondary, true);
//     set_star_continuum::set_star_continuum(&model, &mut star1_fine_grid, &mut star2_fine_grid);
// }