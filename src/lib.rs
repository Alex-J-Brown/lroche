use pyo3::prelude::*;

pub mod pparam;
pub mod ldc;
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
    Ok(())
}
