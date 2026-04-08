use pyo3::prelude::*;

pub mod binary_model;
pub mod comp_gravity;
pub mod comp_light;
pub mod comp_radius;
pub mod constants;
pub mod ginterp;
pub mod ldc;
pub mod model;
pub mod numface;
pub mod pparam;
pub mod set_bright_spot_grid;
pub mod set_disc_continuum;
pub mod set_disc_grid;
pub mod set_star_continuum;
pub mod set_star_grid;

#[pymodule]
fn lroche(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<binary_model::BinaryModel>()?;
    Ok(())
}
