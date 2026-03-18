Python extension of the lroche routine from Tom Marsh's LCURVE (cpp-lcurve) translated into Rust
===================================


.. code-block:: python

    import lroche

    model = lroche.BinaryModel.from_file("WDdM.mod")

load in your data in python and then convert each column to a contiguous array.

.. code-block:: python

    flux = np.ascontiguousarray(time)
    t_exp = np.ascontiguousarray(t_exp)
    flux = np.ascontiguousarray(flux)
    flux_err = np.ascontiguousarray(flux_err)
    weight = np.ascontiguousarray(weight)
    n_div = np.ascontiguousarray(n_div)

then create empty arrays for rust to fill

.. code-block:: python

    model_flux1 = np.ascontiguousarray(np.empty_like(time))
    model_flux2 = np.ascontiguousarray(np.empty_like(time))

And calculate the model

.. code-block:: python

    model.compute_light_curve(time, exp, flux, flux_err, weight, n_div, model_flux1, model_flux2)

    fig, ax = plt.subplots
    ax.plot(time, model_flux1+model_flux2)
    plt.show()