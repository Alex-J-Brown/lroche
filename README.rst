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


Then calculate the model

.. code-block:: python

    lc = model.compute_light_curve(time, exp, flux, flux_err, weight, n_div)

    fig, ax = plt.subplots
    ax.plot(time, lc.star1)
    ax.plot(time, lc.star2)
    ax.plot(time, lc.disc)
    ax.plot(time, lc.disc_edge)
    ax.plot(time, lc.bright_spot)
    ax.plot(time, lc.total)
    plt.show()