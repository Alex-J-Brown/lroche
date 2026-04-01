Python extension of the lroche routine from Tom Marsh's LCURVE (cpp-lcurve) translated into Rust
===================================


.. code-block:: python

    import lroche

    model = lroche.BinaryModel.from_file("WDdM.mod")

load in your data in python and then convert each column to a contiguous array.

.. code-block:: python

    time, t_exp, flux, flux_err, weight, n_div = np.loadtxt("example_data.dat", unpack=True)
    time = np.ascontiguousarray(time)
    t_exp = np.ascontiguousarray(t_exp)
    flux = np.ascontiguousarray(flux)
    flux_err = np.ascontiguousarray(flux_err)
    weight = np.ascontiguousarray(weight)
    n_div = np.ascontiguousarray(n_div)


Then calculate the model. Note that flux, flux_err, and weight are optional. If flux and flux_err are supplied then the returned light curve model will be scaled to match the data. If weights are not given they will be assumed to be 1.

.. code-block:: python

    lc = model.compute_light_curve(time, t_exp, n_div, flux, flux_err, weight)

    fig, ax = plt.subplots()
    ax.plot(time, lc.star1)
    ax.plot(time, lc.star2)
    ax.plot(time, lc.disc)
    ax.plot(time, lc.disc_edge)
    ax.plot(time, lc.bright_spot)
    ax.plot(time, lc.total)
    plt.show()

Volume-averaged radii and flux-weighted log(g) calculated from the model can be accessed as attributes

.. code-block:: python

    print(lc.rva1)
    print(lc.rva2)
    print(lc.logg1)
    print(lc.logg2)

And if arrays of flux and flux_err were given then the chi-squared and log(probability) can also be accessed

.. code-block:: python

    print(lc.chi2)
    print(lc.log_prob)

If flux and flux_err are not given then these will both be None.

Once initialised, the model can be updated from a python dictionary and also accessed as a python dictionary.

.. code-block:: python

    model.update({"iangle": 85.5, "t1": 20000.0})
    current_parameters = model.model.to_dict()

If you are updating a parameter that was not previously 'defined', then defined will automatically  be set to 'true' when updated.
If you want to set it back to 'false' then the whole parameter details must be given e.g.

.. code-block:: python

    model.update({"radius_spot": {"value": 0.326, "range": 0.0, "dstep": 0.0, "vary": False, "defined": False}})

Note that `model.update` checks if the grid needs to be rebuilt depending on what parameters are given. This is useful for fitting multiband data where the grid will remain the same for all bands but the fluxes and wavelengths will change, saving time by only updating the fluxes of the gridpoints.

**Warning!** Currently there are no checks made on the models or updates to make sure these are valid. If an update is sent with an incorrect parameter name then it will not complain but also will not update the intended parameter.

The LCURVE fitting routines, `levmarq`, and `simplex` have not been implemented in lroche, however the relatively simple python interface means that fitting algorithms from scipy can be used easily.

.. code-block:: python

    from scipy.minimize import curve_fit, minimize


    def fit_function_minimize(params):
        model.update({"t0": params[0]})
        return model.compute_light_curve(time, exp, n_div, flux, flux_err, weight).chi2


    res = minimize(fit_function_minimize, [58674.006198], method='Nelder-Mead')
    print(f"t0 = {res.x[0]}")


    def fit_function_curve_fit(time, t0):
        model.update({"t0": t0})
        return model.compute_light_curve(time, exp, n_div, flux, flux_err, weight).total

    popt, pcov = curve_fit(fit_function_curve_fit, xdata=time, ydata=flux, p0=[58674.006198], sigma=flux_err, method='lm', xtol=1e-12)
    print(f"t0 = {popt[0]} +- {np.sqrt(np.diag(pcov))[0]}")
    
