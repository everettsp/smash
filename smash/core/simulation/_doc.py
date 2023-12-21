from __future__ import annotations

from textwrap import dedent

from smash._constant import (
    DEFAULT_SIMULATION_COST_OPTIONS,
    DEFAULT_SIMULATION_COMMON_OPTIONS,
    DEFAULT_SIMULATION_RETURN_OPTIONS,
)
from smash.util._doctools import DocAppender, DocSubstitution

# TODO: mapping and optimizer arguments are duplicated between each docstrings.
# Maybe create intermediate variables to store this with or without ann optimizer.
# Idem for control_info between bayesian_optimize and optimize

OPTIMIZE_OPTIONS_BASE_DOC = {
    "parameters": (
        """
        `str`, `list[str, ...]` or None, default None
        """,
        """
        Name of parameters to optimize. Should be one or a sequence of any key of:

        - `Model.rr_parameters`
        - `Model.rr_initial_states`
        %(parameters_serr_mu_parameters)s
        %(parameters_serr_sigma_parameters)s

        >>> optimize_options = {
            "parameters": "cp",
        }
        >>> optimize_options = {
            "parameters": ["cp", "ct", "kexc", "llr"],
        }

        .. note::
            If not given, all parameters in `Model.rr_parameters`%(parameters_note_serr_parameters)s will be optimized.
        """,
    ),
    "bounds": (
        """
        `dict[str, tuple[float, float]]` or None, default None
        """,
        """
        Bounds on optimized parameters. A dictionary where the keys represent parameter names, and the values are pairs of ``(min, max)`` 
        values (i.e., a list or tuple) with ``min`` lower than ``max``. The keys must be included in **parameters**.

        >>> optimize_options = {
            "bounds": {
                "cp": (1, 2000),
                "ct": (1, 1000),
                "kexc": (-10, 5), 
                "llr": (1, 1000)
            },
        }

        .. note::
            If not given, default bounds will be applied to each parameter. 
            See `Model.get_rr_parameters_bounds`, `Model.get_rr_initial_states_bounds`%(bounds_get_serr_parameters_bounds)s
        """,
    ),
    "control_tfm": (
        """
        `str` or None, default None
        """,
        """
        Transformation method applied to the control vector. Should be one of:
        
        - ``'keep'`` (all **optimizer**)
        - ``'normalize'`` (``'sbs'`` or ``'lbfgsb'`` **optimizer** only)
        - ``'sbs'`` (``'sbs'`` **optimizer** only)

        .. note::
            If not given, a default control vector transformation will be set depending on the optimizer:

            - **optimizer** = ``'sbs'``; **control_tfm** = ``'sbs'``
            - **optimizer** = ``'lbfgsb'``; **control_tfm** = ``'normalize'``
            - **optimizer** = ``'ann'``; **control_tfm** = ``'keep'``

        .. hint::
            See a detailed explanation on the control vector transformation in (TODO FC: link Math/Num) section.
        """,
    ),
    "descriptor": (
        """
        `dict[str, list[str, ...]]` or None, default None
        """,
        """
        Descriptors linked to optimized parameters. A dictionary where the keys represent parameter names, and the values are list of descriptor names.
        The keys must be included in **parameters**.

        >>> optimize_options = {
            "descriptor": {
                "cp": ["slope", "dd"],
                "ct": ["slope"],
                "kexc": ["slope", "dd"],
                "llr": ["dd"],
            },
        }

        .. note::
            If not given, all descriptors will be used for each parameter.
            This option is only be used when **mapping** is ``'multi-linear'`` or ``'multi-polynomial'``. In case of ``'ann'``, all descriptors will be used.
        """,
    ),
    "net": (
        """
        `Net <smash.factory.Net>` or None, default None
        """,
        """
        The neural network used to learn the descriptors-to-parameters mapping.
        
        .. note::
            If not given, a default neural network will be used. This option is only used when **mapping** is ``'ann'``.
            See `Net <smash.factory.Net>` to learn how to create a customized neural network for training.
        """,
    ),
    "learning_rate": (
        """
        `float` or None, default None
        """,
        """
        The learning rate used for weight updates during training.

        .. note::
            If not given, a default learning rate will be used. This option is only used when **mapping** is ``'ann'``.
        """,
    ),
    "random_state": (
        """
        `int` or None, default None
        """,
        """
        A random seed used to initialize neural network weights. 
        
        .. note::
            If not given, the weights will be initialized with a random seed. This options is only used when **mapping** is ``'ann'``.
        """,
    ),
    "termination_crit": (
        """
        `dict[str, Any]` or None, default None
        """,
        """
        Termination criteria. The elements are:

        - ``'maxiter'``: The maximum number of iterations. Only used when **optimizer** is ``'sbs'`` or ``'lbfgsb'``.
        - ``'factr'``: An additional termination criterion based on cost values. Only used when **optimizer** is ``'lbfgsb'``.
        - ``'pgtol'``: An additional termination criterion based on the projected gradient of the cost function. Only used when **optimizer** is ``'lbfgsb'``.
        - ``'epochs'``: The number of training epochs for the neural network. Only used when **mapping** is ``'ann'``.
        - ``'early_stopping'``: A positive number to stop training when the loss function does not decrease below the current optimal value for **early_stopping** consecutive epochs. When set to zero, early stopping is disabled, and the training continues for the full number of epochs. Only used when **mapping** is ``'ann'``.
        
        >>> optimize_options = {
            "termination_crit": {
                "maxiter": 10,
                "factr": 1e6,
            },
        }
        >>> optimize_options = {
            "termination_crit": {
                "epochs": 200,
            },
        }

        .. note::
            If not given, default values are set to each elements.
        """,
    ),
}

COST_OPTIONS_BASE_DOC = {
    "jobs_cmpt": (
        """
        `str` or `list[str, ...]`, default 'nse'
        """,
        """
        Type of observation objective function(s) to be computed. Should be one or a sequence of any of

        - ``'nse'``, ``'nnse'``, ``'kge'``, ``'mae'``, ``'mape'``, ``'mse'``, ``'rmse'``, ``'lgrm'`` (classical evaluation metrics)
        - ``'Crc'``, ``'Crchf'``, ``'Crclf'``, ``'Crch2r'``, ``'Cfp2'``, ``'Cfp10'``, ``'Cfp50'``, ``'Cfp90'`` (continuous signatures-based error metrics)
        - ``'Eff'``, ``'Ebf'``, ``'Erc'``, ``'Erchf'``, ``'Erclf'``, ``'Erch2r'``, ``'Elt'``, ``'Epf'`` (flood event signatures-based error metrics)

        >>> cost_options = {
            "jobs_cmpt": "nse",
        }
        >>> cost_options = {
            "jobs_cmpt": ["nse", "Epf"],
        }

        .. hint::
            See a detailed explanation on the objective function in :ref:`Math / Num Documentation <math_num_documentation.signal_analysis.cost_functions>` section.
        """,
    ),
    "jobs_cmpt_tfm": (
        """
        `str` or `list[str, ...]`, default 'keep'
        """,
        """
        Type of transformation applied to discharge in observation objective function(s). Should be one or a sequence of any of

        - ``'keep'`` : No transformation :math:`f:x \\rightarrow x`
        - ``'sqrt'`` : Square root transformation :math:`f:x \\rightarrow \sqrt{x}`
        - ``'inv'`` : Multiplicative inverse transformation :math:`f:x \\rightarrow \\frac{1}{x}`

        >>> cost_options = {
            "jobs_cmpt_tfm": "inv",
        }
        >>> cost_options = {
            "jobs_cmpt_tfm": ["keep", "inv"],
        }

        .. note::
            If **jobs_cmpt** is a multi-criteria and only one transformation is choosen in **jobs_cmpt_tfm**. The transformation will be applied to each 
            observation objective function.
        """,
    ),
    "wjobs_cmpt": (
        """
        `str` or `list[float, ...]`, default 'mean'
        """,
        """
        The corresponding weighting of observation objective functions in case of multi-criteria 
        (i.e., a sequence of objective functions to compute). There are two ways to specify it:

        - An alias among ``'mean'``
        - A sequence of value whose size must be equal to the number of observation objective function(s) in **jobs_cmpt**

        >>> cost_options = {
            "wjobs_cmpt": "mean",
        }
        >>> cost_options = {
            "wjobs_cmpt": [0.7, 0.3],
        }
        """,
    ),
    "wjreg": (
        """
        `float` or `str`, default 0
        """,
        """
        The weighting of regularization term. There are two ways to specify it:

        - A value greater than or equal to 0
        - An alias among ``'fast'`` or ``'lcurve'``. **wjreg** will be auto-computed by one of these methods.

        >>> cost_options = {
            "wjreg": 1e-4,
        }
        >>> cost_options = {
            "wjreg": "lcurve",
        }

        .. hint::
            See a detailed explanation on the weighting of regularization term in (TODO FC: link Math/Num) section.
        """,
    ),
    "jreg_cmpt": (
        """
        `str` or `list[str, ...]`, default 'prior'
        """,
        """
        Type(s) of regularization function(s) to be minimized when regularization term is set (i.e., **wjreg** > 0). Should be one or a sequence of any of

        - ``'prior'``
        - ``'smoothing'``
        - ``'hard-smoothing'``

        >>> cost_options = {
            "jreg_cmpt": "prior",
        }
        >>> cost_options = {
            "jreg_cmpt": ["prior", "smoothing"],
        }

        .. hint::
            See a detailed explanation on the regularization function in (TODO FC: link Math/Num) section.
        """,
    ),
    "wjreg_cmpt": (
        """
        `str` or `list[float, ...]`, default 'mean'
        """,
        """
        The corresponding weighting of regularization functions in case of multi-regularization
        (i.e., a sequence of regularization functions to compute). There are two ways to specify it:

        - An alias among ``'mean'``
        - A sequence of value whose size must be equal to the number of regularization function(s) in **jreg_cmpt**

        >>> cost_options = {
            "wjreg_cmpt": "mean",
        }
        >>> cost_options = {
            "wjreg_cmpt": [1., 2.],
        }
        """,
    ),
    "gauge": (
        """
        `str` or `list[str, ...]`, default 'dws'
        """,
        """
        Type of gauge to be computed. There are two ways to specify it:

        - An alias among ``'all'`` (all gauge codes) or ``'dws'`` (most downstream gauge code(s))
        - A gauge code or any sequence of gauge codes. The gauge code(s) given must belong to the gauge codes defined in the `Model.mesh`

        >>> cost_options = {
            "gauge": "dws",
        }
        >>> cost_options = {
            "gauge": "V3524010",
        }
        >>> cost_options = {
            "gauge": ["V3524010", "V3515010"],
        }
        """,
    ),
    "wgauge": (
        """
        `str` or `list[float, ...]` default 'mean'        
        """,
        """
        Type of gauge weights. There are two ways to specify it:

        - An alias among ``'mean'``, ``'lquartile'`` (1st quantile or lower quantile), ``'median'``, or ``'uquartile'`` (3rd quantile or upper quantile)
        - A sequence of value whose size must be equal to the number of gauges optimized in **gauge**

        >>> cost_options = {
            "wgauge": "mean",
        }
        >>> cost_options = {
            "wgauge": [0.6, 0.4]",
        }
        """,
    ),
    "control_prior": (
        """
        `dict[str, list[str, list[float]]]` or None, default None
        """,
        """
        Prior applied to the control vector.
        A dictionary containing the type of prior to link to control vector. The keys are any control parameter name (i.e. ``'cp0'``, ``'cp1-1'``, ``'cp-slope-a'``, etc.), 
        see `bayesian_optimize_control_info <smash.bayesian_optimize_control_info>` to retrieve control parameters names.
        The values are list of length 2 containing distribution information (i.e. distribution name and parameters). 
        Below, the set of available distributions and the associated number of parameters:

        - ``'FlatPrior'``,   []                                 (0)
        - ``'Uniform'``,     [lower_bound, higher_bound]        (2)
        - ``'Gaussian'``,    [mean, standard_deviation]         (2)
        - ``'Exponential'``, [threshold, scale]                 (2)
        - ``'LogNormal'``,   [mean_log, standard_deviation_log] (2)
        - ``'Triangle'``,    [peak, lower_bound, higher_bound]  (3)

        >>> cost_options = {
            control_prior: {
                "cp0": ["Gaussian", [200, 100]],
                "kexc0": ["Gaussian", [0, 5]],
            }
        }

        .. note::
            If not given, ``'FlatPrior'`` is applied to each control vector parameter (i.e. equivalent to no prior).

        .. hint::
            See a detailed explanation on the distribution in (TODO BR: add link Math/Num) section.
        """,
    ),
    "event_seg": (
        """
        `dict[str, float]`, default {'peak_quant': 0.995, 'max_duration': 240}
        """,
        """
        A dictionary of event segmentation options when calculating flood event signatures for cost computation (i.e., **jobs_cmpt** includes flood events signatures).

        >>> cost_options = {
            event_seg = {
                "peak_quant": 0.998,
                "max_duration": 120,
            }
        }

        .. hint::
            See `hydrograph_segmentation <smash.hydrograph_segmentation>` for more.

        """,
    ),
    "end_warmup": (
        """
        `str`, `pandas.Timestamp` or None, default None
        """,
        """
        The end of the warm-up period, which must be between the start time and the end time defined in `Model.setup`.

        >>> cost_options = {
            "end_warmup": "1997-12-21",
        }
        >>> cost_options = {
            "end_warmup": pd.Timestamp("19971221"),
        }
        
        .. note::
            If not given, it is set to be equal to the `Model.setup` start time.
        """,
    ),
}

COMMON_OPTIONS_BASE_DOC = {
    "verbose": (
        """
        `bool`, default False
        """,
        """
        Whether to display information about the running method.
        """,
    ),
    "ncpu": (
        """
        `int`, default 1
        """,
        """
        Number of CPU(s) to perform a parallel computation.
        """,
    ),
}

RETURN_OPTIONS_BASE_DOC = {
    "time_step": (
        """
        `str`, `pandas.Timestamp`, `pandas.DatetimeIndex` or `list[str, ...]`, default 'all'
        """,
        """
        Returned time steps. There are five ways to specify it:

        - An alias among ``'all'`` (return all time steps).
        - A date as string which respect `pandas.Timestamp` format
        - A `pandas.Timestamp`.
        - A `pandas.DatetimeIndex`.
        - A sequence of dates as strings.

        >>> return_options = {
            "time_step": "all",
        }
        >>> return_options = {
            "time_step": "1997-12-21",
        }
        >>> return_options = {
            "time_step": pd.Timestamp("19971221"),
        }
        >>> return_options = {
            "time_step": pd.date_range(
                start="1997-12-21",
                end="1998-12-21",
                freq="1D"
            ),
        }
        >>> return_options = {
            "time_step": ["1998-05-23", "1998-05-24", "1998-05-25"],
        }

        .. note::
            It only applies to the following variables: ``'rr_states'`` and ``'q_domain'``
        """,
    ),
    "rr_states": (
        """
        `bool`, default False
        """,
        """
        Whether to return rainfall-runoff states for specific time steps.
        """,
    ),
    "q_domain": (
        """
        `bool`, default False
        """,
        """
        Whether to return simulated discharge on the whole domain for specific time steps.
        """,
    ),
    "iter_cost": (
        """
        `bool`, default False
        """,
        """
        Whether to return cost iteration values.
        """,
    ),
    "iter_projg": (
        """
        `bool`, default False
        """,
        """
        Whether to return infinity norm of the projected gardient iteration values.
        """,
    ),
    "control_vector": (
        """
        `bool`, default False
        """,
        """
        Whether to return control vector at end of optimization. In case of optimization with ``'ann'`` **mapping**, the control vector is represented in `Net.layers <smash.factory.Net.layers>` instead.
        """,
    ),
    "net": (
        """
        `bool`, default False
        """,
        """
        Whether to return the trained neural network `Net <smash.factory.Net>`. Only used with ``'ann'`` **mapping**.
        """,
    ),
    "cost": (
        """
        `bool`, default False
        """,
        """
        Whether to return cost value.
        """,
    ),
    "jobs": (
        """
        `bool`, default False
        """,
        """
        Whether to return jobs (observation component of cost) value.
        """,
    ),
    "jreg": (
        """
        `bool`, default False
        """,
        """
        Whether to return jreg (regularization component of cost) value.
        """,
    ),
    "lcurve_wjreg": (
        """
        `bool`, default False
        """,
        """
        Whether to return the wjreg lcurve. Only used if **wjreg** in cost_options is equal to ``'lcurve'``.
        """,
    ),
    "lcurve_multiset": (
        """
        `bool`, default False
        """,
        """
        Whether to return the multiset estimate L-curve.
        """,
    ),
    "log_lkh": (
        """
        `bool`, default False
        """,
        """
        Whether to return log likelihood component value.
        """,
    ),
    "log_prior": (
        """
        `bool`, default False
        """,
        """
        Whether to return log prior component value.
        """,
    ),
    "log_h": (
        """
        `bool`, default False
        """,
        """
        Whether to return log h component value.
        """,
    ),
    "serr_mu": (
        """
        `bool`, default False
        """,
        """
        Whether to return mu, the mean of structural errors. It can also be returned directly from the Model object using the `Model.get_serr_mu` method.
        """,
    ),
    "serr_sigma": (
        """
        `bool`, default False
        """,
        """
        Whether to return sigma, the standard deviation of structural errors. It can also be returned directly from the Model object using the `Model.get_serr_sigma` method.
        """,
    ),
}


def _gen_docstring_from_base_doc(
    base_doc: dict[str, tuple(str, str)], keys: list[str], nindent: int = 0
) -> str:
    ret = []
    for key in keys:
        typ = dedent(base_doc[key][0].strip())
        description = dedent(base_doc[key][1].strip())
        text = f"{key} : {typ}\n\t{description}\n"

        ret.append(nindent * "    " + text)
    return "\n".join(ret)


_forward_run_doc = (
    """
Run the forward Model.

.. hint::
    See a detailed explanation on the forward run usage in the (TODO FC: link user guide) section.

Parameters
----------
%(model_parameter)s
cost_options : `dict[str, Any]` or None, default None
    Dictionary containing computation cost options for simulated and observed responses. The elements are:

"""
    + _gen_docstring_from_base_doc(
        COST_OPTIONS_BASE_DOC,
        DEFAULT_SIMULATION_COST_OPTIONS["forward_run"].keys(),
        nindent=1,
    )
    + """
common_options : `dict[str, Any]` or None, default None
    Dictionary containing common options with two elements:

"""
    + _gen_docstring_from_base_doc(
        COMMON_OPTIONS_BASE_DOC, DEFAULT_SIMULATION_COMMON_OPTIONS.keys(), nindent=1
    )
    + """
return_options : `dict[str, Any]` or None, default None
    Dictionary containing return options to save intermediate variables. The elements are:

"""
    + _gen_docstring_from_base_doc(
        RETURN_OPTIONS_BASE_DOC,
        DEFAULT_SIMULATION_RETURN_OPTIONS["forward_run"].keys(),
        nindent=1,
    )
    + """
Returns
-------
%(model_return)s
forward_run : `ForwardRun <smash.ForwardRun>` or None, default None
    It returns an object containing the intermediate variables defined in **return_options**.
    If no intermediate variables are defined, it returns None.

See Also
--------
ForwardRun : Represents forward run optional results.

Examples
--------
>>> from smash.factory import load_dataset
>>> setup, mesh = load_dataset("cance")
>>> model = smash.Model(setup, mesh)

Run the direct Model

>>> %(model_example_func)s
</> Forward Run

Get the simulated discharges

>>> %(model_example_response)s.response.q
array([[1.9826430e-03, 1.3466669e-07, 6.7617895e-12, ..., 3.2273201e+01,
        3.2118713e+01, 3.1965160e+01],
       [2.3777038e-04, 7.3761623e-09, 1.7551447e-13, ..., 7.9022121e+00,
        7.8704414e+00, 7.8388391e+00],
       [2.9721676e-05, 5.4272520e-10, 8.4623445e-15, ..., 2.0933011e+00,
        2.0847433e+00, 2.0762112e+00]], dtype=float32)
"""
)

_optimize_doc = (
    """
Model assimilation using numerical optimization algorithms.

.. hint::
    See a detailed explanation on the optimize usage in (TODO FC: link user guide) section. 

Parameters
----------
%(model_parameter)s
mapping : `str`, default 'uniform'
    Type of mapping. Should be one of 
    
    - ``'uniform'`` 
    - ``'distributed'``
    - ``'multi-linear'``
    - ``'multi-polynomial'``
    - ``'ann'``

    .. hint::
        See a detailed explanation on the mapping in (TODO FC: link Math/Num) section.

optimizer : `str` or None, default None
    Name of optimizer. Should be one of 
    
    - ``'sbs'`` (``'uniform'`` **mapping** only)
    - ``'lbfgsb'`` (``'uniform'``, ``'distributed'``, ``'multi-linear'`` or ``'multi-polynomial'`` **mapping** only)
    - ``'sgd'`` (``'ann'`` **mapping** only)
    - ``'adam'`` (``'ann'`` **mapping** only)
    - ``'adagrad'`` (``'ann'`` **mapping** only)
    - ``'rmsprop'`` (``'ann'`` **mapping** only)

    .. note::
        If not given, a default optimizer will be set depending on the optimization mapping:

        - **mapping** = ``'uniform'``; **optimizer** = ``'sbs'``
        - **mapping** = ``'distributed'``, ``'multi-linear'``, or ``'multi-polynomial'``; **optimizer** = ``'lbfgsb'``
        - **mapping** = ``'ann'``; **optimizer** = ``'adam'``

    .. hint::
        See a detailed explanation on the optimizer in (TODO FC: link Math/Num) section.

optimize_options : `dict[str, Any]` or None, default None
    Dictionary containing optimization options for fine-tuning the optimization process. 
    See `%(default_optimize_options_func)s` to retrieve the default optimize options based on the **mapping** and **optimizer**.

"""
    + _gen_docstring_from_base_doc(
        OPTIMIZE_OPTIONS_BASE_DOC,
        [
            "parameters",
            "bounds",
            "control_tfm",
            "descriptor",
            "net",
            "learning_rate",
            "random_state",
            "termination_crit",
        ],
        nindent=1,
    )
    + """ 
cost_options : `dict[str, Any]` or None, default None
    Dictionary containing computation cost options for simulated and observed responses. The elements are:

"""
    + _gen_docstring_from_base_doc(
        COST_OPTIONS_BASE_DOC,
        DEFAULT_SIMULATION_COST_OPTIONS["optimize"].keys(),
        nindent=1,
    )
    + """
common_options : `dict[str, Any]` or None, default None
    Dictionary containing common options with two elements:

"""
    + _gen_docstring_from_base_doc(
        COMMON_OPTIONS_BASE_DOC, DEFAULT_SIMULATION_COMMON_OPTIONS.keys(), nindent=1
    )
    + """
return_options : `dict[str, Any]` or None, default None
    Dictionary containing return options to save intermediate variables. The elements are:

"""
    + _gen_docstring_from_base_doc(
        RETURN_OPTIONS_BASE_DOC,
        DEFAULT_SIMULATION_RETURN_OPTIONS["optimize"].keys(),
        nindent=1,
    )
    + """
Returns
-------
%(model_return)s
optimize : `Optimize <smash.Optimize>` or None, default None
    It returns an object containing the intermediate variables defined in **return_options**.
    If no intermediate variables are defined, it returns None.

See Also
--------
Optimize : Represents optimize optional results.

Examples
--------
>>> from smash.factory import load_dataset
>>> setup, mesh = load_dataset("cance")
>>> model = smash.Model(setup, mesh)

Optimize the Model

>>> %(model_example_func)s
</> Optimize
    At iterate      0    nfg =     1    J =      0.643190    ddx = 0.64
    At iterate      1    nfg =    30    J =      0.097397    ddx = 0.64
    At iterate      2    nfg =    59    J =      0.052158    ddx = 0.32
    At iterate      3    nfg =    88    J =      0.043086    ddx = 0.08
    At iterate      4    nfg =   118    J =      0.040684    ddx = 0.02
    At iterate      5    nfg =   152    J =      0.040604    ddx = 0.01
    CONVERGENCE: DDX < 0.01

Get the simulated discharges

>>> %(model_example_response)s.response.q
array([[1.9826430e-03, 1.3466669e-07, 6.7617895e-12, ..., 3.2273201e+01,
        3.2118713e+01, 3.1965160e+01],
       [2.3777038e-04, 7.3761623e-09, 1.7551447e-13, ..., 7.9022121e+00,
        7.8704414e+00, 7.8388391e+00],
       [2.9721676e-05, 5.4272520e-10, 8.4623445e-15, ..., 2.0933011e+00,
        2.0847433e+00, 2.0762112e+00]], dtype=float32)
"""
)

_multiset_estimate_doc = (
    """
Model assimilation using Bayesian-like estimation on multiple sets of solutions.

.. hint::
    See a detailed explanation on the multiset estimate usage in (TODO FC: link user guide) section.

Parameters
----------
%(model_parameter)s

multiset : `MultipleForwardRun <smash.MultipleForwardRun>` or `MultipleOptimize <smash.MultipleOptimize>`
    The returned object created by `multiple_forward_run <smash.multiple_forward_run>` or `multiple_optimize <smash.multiple_optimize>` 
    method containing information about multiple sets of rainfall-runoff parameters or initial states.

alpha : `float`, `list[float, ...]`, or None, default None
    A regularization parameter that controls the decay rate of the likelihood function.
    If **alpha** is a list, the L-curve approach will be used to find an optimal value for the regularization parameter.

    .. note::
        If not given, a default numeric range will be set for optimization through the L-curve process.

common_options : `dict[str, Any]` or None, default None
    Dictionary containing common options with two elements:

"""
    + _gen_docstring_from_base_doc(
        COMMON_OPTIONS_BASE_DOC,
        DEFAULT_SIMULATION_COMMON_OPTIONS.keys(),
        nindent=1,
    )
    + """
return_options : `dict[str, Any]` or None, default None
    Dictionary containing return options to save intermediate variables. The elements are:

"""
    + _gen_docstring_from_base_doc(
        RETURN_OPTIONS_BASE_DOC,
        DEFAULT_SIMULATION_RETURN_OPTIONS["multiset_estimate"].keys(),
        nindent=1,
    )
    + """
Returns
-------
%(model_return)s

multiset_estimate : `MultisetEstimate <smash.MultisetEstimate>` or None, default None
    It returns an object containing the intermediate variables defined in **return_options**. If no intermediate variables are defined, it returns None.

See Also
--------
MultisetEstimate : Represents multiset estimate optional results.
MultipleForwardRun : Represents multiple forward run computation result.
MultipleOptimize : Represents multiple optimize computation result.

Examples
--------
>>> from smash.factory import load_dataset
>>> from smash.factory import generate_samples
>>> setup, mesh = load_dataset("cance")
>>> model = smash.Model(setup, mesh)

Define sampling problem and generate samples

>>> problem = {
...    'num_vars': 4,
...    'names': ['cp', 'ct', 'kexc', 'llr'],
...    'bounds': [[1, 2000], [1, 1000], [-20, 5], [1, 1000]]
... }
>>> sr = generate_samples(problem, n=100, random_state=11)

Run Model with multiple sets of parameters

>>> mfr = smash.multiple_forward_run(model, samples=sr)
</> Multiple Forward Run
    Forward Run 100/100 (100%(percent)s)

Estimate Model on multiple sets of solutions

>>> %(model_example_func)s
</> Multiple Set Estimate
    L-curve Computing: 100%(percent)s|████████████████████████████████████████| 50/50 [00:02<00:00, 17.75it/s]

Get the simulated discharges

>>> %(model_example_response)s.response.q
array([[9.4456947e-05, 9.3808041e-05, 9.3033530e-05, ..., 1.7851851e+01,
        1.7691467e+01, 1.7528925e+01],
       [2.5728115e-05, 2.4717448e-05, 2.3537395e-05, ..., 3.9668279e+00,
        3.9301457e+00, 3.8938913e+00],
       [5.9592148e-06, 5.0831463e-06, 4.2117308e-06, ..., 9.8274386e-01,
        9.7398454e-01, 9.6534210e-01]], dtype=float32)
"""
)

_bayesian_optimize_doc = (
    """
Model bayesian assimilation using numerical optimization algorithms.

.. hint::
    See a detailed explanation on the bayesian optimize usage in (TODO FC: link user guide) section. 

Parameters
----------
%(model_parameter)s
mapping : `str`, default 'uniform'
    Type of mapping. Should be one of 
    
    - ``'uniform'`` 
    - ``'distributed'``
    - ``'multi-linear'``
    - ``'multi-polynomial'``

    .. hint::
        See a detailed explanation on the mapping in (TODO FC: link Math/Num) section.

optimizer : `str` or None, default None
    Name of optimizer. Should be one of 
    
    - ``'sbs'`` (``'uniform'`` **mapping** only)
    - ``'lbfgsb'`` (``'uniform'``, ``'distributed'``, ``'multi-linear'`` or ``'multi-polynomial'`` **mapping** only)

    .. note::
        If not given, a default optimizer will be set depending on the optimization mapping:

        - **mapping** = ``'uniform'``; **optimizer** = ``'sbs'``
        - **mapping** = ``'distributed'``, ``'multi-linear'``, or ``'multi-polynomial'``; **optimizer** = ``'lbfgsb'``

optimize_options : `dict[str, Any]` or None, default None
    Dictionary containing optimization options for fine-tuning the optimization process. 
    See `%(default_optimize_options_func)s` to retrieve the default optimize options based on the **mapping** and **optimizer**.

"""
    + _gen_docstring_from_base_doc(
        OPTIMIZE_OPTIONS_BASE_DOC,
        [
            "parameters",
            "bounds",
            "control_tfm",
            "descriptor",
            "termination_crit",
        ],
        nindent=1,
    )
    + """ 
cost_options : `dict[str, Any]` or None, default None
    Dictionary containing computation cost options for simulated and observed responses. The elements are:

"""
    + _gen_docstring_from_base_doc(
        COST_OPTIONS_BASE_DOC,
        DEFAULT_SIMULATION_COST_OPTIONS["bayesian_optimize"].keys(),
        nindent=1,
    )
    + """
common_options : `dict[str, Any]` or None, default None
    Dictionary containing common options with two elements:

"""
    + _gen_docstring_from_base_doc(
        COMMON_OPTIONS_BASE_DOC, DEFAULT_SIMULATION_COMMON_OPTIONS.keys(), nindent=1
    )
    + """
return_options : `dict[str, Any]` or None, default None
    Dictionary containing return options to save intermediate variables. The elements are:

"""
    + _gen_docstring_from_base_doc(
        RETURN_OPTIONS_BASE_DOC,
        DEFAULT_SIMULATION_RETURN_OPTIONS["bayesian_optimize"].keys(),
        nindent=1,
    )
    + """
Returns
-------
%(model_return)s
bayesian_optimize : `BayesianOptimize <smash.BayesianOptimize>` or None, default None
    It returns an object containing the intermediate variables defined in **return_options**.
    If no intermediate variables are defined, it returns None.

See Also
--------
BayesianOptimize : Represents bayesian optimize optional results.

Examples
--------
>>> from smash.factory import load_dataset
>>> setup, mesh = load_dataset("cance")
>>> model = smash.Model(setup, mesh)

Optimize the Model

>>> %(model_example_func)s
</> Bayesian Optimize
    At iterate      0    nfg =     1    J =     26.510803    ddx = 0.64
    At iterate      1    nfg =    68    J =      2.536702    ddx = 0.64
    At iterate      2    nfg =   136    J =      2.402311    ddx = 0.16
    At iterate      3    nfg =   202    J =      2.329653    ddx = 0.16
    At iterate      4    nfg =   270    J =      2.277469    ddx = 0.04
    At iterate      5    nfg =   343    J =      2.271495    ddx = 0.02
    At iterate      6    nfg =   416    J =      2.270596    ddx = 0.01
    At iterate      7    nfg =   488    J =      2.269927    ddx = 0.01
    At iterate      8    nfg =   561    J =      2.269505    ddx = 0.01
    CONVERGENCE: DDX < 0.01

Get the simulated discharges:

>>> %(model_example_response)s.response.q
array([[2.8263923e-04, 2.6972688e-04, 2.5154054e-04, ..., 2.2405909e+01,
        2.2121796e+01, 2.1843649e+01],
       [7.0431546e-05, 5.6386390e-05, 4.1599422e-05, ..., 5.3495941e+00,
        5.2908549e+00, 5.2333012e+00],
       [1.3507083e-05, 7.9051424e-06, 4.3033124e-06, ..., 1.3771975e+00,
        1.3625814e+00, 1.3482267e+00]], dtype=float32)
"""
)

_multiple_forward_run_doc = (
    """
Run the forward Model with multiple sets of parameters.

.. hint::
    See a detailed explanation on the multiple forward run usage in (TODO FC: link user guide) section.

Parameters
----------
model : `Model <smash.Model>`
    Primary data structure of the hydrological model `smash`.

samples : `Samples <smash.Samples>`
    Represents the generated samples result.

cost_options : `dict[str, Any]` or None, default None
    Dictionary containing computation cost options for simulated and observed responses. The elements are:

"""
    + _gen_docstring_from_base_doc(
        COST_OPTIONS_BASE_DOC,
        DEFAULT_SIMULATION_COST_OPTIONS["forward_run"].keys(),
        nindent=1,
    )
    + """
common_options : `dict[str, Any]` or None, default None
    Dictionary containing common options with two elements:

"""
    + _gen_docstring_from_base_doc(
        COMMON_OPTIONS_BASE_DOC, DEFAULT_SIMULATION_COMMON_OPTIONS.keys(), nindent=1
    )
    + """

Returns
-------
multiple_forward_run : `MultipleForwardRun <smash.MultipleForwardRun>`
    It returns an object containing the results of the multiple forward run.

See Also
--------
Samples : Represents the generated samples result.
MultipleForwardRun : Represents the multiple forward run result.

Examples
--------
>>> from smash.factory import load_dataset
>>> from smash.factory import generate_samples
>>> setup, mesh = load_dataset("cance")
>>> model = smash.Model(setup, mesh)

Define sampling problem and generate samples

>>> problem = {
...            'num_vars': 4,
...            'names': ['cp', 'ct', 'kexc', 'llr'],
...            'bounds': [[1, 2000], [1, 1000], [-20, 5], [1, 1000]]
... }
>>> sr = generate_samples(problem, n=5, random_state=11)

Run Model with multiple sets of parameters

>>> mfr = smash.multiple_forward_run(model, samples=sr)
</> Multiple Forward Run
    Forward Run 5/5 (100%)

Get the cost values through multiple forward runs

>>> mfr.cost
array([1.17112  , 1.0390087, 1.2135248, 1.2491335, 1.2172333], dtype=float32)
    """
)

_multiple_optimize_doc = (
    """
Run multiple optimization processes with multiple sets of parameters (i.e. starting points), yielding multiple solutions.

.. hint::
    See a detailed explanation on the multiple optimize usage in (TODO FC: link user guide) section.

Parameters
----------
model : `Model <smash.Model>`
    Primary data structure of the hydrological model `smash`.

samples : `Samples <smash.Samples>`
    Represents the generated samples result.

mapping : `str`, default 'uniform'
    Type of mapping. Should be one of 
    
    - ``'uniform'`` 
    - ``'distributed'``
    - ``'multi-linear'``
    - ``'multi-polynomial'``

    .. hint::
        See a detailed explanation on the mapping in (TODO FC: link Math/Num) section.

optimizer : `str` or None, default None
    Name of optimizer. Should be one of 
    
    - ``'sbs'`` (``'uniform'`` **mapping** only)
    - ``'lbfgsb'`` (``'uniform'``, ``'distributed'``, ``'multi-linear'`` or ``'multi-polynomial'`` **mapping** only)

    .. note::
        If not given, a default optimizer will be set depending on the optimization mapping:

        - **mapping** = ``'uniform'``; **optimizer** = ``'sbs'``
        - **mapping** = ``'distributed'``, ``'multi-linear'``, or ``'multi-polynomial'``; **optimizer** = ``'lbfgsb'``

    .. hint::
        See a detailed explanation on the optimizer in (TODO FC: link Math/Num) section.

optimize_options : `dict[str, Any]` or None, default None
    Dictionary containing optimization options for fine-tuning the optimization process. 
    See `%(default_optimize_options_func)s` to retrieve the default optimize options based on the **mapping** and **optimizer**.

"""
    + _gen_docstring_from_base_doc(
        OPTIMIZE_OPTIONS_BASE_DOC,
        [
            "parameters",
            "bounds",
            "control_tfm",
            "descriptor",
            "termination_crit",
        ],
        nindent=1,
    )
    + """ 
cost_options : `dict[str, Any]` or None, default None
    Dictionary containing computation cost options for simulated and observed responses. The elements are:

"""
    + _gen_docstring_from_base_doc(
        COST_OPTIONS_BASE_DOC,
        DEFAULT_SIMULATION_COST_OPTIONS["optimize"].keys(),
        nindent=1,
    )
    + """
common_options : `dict[str, Any]` or None, default None
    Dictionary containing common options with two elements:

"""
    + _gen_docstring_from_base_doc(
        COMMON_OPTIONS_BASE_DOC, DEFAULT_SIMULATION_COMMON_OPTIONS.keys(), nindent=1
    )
    + """

Returns
-------
multiple_optimize : `MultipleOptimize <smash.MultipleOptimize>`
    It returns an object containing the results of the multiple optimize.

See Also
--------
Samples : Represents the generated samples result.
MultipleOptimize : Represents the multiple optimize result.

Examples
--------
>>> from smash.factory import load_dataset
>>> from smash.factory import generate_samples
>>> setup, mesh = load_dataset("cance")
>>> model = smash.Model(setup, mesh)

Define sampling problem and generate samples

>>> problem = {
...            'num_vars': 4,
...            'names': ['cp', 'ct', 'kexc', 'llr'],
...            'bounds': [[1, 2000], [1, 1000], [-20, 5], [1, 1000]]
... }
>>> sr = generate_samples(problem, n=3, random_state=11)

Run multiple optimization processes

>>> mopt = smash.multiple_optimize(
...     model,
...     samples=sr,
...     optimize_options={"termination_crit": {"maxiter": 2}}
... )
</> Multiple Optimize
    Optimize 3/3 (100%(percent)s)

Get the cost values through multiple runs of optimization

>>> mopt.cost
array([0.5622911 , 0.0809496 , 0.16873538], dtype=float32)
"""
)

_optimize_control_info_doc = (
    """
Information on the optimization control vector of Model.

Parameters
----------
model : `Model <smash.Model>`
    Primary data structure of the hydrological model `smash`.

mapping : `str`, default 'uniform'
    Type of mapping. Should be one of

    - ``'uniform'``
    - ``'distributed'``
    - ``'multi-linear'``
    - ``'multi-polynomial'``

    .. hint::
        See a detailed explanation on the mapping in (TODO FC: link Math/Num) section.

optimizer : `str` or None, default None
    Name of optimizer. Should be one of

    - ``'sbs'`` (``'uniform'`` **mapping** only)
    - ``'lbfgsb'`` (``'uniform'``, ``'distributed'``, ``'multi-linear'`` or ``'multi-polynomial'`` **mapping** only)

    .. note::
        If not given, a default optimizer will be set depending on the optimization mapping:

        - **mapping** = ``'uniform'``; **optimizer** = ``'sbs'``
        - **mapping** = ``'distributed'``, ``'multi-linear'``, or ``'multi-polynomial'``; **optimizer** = ``'lbfgsb'``

    .. hint::
        See a detailed explanation on the optimizer in (TODO FC: link Math/Num) section.

optimize_options : `dict[str, Any]` or None, default None
    Dictionary containing optimization options for fine-tuning the optimization process. 
    See `%(default_optimize_options_func)s` to retrieve the default optimize options based on the **mapping** and **optimizer**.

"""
    + _gen_docstring_from_base_doc(
        OPTIMIZE_OPTIONS_BASE_DOC,
        [
            "parameters",
            "bounds",
            "control_tfm",
            "descriptor",
            "termination_crit",
        ],
        nindent=1,
    )
    + """ 
cost_options : `dict[str, Any]` or None, default None
    Dictionary containing computation cost options for simulated and observed responses. The elements are:

"""
    + _gen_docstring_from_base_doc(
        COST_OPTIONS_BASE_DOC,
        DEFAULT_SIMULATION_COST_OPTIONS["optimize"].keys(),
        nindent=1,
    )
    + """
Returns
-------
control_info : `dict[str, Any]`
    A dictionary containing optimize control information of Model. The elements are:

    - n : `int`
        The size of the control vector.

    - nbk : `numpy.ndarray`
        An array of shape *(4,)* containing the number of elements by kind (`Model.rr_parameters`, `Model.rr_initial_states`, `Model.serr_mu_parameters`, `Model.serr_sigma_parameters`) of the control vector (``sum(nbk) = n``).

    - x : `numpy.ndarray`
        An array of shape *(n,)* containing the initial values of the control vector (it can be transformed).

    - l : `numpy.ndarray`
        An array of shape *(n,)* containing the lower bounds of the control vector (it can be transformed).

    - u : `numpy.ndarray`
        An array of shape *(n,)* containing the upper bounds of the control vector (it can be transformed).

    - nbd : `numpy.ndarray`
        An array of shape *(n,)* containing the type of bounds of the control vector. The values are:

        - ``0``: unbounded
        - ``1``: only lower bound
        - ``2``: both lower and upper bounds
        - ``3``: only upper bound

    - name : `numpy.ndarray`
        An array of shape *(n,)* containing the names of the control vector. The naming convention is:

        - ``<key>0``: Spatially uniform parameter or multi-linear/polynomial intercept where ``<key>`` is the name of any rainfall-runoff parameters or initial_states (``'cp0'``, ``'llr0'``, ``'ht0'``, etc).
        - ``<key><row>-<col>``: Spatially distributed parameter where ``<key>`` is the name of any rainfall-runoff parameters or initial_states and ``<row>``, ``<col>``, the corresponding position in the spatial domain (``'cp1-1'``, ``'llr20-2'``, ``'ht3-12'``, etc). It's one based indexing.
        - ``<key>-<desc>-<kind>``: Multi-linear/polynomial descriptor linked parameter where ``<key>`` is the name of any rainfall-runoff parameters or initial_states, ``<desc>`` the corresponding descriptor and ``<kind>``, the kind of parameter (coefficient or exposant) (``'cp-slope-a'``, ``'llr-slope-b'``, ``'ht-dd-a'``).

    - x_bkg : `numpy.ndarray`
        An array of shape *(n,)* containing the background values of the control vector.

    - l_bkg : `numpy.ndarray`
        An array of shape *(n,)* containing the background lower bounds of the control vector.

    - u_bkg : `numpy.ndarray`
        An array of shape *(n,)* containing the background upper bounds of the control vector.

Examples
--------
>>> from smash.factory import load_dataset
>>> setup, mesh = load_dataset("cance")
>>> model = smash.Model(setup, mesh)

Default optimize control vector information

>>> control_info = smash.optimize_control_info(model)
>>> control_info
{
    'l': array([-13.815511 , -13.815511 ,  -4.6052704, -13.815511 ], dtype=float32),
    'l_bkg': array([ 1.e-06,  1.e-06, -5.e+01,  1.e-06], dtype=float32),
    'n': 4,
    'name': array(['cp0', 'ct0', 'kexc0', 'llr0'], dtype='<U5'),
    'nbd': array([2, 2, 2, 2], dtype=int32),
    'nbk': array([4, 0, 0, 0], dtype=int32),
    'u': array([6.9077554, 6.9077554, 4.6052704, 6.9077554], dtype=float32), 
    'u_bkg': array([1000., 1000.,   50., 1000.], dtype=float32), 
    'x': array([5.2983174, 6.214608 , 0.       , 1.609438 ], dtype=float32), 
    'x_bkg': array([200., 500.,   0.,   5.], dtype=float32),
}

This gives a direct indication of what the optimizer takes as input, depending on the optimization configuration set up.
4 rainfall-runoff parameters are uniformly optimized (``'cp0'``, ``'ct0'``, ``'kexc0'`` and ``'llr0'``).
Each parameter has a lower and upper bound (``2`` in ``nbd``) and a transformation was applied to the control (``x`` relative to ``x_bkg``)

With a customize optimize configuration. Here, choosing a ``multi-linear`` mapping and optimizing only ``cp`` and ``kexc`` with different descriptors

>>> control_info = smash.optimize_control_info(
        model,
        mapping="multi-linear",
        optimize_options={
            "parameters": ["cp", "kexc"],
            "descriptor": {"kexc": ["dd"]},
        },
    )
>>> control_info
{
    'l': array([-99., -99., -99., -99., -99.], dtype=float32),
    'l_bkg': array([-99., -99., -99., -99., -99.], dtype=float32),
    'n': 5,
    'name': array(['cp0', 'cp-slope-a', 'cp-dd-a', 'kexc0', 'kexc-dd-a'], dtype='<U10'),
    'nbd': array([0, 0, 0, 0, 0], dtype=int32),
    'nbk': array([5, 0, 0, 0], dtype=int32),
    'u': array([-99., -99., -99., -99., -99.], dtype=float32),
    'u_bkg': array([-99., -99., -99., -99., -99.], dtype=float32),
    'x': array([-1.3862944,  0.       ,  0.       ,  0.       ,  0.       ], dtype=float32),
    'x_bkg': array([-1.3862944,  0.       ,  0.       ,  0.       ,  0.       ], dtype=float32),
}

5 parameters are optimized which are the intercepts (``'cp0'`` and  ``'kexc0'``) 
and the coefficients (``'cp-slope-a'``, ``'cp-dd-a'`` and ``'kexc-dd-a'``) of the regression between the descriptors 
(``slope`` and ``dd``) and the rainfall-runoff parameters (``cp`` and ``kexc``)
"""
)

_bayesian_optimize_control_info_doc = (
    """
Information on the bayesian optimization control vector of Model.

Parameters
----------
model : `Model <smash.Model>`
    Primary data structure of the hydrological model `smash`.

mapping : `str`, default 'uniform'
    Type of mapping. Should be one of

    - ``'uniform'``
    - ``'distributed'``
    - ``'multi-linear'``
    - ``'multi-polynomial'``

    .. hint::
        See a detailed explanation on the mapping in (TODO FC: link Math/Num) section.

optimizer : `str` or None, default None
    Name of optimizer. Should be one of

    - ``'sbs'`` (``'uniform'`` **mapping** only)
    - ``'lbfgsb'`` (``'uniform'``, ``'distributed'``, ``'multi-linear'`` or ``'multi-polynomial'`` **mapping** only)

    .. note::
        If not given, a default optimizer will be set depending on the optimization mapping:

        - **mapping** = ``'uniform'``; **optimizer** = ``'sbs'``
        - **mapping** = ``'distributed'``, ``'multi-linear'``, or ``'multi-polynomial'``; **optimizer** = ``'lbfgsb'``

    .. hint::
        See a detailed explanation on the optimizer in (TODO FC: link Math/Num) section.

optimize_options : `dict[str, Any]` or None, default None
    Dictionary containing optimization options for fine-tuning the optimization process. 
    See `%(default_optimize_options_func)s` to retrieve the default optimize options based on the **mapping** and **optimizer**.

"""
    + _gen_docstring_from_base_doc(
        OPTIMIZE_OPTIONS_BASE_DOC,
        [
            "parameters",
            "bounds",
            "control_tfm",
            "descriptor",
            "termination_crit",
        ],
        nindent=1,
    )
    + """ 
cost_options : `dict[str, Any]` or None, default None
    Dictionary containing computation cost options for simulated and observed responses. The elements are:

"""
    + _gen_docstring_from_base_doc(
        COST_OPTIONS_BASE_DOC,
        DEFAULT_SIMULATION_COST_OPTIONS["bayesian_optimize"].keys(),
        nindent=1,
    )
    + """
Returns
-------
control_info : `dict[str, Any]`
    A dictionary containing optimize control information of Model. The elements are:

    - n : `int`
        The size of the control vector.

    - nbk : `numpy.ndarray`
        An array of shape *(4,)* containing the number of elements by kind (`Model.rr_parameters`, `Model.rr_initial_states`, `Model.serr_mu_parameters`, `Model.serr_sigma_parameters`) of the control vector (``sum(nbk) = n``).

    - x : `numpy.ndarray`
        An array of shape *(n,)* containing the initial values of the control vector (it can be transformed).

    - l : `numpy.ndarray`
        An array of shape *(n,)* containing the lower bounds of the control vector (it can be transformed).

    - u : `numpy.ndarray`
        An array of shape *(n,)* containing the upper bounds of the control vector (it can be transformed).

    - nbd : `numpy.ndarray`
        An array of shape *(n,)* containing the type of bounds of the control vector. The values are:

        - ``0``: unbounded
        - ``1``: only lower bound
        - ``2``: both lower and upper bounds
        - ``3``: only upper bound

    - name : `numpy.ndarray`
        An array of shape *(n,)* containing the names of the control vector. The naming convention is:

        - ``<key>0``: Spatially uniform parameter or multi-linear/polynomial intercept where ``<key>`` is the name of any rainfall-runoff parameters or initial_states (``'cp0'``, ``'llr0'``, ``'ht0'``, etc).
        - ``<key><row>-<col>``: Spatially distributed parameter where ``<key>`` is the name of any rainfall-runoff parameters or initial_states and ``<row>``, ``<col>``, the corresponding position in the spatial domain (``'cp1-1'``, ``'llr20-2'``, ``'ht3-12'``, etc). It's one based indexing.
        - ``<key>-<desc>-<kind>``: Multi-linear/polynomial descriptor linked parameter where ``<key>`` is the name of any rainfall-runoff parameters or initial_states, ``<desc>`` the corresponding descriptor and ``<kind>``, the kind of parameter (coefficient or exposant) (``'cp-slope-a'``, ``'llr-slope-b'``, ``'ht-dd-a'``).
        - ``<key>-<code>``: Structural error parameter where ``<key>`` is the name of any structural error mu or sigma parameters and ``<code>``, the corresponding gauge (``'sg0-V3524010'``, ``'sg1-V3524010'``, etc)

    - x_bkg : `numpy.ndarray`
        An array of shape *(n,)* containing the background values of the control vector.

    - l_bkg : `numpy.ndarray`
        An array of shape *(n,)* containing the background lower bounds of the control vector.

    - u_bkg : `numpy.ndarray`
        An array of shape *(n,)* containing the background upper bounds of the control vector.

Examples
--------
>>> from smash.factory import load_dataset
>>> setup, mesh = load_dataset("cance")
>>> model = smash.Model(setup, mesh)

Default optimize control vector information

>>> control_info = smash.bayesian_optimize_control_info(model)
>>> control_info
{
    'l': array([-1.3815511e+01, -1.3815511e+01, -4.6052704e+00, -1.3815511e+01, 1.0000000e-06, 1.0000000e-06], dtype=float32),
    'l_bkg': array([ 1.e-06,  1.e-06, -5.e+01,  1.e-06,  1.e-06,  1.e-06], dtype=float32),
    'n': 6,
    'name': array(['cp0', 'ct0', 'kexc0', 'llr0', 'sg0-V3524010', 'sg1-V3524010'], dtype='<U12'),
    'nbd': array([2, 2, 2, 2, 2, 2], dtype=int32),
    'nbk': array([4, 0, 0, 2], dtype=int32),
    'u': array([   6.9077554,    6.9077554,    4.6052704,    6.9077554, 1000.       ,   10.       ], dtype=float32),
    'u_bkg': array([1000., 1000.,   50., 1000., 1000.,   10.], dtype=float32),
    'x': array([5.2983174, 6.214608 , 0.       , 1.609438 , 1.       , 0.2      ], dtype=float32),
    'x_bkg': array([2.e+02, 5.e+02, 0.e+00, 5.e+00, 1.e+00, 2.e-01], dtype=float32),
}

This gives a direct indication of what the optimizer takes as input, depending on the optimization configuration set up.
4 rainfall-runoff parameters are uniformly optimized (``'cp0'``, ``'ct0'``, ``'kexc0'`` and ``'llr0'``) and 2 structural error sigma parameters at gauge ``'V3524010'`` (``'sg0-V3524010'``, ``'sg1-V3524010'``)
Each parameter has a lower and upper bound (``2`` in ``nbd``) and a transformation was applied to the control (``x`` relative to ``x_bkg``)

With a customize optimize configuration. Here, choosing a ``multi-linear`` mapping and
optimizing only 2 rainfall-runoff parameters ``cp``, ``kexc`` with different descriptors and 2 structural error sigma parameters ``sg0`` and ``sg1``.

>>> control_info = smash.bayesian_optimize_control_info(
        model,
        mapping="multi-linear",
        optimize_options={
            "parameters": ["cp", "kexc", "sg0", "sg1"],
            "descriptor": {"kexc": ["dd"]},
        },
    )
>>> control_info
{
    'l': array([-99., -99., -99., -99., -99.,   0.,   0.], dtype=float32),
    'l_bkg': array([-9.9e+01, -9.9e+01, -9.9e+01, -9.9e+01, -9.9e+01,  1.0e-06, 1.0e-06], dtype=float32),
    'n': 7,
    'name': array(['cp0', 'cp-slope-a', 'cp-dd-a', 'kexc0', 'kexc-dd-a', 'sg0-V3524010', 'sg1-V3524010'], dtype='<U12'),
    'nbd': array([0, 0, 0, 0, 0, 2, 2], dtype=int32),
    'nbk': array([5, 0, 0, 2], dtype=int32),
    'u': array([-99., -99., -99., -99., -99.,   1.,   1.], dtype=float32),
    'u_bkg': array([ -99.,  -99.,  -99.,  -99.,  -99., 1000.,   10.], dtype=float32),
    'x': array([-1.3862944e+00,  0.0000000e+00,  0.0000000e+00,  0.0000000e+00, 0.0000000e+00,  9.9999900e-04,  1.9999903e-02], dtype=float32),
    'x_bkg': array([-1.3862944,  0.       ,  0.       ,  0.       ,  0.       , 1.       ,  0.2      ], dtype=float32),
}

7 parameters are optimized which are the intercepts (``'cp0'`` and  ``'kexc0'``) 
and the coefficients (``'cp-slope-a'``, ``'cp-dd-a'`` and ``'kexc-dd-a'``) of the regression between the descriptors 
(``slope`` and ``dd``) and the rainfall-runoff parameters (``cp`` and ``kexc``) and the 2 structural error sigma parameters
(``'sg0'`` and ``'sg1'``) associated to the gauge ``'V3524010'``.

Retrieving information from the control vector is particularly useful for defining priors on the parameters.
During a bayesian optimization, it is possible to define these priors in the **cost_options** argument within
the ``'control_prior'`` key. The problem is that we don't know the control vector in advance until we've filled in all
the optimization options. This is why we can define all the optimization options in the 
`bayesian_optimize_control_info <smash.bayesian_optimize_control_info>` method, retrieve the names of the parameters 
that make up the control vector and then call the optimization function, assigning the priors we want to.

Assign Gaussian priors to the two rainfall-runoff parameters ``'cp0'`` and ``'kexc0'`` and perform a spatially uniform optimization

>>> model.bayesian_optimize(
        cost_options={
            "control_prior": {"cp0": ["Gaussian", [200, 100]], "kexc0": ["Gaussian", [0, 5]]}
        },
    )
</> Bayesian Optimize
    At iterate      0    nfg =     1    J =     29.988619    ddx = 0.64
    At iterate      1    nfg =    68    J =      2.749650    ddx = 0.64
    At iterate      2    nfg =   135    J =      2.580721    ddx = 0.32
    At iterate      3    nfg =   202    J =      2.473113    ddx = 0.16
    At iterate      4    nfg =   269    J =      2.432616    ddx = 0.08
    At iterate      5    nfg =   342    J =      2.430271    ddx = 0.02
    At iterate      6    nfg =   416    J =      2.429243    ddx = 0.01
    CONVERGENCE: DDX < 0.01
"""
)

_forward_run_doc_appender = DocAppender(_forward_run_doc, indents=0)
_smash_forward_run_doc_substitution = DocSubstitution(
    model_parameter="model : `Model <smash.Model>`\n\tPrimary data structure of the hydrological model `smash`.",
    model_return="model : `Model <smash.Model>`\n\t It returns an updated copy of the initial Model object.",
    model_example_func="model_fwd = smash.forward_run()",
    model_example_response="model_fwd",
    percent="%",
)
_model_forward_run_doc_substitution = DocSubstitution(
    model_parameter="",
    model_return="",
    model_example_func="model.forward_run()",
    model_example_response="model",
    percent="%",
)

_optimize_doc_appender = DocAppender(_optimize_doc, indents=0)
_smash_optimize_doc_substitution = DocSubstitution(
    model_parameter="model : `Model <smash.Model>`\n\tPrimary data structure of the hydrological model `smash`.",
    default_optimize_options_func="default_optimize_options <smash.default_optimize_options>",
    parameters_serr_mu_parameters="",
    parameters_serr_sigma_parameters="",
    parameters_note_serr_parameters="",
    bounds_get_serr_parameters_bounds="",
    model_return="model : `Model <smash.Model>`\n\t It returns an updated copy of the initial Model object.",
    model_example_func="model_opt = smash.optimize()",
    model_example_response="model_opt",
    percent="%",
)
_model_optimize_doc_substitution = DocSubstitution(
    model_parameter="",
    default_optimize_options_func="default_optimize_options <smash.default_optimize_options>",
    parameters_serr_mu_parameters="",
    parameters_serr_sigma_parameters="",
    parameters_note_serr_parameters="",
    bounds_get_serr_parameters_bounds="",
    model_return="",
    model_example_func="model.optimize()",
    model_example_response="model",
    percent="%",
)

_multiset_estimate_doc_appender = DocAppender(_multiset_estimate_doc, indents=0)
_smash_multiset_estimate_doc_substitution = DocSubstitution(
    model_parameter="model : `Model <smash.Model>`\n\tPrimary data structure of the hydrological model `smash`.",
    model_return="model : `Model <smash.Model>`\n\t It returns an updated copy of the initial Model object.",
    model_example_func="model_estim = smash.multiset_estimate(model, multiset=mfr)",
    model_example_response="model_estim",
    percent="%",
)
_model_multiset_estimate_doc_substitution = DocSubstitution(
    model_parameter="",
    model_return="",
    model_example_func="model.multiset_estimate(multiset=mfr)",
    model_example_response="model",
    percent="%",
)

_bayesian_optimize_doc_appender = DocAppender(_bayesian_optimize_doc, indents=0)
_smash_bayesian_optimize_doc_substitution = DocSubstitution(
    model_parameter="model : `Model <smash.Model>`\n\tPrimary data structure of the hydrological model `smash`.",
    default_optimize_options_func="default_bayesian_optimize_options <smash.default_bayesian_optimize_options>",
    parameters_serr_mu_parameters="- `Model.serr_mu_parameters`",
    parameters_serr_sigma_parameters="- `Model.serr_sigma_parameters`",
    parameters_note_serr_parameters=", `Model.serr_mu_parameters` and `Model.serr_sigma_parameters`",
    bounds_get_serr_parameters_bounds=", `Model.get_serr_mu_parameters_bounds` and `Model.get_serr_sigma_parameters_bounds`",
    model_return="model : `Model <smash.Model>`\n\t It returns an updated copy of the initial Model object.",
    model_example_func="model_bayes_opt = smash.bayesian_optimize()",
    model_example_response="model_bayes_opt",
    percent="%",
)
_model_bayesian_optimize_doc_substitution = DocSubstitution(
    model_parameter="",
    default_optimize_options_func="default_bayesian_optimize_options <smash.default_bayesian_optimize_options>",
    parameters_serr_mu_parameters="- `Model.serr_mu_parameters`",
    parameters_serr_sigma_parameters="- `Model.serr_sigma_parameters`",
    parameters_note_serr_parameters=", `Model.serr_mu_parameters` and `Model.serr_sigma_parameters`",
    bounds_get_serr_parameters_bounds=", `Model.get_serr_mu_parameters_bounds` and `Model.get_serr_sigma_parameters_bounds`",
    model_return="",
    model_example_func="model.bayesian_optimize()",
    model_example_response="model",
    percent="%",
)

_multiple_forward_run_doc_appender = DocAppender(_multiple_forward_run_doc, indents=0)

_multiple_optimize_doc_appender = DocAppender(_multiple_optimize_doc, indents=0)
_smash_multiple_optimize_doc_substitution = DocSubstitution(
    default_optimize_options_func="default_optimize_options <smash.default_optimize_options>",
    parameters_serr_mu_parameters="",
    parameters_serr_sigma_parameters="",
    parameters_note_serr_parameters="",
    bounds_get_serr_parameters_bounds="",
    percent="%",
)

_optimize_control_info_doc_appender = DocAppender(_optimize_control_info_doc, indents=0)
_smash_optimize_control_info_doc_substitution = DocSubstitution(
    default_optimize_options_func="default_optimize_options <smash.default_optimize_options>",
    parameters_serr_mu_parameters="",
    parameters_serr_sigma_parameters="",
    parameters_note_serr_parameters="",
    bounds_get_serr_parameters_bounds="",
)

_bayesian_optimize_control_info_doc_appender = DocAppender(
    _bayesian_optimize_control_info_doc, indents=0
)
_smash_bayesian_optimize_control_info_doc_substitution = DocSubstitution(
    default_optimize_options_func="default_bayesian_optimize_options <smash.default_bayesian_optimize_options>",
    parameters_serr_mu_parameters="- `Model.serr_mu_parameters`",
    parameters_serr_sigma_parameters="- `Model.serr_sigma_parameters`",
    parameters_note_serr_parameters=", `Model.serr_mu_parameters` and `Model.serr_sigma_parameters`",
    bounds_get_serr_parameters_bounds=", `Model.get_serr_mu_parameters_bounds` and `Model.get_serr_sigma_parameters_bounds`",
)