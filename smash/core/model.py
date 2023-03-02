from __future__ import annotations

from smash.solver._mwd_setup import SetupDT
from smash.solver._mwd_mesh import MeshDT
from smash.solver._mwd_input_data import Input_DataDT
from smash.solver._mwd_parameters import ParametersDT
from smash.solver._mwd_states import StatesDT
from smash.solver._mwd_output import OutputDT

from smash.solver._mw_forward import forward

from smash.solver._mw_derived_type_update import update_optimize_setup

from smash.core._constant import OPTIM_FUNC

from smash.core._build_model import (
    _parse_derived_type,
    _build_setup,
    _build_mesh,
    _build_input_data,
    _build_parameters,
)


from smash.core.optimize._optimize import (
    _optimize_sbs,
    _optimize_nelder_mead,
    _optimize_lbfgsb,
)

from smash.core.optimize._ann_optimize import _ann_optimize

from smash.core.optimize.bayes_optimize import _bayes_computation

from smash.core.optimize._standardize import (
    _standardize_optimize_args,
    _standardize_optimize_options,
    _standardize_bayes_estimate_args,
    _standardize_bayes_optimize_args,
    _standardize_ann_optimize_args,
)

from smash.core._event_segmentation import (
    _event_segmentation,
    _standardize_event_seg_options,
)

from smash.core.signatures import (
    _standardize_signatures,
    _signatures,
    _signatures_sensitivity,
)

from smash.core.prcp_indices import _prcp_indices

from smash.core.generate_samples import _get_bound_constraints, _standardize_problem

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import pandas as pd
    from smash.core.net import Net

import numpy as np


__all__ = ["Model"]


class Model(object):

    """
    Primary data structure of the hydrological model `smash`.

    Parameters
    ----------
    setup : dict
        Model initialization setup dictionary (see: :ref:`setup arguments <user_guide.model_initialization.setup>`).

    mesh : dict
        Model initialization mesh dictionary. (see: :ref:`mesh arguments <user_guide.model_initialization.mesh>`).

    See Also
    --------
    save_setup: Save Model initialization setup dictionary.
    read_setup: Read Model initialization setup dictionary.
    save_mesh: Save Model initialization mesh dictionary.
    read_mesh: Read Model initialization mesh dictionary.
    generate_mesh: Automatic mesh generation.
    save_model: Save Model object.
    read_model: Read Model object.

    Examples
    --------
    >>> setup, mesh = smash.load_dataset("cance")
    >>> model = smash.Model(setup, mesh)
    >>> model
    Structure: 'gr-a'
    Spatio-Temporal dimension: (x: 28, y: 28, time: 1440)
    Last update: Initialization
    """

    def __init__(self, setup: dict, mesh: dict):
        if setup or mesh:
            if isinstance(setup, dict):
                descriptor_name = setup.get("descriptor_name", [])

                nd = 1 if isinstance(descriptor_name, str) else len(descriptor_name)

                self.setup = SetupDT(nd, mesh["ng"])

                _parse_derived_type(self.setup, setup)

            else:
                raise TypeError(
                    f"'setup' argument must be dictionary, not {type(setup)}"
                )

            _build_setup(self.setup)

            if isinstance(mesh, dict):
                self.mesh = MeshDT(self.setup, mesh["nrow"], mesh["ncol"], mesh["ng"])

                _parse_derived_type(self.mesh, mesh)

            else:
                raise TypeError(f"'mesh' argument must be dictionary, not {type(mesh)}")

            _build_mesh(self.setup, self.mesh)

            self.input_data = Input_DataDT(self.setup, self.mesh)

            _build_input_data(self.setup, self.mesh, self.input_data)

            self.parameters = ParametersDT(self.mesh)

            _build_parameters(self.setup, self.mesh, self.input_data, self.parameters)

            self.states = StatesDT(self.mesh)

            self.output = OutputDT(self.setup, self.mesh)

            self._last_update = "Initialization"

    def __repr__(self):
        structure = f"Structure: '{self.setup.structure}'"
        dim = f"Spatio-Temporal dimension: (x: {self.mesh.ncol}, y: {self.mesh.nrow}, time: {self.setup._ntime_step})"
        last_update = f"Last update: {self._last_update}"

        return f"{structure}\n{dim}\n{last_update}"

    @property
    def setup(self):
        """
        The setup of the Model.

        The model setup is represented as a SetupDT object. See `SetupDT <smash.solver._mwd_setup.SetupDT>`.

        Examples
        --------
        >>> setup, mesh = smash.load_dataset("cance")
        >>> model = smash.Model(setup, mesh)

        If you are using IPython, tab completion allows you to visualize all the attributes and methods:

        >>> model.setup.<TAB>
        model.setup.copy(                   model.setup.prcp_directory
        model.setup.daily_interannual_pet   model.setup.prcp_format
        model.setup.descriptor_directory    model.setup.qobs_directory
        model.setup.descriptor_format       model.setup.read_descriptor
        model.setup.descriptor_name         model.setup.read_pet
        model.setup.dt                      model.setup.read_prcp
        model.setup.end_time                model.setup.read_qobs
        model.setup.from_handle(            model.setup.save_net_prcp_domain
        model.setup.mean_forcing            model.setup.save_qsim_domain
        model.setup.pet_conversion_factor   model.setup.sparse_storage
        model.setup.pet_directory           model.setup.start_time
        model.setup.pet_format              model.setup.structure
        model.setup.prcp_conversion_factor

        Notes
        -----
        This object is a wrapped derived type from `f90wrap <https://github.com/jameskermode/f90wrap>`__.
        """

        return self._setup

    @setup.setter
    def setup(self, value):
        if isinstance(value, SetupDT):
            self._setup = value

        else:
            raise TypeError(
                f"'setup' attribute must be set with {type(SetupDT())}, not {type(value)}"
            )

    @property
    def mesh(self):
        """
        The mesh of the Model.

        The model mesh is represented as a MeshDT object. See `MeshDT <smash.solver._mwd_mesh.MeshDT>`.

        Examples
        --------
        >>> setup, mesh = smash.load_dataset("cance")
        >>> model = smash.Model(setup, mesh)

        If you are using IPython, tab completion allows you to visualize all the attributes and methods:

        >>> model.mesh.<TAB>
        model.mesh.active_cell   model.mesh.gauge_pos
        model.mesh.area          model.mesh.nac
        model.mesh.code          model.mesh.ncol
        model.mesh.copy(         model.mesh.ng
        model.mesh.dx            model.mesh.nrow
        model.mesh.flwacc        model.mesh.path
        model.mesh.flwdir        model.mesh.xmin
        model.mesh.flwdst        model.mesh.ymax
        model.mesh.from_handle(


        Notes
        -----
        This object is a wrapped derived type from `f90wrap <https://github.com/jameskermode/f90wrap>`__.
        """

        return self._mesh

    @mesh.setter
    def mesh(self, value):
        if isinstance(value, MeshDT):
            self._mesh = value

        else:
            raise TypeError(
                f"'mesh' attribute must be set with {type(MeshDT())}, not {type(value)}"
            )

    @property
    def input_data(self):
        """
        The input data of the Model.

        The model input data is represented as a Input_DataDT object. See `Input_DataDT <smash.solver._mwd_input_data.Input_DataDT>`.

        Examples
        --------
        >>> setup, mesh = smash.load_dataset("cance")
        >>> model = smash.Model(setup, mesh)

        If you are using IPython, tab completion allows you to visualize all the attributes and methods:

        >>> model.input_data.<TAB>
        model.input_data.copy(         model.input_data.pet
        model.input_data.descriptor    model.input_data.prcp
        model.input_data.from_handle(  model.input_data.qobs
        model.input_data.mean_pet      model.input_data.sparse_pet
        model.input_data.mean_prcp     model.input_data.sparse_prcp

        Notes
        -----
        This object is a wrapped derived type from `f90wrap <https://github.com/jameskermode/f90wrap>`__.
        """
        return self._input_data

    @input_data.setter
    def input_data(self, value):
        if isinstance(value, Input_DataDT):
            self._input_data = value

        else:
            raise TypeError(
                f"'input_data' attribute must be set with {type(Input_DataDT())}, not {type(value)}"
            )

    @property
    def parameters(self):
        """
        The parameters of the Model.

        The model parameters is represented as a ParametersDT object. See `ParametersDT <smash.solver._mwd_parameters.ParametersDT>`.

        Examples
        --------
        >>> setup, mesh = smash.load_dataset("cance")
        >>> model = smash.Model(setup, mesh)

        If you are using IPython, tab completion allows you to visualize all the attributes and methods:

        >>> model.parameters.<TAB>
        model.parameters.alpha         model.parameters.cusl1
        model.parameters.b             model.parameters.cusl2
        model.parameters.beta          model.parameters.ds
        model.parameters.cft           model.parameters.dsm
        model.parameters.ci            model.parameters.exc
        model.parameters.clsl          model.parameters.from_handle(
        model.parameters.copy(         model.parameters.ks
        model.parameters.cp            model.parameters.lr
        model.parameters.cst           model.parameters.ws

        Notes
        -----
        This object is a wrapped derived type from `f90wrap <https://github.com/jameskermode/f90wrap>`__.
        """

        return self._parameters

    @parameters.setter
    def parameters(self, value):
        if isinstance(value, ParametersDT):
            self._parameters = value

        else:
            raise TypeError(
                f"'parameters' attribute must be set with {type(ParametersDT())}, not {type(value)}"
            )

    @property
    def states(self):
        """
        The states of the Model.

        The model states is represented as a StatesDT object. See `StatesDT <smash.solver._mwd_states.StatesDT>`.

        Examples
        --------
        >>> setup, mesh = smash.load_dataset("cance")
        >>> model = smash.Model(setup, mesh)

        If you are using IPython, tab completion allows you to visualize all the attributes and methods:

        >>> model.states.<TAB>
        model.states.copy(         model.states.hlsl
        model.states.from_handle(  model.states.hp
        model.states.hft           model.states.hst
        model.states.hi            model.states.husl1
        model.states.hlr           model.states.husl2

        Notes
        -----
        This object is a wrapped derived type from `f90wrap <https://github.com/jameskermode/f90wrap>`__.
        """

        return self._states

    @states.setter
    def states(self, value):
        if isinstance(value, StatesDT):
            self._states = value

        else:
            raise TypeError(
                f"'states' attribute must be set with {type(StatesDT())}, not {type(value)}"
            )

    @property
    def output(self):
        """
        The output of the Model.

        The model output is represented as a OutputDT object. See `OutputDT <smash.solver._mwd_output.OutputDT>`.

        Examples
        --------
        >>> setup, mesh = smash.load_dataset("cance")
        >>> model = smash.Model(setup, mesh)

        If you are using IPython, tab completion allows you to visualize all the attributes and methods:

        >>> model.output.<TAB>
        model.output.an                   model.output.parameters_gradient
        model.output.copy(                model.output.qsim
        model.output.cost                 model.output.qsim_domain
        model.output.from_handle(         model.output.sp1
        model.output.fstates              model.output.sp2
        model.output.ian                  model.output.sparse_qsim_domain

        Notes
        -----
        This object is a wrapped derived type from `f90wrap <https://github.com/jameskermode/f90wrap>`__.
        """

        return self._output

    @output.setter
    def output(self, value):
        if isinstance(value, OutputDT):
            self._output = value

        else:
            raise TypeError(
                f"'output' attribute must be set with {type(OutputDT())}, not {type(value)}"
            )

    @property
    def _last_update(self):
        return self.__last_update

    @_last_update.setter
    def _last_update(self, value):
        if isinstance(value, str):
            self.__last_update = value

        else:
            raise TypeError(
                f"'_last_update' attribute must be set with {str}, not {type(value)}"
            )

    def copy(self):
        """
        Make a deepcopy of the Model.

        Returns
        -------
        Model
            A copy of Model.

        Examples
        --------
        >>> setup, mesh = smash.load_dataset("cance")
        >>> model = smash.Model(setup, mesh)

        Create a pointer towards Model

        >>> model_ptr = model
        >>> model_ptr.parameters.cp = 1
        >>> model_ptr.parameters.cp[0,0], model.parameters.cp[0,0]
        (1.0, 1.0)

        Create a deepcopy of Model

        >>> model_dc = model.copy()
        >>> model_dc.parameters.cp = 200
        >>> model_dc.parameters.cp[0,0], model.parameters.cp[0,0]
        (200.0, 1.0)
        """

        copy = Model(None, None)
        copy.setup = self.setup.copy()
        copy.mesh = self.mesh.copy()
        copy.input_data = self.input_data.copy()
        copy.parameters = self.parameters.copy()
        copy.states = self.states.copy()
        copy.output = self.output.copy()
        copy._last_update = self._last_update

        return copy

    def run(self, inplace: bool = False):
        """
        Run the Model.

        Parameters
        ----------
        inplace : bool, default False
            if True, perform operation in-place.

        Returns
        -------
        Model : Model or None
            Model with run outputs or None if inplace.

        Notes
        -----
        This method is directly calling the forward model :math:`Y = M(k)`.

        Examples
        --------
        >>> setup, mesh = smash.load_dataset("cance")
        >>> model = smash.Model(setup, mesh)
        >>> model.run(inplace=True)
        >>> model
        Structure: 'gr-a'
        Spatio-Temporal dimension: (x: 28, y: 28, time: 1440)
        Last update: Forward Run

        Access to simulated discharge

        >>> model.output.qsim[0,:]
        array([1.9826449e-03, 1.3466686e-07, 6.7618025e-12, ..., 2.0916510e+01,
               2.0762346e+01, 2.0610489e+01], dtype=float32)
        """

        if inplace:
            instance = self

        else:
            instance = self.copy()

        print("</> Run Model")

        cost = np.float32(0)

        forward(
            instance.setup,
            instance.mesh,
            instance.input_data,
            instance.parameters,
            instance.parameters.copy(),
            instance.states,
            instance.states.copy(),
            instance.output,
            cost,
        )

        instance._last_update = "Forward Run"

        if not inplace:
            return instance

    def optimize(
        self,
        mapping: str = "uniform",
        algorithm: str | None = None,
        control_vector: str | list | tuple | set | None = None,
        bounds: dict | None = None,
        jobs_fun: str | list | tuple | set = "nse",
        wjobs_fun: list | tuple | set | None = None,
        event_seg: dict | None = None,
        gauge: str | list | tuple | set = "downstream",
        wgauge: str | list | tuple | set = "mean",
        ost: str | pd.Timestamp | None = None,
        options: dict | None = None,
        verbose: bool = True,
        inplace: bool = False,
    ):
        """
        Optimize the Model.

        .. hint::
            See the :ref:`user_guide` and :ref:`math_num_documentation` for more.

        Parameters
        ----------
        mapping : str, default 'uniform'
            Type of mapping. Should be one of

            - 'uniform'
            - 'distributed'
            - 'hyper-linear'
            - 'hyper-polynomial'

        algorithm : str or None, default None
            Type of algorithm. Should be one of

            - 'sbs'
            - 'nelder-mead'
            - 'l-bfgs-b'

            .. note::
                If not given, chosen to be one of ``sbs`` or ``l-bfgs-b`` depending on the optimization mapping.

        control_vector : str, sequence or None, default None
            Parameters and/or states to be optimized. The control vector argument
            can be any parameter or state name or any sequence of parameter and/or state names.

            .. note::
                If not given, the control vector will be composed of the parameters of the structure defined in the Model setup.

        bounds : dict or None, default None
            Bounds on control vector. The bounds argument is a dictionary where keys are the name of the
            parameters and/or states in the control vector (can be a subset of control vector sequence)
            and the values are pairs of ``(min, max)`` values (i.e. list, set or tuple) with ``min`` lower than ``max``.
            None value inside the dictionary will be filled in with default bound values.

            .. note::
                If not given, the bounds will be filled in with default bound values.

        jobs_fun : str or sequence, default 'nse'
            Type of objective function(s) to be minimized. Should be one or a sequence of any

            - ``Classical Objective Function``
                'nse', 'kge', 'kge2', 'se', 'rmse', 'logarithmic'
            - ``Continuous Signature``
                'Crc', 'Cfp2', 'Cfp10', 'Cfp50', 'Cfp90'
            - ``Flood Event Signature``
                'Epf', 'Elt', 'Erc'

            .. hint::
                See a detailed explanation on the cost function in :ref:`Math / Num Documentation <math_num_documentation.cost_functions>` section.

        wjobs_fun : sequence or None, default None
            Objective function(s) weights in case of multi-criteria optimization (i.e. a sequence of objective functions to minimize).

            .. note::
                If not given, the weights will correspond to the mean of the objective functions.


        event_seg : dict or None, default None
            A dictionary of event segmentation options when calculating flood event signatures for cost computation. The keys are

            - 'peak_quant'
            - 'max_duration'

            See `smash.Model.event_segmentation` for more.

            .. note::
                If not given in case flood signatures are computed, the default values will be set for these parameters.

        gauge : str, sequence, default 'downstream'
            Type of gauge to be optimized. There are two ways to specify it:

            1. A gauge code or any sequence of gauge codes.
               The gauge code(s) given must belong to the gauge codes defined in the Model mesh.
            2. An alias among ``all`` and ``downstream``. ``all`` is equivalent to a sequence of all gauge codes.
               ``downstream`` is equivalent to the gauge code of the most downstream gauge.

        wgauge : str, sequence, default 'mean'
            Type of gauge weights. There are two ways to specify it:

            1. A sequence of value whose size must be equal to the number of gauges optimized.
            2. An alias among ``mean``, ``median``, ``area`` or ``minv_area``.

        ost : str, pandas.Timestamp or None, default None
            The optimization start time. The optimization will only be performed between the
            optimization start time ``ost`` and the end time. The value can be a str which can be interpreted by
            pandas.Timestamp `(see here) <https://pandas.pydata.org/docs/reference/api/pandas.Timestamp.html>`__.
            The ``ost`` date value must be between the start time and the end time defined in the Model setup.

            .. note::
                If not given, the optimization start time will be equal to the start time.

        verbose : bool, default True
            Display information while optimizing.

        inplace : bool, default False
            if True, perform operation in-place.

        options : dict or None, default None
            A dictionary of algorithm options. The following options are available:
        
        maxiter : int, default is 40
            Maximum number of iterations for the optimization
        
        jreg_fun : str, sequence or None, default is None
            Type of regularization function(s) to be minimized. Should be one or a sequence of any

        - ``Regularization Function``
            'prior', 'smoothing'

        .. hint::
            See a detailed explanation on the cost function in :ref:`Math / Num Documentation <math_num_documentation.cost_functions>` section.

        wjreg_fun : sequence or None, default None
            Regularization function(s) weights in case of multi-regularization (i.e. a sequence of regularization functions to minimize).
            
            .. note::
            If not given, the weights is set to 1.
        
        wjreg: float, default is 0.
            Global regularization weith
        
        auto_regul: str {lcurve, fast} | None
            Methods to automatically compute the wjreg weith. the "fast" method consist of 2 optimizations cycle. The lcurve methods is more accurate and consist to nb_wjreg_lcurve optimization cycle.
        
        nb_wjreg_lcurve: int, default is 6
            Number of optimization cycle during the lcurve process. 6 is the minimum required.
            

        Returns
        -------
        Model : Model or None
            Model with optimize outputs or None if inplace.

        Notes
        -----
        This method is directly calling the forward model :math:`Y = M(k)` and the adjoint model
        :math:`\\delta k^* = \\left( \\frac{\\delta M}{\\delta k} \\right)^* . \\delta Y^*`
        if the algorithm ``l-bfgs-b`` is choosen to retrieve the gradient of the cost function wrt the control vector.

        Examples
        --------
        >>> setup, mesh = smash.load_dataset("cance")
        >>> model = smash.Model(setup, mesh)
        >>> model.optimize(inplace=True)
        >>> model
        Structure: 'gr-a'
        Spatio-Temporal dimension: (x: 28, y: 28, time: 1440)
        Last update: SBS Optimization

        Access to simulated discharge

        >>> model.output.qsim[0,:]
        array([5.7140866e-04, 4.7018618e-04, 3.5345653e-04, ..., 1.9017689e+01,
               1.8781073e+01, 1.8549627e+01], dtype=float32)

        Access to optimized parameters

        >>> ind = tuple(model.mesh.gauge_pos[0,:])
        >>> ind
        (20, 27)
        >>> (
        ... "cp", model.parameters.cp[ind],
        ... "cft", model.parameters.cft[ind],
        ... "exc", model.parameters.exc[ind],
        ... "lr", model.parameters.lr[ind],
        ... )
        ('cp', 76.57858, 'cft', 263.64627, 'exc', -1.455813, 'lr', 30.859276)
        """

        if inplace:
            instance = self

        else:
            instance = self.copy()

        print("</> Optimize Model")
        

        # % standardize args
        (
            mapping,
            algorithm,
            control_vector,
            jobs_fun,
            wjobs_fun,
            event_seg,
            bounds,
            wgauge,
            ost,
        ) = _standardize_optimize_args(
            mapping,
            algorithm,
            control_vector,
            jobs_fun,
            wjobs_fun,
            event_seg,
            bounds,
            gauge,
            wgauge,
            ost,
            instance.setup,
            instance.mesh,
            instance.input_data,
        )

        options = _standardize_optimize_options(options,instance.setup)
        
        if not "jreg_fun" in options:
            njr=0
        else:
            njr=options["jreg_fun"].size
        
        # % Update optimize setup derived type according to new optimize args and options !
        # % This Fortran subroutine reset optimize_setup values and realloc arrays.
        update_optimize_setup(
            instance.setup._optimize,
            instance.setup._ntime_step,
            instance.setup._nd,
            instance.mesh.ng,
            mapping,
            jobs_fun.size,
            njr, 
        )

        OPTIM_FUNC[algorithm](
        instance,
        control_vector,
        mapping,
        jobs_fun,
        wjobs_fun,
        event_seg,
        bounds,
        wgauge,
        ost,
        verbose,
        **options,
        )

        instance._last_update = f"{algorithm.upper()} Optimization"

        if not inplace:
            return instance

    def bayes_estimate(
        self,
        k: int | float | range | list | tuple | set | np.ndarray = 4,
        generator: str = "uniform",
        n: int = 1000,
        random_state: int | None = None,
        backg_sol: np.ndarray | None = None,
        coef_std: float | None = None,
        control_vector: str | list | tuple | set | None = None,
        bounds: list | tuple | set | None = None,
        jobs_fun: str | list | tuple | set = "nse",
        wjobs_fun: list | tuple | set | None = None,
        event_seg: dict | None = None,
        gauge: str | list | tuple | set = "downstream",
        wgauge: str | list | tuple | set = "mean",
        ost: str | pd.Timestamp | None = None,
        ncpu: int = 1,
        verbose: bool = True,
        inplace: bool = False,
        return_br: bool = False,
    ):
        """
        Estimate prior Model parameters/states using Bayesian approach.

        .. hint::
            See the :ref:`user_guide` and :ref:`math_num_documentation` for more.

        Parameters
        ----------
        k : int, float or sequence, default 4
            A regularisation parameter that controls the decay rate of the likelihood function.

            .. note::
                If k is a sequence, then the L-curve approach will be used to find an optimal value of k.

        generator, n, random_state, backg_sol, coef_std : multiple types
            Multiple arguments to generate spatially uniform Model parameters/states.
            See `smash.generate_samples` for more.

            .. hint::
                The generating samples problem can be redefined by using control_vector and bounds arguments.

        control_vector, bounds, jobs_fun, wjobs_fun, event_seg, gauge, wgauge, ost : multiple types
                Optimization setting to run the forward hydrological model and compute the cost values.
                See `smash.Model.optimize` for more.

        ncpu : integer, default 1
                If ncpu > 1, perform a parallel computation for all parameter sets.

        verbose : bool, default True
            Display information while estimating.

        inplace : bool, default False
            if True, perform operation in-place.

        return_br : bool, default False
            If True, also return the Bayesian estimation results ``BayesResult``.


        Returns
        -------
        Model : Model or None
            Model with optimize outputs if not inplace.

        res : BayesResult
            The Bayesian estimation results represented as a ``BayesResult`` object if return_br.

        See Also
        --------
        BayesResult: Represents the Bayesian estimation or optimization result.

        Examples
        --------
        >>> setup, mesh = smash.load_dataset("cance")
        >>> model = smash.Model(setup, mesh)
        >>> br = model.bayes_estimate(n=200, inplace=True, return_br=True, random_state=99)

        Access to cost values of the direct simulations

        >>> cost = br.data["cost"]
        >>> cost.sort()  # sort the values by ascending order
        >>> cost
        array([4.88620669e-01, 5.70850313e-01, 7.37333179e-01, 7.59980202e-01,
            ...
            9.26341797e+03, 1.12409111e+04, 1.13607480e+04, 1.18410371e+04])

        Compare to the cost value of the Model with the estimated parameters using Bayesian apporach

        >>> model.output.cost
        0.41908782720565796

        """

        if inplace:
            instance = self

        else:
            instance = self.copy()

        print("</> Bayes Estimate Model")

        # % standardize args
        (
            control_vector,
            jobs_fun,
            wjobs_fun,
            event_seg,
            bounds,
            wgauge,
            ost,
            k,
        ) = _standardize_bayes_estimate_args(
            control_vector,
            jobs_fun,
            wjobs_fun,
            event_seg,
            bounds,
            gauge,
            wgauge,
            ost,
            instance.setup,
            instance.mesh,
            instance.input_data,
            k,
        )

        res = _bayes_computation(
            instance,
            generator,
            n,
            random_state,
            backg_sol,
            coef_std,
            k,
            None,
            None,
            None,
            None,
            control_vector,
            None,
            jobs_fun,
            wjobs_fun,
            event_seg,
            bounds,
            wgauge,
            ost,
            verbose,
            None,
            ncpu,
        )

        instance._last_update = "Bayesian Estimation"

        if return_br:
            if not inplace:
                return instance, res

            else:
                return res

        else:
            if not inplace:
                return instance

    def bayes_optimize(
        self,
        k: int | float | range | list | tuple | set | np.ndarray = 4,
        density_estimate: bool = True,
        de_bw_method: str | None = None,
        de_weights: np.ndarray | None = None,
        generator: str = "uniform",
        n: int = 1000,
        random_state: int | None = None,
        backg_sol: np.ndarray | None = None,
        coef_std: float | None = None,
        mapping: str = "uniform",
        algorithm: str | None = None,
        control_vector: str | list | tuple | set | None = None,
        bounds: list | tuple | set | None = None,
        jobs_fun: str | list | tuple | set = "nse",
        wjobs_fun: list | tuple | set | None = None,
        event_seg: dict | None = None,
        gauge: str | list | tuple | set = "downstream",
        wgauge: str | list | tuple | set = "mean",
        ost: str | pd.Timestamp | None = None,
        options: dict | None = None,
        ncpu: int = 1,
        verbose: bool = True,
        inplace: bool = False,
        return_br: bool = False,
    ):
        """
        Optimize the Model using Bayesian approach.

        .. hint::
            See the :ref:`user_guide` and :ref:`math_num_documentation` for more.

        Parameters
        ----------
        k : int, float or sequence, default 4
            A regularisation parameter that controls the decay rate of the likelihood function.

            .. note::
                If k is a sequence, then the L-curve approach will be used to find an optimal value of k.

        density_estimate : bool, default True
            Take into account the density function estimated using Gaussian kernel.

        de_bw_method : str, scalar, callable or None, default None
            The method used to calculate the estimator bandwidth if density_estimate.
            This can be 'scott', 'silverman', a scalar constant or a callable.

            .. note::
                If not given and density_estimate=True, 'scott' is used as default.

            See `here <https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.gaussian_kde.html>`__ for more details.

        de_weights : array-like or None, default None
            A parameter related to weights of datapoints when estimating the density distribution.

            .. note::
                If not given and density_estimate=True, the samples are assumed to be equally weighted.

            See `here <https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.gaussian_kde.html>`__ for more details.

        generator, n, random_state, backg_sol, coef_std : multiple types
            Multiple arguments to generate spatially uniform Model parameters/states.
            See `smash.generate_samples` for more.

            .. hint::
                The generating samples problem can be redefined by using control_vector and bounds arguments.

        mapping, algorithm, control_vector, bounds, jobs_fun, wjobs_fun, event_seg, gauge, wgauge, ost, options : multiple types
                Optimization setting to optimize the Model using each generated spatially uniform parameters/states set as a first guess.
                See `smash.Model.optimize` for more.

        ncpu : integer, default 1
                If ncpu > 1, perform a parallel computation for all parameter sets.

        verbose : bool, default True
            Display information while optimizing.

        inplace : bool, default False
            if True, perform operation in-place.

        return_br : bool, default False
            If True, also return the Bayesian optimization results ``BayesResult``.

        Returns
        -------
        Model : Model or None
            Model with optimize outputs if not inplace.

        res : BayesResult
            The Bayesian optimization results represented as a ``BayesResult`` object if return_br.

        See Also
        --------
        BayesResult: Represents the Bayesian estimation or optimization result.

        Examples
        --------
        >>> setup, mesh = smash.load_dataset("cance")
        >>> model = smash.Model(setup, mesh)
        >>> br = model.bayes_optimize(k=1.75, n=100, inplace=True, options={"maxiter": 2}, return_br=True, random_state=99)

        Access to cost values of the optimizations with different set of Model parameters

        >>> cost = br.data["cost"]
        >>> cost.sort()  # sort the values by ascending order
        >>> cost
        array([0.04417887, 0.04713613, 0.05019313, 0.0512022 , 0.0563822 ,
            ...
            1.14625084, 1.15444124, 1.17005932, 1.19171917, 1.19486201])

        Compare to the cost value of the Model with the optimized parameters using Bayesian apporach

        >>> model.output.cost
        0.04369253292679787
        """

        if inplace:
            instance = self

        else:
            instance = self.copy()

        print("</> Bayes Optimize Model")

        # % standardize args
        (
            mapping,
            algorithm,
            control_vector,
            jobs_fun,
            wjobs_fun,
            event_seg,
            bounds,
            wgauge,
            ost,
            k,
        ) = _standardize_bayes_optimize_args(
            mapping,
            algorithm,
            control_vector,
            jobs_fun,
            wjobs_fun,
            event_seg,
            bounds,
            gauge,
            wgauge,
            ost,
            instance.setup,
            instance.mesh,
            instance.input_data,
            k,
        )

        options = _standardize_optimize_options(options,instance.setup)

        res = _bayes_computation(
            instance,
            generator,
            n,
            random_state,
            backg_sol,
            coef_std,
            k,
            density_estimate,
            de_bw_method,
            de_weights,
            algorithm,
            control_vector,
            mapping,
            jobs_fun,
            wjobs_fun,
            event_seg,
            bounds,
            wgauge,
            ost,
            verbose,
            options,
            ncpu,
        )

        instance._last_update = "Bayesian Optimization"

        if return_br:
            if not inplace:
                return instance, res

            else:
                return res

        else:
            if not inplace:
                return instance

    def ann_optimize(
        self,
        net: Net | None = None,
        control_vector: str | list | tuple | set | None = None,
        bounds: list | tuple | set | None = None,
        jobs_fun: str | list | tuple | set = "nse",
        wjobs_fun: list | tuple | set | None = None,
        event_seg: dict | None = None,
        gauge: str | list | tuple | set = "downstream",
        wgauge: str | list | tuple | set = "mean",
        ost: str | pd.Timestamp | None = None,
        epochs: int = 500,
        early_stopping: bool = False,
        verbose: bool = True,
        inplace: bool = False,
        return_net: bool = False,
    ):
        """
        Optimize the Model using Artificial Neural Network.

        .. hint::
            See the :ref:`user_guide` and :ref:`math_num_documentation` for more.

        Parameters
        ----------
        net : Net or None, default None
            The neural network Net will be trained to learn the descriptors to parameters mapping.

            .. note::
                If not given, a default network will be used. Otherwise, perform operation in-place on this Net.

        control_vector, bounds, jobs_fun, wjobs_fun, event_seg, gauge, wgauge, ost : multiple types
                Optimization setting to run the forward hydrological model and compute the cost values.
                See `smash.Model.optimize` for more.

        epochs : int, default 500
            The number of epochs to train the network.

        early_stopping : bool, default False
            Stop updating weights and biases when the loss function stops decreasing.

        verbose : bool, default True
            Display information while training.

        inplace : bool, default False
            if True, perform operation in-place.

        return_net : bool, default False
            If True and the default graph is used (net is None), also return the trained neural network.

        Returns
        -------
        Model : Model or None
            Model with optimize outputs if not inplace.

        Net : Net or None
            Net with trained weights and biases if return_net and the default graph is used.

        See Also
        --------
        Net : Artificial Neural Network initialization.

        Examples
        --------
        >>> setup, mesh = smash.load_dataset("cance")
        >>> model = smash.Model(setup, mesh)
        >>> net = model.ann_optimize(inplace=True, return_net=True)
        >>> model
        Structure: 'gr-a'
        Spatio-Temporal dimension: (x: 28, y: 28, time: 1440)
        Last update: ANN Optimization

        Display a summary of the neural network

        >>> net
        +-------------+
        | Net summary |
        +-------------+
        Input Shape: (2,)
        +----------------------+--------------+---------+
        | Layer (type)         | Output Shape | Param # |
        +----------------------+--------------+---------+
        | Dense                | (18,)        | 54      |
        | Activation (ReLU)    | (18,)        | 0       |
        | Dense                | (9,)         | 171     |
        | Activation (ReLU)    | (9,)         | 0       |
        | Dense                | (4,)         | 40      |
        | Activation (Sigmoid) | (4,)         | 0       |
        | Scale (MinMaxScale)  | (4,)         | 0       |
        +----------------------+--------------+---------+
        Total params: 265
        Trainable params: 265
        Non-trainable params: 0

        Access to some training information

        >>> net.history['loss_train']  # training loss
        [1.2267546653747559, ..., 0.03432881459593773]
        >>> net.layers[0].weight  # trained weights of the first layer
        array([[ 0.07801535,  0.10680847, -0.33354243, -0.17218271,  0.09706582,
                 0.63727553,  0.00399343,  0.00982828,  0.13701385, -0.00708624,
                -0.25468548, -0.29651733, -0.02438137, -0.10962573,  0.19415941,
                -0.2962292 ,  0.54098361,  0.70156156],
               [ 0.55410658,  0.12130747,  0.12366326,  0.53987803,  0.16916006,
                 0.30923786, -0.20702523,  0.48505195,  0.24157872,  0.38465208,
                -0.45171199, -0.29186438,  0.73921878, -0.00474557, -0.16782353,
                 0.27061135,  0.55383291,  0.71541684]])
        """

        if inplace:
            instance = self

        else:
            instance = self.copy()

        print("</> ANN Optimize Model")

        if net is None:
            use_default_graph = True

        else:
            use_default_graph = False

        # % standardize args
        (
            control_vector,
            jobs_fun,
            wjobs_fun,
            event_seg,
            bounds,
            wgauge,
            ost,
        ) = _standardize_ann_optimize_args(
            control_vector,
            jobs_fun,
            wjobs_fun,
            event_seg,
            bounds,
            gauge,
            wgauge,
            ost,
            instance.setup,
            instance.mesh,
            instance.input_data,
        )

        net = _ann_optimize(
            instance,
            control_vector,
            jobs_fun,
            wjobs_fun,
            event_seg,
            bounds,
            wgauge,
            ost,
            net,
            epochs,
            early_stopping,
            verbose,
        )

        instance._last_update = "ANN Optimization"

        if return_net and use_default_graph:
            if not inplace:
                return instance, net

            else:
                return net

        else:
            if not inplace:
                return instance

    def event_segmentation(self, peak_quant: float = 0.999, max_duration: float = 240):
        """
        Compute segmentation information of flood events over all catchments of the Model.

        .. hint::
            See the :ref:`User Guide <user_guide.event_segmentation>` and :ref:`Math / Num Documentation <math_num_documentation.hydrograph_segmentation>` for more.

        Parameters
        ----------
        peak_quant: float, default 0.999
            An event will be selected if its discharge exceed this quantile of the observed discharge timeseries.

        max_duration: float, default 240
            The maximum duration of an event (in hour).

        Returns
        -------
        res : pandas.DataFrame
            Flood events information obtained from segmentation algorithm.
            The dataframe has 6 columns which are

            - 'code' : the catchment code.
            - 'start' : the beginning of event.
            - 'end' : the end of event.
            - 'maxrainfall' : the moment that the maximum precipation is observed.
            - 'flood' : the moment that the maximum discharge is observed.
            - 'season' : the season that event occurrs.

        Examples
        --------
        >>> setup, mesh = smash.load_dataset("cance")
        >>> model = smash.Model(setup, mesh)

        Perform segmentation algorithm and display flood events infomation:

        >>> res = model.event_segmentation()
        >>> res
               code               start                   flood  season
        0  V3524010 2014-11-03 03:00:00 ... 2014-11-04 19:00:00  autumn
        1  V3515010 2014-11-03 10:00:00 ... 2014-11-04 20:00:00  autumn
        2  V3517010 2014-11-03 08:00:00 ... 2014-11-04 16:00:00  autumn

        [3 rows x 6 columns]

        """
        print("</> Model Event Segmentation")

        return _event_segmentation(self, peak_quant, max_duration)

    def signatures(
        self,
        sign: str | list | None = None,
        obs_comp: bool = True,
        event_seg: dict | None = None,
    ):
        """
        Compute continuous or/and flood event signatures of the Model.

        .. hint::
            See the :ref:`User Guide <user_guide.signatures.computation>` and :ref:`Math / Num Documentation <math_num_documentation.hydrological_signature>` for more.

        Parameters
        ----------
        sign : str, list of str or None, default None
            Define signature(s) to compute. Should be one of

            - 'Crc', 'Crchf', 'Crclf', 'Crch2r', 'Cfp2', 'Cfp10', 'Cfp50', 'Cfp90' (continuous signatures)
            - 'Eff', 'Ebf', 'Erc', 'Erchf', 'Erclf', 'Erch2r', 'Elt', 'Epf' (flood event signatures)

            .. note::
                If not given, all of continuous and flood event signatures will be computed.

        obs_comp : bool, default True
            If True, compute also the signatures from observed discharge in addition to simulated discharge.

        event_seg : dict or None, default None
            A dictionary of event segmentation options when calculating flood event signatures. The keys are

            - 'peak_quant'
            - 'max_duration'

            See `smash.Model.event_segmentation` for more.

            .. note::
                If not given in case flood signatures are computed, the default values will be set for these parameters.

        Returns
        -------
        res : SignResult
            The signatures computation results represented as a ``SignResult`` object.

        See Also
        --------
        SignResult: Represents signatures computation result.

        Examples
        --------
        >>> setup, mesh = smash.load_dataset("cance")
        >>> model = smash.Model(setup, mesh)
        >>> model.run(inplace=True)

        Compute all observed and simulated signatures:

        >>> res = model.signatures()
        >>> res.cont["obs"]  # observed continuous signatures
               code       Crc     Crchf  ...   Cfp50      Cfp90
        0  V3524010  0.516207  0.191349  ...  3.3225  42.631802
        1  V3515010  0.509180  0.147217  ...  1.5755  10.628400
        2  V3517010  0.514302  0.148364  ...  0.3235   2.776700

        [3 rows x 9 columns]

        >>> res.event["sim"]  # simulated flood event signatures
               code  season               start  ...  Elt         Epf
        0  V3524010  autumn 2014-11-03 03:00:00  ...    3  106.190651
        1  V3515010  autumn 2014-11-03 10:00:00  ...    0   21.160324
        2  V3517010  autumn 2014-11-03 08:00:00  ...    1   5.613392

        [3 rows x 12 columns]

        """

        print("</> Model Signatures")

        cs, es = _standardize_signatures(sign)

        event_seg = _standardize_event_seg_options(event_seg)

        return _signatures(self, cs, es, obs_comp, **event_seg)

    def signatures_sensitivity(
        self,
        problem: dict | None = None,
        n: int = 64,
        sign: str | list[str] | None = None,
        event_seg: dict | None = None,
        random_state: int | None = None,
    ):
        """
        Compute the first- and total-order variance-based sensitivity (Sobol indices) of spatially uniform hydrological model parameters on the output signatures.

        .. hint::
            See the :ref:`User Guide <user_guide.signatures.sensitivity>` and :ref:`Math / Num Documentation <math_num_documentation.hydrological_signature>` for more.

        Parameters
        ----------
        problem : dict or None, default None
            Problem definition to generate a multiple set of spatially uniform Model parameters. The keys are

            - 'num_vars' : the number of Model parameters.
            - 'names' : the name of Model parameters.
            - 'bounds' : the upper and lower bounds of each Model parameters (a sequence of (min, max)).

            .. note::
                If not given, the problem will be set automatically using `smash.Model.get_bound_constraints` method.

        n : int, default 64
            Number of trajectories to generate for each model parameter (ideally a power of 2).
            Then the number of sample to generate is equal to :math:`N(D+2)`
            where :math:`D` is the number of model parameters.

            See `here <https://salib.readthedocs.io/en/latest/api/SALib.sample.html#SALib.sample.sobol.sample>`__ for more details.

        sign : str, list or None, default None
            Define signature(s) to compute. Should be one of

            - 'Crc', 'Crchf', 'Crclf', 'Crch2r', 'Cfp2', 'Cfp10', 'Cfp50', 'Cfp90' (continuous signatures)
            - 'Eff', 'Ebf', 'Erc', 'Erchf', 'Erclf', 'Erch2r', 'Elt', 'Epf' (flood event signatures)

            .. note::
                If not given, all of continuous and flood event signatures will be computed.

        event_seg : dict or None, default None
            A dictionary of event segmentation options when calculating flood event signatures. The keys are

            - 'peak_quant'
            - 'max_duration'

            See `smash.Model.event_segmentation` for more.

            .. note::
                If not given in case flood signatures are computed, the default values will be set for these parameters.

        random_state : int or None, default None
            Random seed used to generate samples for sensitivity computation.

            .. note::
                If not given, generates parameters sets with a random seed.

        Returns
        -------
        res : SignSensResult
            The signatures sensitivity computation results represented as a ``SignSensResult`` object.

        See Also
        --------
        SignSensResult: Represents signatures sensitivity computation result.

        Examples
        --------
        >>> setup, mesh = smash.load_dataset("cance")
        >>> model = smash.Model(setup, mesh)
        >>> res = model.signatures_sensitivity(random_state=99)

        Total-order sensitivity indices of production parameter `cp` on continuous signatures:

        >>> res.cont["total_si"]["cp"]
               code       Crc     Crchf  ...     Cfp50     Cfp90
        0  V3524010  0.089630  0.023528  ...  1.092226  0.009049
        1  V3515010  0.064792  0.023358  ...  0.353862  0.011404
        2  V3517010  0.030655  0.028925  ...  0.113354  0.010977

        [3 rows x 9 columns]

        First-order sensitivity indices of linear routing parameter `lr` on flood event signatures:

        >>> res.event["first_si"]["lr"]
               code  season               start  ...       Elt       Epf
        0  V3524010  autumn 2014-11-03 03:00:00  ...  0.358599  0.015600
        1  V3515010  autumn 2014-11-03 10:00:00  ...  0.344644  0.010241
        2  V3517010  autumn 2014-11-03 08:00:00  ...  0.063007  0.010851

        [3 rows x 12 columns]

        """

        print("</> Model Signatures Sensitivity")

        instance = self.copy()

        cs, es = _standardize_signatures(sign)

        problem = _standardize_problem(problem, instance.setup, False)

        event_seg = _standardize_event_seg_options(event_seg)

        res = _signatures_sensitivity(
            instance, problem, n, cs, es, random_state, **event_seg
        )

        return res

    def prcp_indices(self):
        """
        Compute precipitations indices of the Model.

        4 precipitation indices are calculated for each gauge and each time step:

        - ``std`` : The precipitation spatial standard deviation,
        - ``d1`` : The first scaled moment, :cite:p:`zocatelli_2011`
        - ``d2`` : The second scaled moment, :cite:p:`zocatelli_2011`
        - ``vg`` : The vertical gap. :cite:p:`emmanuel_2015`

        .. hint::
            See the :ref:`User Guide <user_guide.prcp_indices>` for more.

        Returns
        -------
        res : PrcpIndicesResult
            The precipitation indices results represented as a ``PrcpIndicesResult`` object.

        See Also
        --------
        PrcpIndicesResult: Represents the precipitation indices result.

        Examples
        --------
        >>> setup, mesh = smash.load_dataset("cance")
        >>> model = smash.Model(setup, mesh)
        >>> res = model.prcp_indices()
        >>> res
        d1: array([[nan, nan, nan, ..., nan, nan, nan],
           [nan, nan, nan, ..., nan, nan, nan],
           [nan, nan, nan, ..., nan, nan, nan]], dtype=float32)
        ...
        vg: array([[nan, nan, nan, ..., nan, nan, nan],
           [nan, nan, nan, ..., nan, nan, nan],
           [nan, nan, nan, ..., nan, nan, nan]], dtype=float32)

        Each attribute is a numpy.ndarray of shape (number of gauge, number of time step)

        >>> res.d1.shape
        (3, 1440)

        NaN value means that there is no precipitation at this specific gauge and time step.
        Using numpy.where to find the index where precipitation indices were calculated on the most downstream gauge for the first scaled moment.

        >>> ind = np.argwhere(~np.isnan(res.d1[0,:])).squeeze()

        Viewing the first scaled moment on the first time step where rainfall occured on the most downstream gauge.

        >>> res.d1[0, ind[0]]
        1.209175
        """

        print("</> Model Precipitation Indices")

        return _prcp_indices(self)

    def get_bound_constraints(self, states: bool = False):
        """
        Get the boundary constraints of the Model parameters/states.

        Parameters
        ----------
        states : bool, default True
            If True, return boundary constraints of the Model states instead of Model parameters.

        Returns
        -------
        problem : dict
            The boundary constraint problem of the Model parameters/states. The keys are

            - 'num_vars': The number of Model parameters/states.
            - 'names': The name of Model parameters/states.
            - 'bounds': The upper and lower bounds of each Model parameters/states (a sequence of (min, max)).

        Examples
        --------
        >>> setup, mesh = smash.load_dataset("cance")
        >>> model = smash.Model(setup, mesh)
        >>> problem = model.get_bound_constraints()
        >>> problem
        {
            'num_vars': 4,
            'names': ['cp', 'cft', 'exc', 'lr'],
            'bounds': [[1e-06, 1000.0], [1e-06, 1000.0], [-50.0, 50.0], [1e-06, 1000.0]]
        }

        """

        return _get_bound_constraints(self.setup, states)
