from __future__ import annotations

from smash._constant import (
    STRUCTURE_OPR_PARAMETERS,
    STRUCTURE_OPR_STATES,
    STRUCTURE_COMPUTE_CI,
    DEFAULT_OPR_PARAMETERS,
    DEFAULT_OPR_INITIAL_STATES,
    SERR_MU_MAPPING_PARAMETERS,
    SERR_SIGMA_MAPPING_PARAMETERS,
    DEFAULT_SERR_MU_PARAMETERS,
    DEFAULT_SERR_SIGMA_PARAMETERS,
)

from smash.core.model._read_input_data import (
    _read_qobs,
    _read_prcp,
    _read_pet,
    _read_descriptor,
)
from smash.core.model._standardize import _standardize_setup

from smash.fcore._mwd_sparse_matrix_manipulation import (
    compute_rowcol_to_ind_sparse as wrap_compute_rowcol_to_ind_sparse,
)
from smash.fcore._mw_atmos_statistic import (
    compute_mean_atmos as wrap_compute_mean_atmos,
)
from smash.fcore._mw_interception_capacity import (
    compute_interception_capacity as wrap_compute_interception_capacity,
)

import pandas as pd
import numpy as np

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from smash.fcore._mwd_setup import SetupDT
    from smash.fcore._mwd_mesh import MeshDT
    from smash.fcore._mwd_input_data import Input_DataDT
    from smash.fcore._mwd_parameters import ParametersDT
    from smash.fcore._mwd_output import OutputDT


# TODO: Move this function to a generic common function file
def _map_dict_to_object(dct: dict, obj: object, skip: list = []):
    for key, value in dct.items():
        if key in skip:
            continue
        if hasattr(obj, key):
            setattr(obj, key, value)
        # % Apply to the same object and not sub-object
        elif isinstance(value, dict):
            _map_dict_to_object(value, obj)


def _build_setup(setup: SetupDT):
    _standardize_setup(setup)

    st = pd.Timestamp(setup.start_time)

    et = pd.Timestamp(setup.end_time)

    setup.ntime_step = (et - st).total_seconds() / setup.dt

    setup.nop = len(STRUCTURE_OPR_PARAMETERS[setup.structure])
    setup.nos = len(STRUCTURE_OPR_STATES[setup.structure])
    setup.nsep_mu = len(SERR_MU_MAPPING_PARAMETERS[setup.serr_mu_mapping])
    setup.nsep_sigma = len(SERR_SIGMA_MAPPING_PARAMETERS[setup.serr_sigma_mapping])


def _build_mesh(setup: SetupDT, mesh: MeshDT):
    if setup.sparse_storage:
        wrap_compute_rowcol_to_ind_sparse(mesh)  # % Fortran subroutine
    mesh.local_active_cell = mesh.active_cell.copy()


def _build_input_data(setup: SetupDT, mesh: MeshDT, input_data: Input_DataDT):
    if setup.read_qobs:
        _read_qobs(setup, mesh, input_data)

    if setup.read_prcp:
        _read_prcp(setup, mesh, input_data)

    if setup.read_pet:
        _read_pet(setup, mesh, input_data)

    if setup.read_descriptor:
        _read_descriptor(setup, mesh, input_data)

    wrap_compute_mean_atmos(setup, mesh, input_data)  # % Fortran subroutine


def _build_parameters(
    setup: SetupDT,
    mesh: MeshDT,
    input_data: Input_DataDT,
    parameters: ParametersDT,
):
    # % Build parameters
    parameters.opr_parameters.keys = STRUCTURE_OPR_PARAMETERS[setup.structure]

    for i, key in enumerate(parameters.opr_parameters.keys):
        value = DEFAULT_OPR_PARAMETERS[key]
        if key == "llr":
            value *= setup.dt / 3_600
        parameters.opr_parameters.values[..., i] = value

    # % Build initial states
    parameters.opr_initial_states.keys = STRUCTURE_OPR_STATES[setup.structure]

    for i, key in enumerate(parameters.opr_initial_states.keys):
        value = DEFAULT_OPR_INITIAL_STATES[key]
        parameters.opr_initial_states.values[..., i] = value

    # TODO: Add a key in setup to trigger the auto computation of CI (True by default)
    if STRUCTURE_COMPUTE_CI[setup.structure] and setup.dt < 86_400:
        # % Date
        day_index = pd.date_range(
            start=setup.start_time, end=setup.end_time, freq=f"{int(setup.dt)}s"
        )[1:].to_series()

        # % Date to proleptic Gregorian ordinal
        day_index = day_index.apply(lambda x: x.toordinal()).to_numpy()

        # % Scale to 1 (Fortran indexing)
        day_index = day_index - day_index[0] + 1

        ind = np.argwhere(parameters.opr_parameters.keys == "ci").item()

        wrap_compute_interception_capacity(
            setup,
            mesh,
            input_data,
            day_index,
            day_index[-1],
            parameters.opr_parameters.values[..., ind],
        )  # % Fortran subroutine

    # % Build structural error mu parameters
    parameters.serr_mu_parameters.keys = SERR_MU_MAPPING_PARAMETERS[
        setup.serr_mu_mapping
    ]

    for i, key in enumerate(parameters.serr_mu_parameters.keys):
        value = DEFAULT_SERR_MU_PARAMETERS[key]
        parameters.serr_mu_parameters.values[..., i] = value

    # % Build structural error sigma parameters
    parameters.serr_sigma_parameters.keys = SERR_SIGMA_MAPPING_PARAMETERS[
        setup.serr_sigma_mapping
    ]

    for i, key in enumerate(parameters.serr_sigma_parameters.keys):
        value = DEFAULT_SERR_SIGMA_PARAMETERS[key]
        parameters.serr_sigma_parameters.values[..., i] = value


def _build_output(
    setup: SetupDT,
    output: OutputDT,
):
    # % Build final states
    output.opr_final_states.keys = STRUCTURE_OPR_STATES[setup.structure]
