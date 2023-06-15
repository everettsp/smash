from __future__ import annotations

from smash.core._build_model_dev import (
    _map_dict_to_object_dev,
    _build_setup_dev,
    _build_mesh_dev,
    _build_input_data_dev,
    _build_parameters_dev,
)

from smash.solver._mwd_setup_dev import SetupDT_dev
from smash.solver._mwd_mesh_dev import MeshDT_dev
from smash.solver._mwd_input_data_dev import Input_DataDT_dev
from smash.solver._mwd_parameters_dev import ParametersDT_dev
from smash.solver._mwd_output_dev import OutputDT_dev

import numpy as np

__all__ = ["Model_dev"]


class Model_dev(object):
    def __init__(self, setup: dict | None, mesh: dict | None):
        if setup and mesh:
            if isinstance(setup, dict):
                nd = np.array(setup.get("descriptor_name", [])).size

                self.setup = SetupDT_dev(nd)

                _map_dict_to_object_dev(setup, self.setup)

                _build_setup_dev(self.setup)

            else:
                raise TypeError(f"setup argument must be dict")

            if isinstance(mesh, dict):
                self.mesh = MeshDT_dev(
                    self.setup, mesh["nrow"], mesh["ncol"], mesh["ng"]
                )

                _map_dict_to_object_dev(mesh, self.mesh)

                _build_mesh_dev(self.setup, self.mesh)

            else:
                raise TypeError(f"mesh argument must be dict")

            self._input_data = Input_DataDT_dev(self.setup, self.mesh)
            self.obs_response = self._input_data.obs_response
            self.physio_data = self._input_data.physio_data
            self.atmos_data = self._input_data.atmos_data

            _build_input_data_dev(self.setup, self.mesh, self._input_data)

            self._parameters = ParametersDT_dev(self.setup, self.mesh)
            self.opr_parameters = getattr(
                self._parameters.opr_parameters, self.setup.structure
            )
            self.opr_initial_states = getattr(
                self._parameters.opr_initial_states, self.setup.structure
            )

            _build_parameters_dev(
                self.setup, self.mesh, self._input_data, self._parameters
            )

            self._output = OutputDT_dev(self.setup, self.mesh)
            self.sim_response = self._output.sim_response
            self.opr_final_states = getattr(
                self._output.opr_final_states, self.setup.structure
            )

    def __copy__(self):
        copy = Model_dev(None, None)
        copy.setup = self.setup.copy()
        copy.mesh = self.mesh.copy()
        copy._input_data = self._input_data.copy()
        copy._parameters = self._parameters.copy()
        copy._output = self._output.copy()

        return copy

    @property
    def setup(self):
        return self._setup

    @setup.setter
    def setup(self, value):
        self._setup = value

    @property
    def mesh(self):
        return self._mesh

    @mesh.setter
    def mesh(self, value):
        self._mesh = value

    @property
    def obs_response(self):
        return self._obs_response

    @obs_response.setter
    def obs_response(self, value):
        self._obs_response = value

    @property
    def physio_data(self):
        return self._physio_data

    @physio_data.setter
    def physio_data(self, value):
        self._physio_data = value

    @property
    def atmos_data(self):
        return self._atmos_data

    @atmos_data.setter
    def atmos_data(self, value):
        self._atmos_data = value

    @property
    def opr_parameters(self):
        return self._opr_parameters

    @opr_parameters.setter
    def opr_parameters(self, value):
        self._opr_parameters = value

    @property
    def opr_initial_states(self):
        return self._opr_initial_states

    @opr_initial_states.setter
    def opr_initial_states(self, value):
        self._opr_initial_states = value

    @property
    def sim_response(self):
        return self._sim_response

    @sim_response.setter
    def sim_response(self, value):
        self._sim_response = value

    @property
    def opr_final_states(self):
        return self._opr_final_states

    @opr_final_states.setter
    def opr_final_states(self, value):
        self._opr_final_states = value

    def copy(self):
        return self.__copy__()