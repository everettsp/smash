# smash main-package
from smash.core.model.model import Model

from smash.core.signal_analysis.metrics.metrics import metrics
from smash.core.signal_analysis.segmentation.segmentation import hydrograph_segmentation
from smash.core.signal_analysis.signatures.signatures import signatures
from smash.core.signal_analysis.prcp_indices.prcp_indices import precipitation_indices

from smash.core.simulation.run.run import forward_run, multiple_forward_run
from smash.core.simulation.optimize.optimize import (
    optimize,
    multiple_optimize,
    bayesian_optimize,
)
from smash.core.simulation.options import (
    default_optimize_options,
    default_bayesian_optimize_options,
)
from smash.core.simulation.control import (
    optimize_control_info,
    bayesian_optimize_control_info,
)
from smash.core.simulation.estimate.estimate import multiset_estimate

# smash sub-packages
from smash import io, factory

# smash objects
from smash.core.simulation.run.run import ForwardRun, MultipleForwardRun
from smash.core.simulation.optimize.optimize import (
    Optimize,
    MultipleOptimize,
    BayesianOptimize,
)
from smash.core.simulation.estimate.estimate import MultisetEstimate

from smash.core.signal_analysis.signatures.signatures import Signatures
from smash.core.signal_analysis.prcp_indices.prcp_indices import PrecipitationIndices

from smash.factory.samples.samples import Samples

from . import _version


def __getattr__(name):
    raise AttributeError(f"module 'smash' has no attribute '{name}'")


__all__ = [
    "Model",
    "metrics",
    "hydrograph_segmentation",
    "signatures",
    "precipitation_indices",
    "forward_run",
    "multiple_forward_run",
    "optimize",
    "multiple_optimize",
    "bayesian_optimize",
    "default_optimize_options",
    "default_bayesian_optimize_options",
    "optimize_control_info",
    "bayesian_optimize_control_info",
    "multiset_estimate",
    "io",
    "factory",
    "ForwardRun",
    "MultipleForwardRun",
    "Optimize",
    "MultipleOptimize",
    "BayesianOptimize",
    "MultisetEstimate",
    "Signatures",
    "PrecipitationIndices",
    "Samples",
]

__version__ = _version.get_versions()["version"]

__all__.extend(["__version__"])
