from __future__ import annotations

from smash._constant import PEAK_QUANT, MAX_DURATION

from smash.core.signal_analysis.segmentation._tools import _get_season, _events_grad

from smash.core.signal_analysis.segmentation._standardize import (
    _standardize_hydrograph_segmentation_args,
)

import numpy as np
import pandas as pd
import warnings

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from smash.core.model.model import Model
    from smash._typing import Numeric

__all__ = ["hydrograph_segmentation"]


def hydrograph_segmentation(
    model: Model,
    peak_quant: float = PEAK_QUANT,
    max_duration: Numeric = MAX_DURATION,
    by: str = "obs",
):
    """
    Compute segmentation information of flood events over all catchments of the Model.

    .. hint::
        See the (TODO: Fill) for more.

    Parameters
    ----------
    model : Model
        Model object.

    peak_quant : float, default 0.995
        Events will be selected if their discharge peaks exceed the **peak_quant**-quantile of the observed discharge timeseries.

    max_duration : int or float, default 240
        The expected maximum duration of an event (in hours). If multiple events are detected, their duration may exceed this value.

    by : str, default 'obs'
        Compute segmentation information based on observed (obs) or simulated (sim) discharges.
        A simulation (forward run or optimization) is required to obtain the simulated discharge when **by** is 'sim'.

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
    >>> import smash
    >>> from smash.factory import load_dataset
    >>> setup, mesh = load_dataset("cance")
    >>> model = smash.Model(setup, mesh)

    Perform segmentation algorithm and display flood events information:

    >>> hydro_seg = smash.hydrograph_segmentation(model)
    >>> hydro_seg
           code               start                   flood  season
    0  V3524010 2014-11-03 03:00:00 ... 2014-11-04 19:00:00  autumn
    1  V3515010 2014-11-03 10:00:00 ... 2014-11-04 20:00:00  autumn
    2  V3517010 2014-11-03 08:00:00 ... 2014-11-04 16:00:00  autumn

    [3 rows x 6 columns]
    """

    peak_quant, max_duration, by = _standardize_hydrograph_segmentation_args(
        peak_quant, max_duration, by
    )

    return _hydrograph_segmentation(model, peak_quant, max_duration, by)


def _hydrograph_segmentation(
    instance: Model, peak_quant: float, max_duration: Numeric, by: str
):
    date_range = pd.date_range(
        start=instance.setup.start_time,
        periods=instance.atmos_data.mean_prcp.shape[1],
        freq=f"{int(instance.setup.dt)}s",
    )

    col_name = ["code", "start", "end", "maxrainfall", "flood", "season"]

    df = pd.DataFrame(columns=col_name)

    for i, catchment in enumerate(instance.mesh.code):
        prcp = instance.atmos_data.mean_prcp[i, :].copy()

        suffix = "_data" if by == "obs" else ""
        q = getattr(instance, f"response{suffix}").q[i, :].copy()

        if (prcp < 0).all() or (q < 0).all():
            warnings.warn(
                f"Catchment {catchment} has no precipitation or/and discharge data"
            )

            pdrow = pd.DataFrame(
                [[catchment] + [np.nan] * (len(col_name) - 1)], columns=col_name
            )
            df = pdrow.copy() if df.empty else pd.concat([df, pdrow], ignore_index=True)

        else:
            list_events = _events_grad(
                prcp, q, peak_quant, max_duration, instance.setup.dt
            )

            for t in list_events:
                ts = date_range[t["start"]]
                te = date_range[t["end"]]
                peakq = date_range[t["peakQ"]]
                peakp = date_range[t["peakP"]]
                season = _get_season(ts)

                pdrow = pd.DataFrame(
                    [[catchment, ts, te, peakp, peakq, season]], columns=col_name
                )
                df = (
                    pdrow.copy()
                    if df.empty
                    else pd.concat([df, pdrow], ignore_index=True)
                )

    return df