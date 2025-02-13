# coding: utf-8

"""
Histogram hooks.
"""

from __future__ import annotations

from collections import defaultdict

import law
import order as od
import scinum as sn
import json
import pickle
import numpy as np
from collections import defaultdict

from columnflow.util import maybe_import, DotDict

np = maybe_import("numpy")
hist = maybe_import("hist")



logger = law.logger.get_logger(__name__)


def add_hist_hooks(config: od.Config) -> None:
    """
    Add histogram hooks to a configuration.
    """
    # helper to convert a histogram to a number object containing bin values and uncertainties
    # from variances stored in an array of values
    def hist_to_num(h: hist.Histogram, unc_name=str(sn.DEFAULT)) -> sn.Number:
        return sn.Number(h.values(), {unc_name: h.variances()**0.5})

    # helper to integrate values stored in an array based number object
    def integrate_num(num: sn.Number, axis=None) -> sn.Number:
        return sn.Number(
            nominal=num.nominal.sum(axis=axis),
            uncertainties={
                unc_name: (
                    (unc_values_up**2).sum(axis=axis)**0.5,
                    (unc_values_down**2).sum(axis=axis)**0.5,
                )
                for unc_name, (unc_values_up, unc_values_down) in num.uncertainties.items()
            },
        )

    def qcd_estimation(task, hists):
        if not hists:
            return hists

        # get the qcd process
        qcd_proc = config.get_process("qcd", default=None)
        if not qcd_proc:
            return hists

        # extract all unique category ids and verify that the axis order is exactly
        # "category -> shift -> variable" which is needed to insert values at the end
        CAT_AXIS, SHIFT_AXIS, VAR_AXIS = range(3)
        category_ids = set()
        for proc, h in hists.items():
            # validate axes
            assert len(h.axes) == 3
            assert h.axes[CAT_AXIS].name == "category"
            assert h.axes[SHIFT_AXIS].name == "shift"
            # get the category axis
            cat_ax = h.axes["category"]
            for cat_index in range(cat_ax.size):
                category_ids.add(cat_ax.value(cat_index))

        # create qcd groups
        qcd_groups: dict[str, dict[str, od.Category]] = defaultdict(DotDict)
        for cat_id in category_ids:
            cat_inst = config.get_category(cat_id)
            if cat_inst.has_tag({"os", "iso"}, mode=all):
                qcd_groups[cat_inst.x.qcd_group].os_iso = cat_inst
            elif cat_inst.has_tag({"os", "noniso"}, mode=all):
                qcd_groups[cat_inst.x.qcd_group].os_noniso = cat_inst
            elif cat_inst.has_tag({"ss", "iso"}, mode=all):
                qcd_groups[cat_inst.x.qcd_group].ss_iso = cat_inst
            elif cat_inst.has_tag({"ss", "noniso"}, mode=all):
                qcd_groups[cat_inst.x.qcd_group].ss_noniso = cat_inst

        # get complete qcd groups
        complete_groups = [name for name, cats in qcd_groups.items() if len(cats) == 4]

        # nothing to do if there are no complete groups
        if not complete_groups:
            return hists

        # sum up mc and data histograms, stop early when empty
        mc_hists = [h for p, h in hists.items() if p.is_mc and not p.has_tag("signal")]
        data_hists = [h for p, h in hists.items() if p.is_data]
        if not mc_hists or not data_hists:
            return hists
        mc_hist = sum(mc_hists[1:], mc_hists[0].copy())
        data_hist = sum(data_hists[1:], data_hists[0].copy())

        # start by copying the mc hist and reset it, then fill it at specific category slices
        hists[qcd_proc] = qcd_hist = mc_hist.copy().reset()
        for group_name in complete_groups:
            group = qcd_groups[group_name]

            # get the corresponding histograms and convert them to number objects,
            # each one storing an array of values with uncertainties
            # shapes: (SHIFT, VAR)
            get_hist = lambda h, region_name: h[{"category": hist.loc(group[region_name].id)}]
            os_noniso_mc = hist_to_num(get_hist(mc_hist, "os_noniso"), "os_noniso_mc")
            ss_noniso_mc = hist_to_num(get_hist(mc_hist, "ss_noniso"), "ss_noniso_mc")
            ss_iso_mc = hist_to_num(get_hist(mc_hist, "ss_iso"), "ss_iso_mc")
            os_noniso_data = hist_to_num(get_hist(data_hist, "os_noniso"), "os_noniso_data")
            ss_noniso_data = hist_to_num(get_hist(data_hist, "ss_noniso"), "ss_noniso_data")
            ss_iso_data = hist_to_num(get_hist(data_hist, "ss_iso"), "ss_iso_data")

            # estimate qcd shapes in the three sideband regions
            # shapes: (SHIFT, VAR)
            os_noniso_qcd = os_noniso_data - os_noniso_mc
            ss_iso_qcd = ss_iso_data - ss_iso_mc
            ss_noniso_qcd = ss_noniso_data - ss_noniso_mc

            # get integrals in ss regions for the transfer factor
            # shapes: (SHIFT,)
            int_ss_iso = integrate_num(ss_iso_qcd, axis=1)
            int_ss_noniso = integrate_num(ss_noniso_qcd, axis=1)

            # complain about negative integrals
            int_ss_iso_neg = int_ss_iso <= 0
            int_ss_noniso_neg = int_ss_noniso <= 0
            if int_ss_iso_neg.any():
                shift_ids = list(map(mc_hist.axes["shift"].value, np.where(int_ss_iso_neg)[0]))
                shifts = list(map(config.get_shift, shift_ids))
                logger.warning(
                    f"negative QCD integral in ss_iso region for group {group_name} and shifts: "
                    f"{', '.join(map(str, shifts))}",
                )
            if int_ss_noniso_neg.any():
                shift_ids = list(map(mc_hist.axes["shift"].value, np.where(int_ss_noniso_neg)[0]))
                shifts = list(map(config.get_shift, shift_ids))
                logger.warning(
                    f"negative QCD integral in ss_noniso region for group {group_name} and shifts: "
                    f"{', '.join(map(str, shifts))}",
                )

            # ABCD method
            # shape: (SHIFT, VAR)
            os_iso_qcd = os_noniso_qcd * ((int_ss_iso / int_ss_noniso)[:, None])

            # combine uncertainties and store values in bare arrays
            os_iso_qcd_values = os_iso_qcd()
            os_iso_qcd_variances = os_iso_qcd(sn.UP, sn.ALL, unc=True)**2

            # define uncertainties
            unc_data = os_iso_qcd(sn.UP, ["os_noniso_data", "ss_iso_data", "ss_noniso_data"], unc=True)
            unc_mc = os_iso_qcd(sn.UP, ["os_noniso_mc", "ss_iso_mc", "ss_noniso_mc"], unc=True)
            unc_data_rel = abs(unc_data / os_iso_qcd_values)
            unc_mc_rel = abs(unc_mc / os_iso_qcd_values)

            # only keep the MC uncertainty if it is larger than the data uncertainty and larger than 15%
            keep_variance_mask = (
                np.isfinite(unc_mc_rel) &
                (unc_mc_rel > unc_data_rel) &
                (unc_mc_rel > 0.15)
            )
            os_iso_qcd_variances[keep_variance_mask] = unc_mc[keep_variance_mask]**2
            os_iso_qcd_variances[~keep_variance_mask] = 0

            # retro-actively set values to zero for shifts that had negative integrals
            neg_int_mask = int_ss_iso_neg | int_ss_noniso_neg
            os_iso_qcd_values[neg_int_mask] = 1e-5
            os_iso_qcd_variances[neg_int_mask] = 0

            # residual zero filling
            zero_mask = os_iso_qcd_values <= 0
            os_iso_qcd_values[zero_mask] = 1e-5
            os_iso_qcd_variances[zero_mask] = 0

            # insert values into the qcd histogram
            cat_axis = qcd_hist.axes["category"]
            for cat_index in range(cat_axis.size):
                if cat_axis.value(cat_index) == group.os_iso.id:
                    qcd_hist.view().value[cat_index, ...] = os_iso_qcd_values
                    qcd_hist.view().variance[cat_index, ...] = os_iso_qcd_variances
                    break
            else:
                raise RuntimeError(
                    f"could not find index of bin on 'category' axis of qcd histogram {qcd_hist} "
                    f"for category {group.os_iso}",
                )

        return hists

    def qcd_inverted(task, hists):
        if not hists:
            return hists

        # get the qcd process
        qcd_proc = config.get_process("qcd", default=None)
        if not qcd_proc:
            return hists

        # extract all unique category ids and verify that the axis order is exactly
        # "category -> shift -> variable" which is needed to insert values at the end
        CAT_AXIS, SHIFT_AXIS, VAR_AXIS = range(3)
        category_ids = set()
        for proc, h in hists.items():
            # validate axes
            assert len(h.axes) == 3
            assert h.axes[CAT_AXIS].name == "category"
            assert h.axes[SHIFT_AXIS].name == "shift"
            # get the category axis
            cat_ax = h.axes["category"]
            for cat_index in range(cat_ax.size):
                category_ids.add(cat_ax.value(cat_index))

        # create qcd groups
        qcd_groups: dict[str, dict[str, od.Category]] = defaultdict(DotDict)
        for cat_id in category_ids:
            cat_inst = config.get_category(cat_id)
            if cat_inst.has_tag({"os", "iso"}, mode=all):
                qcd_groups[cat_inst.x.qcd_group].os_iso = cat_inst
            elif cat_inst.has_tag({"os", "noniso"}, mode=all):
                qcd_groups[cat_inst.x.qcd_group].os_noniso = cat_inst
            elif cat_inst.has_tag({"ss", "iso"}, mode=all):
                qcd_groups[cat_inst.x.qcd_group].ss_iso = cat_inst
            elif cat_inst.has_tag({"ss", "noniso"}, mode=all):
                qcd_groups[cat_inst.x.qcd_group].ss_noniso = cat_inst

        # get complete qcd groups
        complete_groups = [name for name, cats in qcd_groups.items() if len(cats) == 4]

        # nothing to do if there are no complete groups
        if not complete_groups:
            return hists

        # sum up mc and data histograms, stop early when empty
        mc_hists = [h for p, h in hists.items() if p.is_mc and not p.has_tag("signal")]
        data_hists = [h for p, h in hists.items() if p.is_data]
        if not mc_hists or not data_hists:
            return hists
        mc_hist = sum(mc_hists[1:], mc_hists[0].copy())
        data_hist = sum(data_hists[1:], data_hists[0].copy())

        # start by copying the mc hist and reset it, then fill it at specific category slices
        hists[qcd_proc] = qcd_hist = mc_hist.copy().reset()
        for group_name in complete_groups:
            group = qcd_groups[group_name]

            # get the corresponding histograms and convert them to number objects,
            # each one storing an array of values with uncertainties
            # shapes: (SHIFT, VAR)
            get_hist = lambda h, region_name: h[{"category": hist.loc(group[region_name].id)}]
            os_noniso_mc = hist_to_num(get_hist(mc_hist, "os_noniso"), "os_noniso_mc")
            ss_noniso_mc = hist_to_num(get_hist(mc_hist, "ss_noniso"), "ss_noniso_mc")
            ss_iso_mc = hist_to_num(get_hist(mc_hist, "ss_iso"), "ss_iso_mc")
            os_noniso_data = hist_to_num(get_hist(data_hist, "os_noniso"), "os_noniso_data")
            ss_noniso_data = hist_to_num(get_hist(data_hist, "ss_noniso"), "ss_noniso_data")
            ss_iso_data = hist_to_num(get_hist(data_hist, "ss_iso"), "ss_iso_data")

            # estimate qcd shapes in the three sideband regions
            # shapes: (SHIFT, VAR)
            os_noniso_qcd = os_noniso_data - os_noniso_mc
            ss_iso_qcd = ss_iso_data - ss_iso_mc
            ss_noniso_qcd = ss_noniso_data - ss_noniso_mc

            # get integrals in ss regions for the transfer factor
            # shapes: (SHIFT,)
            int_ss_noniso = integrate_num(ss_noniso_qcd, axis=1)
            int_os_noniso = integrate_num(os_noniso_qcd, axis=1)

            # complain about negative integrals
            int_ss_noniso_neg = int_ss_noniso <= 0
            int_os_noniso_neg = int_os_noniso <= 0
            if int_ss_noniso_neg.any():
                shift_ids = list(map(mc_hist.axes["shift"].value, np.where(int_ss_noniso_neg)[0]))
                shifts = list(map(config.get_shift, shift_ids))
                logger.warning(
                    f"negative QCD integral in ss_iso region for group {group_name} and shifts: "
                    f"{', '.join(map(str, shifts))}",
                )
            if int_os_noniso_neg.any():
                shift_ids = list(map(mc_hist.axes["shift"].value, np.where(int_os_noniso_neg)[0]))
                shifts = list(map(config.get_shift, shift_ids))
                logger.warning(
                    f"negative QCD integral in ss_noniso region for group {group_name} and shifts: "
                    f"{', '.join(map(str, shifts))}",
                )

            # ABCD method
            # shape: (SHIFT, VAR)
            os_iso_qcd = ss_iso_qcd * ((int_os_noniso / int_ss_noniso)[:, None])

            # combine uncertainties and store values in bare arrays
            os_iso_qcd_values = os_iso_qcd()
            os_iso_qcd_variances = os_iso_qcd(sn.UP, sn.ALL, unc=True)**2

            # define uncertainties
            unc_data = os_iso_qcd(sn.UP, ["ss_iso_data", "os_noniso_data", "ss_noniso_data"], unc=True)
            unc_mc = os_iso_qcd(sn.UP, ["ss_iso_mc", "os_noniso_mc", "ss_noniso_mc"], unc=True)
            unc_data_rel = abs(unc_data / os_iso_qcd_values)
            unc_mc_rel = abs(unc_mc / os_iso_qcd_values)

            # only keep the MC uncertainty if it is larger than the data uncertainty and larger than 15%
            keep_variance_mask = (
                np.isfinite(unc_mc_rel) &
                (unc_mc_rel > unc_data_rel) &
                (unc_mc_rel > 0.15)
            )
            os_iso_qcd_variances[keep_variance_mask] = unc_mc[keep_variance_mask]**2
            os_iso_qcd_variances[~keep_variance_mask] = 0

            # retro-actively set values to zero for shifts that had negative integrals
            neg_int_mask = int_os_noniso_neg | int_ss_noniso_neg
            os_iso_qcd_values[neg_int_mask] = 1e-5
            os_iso_qcd_variances[neg_int_mask] = 0

            # residual zero filling
            zero_mask = os_iso_qcd_values <= 0
            os_iso_qcd_values[zero_mask] = 1e-5
            os_iso_qcd_variances[zero_mask] = 0

            # insert values into the qcd histogram
            cat_axis = qcd_hist.axes["category"]
            for cat_index in range(cat_axis.size):
                if cat_axis.value(cat_index) == group.os_iso.id:
                    qcd_hist.view().value[cat_index, ...] = os_iso_qcd_values
                    qcd_hist.view().variance[cat_index, ...] = os_iso_qcd_variances
                    break
            else:
                raise RuntimeError(
                    f"could not find index of bin on 'category' axis of qcd histogram {qcd_hist} "
                    f"for category {group.os_iso}",
                )

        return hists

    def closure_test(task, hists):

        print("------------------------------------------")
        print("------ Entering closure test hook --------")
        print("------------------------------------------")

        if not hists:
            return hists

        # get the qcd process
        qcd_proc = config.get_process("qcd", default=None)
        if not qcd_proc:
            return hists

        # extract all unique category ids and verify that the axis order is exactly
        # "category -> shift -> variable" which is needed to insert values at the end
        CAT_AXIS, SHIFT_AXIS, VAR_AXIS = range(3)
        category_ids = set()
        for proc, h in hists.items():
            # validate axes
            assert len(h.axes) == 3
            assert h.axes[CAT_AXIS].name == "category"
            assert h.axes[SHIFT_AXIS].name == "shift"
            # get the category axis
            cat_ax = h.axes["category"]
            for cat_index in range(cat_ax.size):
                category_ids.add(cat_ax.value(cat_index))

        # create qcd groups
        qcd_groups: dict[str, dict[str, od.Category]] = defaultdict(DotDict)
        for cat_id in category_ids:
            cat_inst = config.get_category(cat_id)
            if cat_inst.has_tag({"os", "iso"}, mode=all):
                qcd_groups[cat_inst.x.qcd_group].os_iso = cat_inst
            elif cat_inst.has_tag({"os", "noniso"}, mode=all):
                qcd_groups[cat_inst.x.qcd_group].os_noniso = cat_inst
            elif cat_inst.has_tag({"ss", "iso"}, mode=all):
                qcd_groups[cat_inst.x.qcd_group].ss_iso = cat_inst
            elif cat_inst.has_tag({"ss", "noniso"}, mode=all):
                qcd_groups[cat_inst.x.qcd_group].ss_noniso = cat_inst

        # get complete qcd groups
        complete_groups = [name for name, cats in qcd_groups.items() if len(cats) == 4]

        # nothing to do if there are no complete groups
        if not complete_groups:
            return hists

        # sum up mc, use MC qcd as pseudo data and stop early when either is empty
        mc_hists = [h for p, h in hists.items() if p.is_mc and not p.has_tag("signal") and not p.has_tag("qcd")]
        data_qcd_hists = [h for p, h in hists.items() if p.has_tag("qcd")]
        if not mc_hists or not data_qcd_hists:
            return hists
        mc_hist = sum(mc_hists[1:], mc_hists[0].copy())
        data_qcd_hists = sum(data_qcd_hists[1:], mc_hists[0].copy())
        data_hist_tmp = data_qcd_hists + mc_hist

        data_hists = [h for p, h in hists.items() if p.is_data]
        data_hist = sum(data_hists[1:], data_hists[0].copy())
        data_qcd_hist = data_hist.copy().reset()
        data_qcd_hist.fill(data_hist_tmp)

        print("")
        print("mc_hist len: ", len(mc_hist))
        print("data_qcd_hists len: ", len(data_qcd_hists))
        print("data_qcd_hist len: ", len(data_qcd_hist))
        print("")

        # start by copying the mc hist and reset it, then fill it at specific category slices
        hists[qcd_proc] = qcd_hist = mc_hist.copy().reset()
        for group_name in complete_groups:
            group = qcd_groups[group_name]

            # get the corresponding histograms and convert them to number objects,
            # each one storing an array of values with uncertainties
            # shapes: (SHIFT, VAR)
            get_hist = lambda h, region_name: h[{"category": hist.loc(group[region_name].id)}]
            os_noniso_mc = hist_to_num(get_hist(mc_hist, "os_noniso"), "os_noniso_mc")
            ss_noniso_mc = hist_to_num(get_hist(mc_hist, "ss_noniso"), "ss_noniso_mc")
            ss_iso_mc = hist_to_num(get_hist(mc_hist, "ss_iso"), "ss_iso_mc")
            os_noniso_data = hist_to_num(get_hist(data_qcd_hist, "os_noniso"), "os_noniso_data")
            ss_noniso_data = hist_to_num(get_hist(data_qcd_hist, "ss_noniso"), "ss_noniso_data")
            ss_iso_data = hist_to_num(get_hist(data_qcd_hist, "ss_iso"), "ss_iso_data")

            # estimate qcd shapes in the three sideband regions
            # shapes: (SHIFT, VAR)
            os_noniso_qcd = os_noniso_data - os_noniso_mc
            ss_iso_qcd = ss_iso_data - ss_iso_mc
            ss_noniso_qcd = ss_noniso_data - ss_noniso_mc

            # get integrals in ss regions for the transfer factor
            # shapes: (SHIFT,)
            int_ss_iso = integrate_num(ss_iso_qcd, axis=1)
            int_ss_noniso = integrate_num(ss_noniso_qcd, axis=1)

            # complain about negative integrals
            int_ss_iso_neg = int_ss_iso <= 0
            int_ss_noniso_neg = int_ss_noniso <= 0
            if int_ss_iso_neg.any():
                shift_ids = list(map(mc_hist.axes["shift"].value, np.where(int_ss_iso_neg)[0]))
                shifts = list(map(config.get_shift, shift_ids))
                logger.warning(
                    f"negative QCD integral in ss_iso region for group {group_name} and shifts: "
                    f"{', '.join(map(str, shifts))}",
                )
            if int_ss_noniso_neg.any():
                shift_ids = list(map(mc_hist.axes["shift"].value, np.where(int_ss_noniso_neg)[0]))
                shifts = list(map(config.get_shift, shift_ids))
                logger.warning(
                    f"negative QCD integral in ss_noniso region for group {group_name} and shifts: "
                    f"{', '.join(map(str, shifts))}",
                )

            # ABCD method
            # shape: (SHIFT, VAR)
            os_iso_qcd = os_noniso_qcd * ((int_ss_iso / int_ss_noniso)[:, None])

            # combine uncertainties and store values in bare arrays
            os_iso_qcd_values = os_iso_qcd()
            os_iso_qcd_variances = os_iso_qcd(sn.UP, sn.ALL, unc=True)**2

            # define uncertainties
            unc_data = os_iso_qcd(sn.UP, ["os_noniso_data", "ss_iso_data", "ss_noniso_data"], unc=True)
            unc_mc = os_iso_qcd(sn.UP, ["os_noniso_mc", "ss_iso_mc", "ss_noniso_mc"], unc=True)
            unc_data_rel = abs(unc_data / os_iso_qcd_values)
            unc_mc_rel = abs(unc_mc / os_iso_qcd_values)

            # only keep the MC uncertainty if it is larger than the data uncertainty and larger than 15%
            keep_variance_mask = (
                np.isfinite(unc_mc_rel) &
                (unc_mc_rel > unc_data_rel) &
                (unc_mc_rel > 0.15)
            )
            os_iso_qcd_variances[keep_variance_mask] = unc_mc[keep_variance_mask]**2
            os_iso_qcd_variances[~keep_variance_mask] = 0

            # retro-actively set values to zero for shifts that had negative integrals
            neg_int_mask = int_ss_iso_neg | int_ss_noniso_neg
            os_iso_qcd_values[neg_int_mask] = 1e-5
            os_iso_qcd_variances[neg_int_mask] = 0

            # residual zero filling
            zero_mask = os_iso_qcd_values <= 0
            os_iso_qcd_values[zero_mask] = 1e-5
            os_iso_qcd_variances[zero_mask] = 0

            # insert values into the qcd histogram
            cat_axis = qcd_hist.axes["category"]
            for cat_index in range(cat_axis.size):
                if cat_axis.value(cat_index) == group.os_iso.id:
                    qcd_hist.view().value[cat_index, ...] = os_iso_qcd_values
                    qcd_hist.view().variance[cat_index, ...] = os_iso_qcd_variances
                    break
            else:
                raise RuntimeError(
                    f"could not find index of bin on 'category' axis of qcd histogram {qcd_hist} "
                    f"for category {group.os_iso}",
                )

        return hists

    # calculate the transfer factor for a chosen category and single decay channel (etau OR mutau OR tautau)
    def fake_factor(task, hists):
        if not hists:
            return hists

        # get dummy processes
        factor_bin = config.get_process("qcd", default=None)
        if not factor_bin:
            return hists

        factor_int = config.get_process("dy", default=None)
        if not factor_int:
            return hists

        # extract all unique category ids and verify that the axis order is exactly
        # "category -> shift -> variable" which is needed to insert values at the end
        CAT_AXIS, SHIFT_AXIS, VAR_AXIS = range(3)
        category_ids = set()
        for proc, h in hists.items():
            # validate axes
            assert len(h.axes) == 3
            assert h.axes[CAT_AXIS].name == "category"
            assert h.axes[SHIFT_AXIS].name == "shift"
            # get the category axis
            cat_ax = h.axes["category"]
            for cat_index in range(cat_ax.size):
                category_ids.add(cat_ax.value(cat_index))

        # create qcd groups
        qcd_groups: dict[str, dict[str, od.Category]] = defaultdict(DotDict)
        for cat_id in category_ids:
            cat_inst = config.get_category(cat_id)
            if cat_inst.has_tag({"os", "iso1"}, mode=all):
                qcd_groups[cat_inst.x.qcd_group].os_iso = cat_inst
            elif cat_inst.has_tag({"os", "noniso1"}, mode=all):
                qcd_groups[cat_inst.x.qcd_group].os_noniso = cat_inst
            elif cat_inst.has_tag({"ss", "iso1"}, mode=all):
                qcd_groups[cat_inst.x.qcd_group].ss_iso = cat_inst
            elif cat_inst.has_tag({"ss", "noniso1"}, mode=all):
                qcd_groups[cat_inst.x.qcd_group].ss_noniso = cat_inst

        # get complete qcd groups
        complete_groups = [name for name, cats in qcd_groups.items() if len(cats) == 4]

        # nothing to do if there are no complete groups
        if not complete_groups:
            return hists

        # sum up mc and data histograms, stop early when empty
        mc_hists = [h for p, h in hists.items() if p.is_mc and not p.has_tag("signal")]
        data_hists = [h for p, h in hists.items() if p.is_data]
        if not mc_hists or not data_hists:
            return hists
        mc_hist = sum(mc_hists[1:], mc_hists[0].copy())
        data_hist = sum(data_hists[1:], data_hists[0].copy())

        # start by copying the mc hist and reset it, then fill it at specific category slices
        hists = {}
        hists[factor_bin] = factor_hist = mc_hist.copy().reset()
        hists[factor_int] = factor_hist_int = mc_hist.copy().reset()
        for group_name in complete_groups:
            group = qcd_groups[group_name]
            # get the corresponding histograms and convert them to number objects,
            # each one storing an array of values with uncertainties
            # shapes: (SHIFT, VAR)
            get_hist = lambda h, region_name: h[{"category": hist.loc(group[region_name].id)}]
            ss_noniso_mc = hist_to_num(get_hist(mc_hist, "ss_noniso1"), "ss_noniso1_mc")
            ss_iso_mc = hist_to_num(get_hist(mc_hist, "ss_iso1"), "ss_iso1_mc")
            ss_noniso_data = hist_to_num(get_hist(data_hist, "ss_noniso1"), "ss_noniso1_data")
            ss_iso_data = hist_to_num(get_hist(data_hist, "ss_iso1"), "ss_iso1_data")

            # take the difference between data and MC in the control regions
            ss_iso_qcd = ss_iso_data - ss_iso_mc
            ss_noniso_qcd = ss_noniso_data - ss_noniso_mc

            # calculate the pt-independent fake factor
            int_ss_iso = integrate_num(ss_iso_qcd, axis=1)
            int_ss_noniso = integrate_num(ss_noniso_qcd, axis=1)
            fake_factor_int = (int_ss_iso / int_ss_noniso)[0, None]

            # calculate the pt-dependent fake factor
            fake_factor = (ss_iso_qcd / ss_noniso_qcd)[:, None]
            fake_factor_values = np.squeeze(np.nan_to_num(fake_factor()), axis=0)
            fake_factor_variances = fake_factor(sn.UP, sn.ALL, unc=True)**2

            # change shape of fake_factor_int for plotting
            fake_factor_int_values = fake_factor_values.copy()
            fake_factor_int_values.fill(fake_factor_int()[0])

            # insert values into the qcd histogram
            cat_axis = factor_hist.axes["category"]
            for cat_index in range(cat_axis.size):
                if cat_axis.value(cat_index) == group.os_iso.id:
                    factor_hist.view().value[cat_index, ...] = fake_factor_values
                    factor_hist.view().variance[cat_index, ...] = fake_factor_variances
                    factor_hist_int.view().value[cat_index, ...] = fake_factor_int_values
                    break
            else:
                raise RuntimeError(
                    f"could not find index of bin on 'category' axis of qcd histogram {factor_hist} "
                    f"for category {group.os_iso}",
                )
        return hists

    # calculate the transfer factor for a chosen category for all decay channels (etau AND muta AND tautau)
    def fake_factor_incl(task, hists):
        if not hists:
            return hists

        # get dummy processes
        factor_bin = config.get_process("qcd", default=None)
        if not factor_bin:
            return hists

        factor_int = config.get_process("dy", default=None)
        if not factor_int:
            return hists

        # extract all unique category ids and verify that the axis order is exactly
        # "category -> shift -> variable" which is needed to insert values at the end
        CAT_AXIS, SHIFT_AXIS, VAR_AXIS = range(3)
        category_ids = set()
        for proc, h in hists.items():
            # validate axes
            assert len(h.axes) == 3
            assert h.axes[CAT_AXIS].name == "category"
            assert h.axes[SHIFT_AXIS].name == "shift"
            # get the category axis
            cat_ax = h.axes["category"]
            for cat_index in range(cat_ax.size):
                category_ids.add(cat_ax.value(cat_index))

        # create qcd groups
        qcd_groups: dict[str, dict[str, od.Category]] = defaultdict(DotDict)
        for cat_id in category_ids:
            cat_inst = config.get_category(cat_id)
            if cat_inst.has_tag({"os", "iso"}, mode=all):
                qcd_groups[cat_inst.x.qcd_group].os_iso = cat_inst
            elif cat_inst.has_tag({"os", "noniso"}, mode=all):
                qcd_groups[cat_inst.x.qcd_group].os_noniso = cat_inst
            elif cat_inst.has_tag({"ss", "iso"}, mode=all):
                qcd_groups[cat_inst.x.qcd_group].ss_iso = cat_inst
            elif cat_inst.has_tag({"ss", "noniso"}, mode=all):
                qcd_groups[cat_inst.x.qcd_group].ss_noniso = cat_inst

        # get complete qcd groups
        complete_groups = [name for name, cats in qcd_groups.items() if len(cats) == 4]

        # nothing to do if there are no complete groups
        if not complete_groups:
            return hists

        # sum up mc and data histograms, stop early when empty
        mc_hists = [h for p, h in hists.items() if p.is_mc and not p.has_tag("signal")]
        data_hists = [h for p, h in hists.items() if p.is_data]
        if not mc_hists or not data_hists:
            return hists
        mc_hist = sum(mc_hists[1:], mc_hists[0].copy())
        data_hist = sum(data_hists[1:], data_hists[0].copy())

        # start by copying the mc hist and reset it, then fill it at specific category slices
        hists = {}
        channels = {}
        hists[factor_bin] = factor_hist = mc_hist.copy().reset()
        hists[factor_int] = factor_hist_int = mc_hist.copy().reset()

        for group_name in complete_groups:
            group = qcd_groups[group_name]
            # get the corresponding histograms and convert them to number objects,
            # each one storing an array of values with uncertainties
            # shapes: (SHIFT, VAR)
            get_hist = lambda h, region_name: h[{"category": hist.loc(group[region_name].id)}]
            ss_noniso_mc = hist_to_num(get_hist(mc_hist, "ss_noniso"), "ss_noniso_mc")
            ss_iso_mc = hist_to_num(get_hist(mc_hist, "ss_iso"), "ss_iso_mc")
            ss_noniso_data = hist_to_num(get_hist(data_hist, "ss_noniso"), "ss_noniso_data")
            ss_iso_data = hist_to_num(get_hist(data_hist, "ss_iso"), "ss_iso_data")

            channels[group_name] = {}
            channels[group_name]["ss_iso_mc"] = ss_iso_mc
            channels[group_name]["ss_noniso_mc"] = ss_noniso_mc
            channels[group_name]["ss_iso_data"] = ss_iso_data
            channels[group_name]["ss_noniso_data"] = ss_noniso_data

        for group_name in complete_groups:
            for k in ["incl"]:  # INDICATE WHICH CATEGORY TO CALCULATE THE FACTOR FOR ! e.g. "incl", "2j" ...
                if group_name == f"etau__{k}":  # DUMMY CATEGORY! mutau and tautau categories will have no factor calculated and etau is a proxy for the incl/2j/... category chosen
                    group = qcd_groups[group_name]

                    ss_iso_mc = channels[f"etau__{k}"]["ss_iso_mc"] + channels[f"mutau__{k}"]["ss_iso_mc"] + channels[f"tautau__{k}"]["ss_iso_mc"]
                    ss_noniso_mc = channels[f"etau__{k}"]["ss_noniso_mc"] + channels[f"mutau__{k}"]["ss_noniso_mc"] + channels[f"tautau__{k}"]["ss_noniso_mc"]
                    ss_iso_data = channels[f"etau__{k}"]["ss_iso_data"] + channels[f"mutau__{k}"]["ss_iso_data"] + channels[f"tautau__{k}"]["ss_iso_data"]
                    ss_noniso_data = channels[f"etau__{k}"]["ss_noniso_data"] + channels[f"mutau__{k}"]["ss_noniso_data"] + channels[f"tautau__{k}"]["ss_noniso_data"]

                    # take the difference between data and MC in the control regions
                    ss_iso_qcd = ss_iso_data - ss_iso_mc
                    ss_noniso_qcd = ss_noniso_data - ss_noniso_mc

                    # calculate the pt-independent fake factor
                    int_ss_iso = integrate_num(ss_iso_qcd, axis=1)
                    int_ss_noniso = integrate_num(ss_noniso_qcd, axis=1)
                    fake_factor_int = (int_ss_iso / int_ss_noniso)[0, None]

                    # calculate the pt-dependent fake factor
                    fake_factor = (ss_iso_qcd / ss_noniso_qcd)[:, None]
                    fake_factor_values = np.squeeze(np.nan_to_num(fake_factor()), axis=0)
                    fake_factor_variances = fake_factor(sn.UP, sn.ALL, unc=True)**2

                    # change shape of fake_factor_int for plotting
                    fake_factor_int_values = fake_factor_values.copy()
                    fake_factor_int_values.fill(fake_factor_int()[0])

                    # insert values into the qcd histogram
                    cat_axis = factor_hist.axes["category"]
                    for cat_index in range(cat_axis.size):
                        if cat_axis.value(cat_index) == group.os_iso.id:
                            factor_hist.view().value[cat_index, ...] = fake_factor_values
                            factor_hist.view().variance[cat_index, ...] = fake_factor_variances
                            factor_hist_int.view().value[cat_index, ...] = fake_factor_int_values
                            break
                    else:
                        raise RuntimeError(
                            f"could not find index of bin on 'category' axis of qcd histogram {factor_hist} "
                            f"for category {group.os_iso}",
                        )
        from IPython import embed
        embed()
        return hists

   

    def apply_fake_factor(task, hists):
        # Check if histograms are available
        if not hists:
            print("no hists")
            return hists

        # Get the qcd process
        qcd_proc = config.get_process("qcd", default=None)
        if not qcd_proc:
            print("no fake") 
            return hists

        # extract all unique category ids and verify that the axis order is exactly
        # "category -> shift -> variable" which is needed to insert values at the end
        CAT_AXIS, SHIFT_AXIS, VAR_AXIS = range(3)
        category_ids = set()
        for proc, h in hists.items():
            # validate axes
            assert len(h.axes) == 3
            assert h.axes[CAT_AXIS].name == "category"
            assert h.axes[SHIFT_AXIS].name == "shift"
            # get the category axis
            cat_ax = h.axes["category"]
            for cat_index in range(cat_ax.size):
                category_ids.add(cat_ax.value(cat_index))

        # create qcd groups
        qcd_groups: dict[str, dict[str, od.Category]] = defaultdict(DotDict)

        ### Create a QCD group ofr each DM and Njet category

        dms = ["tau1a1DM11", "tau1a1DM10", "tau1a1DM2", "tau1pi", "tau1rho"]  # Decay modes
        njets = ["has0j", "has1j", "has2j"]  # Jet multiplicity

        # dms = ["tau1rho"]  # Decay modes
        # njets = ["has0j"]  # Jet multiplicity

        # Loop over all categories and create a QCD group for each DM and Njet category
        for dm in dms:
            for njet in njets:
                for cat_id in category_ids:
                    cat_inst = config.get_category(cat_id)
                    if cat_inst.has_tag({"os", "iso1", njet, dm}, mode=all) and not cat_inst.has_tag("noniso2"): # cat D 
                        qcd_groups[f"dm_{dm}_njet_{njet}"].os_iso = cat_inst
                    elif cat_inst.has_tag({"os", "noniso1", njet, dm}, mode=all) and not cat_inst.has_tag("noniso2"): # cat C
                        qcd_groups[f"dm_{dm}_njet_{njet}"].os_noniso = cat_inst
                    #if cat_inst.has_tag({"os", "iso1", njet, dm}, mode=all) and not cat_inst.has_tag("iso2"): # cat D0 
                    #    qcd_groups[f"dm_{dm}_njet_{njet}"].os_iso = cat_inst
                    #elif cat_inst.has_tag({"os", "noniso1", njet, dm}, mode=all) and not cat_inst.has_tag("iso2"): # cat C0
                    #    qcd_groups[f"dm_{dm}_njet_{njet}"].os_noniso = cat_inst
   

        # Get complete qcd groups
        complete_groups = [name for name, cats in qcd_groups.items() if len(cats) == 2]
    
        # Nothing to do if there are no complete groups, you need C to apply Fake to D 
        if not complete_groups:
            print("no complete groups")
            return hists

        # Sum up mc and data histograms, stop early when empty
        mc_hists = [h for p, h in hists.items() if p.is_mc and not p.has_tag("signal")]
        data_hists = [h for p, h in hists.items() if p.is_data]
        if not mc_hists or not data_hists:
            return hists
        mc_hist = sum(mc_hists[1:], mc_hists[0].copy())
        data_hist = sum(data_hists[1:], data_hists[0].copy())
        
        # Start by copying the data hist and reset it, then fill it at specific category slices
        hists[qcd_proc] = qcd_hist = data_hist.copy().reset()
        for group_name in complete_groups:
            group = qcd_groups[group_name]
   
            # Get the corresponding histograms of the id, if not present, create a zeroed histogram
            get_hist = lambda h, region_name: (
                h[{"category": hist.loc(group[region_name].id)}]
                if group[region_name].id in h.axes["category"]
                else hist.Hist(*[axis for axis in (h[{"category": [0]}] * 0).axes if axis.name != 'category'])
            ) 

            # Get the corresponding histograms and convert them to number objects,
            os_noniso_mc  = hist_to_num(get_hist(mc_hist, "os_noniso"), "os_noniso_mc")
            os_noniso_data = hist_to_num(get_hist(data_hist, "os_noniso"), "os_noniso_data")

            ## DATA - MC of region C (FF are already apply to them)
            fake_hist = os_noniso_data - os_noniso_mc

            # combine uncertainties and store values in bare arrays
            fake_hist_values = fake_hist()
            fake_hist_variances = fake_hist(sn.UP, sn.ALL, unc=True)**2

            # Guaranty positive values of fake_hist
            neg_int_mask = fake_hist_values <= 0
            fake_hist_values[neg_int_mask] = 1e-5
            fake_hist_variances[neg_int_mask] = 0

            ## Use fake_hist as qcd histogram for category D (os_iso)
            cat_axis = qcd_hist.axes["category"]
            for cat_index in range(cat_axis.size):
                if cat_axis.value(cat_index) == group.os_iso.id:
                    qcd_hist.view().value[cat_index, ...] = fake_hist_values
                    qcd_hist.view().variance[cat_index, ...] = fake_hist_variances
                    break
            else:
                raise RuntimeError(
                    f"could not find index of bin on 'category' axis of qcd histogram {mc_hist} "
                    f"for category {group.os_iso}",
                )

        return hists




    # ABB
    def apply_fake_factor_AB(task, hists):
        # Check if histograms are available
        if not hists:
            print("no hists")
            return hists

        # Get the qcd process
        qcd_proc = config.get_process("qcd", default=None)
        if not qcd_proc:
            print("no fake") 
            return hists

        # extract all unique category ids and verify that the axis order is exactly
        # "category -> shift -> variable" which is needed to insert values at the end
        CAT_AXIS, SHIFT_AXIS, VAR_AXIS = range(3)
        category_ids = set()
        for proc, h in hists.items():
            # validate axes
            assert len(h.axes) == 3
            assert h.axes[CAT_AXIS].name == "category"
            assert h.axes[SHIFT_AXIS].name == "shift"
            # get the category axis
            cat_ax = h.axes["category"]
            for cat_index in range(cat_ax.size):
                category_ids.add(cat_ax.value(cat_index))

        # create qcd groups
        qcd_groups: dict[str, dict[str, od.Category]] = defaultdict(DotDict)

        ### Create a QCD group ofr each DM and Njet category

        dms = ["tau1a1DM11", "tau1a1DM10", "tau1a1DM2", "tau1pi", "tau1rho"]  # Decay modes
        njets = ["has0j", "has1j", "has2j"]  # Jet multiplicity

        # dms = ["tau1rho"]  # Decay modes
        # njets = ["has0j"]  # Jet multiplicity

        # Loop over all categories and create a QCD group for each DM and Njet category
        for dm in dms:
            for njet in njets:
                for cat_id in category_ids:
                    cat_inst = config.get_category(cat_id)
                    if cat_inst.has_tag({"ss", "iso1", njet, dm}, mode=all) and not cat_inst.has_tag("noniso2"): # cat A 
                        qcd_groups[f"dm_{dm}_njet_{njet}"].os_iso = cat_inst
                    elif cat_inst.has_tag({"ss", "noniso1", njet, dm}, mode=all) and not cat_inst.has_tag("noniso2"): # cat B
                        qcd_groups[f"dm_{dm}_njet_{njet}"].os_noniso = cat_inst
                    #if cat_inst.has_tag({"os", "iso1", njet, dm}, mode=all) and not cat_inst.has_tag("iso2"): # cat D0 
                    #    qcd_groups[f"dm_{dm}_njet_{njet}"].os_iso = cat_inst
                    #elif cat_inst.has_tag({"os", "noniso1", njet, dm}, mode=all) and not cat_inst.has_tag("iso2"): # cat C0
                    #    qcd_groups[f"dm_{dm}_njet_{njet}"].os_noniso = cat_inst
                    # Ignore AB categories
                    # elif cat_inst.has_tag({"ss", "iso1", njet, dm}, mode=all): # cat A
                    #     qcd_groups[f"dm_{dm}_njet_{njet}"].ss_iso = cat_inst
                    # elif cat_inst.has_tag({"ss", "noniso1", njet, dm}, mode=all):  # cat B
                    #     qcd_groups[f"dm_{dm}_njet_{njet}"].ss_noniso = cat_inst 
   

        # Get complete qcd groups
        complete_groups = [name for name, cats in qcd_groups.items() if len(cats) == 2]
    
        # Nothing to do if there are no complete groups, you need C to apply Fake to D 
        if not complete_groups:
            print("no complete groups")
            return hists

        # Sum up mc and data histograms, stop early when empty
        mc_hists = [h for p, h in hists.items() if p.is_mc and not p.has_tag("signal")]
        data_hists = [h for p, h in hists.items() if p.is_data]
        if not mc_hists or not data_hists:
            return hists
        mc_hist = sum(mc_hists[1:], mc_hists[0].copy())
        data_hist = sum(data_hists[1:], data_hists[0].copy())
        
        # Start by copying the data hist and reset it, then fill it at specific category slices
        hists[qcd_proc] = qcd_hist = data_hist.copy().reset()
        mc_hist_incl = mc_hist.copy().reset()
        data_hist_incl = data_hist.copy().reset()
        os_iso_mc_incl = None
        os_iso_data_incl = None
        for gidx, group_name in enumerate(complete_groups):
            group = qcd_groups[group_name]
   
            # Get the corresponding histograms of the id, if not present, create a zeroed histogram
            get_hist = lambda h, region_name: (
                h[{"category": hist.loc(group[region_name].id)}]
                if group[region_name].id in h.axes["category"]
                else hist.Hist(*[axis for axis in (h[{"category": [0]}] * 0).axes if axis.name != 'category'])
            ) 

            # Get the corresponding histograms and convert them to number objects,
            os_noniso_mc  = hist_to_num(get_hist(mc_hist, "os_noniso"), "os_noniso_mc")
            os_noniso_data = hist_to_num(get_hist(data_hist, "os_noniso"), "os_noniso_data")

            ## DATA - MC of region C (FF are already apply to them)
            fake_hist = os_noniso_data - os_noniso_mc

            # combine uncertainties and store values in bare arrays
            fake_hist_values = fake_hist()
            fake_hist_variances = fake_hist(sn.UP, sn.ALL, unc=True)**2

            # Guaranty positive values of fake_hist
            neg_int_mask = fake_hist_values <= 0
            fake_hist_values[neg_int_mask] = 1e-5
            fake_hist_variances[neg_int_mask] = 0

            ## Use fake_hist as qcd histogram for category D (os_iso)
            cat_axis = qcd_hist.axes["category"]
            for cat_index in range(cat_axis.size):
                if cat_axis.value(cat_index) == group.os_iso.id:
                    qcd_hist.view().value[cat_index, ...] = fake_hist_values
                    qcd_hist.view().variance[cat_index, ...] = fake_hist_variances
                    break
            else:
                raise RuntimeError(
                    f"could not find index of bin on 'category' axis of qcd histogram {mc_hist} "
                    f"for category {group.os_iso}",
                )

            hname = qcd_hist.axes[2].name
            path = "/eos/user/g/gsaha/CPinHToTauTauOutput/cf_store/analysis_httcp/cf.PlotVariables1D/QCD0"
            with open(f"{path}/qcd_{hname}_{group_name}.pkl", "wb") as f:
                pickle.dump(qcd_hist, f)

            # create a hist clone of the data_hist
            ratio_hist = data_hist.copy()

            # calultate sum_mc_hist
            #for cat_index in range(cat_axis.size):
            #    if cat_axis.value(cat_index) == group.os_iso.id:
            mc_hists = [h for p, h in hists.items() if p.is_mc and not p.has_tag("signal")]
            mc_hist_sum = sum(mc_hists[1:], mc_hists[0].copy())

            mc_hist_incl = mc_hist_incl + mc_hist_sum.copy()
            data_hist_incl = data_hist_incl + data_hist.copy()
            
            os_iso_mc  = hist_to_num(get_hist(mc_hist_sum, "os_iso"), "os_iso_mc")
            os_iso_data = hist_to_num(get_hist(data_hist, "os_iso"), "os_iso_data")
            #fake_iso_mc = hist_to_num(get_hist(qcd_hist, "os_iso"), "os_iso_mc")
            #os_iso_mc_sum = sum(os_iso_mc, fake_iso_mc)

            #from IPython import embed; embed()

            #ratio = os_iso_data/os_iso_mc_sum
            ratio = os_iso_data/os_iso_mc

            #fake_hist = os_noniso_data - os_noniso_mc

            # total MC
            os_iso_mc_incl   = (os_iso_mc + os_iso_mc_incl) if gidx > 0 else os_iso_mc
            os_iso_data_incl = (os_iso_data + os_iso_data_incl) if gidx > 0 else os_iso_data
            
            # combine uncertainties and store values in bare arrays
            ratio_hist_values = ratio()
            ratio_hist_variances = ratio(sn.UP, sn.ALL, unc=True)**2


            # Guaranty positive values of fake_hist
            neg_int_mask = ratio_hist_values <= 0
            ratio_hist_values[neg_int_mask] = 1e-5
            ratio_hist_variances[neg_int_mask] = 0

            ## Use fake_hist as qcd histogram for category D (os_iso)
            cat_axis = qcd_hist.axes["category"]
            for cat_index in range(cat_axis.size):
                if cat_axis.value(cat_index) == group.os_iso.id:
                    ratio_hist.view().value[cat_index, ...] = ratio_hist_values
                    ratio_hist.view().variance[cat_index, ...] = ratio_hist_variances
                    break
            else:
                raise RuntimeError(
                    f"could not find index of bin on 'category' axis of qcd histogram {mc_hist} "
                    f"for category {group.os_iso}",
                )
            
            path = f"{path}/Ratio"
            # save the ratio in a pickle file
            with open(f"{path}/ratio_{hname}_{group_name}_{group.os_iso.id}.pkl", "wb") as f:
                pickle.dump(ratio_hist, f)

        #from IPython import embed; embed()
        #mc_hist_incl_values = mc_hist_incl()
        #mc_hist_incl_variances = mc_hist_incl(sn.UP, sn.ALL, unc=True)**2
        #data_hist_incl_values = data_hist_incl()
        #data_hist_incl_variances = data_hist_incl(sn.UP, sn.ALL, unc=True)**2


        incl_ratio = os_iso_data_incl/os_iso_mc_incl

        incl_ratio_hist_values = incl_ratio()
        incl_ratio_hist_variances = incl_ratio(sn.UP, sn.ALL, unc=True)**2

        incl_ratio_hist = mc_hist_incl.copy().reset()
        incl_ratio_hist.view().value[0, ...] = incl_ratio_hist_values
        incl_ratio_hist.view().variance[0, ...] = incl_ratio_hist_variances

        
        #os_iso_mc_incl_hist_values = os_iso_mc_incl()
        #os_iso_mc_incl_hist_variances = os_iso_mc_incl(sn.UP, sn.ALL, unc=True)**2
        #os_iso_data_incl_hist_values = os_iso_data_incl()
        #os_iso_data_incl_hist_variances = os_iso_data_incl(sn.UP, sn.ALL, unc=True)**2

        #os_iso_mc_incl_hist = mc_hist_incl.copy().reset()
        #os_iso_data_incl_hist = mc_hist_incl.copy().reset()

        #os_iso_mc_incl_hist.view().value[0, ...] = os_iso_mc_incl_hist_values
        #os_iso_mc_incl_hist.view().variance[0, ...] = os_iso_mc_incl_hist_variances
        #os_iso_data_incl_hist.view().value[0, ...] = os_iso_data_incl_hist_values
        #os_iso_data_incl_hist.view().variance[0, ...] = os_iso_data_incl_hist_variances
        

        #with open(f"{path}/MC_{hname}_inclusive.pkl", "wb") as f:
        #    pickle.dump(os_iso_mc_incl_hist, f) #mc_hist_incl, f)
        with open(f"{path}/RATIO_{hname}_inclusive.pkl", "wb") as f:
            pickle.dump(incl_ratio_hist, f) #data_hist_incl, f)

        
                
        return hists
    
    

    config.x.hist_hooks = {
        "qcd": qcd_estimation,
        "qcd_inverted": qcd_inverted,
        "fake_factor": fake_factor,
        "fake_factor_incl": fake_factor_incl,
        "closure": closure_test,
        "apply_fake_factor": apply_fake_factor,
        "apply_fake_factor_AB": apply_fake_factor_AB
    }
