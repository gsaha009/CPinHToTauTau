# coding: utf-8

"""
Histogram hooks.
"""

from __future__ import annotations

from collections import defaultdict

import os
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


def add_axis_with_value(hist_orig, value_for_new_axis):
    new_axis = hist.axis.IntCategory([value_for_new_axis], growth=True, name="process")

    # Step 1: Build new axes (insert new_axis after the first axis or at any position)
    old_axes = list(hist_orig.axes)
    
    # Insert new axis at desired position — here: second
    new_axes = old_axes[:1] + [new_axis] + old_axes[1:]

    # Step 2: Create new histogram with same storage
    hist_new = hist.Hist(*new_axes, storage=hist_orig._storage_type())

    # Fill using direct view (low-level assignment)
    # Add a new axis at the insert position using np.newaxis broadcasting
    hist_new.view().value[(slice(None),) * 1 + (0,)] = hist_orig.view().value
    hist_new.view().variance[(slice(None),) * 1 + (0,)] = hist_orig.view().variance
    """
    # Step 3: Build index mapping and fill
    for idx in hist_orig:
        # idx is a tuple of bin indices (e.g., (i_cat, i_shift, i_phi))
        value = hist_orig[idx].value
        if value == 0:
            continue

        # Get bin centers for all original axes
        centers = [axis.centers[i] for axis, i in zip(hist_orig.axes, idx)]

        # Insert the fixed value in the right position for the new axis
        # Match axis order: (cat, process, shift, phi, ...) for example
        fill_kwargs = {
            axis.name: center for axis, center in zip(old_axes, centers)
        }
        fill_kwargs[new_axis.name] = value_for_new_axis

        hist_new.fill(**fill_kwargs, weight=value)
    """
    return hist_new


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


    ######################################################
    ###              FAKE FACTOR HISTOGRAMS            ###
    ######################################################
    def fake_factors_onthefly(task, hists):

        # Check if histograms are available
        if not hists:
            print("no hists")
            return hists

        # Get the qcd process, it will be used to store the fake factor histograms
        qcd_proc = config.get_process("qcd", default=None)
        if not qcd_proc:
            print("no qcd process")
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

        # create qcd groups: A, B, C, D of A0, B0, C0, D0 for each DM and Njet category
        qcd_groups: dict[str, dict[str, od.Category]] = defaultdict(DotDict)


        #dms = ["tau2a1DM11", "tau2a1DM10", "tau2a1DM2", "tau2pi", "tau2rho"]  # Decay modes
        dms = ["tau2a1DM10", "tau2a1DM2", "tau2pi", "tau2rho"]  # Decay modes

        # Loop over all categories and create a QCD group for each DM and Njet category
        for dm in dms:
            for cat_id in category_ids:
                cat_inst = config.get_category(cat_id)
                if cat_inst.has_tag({"os", "noniso1", "iso2", "lowmt", dm}, mode=all): # DR_Num
                    qcd_groups[f"dm_{dm}"].dr_num = cat_inst
                elif cat_inst.has_tag({"ss", "noniso1", "iso2", "lowmt", dm}, mode=all): # DR_Den
                    qcd_groups[f"dm_{dm}"].dr_den = cat_inst
                elif cat_inst.has_tag({"ss", "iso1", "iso2", "lowmt", dm}, mode=all): # AR
                    qcd_groups[f"dm_{dm}"].ar = cat_inst
                elif cat_inst.has_tag({"os", "iso1", "iso2", "lowmt", dm}, mode=all): # SR
                    qcd_groups[f"dm_{dm}"].sr = cat_inst

                    
        # Get complete qcd groups
        complete_groups = [name for name, cats in qcd_groups.items() if len(cats) == 4]

        #from IPython import embed; embed()
        
        # Nothing to do if there are no complete groups, you need A and B to estimate FF
        if not complete_groups:
            print("no complete groups")
            return hists
        
        # Sum up mc and data histograms, stop early when empty, this is done for all categories
        mc_hists = [h for p, h in hists.items() if p.is_mc and not p.has_tag("signal")]
        data_hists = [h for p, h in hists.items() if p.is_data]
        if not mc_hists or not data_hists:
            return hists
        mc_hist = sum(mc_hists[1:], mc_hists[0].copy()) # sum all MC histograms, the hist object here contains all categories
        data_hist = sum(data_hists[1:], data_hists[0].copy()) # sum all data histograms, the hist object here contains all categories

        hists[qcd_proc] = qcd_hist = mc_hist.copy().reset()
        for gidx, group_name in enumerate(complete_groups):

            group = qcd_groups[group_name]
            # Get the corresponding histograms of the id, if not present, create a zeroed histogram
            get_hist = lambda h, region_name: (
                h[{"category": hist.loc(group[region_name].id)}]
                if group[region_name].id in h.axes["category"]
                else hist.Hist(*[axis for axis in (h[{"category": [0]}] * 0).axes if axis.name != 'category'])
            )

            # # Get the corresponding histograms and convert them to number objects,
            _dr_num_mc   = hist_to_num(get_hist(mc_hist,   "dr_num"), "dr_num_mc")    # MC in region DR_Num
            _dr_num_data = hist_to_num(get_hist(data_hist, "dr_num"), "dr_num_data")  # Data in region DR_num
            _dr_den_mc   = hist_to_num(get_hist(mc_hist,   "dr_den"), "dr_den_mc")    # MC in region DR_Den
            _dr_den_data = hist_to_num(get_hist(data_hist, "dr_den"), "dr_den_data")  # Data in region DR_Den
            _ar_mc       = hist_to_num(get_hist(mc_hist,   "ar"), "ar_mc")    # MC in region AR
            _ar_data     = hist_to_num(get_hist(data_hist, "ar"), "ar_data")  # Data in region AR
            
            # data will always have a single shift whereas mc might have multiple,
            # broadcast numbers in-place manually if necessary
            if (n_shifts := mc_hist.axes["shift"].size) > 1:
                def broadcast_data_num(num: sn.Number) -> None:
                    num._nominal = np.repeat(num.nominal, n_shifts, axis=0)
                    for name, (unc_up, unc_down) in num._uncertainties.items():
                        num._uncertainties[name] = (
                            np.repeat(unc_up, n_shifts, axis=0),
                            np.repeat(unc_down, n_shifts, axis=0),
                        )
                broadcast_data_num(os_noniso_data)
                broadcast_data_num(ss_noniso_data)
                broadcast_data_num(ss_iso_data)

            dr_num_qcd  = (_dr_num_data - _dr_num_mc) # Data - MC in region DR_Num
            dr_den_qcd  = (_dr_den_data - _dr_den_mc) # Data - MC in region DR_Den
            ar_qcd = (_ar_data - _ar_mc)


            # get integrals in ss regions for the transfer factor
            # shapes: (SHIFT,)
            int_dr_num = integrate_num(dr_num_qcd, axis=1)
            int_dr_den = integrate_num(dr_den_qcd, axis=1)

            # complain about negative integrals
            int_dr_num_neg = int_dr_num <= 0
            int_dr_den_neg = int_dr_den <= 0
            if int_dr_num_neg.any():
                shift_ids = list(map(mc_hist.axes["shift"].value, np.where(int_dr_num_neg)[0]))
                shifts = list(map(config.get_shift, shift_ids))
                logger.warning(
                    f"negative QCD integral in DR Num region for group {group_name} and shifts: "
                    f"{', '.join(map(str, shifts))}",
                )
            if int_dr_den_neg.any():
                shift_ids = list(map(mc_hist.axes["shift"].value, np.where(int_dr_den_neg)[0]))
                shifts = list(map(config.get_shift, shift_ids))
                logger.warning(
                    f"negative QCD integral in DR Den region for group {group_name} and shifts: "
                    f"{', '.join(map(str, shifts))}",
                )

            # ABCD method
            # shape: (SHIFT, VAR)
            sr_qcd = ar_qcd * ((int_dr_num / int_dr_den)[:, None])
                
            # combine uncertainties and store values in bare arrays
            sr_qcd_values = sr_qcd()
            sr_qcd_variances = sr_qcd(sn.UP, sn.ALL, unc=True)**2

            """
            # define uncertainties
            unc_data = sr_qcd(sn.UP, ["_dr_num_data", "_dr_den_data", "_ar_data"], unc=True)
            unc_mc = os_iso_qcd(sn.UP, ["_dr_num_mc", "_dr_den_mc", "_ar_mc"], unc=True)
            unc_data_rel = abs(unc_data / sr_qcd_values)
            unc_mc_rel = abs(unc_mc / sr_qcd_values)

            # only keep the MC uncertainty if it is larger than the data uncertainty and larger than 15%
            keep_variance_mask = (
                np.isfinite(unc_mc_rel) &
                (unc_mc_rel > unc_data_rel) &
                (unc_mc_rel > 0.15)
            )
            sr_qcd_variances[keep_variance_mask] = unc_mc[keep_variance_mask]**2
            sr_qcd_variances[~keep_variance_mask] = 0
            """
            # retro-actively set values to zero for shifts that had negative integrals
            neg_int_mask = int_dr_num_neg | int_dr_den_neg
            sr_qcd_values[neg_int_mask] = 1e-5
            sr_qcd_variances[neg_int_mask] = 0

            # residual zero filling
            zero_mask = sr_qcd_values <= 0
            sr_qcd_values[zero_mask] = 1e-5
            sr_qcd_variances[zero_mask] = 0

            # insert values into the qcd histogram
            cat_axis = qcd_hist.axes["category"]
            for cat_index in range(cat_axis.size):
                if cat_axis.value(cat_index) == group.sr.id:
                    qcd_hist.view().value[cat_index, ...] = sr_qcd_values
                    qcd_hist.view().variance[cat_index, ...] = sr_qcd_variances
                    break
            else:
                raise RuntimeError(
                    f"could not find index of bin on 'category' axis of qcd histogram {qcd_hist} "
                    f"for category {group.sr}",
                )

            
        return hists

    
    ######################################################
    ###         PRODUCE FAKE FACTOR HISTOGRAMS         ###
    ######################################################
    def produce_fake_factor(task, hists):
        '''
        Produce fake factor histograms for ABCD or A0B0C0D0 categories
        FF = (data-mc)A/(data-mc)B
        FF0 = (data-mc)A0/(data-mc)B0
        The resulting fake factors histograms are stroe in a pickle file.

        *** N.B. 
        * This will not save the distribution of the varaible in pdf format, but the FF pkl files will be saved. 
        * One can use one histogram only, e.g. hcand_1_pt
        '''
        # Define if we are calculating the fake factor for A0B0C0D0 or ABCD categories
        iszero_category = False # True for A0B0C0D0 categories

        # Check if histograms are available
        if not hists:
            print("no hists")
            return hists

        # Get the qcd process, it will be used to store the fake factor histograms
        qcd_proc = config.get_process("qcd", default=None)
        if not qcd_proc:
            print("no qcd process")
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

        # create qcd groups: A, B, C, D of A0, B0, C0, D0 for each DM and Njet category
        qcd_groups: dict[str, dict[str, od.Category]] = defaultdict(DotDict)


        dms = ["tau1a1DM11", "tau1a1DM10", "tau1a1DM2", "tau1pi", "tau1rho"]  # Decay modes
        njets = ["has0j", "has1j", "has2j"]  # Jet multiplicity

        # Loop over all categories and create a QCD group for each DM and Njet category
        for dm in dms:
            for njet in njets:
                for cat_id in category_ids:
                    cat_inst = config.get_category(cat_id)

                    if iszero_category: # A0B0C0D0 categories
                        if cat_inst.has_tag({"ss", "iso1", "noniso2", njet, dm}, mode=all): # cat A0
                            qcd_groups[f"dm_{dm}_njet_{njet}"].ss_iso = cat_inst
                        elif cat_inst.has_tag({"ss", "noniso1", "noniso2", njet, dm}, mode=all): # cat BB
                            qcd_groups[f"dm_{dm}_njet_{njet}"].ss_noniso = cat_inst
                    else: # ABCD categories
                        if cat_inst.has_tag({"ss", "iso1", njet, dm}, mode=all) and not cat_inst.has_tag("noniso2"): # cat A 
                            qcd_groups[f"dm_{dm}_njet_{njet}"].ss_iso = cat_inst
                        elif cat_inst.has_tag({"ss", "noniso1", njet, dm}, mode=all) and not cat_inst.has_tag("noniso2"): # cat B
                            qcd_groups[f"dm_{dm}_njet_{njet}"].ss_noniso = cat_inst

        # Get complete qcd groups
        complete_groups = [name for name, cats in qcd_groups.items() if len(cats) == 2]

        # Nothing to do if there are no complete groups, you need A and B to estimate FF
        if not complete_groups:
            print("no complete groups")
            return hists
        
        # Sum up mc and data histograms, stop early when empty, this is done for all categories
        mc_hists = [h for p, h in hists.items() if p.is_mc and not p.has_tag("signal")]
        data_hists = [h for p, h in hists.items() if p.is_data]
        if not mc_hists or not data_hists:
            return hists
        mc_hist = sum(mc_hists[1:], mc_hists[0].copy()) # sum all MC histograms, the hist object here contains all categories
        data_hist = sum(data_hists[1:], data_hists[0].copy()) # sum all data histograms, the hist object here contains all categories

        for gidx, group_name in enumerate(complete_groups):

            group = qcd_groups[group_name]
            # Get the corresponding histograms of the id, if not present, create a zeroed histogram
            get_hist = lambda h, region_name: (
                h[{"category": hist.loc(group[region_name].id)}]
                if group[region_name].id in h.axes["category"]
                else hist.Hist(*[axis for axis in (h[{"category": [0]}] * 0).axes if axis.name != 'category'])
            )

            # # Get the corresponding histograms and convert them to number objects,
            ss_iso_mc = hist_to_num(get_hist(mc_hist, "ss_iso"), "ss_iso_mc") # MC in region A
            ss_iso_data = hist_to_num(get_hist(data_hist, "ss_iso"), "ss_iso_data") # Data in region A
            ss_noniso_mc  = hist_to_num(get_hist(mc_hist, "ss_noniso"), "ss_noniso_mc") # MC in region B
            ss_noniso_data = hist_to_num(get_hist(data_hist, "ss_noniso"), "ss_noniso_data") # Data in region B


            ss_iso_data_minus_mc = (ss_iso_data - ss_iso_mc)[:, None] # Data - MC in region A
            ss_noniso_data_minus_mc = (ss_noniso_data - ss_noniso_mc)[:, None] # Data - MC in region B

            # create histo for them
            ss_iso_data_minus_mc_values = np.squeeze(np.nan_to_num(ss_iso_data_minus_mc()), axis=0) # get the values of the cat A
            ss_iso_data_minus_mc_variances = ss_iso_data_minus_mc(sn.UP, sn.ALL, unc=True)**2
            ss_iso_data_minus_mc_variances = ss_iso_data_minus_mc_variances[0]

            ss_noniso_data_minus_mc_values = np.squeeze(np.nan_to_num(ss_noniso_data_minus_mc()), axis=0) # get the values of the cat B
            ss_noniso_data_minus_mc_variances = ss_noniso_data_minus_mc(sn.UP, sn.ALL, unc=True)**2
            ss_noniso_data_minus_mc_variances = ss_noniso_data_minus_mc_variances[0]

            # guranty positive values 
            neg_int_mask = ss_iso_data_minus_mc_values <= 0
            ss_iso_data_minus_mc_values[neg_int_mask] = 1e-5
            ss_iso_data_minus_mc_variances[neg_int_mask] = 0

            neg_int_mask = ss_noniso_data_minus_mc_values <= 0
            ss_noniso_data_minus_mc_values[neg_int_mask] = 1e-5
            ss_noniso_data_minus_mc_variances[neg_int_mask] = 0

            # create a hist clone of the data_hist
            ss_noniso_data_minus_mc_hist = data_hist.copy()
            ss_iso_data_minus_mc_hist = data_hist.copy()

            # fill the ratio histogram with the values 
            ss_noniso_data_minus_mc_hist.view().value[0, ...] = ss_noniso_data_minus_mc_values
            ss_noniso_data_minus_mc_hist.view().variance[0, ...] = ss_noniso_data_minus_mc_variances
            ss_iso_data_minus_mc_hist.view().value[0, ...] = ss_iso_data_minus_mc_values
            ss_iso_data_minus_mc_hist.view().variance[0, ...] = ss_iso_data_minus_mc_variances

            # Save the fake factor histogram in a pickle file
            #path = "/eos/user/o/oponcet2/analysis/CP_dev/analysis_httcp/cf.PlotVariables1D/FF"
            path = f"{law.wlcg.WLCGFileSystem().base[0].split('root://eosuser.cern.ch')[-1]}/analysis_httcp/cf.PlotVariables1D/{config.campaign.name}/FF"

            hname_noniso = "ss_noniso_data_minus_mc_hist"
            hname_iso = "ss_iso_data_minus_mc_hist"
            # Ensure the folder exists
            if not os.path.exists(path):
                os.makedirs(path)

            with open(f"{path}/fake_factors_{hname_noniso}_{group_name}.pkl", "wb") as f:
                pickle.dump(ss_noniso_data_minus_mc_hist, f)
            
            with open(f"{path}/fake_factors_{hname_iso}_{group_name}.pkl", "wb") as f:
                pickle.dump(ss_iso_data_minus_mc_hist, f)
            

            # calculate the pt-dependent fake factor
            fake_factor = ((ss_iso_data - ss_iso_mc) / (ss_noniso_data - ss_noniso_mc))[:, None] # FF = (data-mc)A/(data-mc)B
            fake_factor_values = np.squeeze(np.nan_to_num(fake_factor()), axis=0) # get the values of the fake factor
            fake_factor_variances = fake_factor(sn.UP, sn.ALL, unc=True)**2 # get the uncertainties of the fake factor

            fake_factor_variances = fake_factor_variances[0]

            # Guaranty positive values of fake_factor
            neg_int_mask = fake_factor_values <= 0
            fake_factor_values[neg_int_mask] = 1e-5
            fake_factor_variances[neg_int_mask] = 0
        

            # create a hist clone of the data_hist
            ratio_hist = data_hist.copy()

            # fill the ratio histogram with the fake factor values in the first category
            ratio_hist.view().value[0, ...] = fake_factor_values
            ratio_hist.view().variance[0, ...] = fake_factor_variances

            # Save the fake factor histogram in a pickle file
            #path = "/eos/user/o/oponcet2/analysis/CP_dev/analysis_httcp/cf.PlotVariables1D/FF"
            #path = f"{law.wlcg.WLCGFileSystem().base[0].split('root://eosuser.cern.ch')[-1]}/analysis_httcp/cf.PlotVariables1D/{config.campaign.name}/FF"
            #hname = ratio_hist.axes[2].name
            
            # Ensure the folder exists
            #if not os.path.exists(path):
            #    os.makedirs(path)
            #with open(f"{path}/fake_factors_{hname}_{group_name}.pkl", "wb") as f:
            #    pickle.dump(ratio_hist, f)
        
        return hists

    def produce_fake_factor_DM0(task, hists):
        '''
        Produce fake factor histograms for ABCD or A0B0C0D0 categories
        FF = (data-mc)A/(data-mc)B
        FF0 = (data-mc)A0/(data-mc)B0
        The resulting fake factors histograms are stroe in a pickle file.

        *** N.B. 
        * This will not save the distribution of the varaible in pdf format, but the FF pkl files will be saved. 
        * One can use one histogram only, e.g. hcand_1_pt
        '''
        # Define if we are calculating the fake factor for A0B0C0D0 or ABCD categories
        iszero_category = False # True for A0B0C0D0 categories

        # Check if histograms are available
        if not hists:
            print("no hists")
            return hists

        # Get the qcd process, it will be used to store the fake factor histograms
        qcd_proc = config.get_process("qcd", default=None)
        if not qcd_proc:
            print("no qcd process")
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

        # create qcd groups: A, B, C, D of A0, B0, C0, D0 for each DM and Njet category
        qcd_groups: dict[str, dict[str, od.Category]] = defaultdict(DotDict)


        dms = ["tau1a1DM11", "tau1a1DM10", "tau1a1DM2", "tau1pi", "tau1rho"]  # Decay modes
        njets = ["has0j", "has1j", "has2j"]  # Jet multiplicity

        # Loop over all categories and create a QCD group for each DM and Njet category
        for dm in dms:
            for njet in njets:
                for cat_id in category_ids:
                    cat_inst = config.get_category(cat_id)

                    if iszero_category: # A0B0C0D0 categories
                        if cat_inst.has_tag({"ss", "iso1", "noniso2", njet, dm}, mode=all): # cat A0
                            qcd_groups[f"dm_{dm}_njet_{njet}"].ss_iso = cat_inst
                        elif cat_inst.has_tag({"ss", "noniso1", "noniso2", njet, dm}, mode=all): # cat BB
                            qcd_groups[f"dm_{dm}_njet_{njet}"].ss_noniso = cat_inst
                    else: # ABCD categories
                        if cat_inst.has_tag({"ss", "iso1", njet, dm}, mode=all) and not cat_inst.has_tag("noniso2"): # cat A 
                            qcd_groups[f"dm_{dm}_njet_{njet}"].ss_iso = cat_inst
                        elif cat_inst.has_tag({"ss", "noniso1", njet, dm}, mode=all) and not cat_inst.has_tag("noniso2"): # cat B
                            qcd_groups[f"dm_{dm}_njet_{njet}"].ss_noniso = cat_inst

        # Get complete qcd groups
        complete_groups = [name for name, cats in qcd_groups.items() if len(cats) == 2]

        # Nothing to do if there are no complete groups, you need A and B to estimate FF
        if not complete_groups:
            print("no complete groups")
            return hists
        
        # Sum up mc and data histograms, stop early when empty, this is done for all categories
        mc_hists = [h for p, h in hists.items() if p.is_mc and not p.has_tag("signal")]
        data_hists = [h for p, h in hists.items() if p.is_data]
        if not mc_hists or not data_hists:
            return hists
        mc_hist = sum(mc_hists[1:], mc_hists[0].copy()) # sum all MC histograms, the hist object here contains all categories
        data_hist = sum(data_hists[1:], data_hists[0].copy()) # sum all data histograms, the hist object here contains all categories

        ss_iso_mc_inclusif_njet = None
        ss_iso_data_inclusif_njet = None
        ss_noniso_mc_inclusif_njet = None
        ss_noniso_data_inclusif_njet = None


        for gidx, group_name in enumerate(complete_groups):

            group = qcd_groups[group_name]
            # Get the corresponding histograms of the id, if not present, create a zeroed histogram
            get_hist = lambda h, region_name: (
                h[{"category": hist.loc(group[region_name].id)}]
                if group[region_name].id in h.axes["category"]
                else hist.Hist(*[axis for axis in (h[{"category": [0]}] * 0).axes if axis.name != 'category'])
            )

            # # Get the corresponding histograms and convert them to number objects,
            ss_iso_mc = hist_to_num(get_hist(mc_hist, "ss_iso"), "ss_iso_mc") # MC in region A
            ss_iso_data = hist_to_num(get_hist(data_hist, "ss_iso"), "ss_iso_data") # Data in region A
            ss_noniso_mc  = hist_to_num(get_hist(mc_hist, "ss_noniso"), "ss_noniso_mc") # MC in region B
            ss_noniso_data = hist_to_num(get_hist(data_hist, "ss_noniso"), "ss_noniso_data") # Data in region B

            print("group_name", group_name)
            ## Add ss_iso_mc if dm is tau1pi
            if group_name.split("_")[1] == "tau1pi":
                print("choosen group", group_name)
                # if no already found one category, initialize the sum
                if ss_iso_mc_inclusif_njet is None:
                    ss_iso_mc_inclusif_njet = ss_iso_mc
                    ss_noniso_mc_inclusif_njet = ss_noniso_mc
                    ss_iso_data_inclusif_njet = ss_iso_data
                    ss_noniso_data_inclusif_njet = ss_noniso_data
                    print("ss_iso_data_inclusif_njet", ss_iso_data_inclusif_njet)
                else:
                    # add the current category to the sum
                    ss_iso_mc_inclusif_njet += ss_iso_mc
                    ss_noniso_mc_inclusif_njet += ss_noniso_mc
                    ss_iso_data_inclusif_njet += ss_iso_data
                    ss_noniso_data_inclusif_njet += ss_noniso_data


                print("ss_iso_data_inclusif_njet", ss_iso_data_inclusif_njet)



            else:
                # skip the category if it is not tau1pi
                continue

        
        # calculate the pt-dependent fake factor
        fake_factor = ((ss_iso_data_inclusif_njet - ss_iso_mc_inclusif_njet) / (ss_noniso_data_inclusif_njet - ss_noniso_mc_inclusif_njet))[:, None] # FF = (data-mc)A/(data-mc)B
        fake_factor_values = np.squeeze(np.nan_to_num(fake_factor()), axis=0) # get the values of the fake factor
        fake_factor_variances = fake_factor(sn.UP, sn.ALL, unc=True)**2 # get the uncertainties of the fake factor

        fake_factor_variances = fake_factor_variances[0]

        # Guaranty positive values of fake_factor
        neg_int_mask = fake_factor_values <= 0
        fake_factor_values[neg_int_mask] = 1e-5
        fake_factor_variances[neg_int_mask] = 0
    

        # create a hist clone of the data_hist
        ratio_hist = data_hist.copy()

        # fill the ratio histogram with the fake factor values in the first category
        ratio_hist.view().value[0, ...] = fake_factor_values
        ratio_hist.view().variance[0, ...] = fake_factor_variances

        # Save the fake factor histogram in a pickle file
        #path = "/eos/user/o/oponcet2/analysis/CP_dev/analysis_httcp/cf.PlotVariables1D/FF"
        path = f"{law.wlcg.WLCGFileSystem().base[0].split('root://eosuser.cern.ch')[-1]}/analysis_httcp/cf.PlotVariables1D/{config.campaign.name}/FF"
        hname = ratio_hist.axes[2].name
        
        # Ensure the folder exists
        if not os.path.exists(path):
            os.makedirs(path)
        # with open(f"{path}/fake_factors_{hname}_{group_name}.pkl", "wb") as f:
        #     pickle.dump(ratio_hist, f)

        with open(f"{path}/fake_factors_{hname}_DM0.pkl", "wb") as f:
            pickle.dump(ratio_hist, f)
        
        return hists


 
    ######################################################
    ###           EXTRAPOLATE FAKE PROCESS             ###
    ######################################################
    def extrapolate_fake(task, hists):
        '''
        This is a multi task function that can create extrapolate the fake process of one region to another.
        - FF x B -> A : type_extrapolation = "AB"; also create ratio plot of DATA/MC of A region for closure correction
        - FF x C -> D : type_extrapolation = "CD" it's control plots
        - FF0 x C0 -> D0 : type_extrapolation = "C0D0"; also create ratio plot of DATA/MC oof D0 region for extrapolation correction
        '''
        # Choose the type of extrapolation
        #type_extrapolation = "CD" # "AB" or "CD" or "C0D0"
        type_extrapolation = config.x.regions_to_extrapolate_fake

        
        # Check if histograms are available
        if not hists:
            print("no hists")
            return hists

        # Get the qcd proces, this will be used as the fake process
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

        ### Create a QCD group ofr each DM and Njet category
        qcd_groups: dict[str, dict[str, od.Category]] = defaultdict(DotDict)

        #dms = ["tau1a1DM11", "tau1a1DM10", "tau1a1DM2", "tau1pi", "tau1rho"]  # Decay modes
        dms = ["tau1a1DM10", "tau1a1DM2", "tau1pi", "tau1rho"]  # Decay modes
        njets = ["has0j", "has1j", "has2j"]  # Jet multiplicity


        # Loop over all categories and create a QCD group for each DM and Njet category
        for dm in dms:
            for njet in njets:
                for cat_id in category_ids:
                    cat_inst = config.get_category(cat_id)

                    # CASE OF AB CLOSURE PLOTS
                    if type_extrapolation == "AB":
                        if cat_inst.has_tag({"ss", "iso1", njet, dm}, mode=all) and not cat_inst.has_tag("noniso2"): # cat A 
                            qcd_groups[f"dm_{dm}_njet_{njet}"].os_iso = cat_inst
                        elif cat_inst.has_tag({"ss", "noniso1", njet, dm}, mode=all) and not cat_inst.has_tag("noniso2"): # cat B
                            qcd_groups[f"dm_{dm}_njet_{njet}"].os_noniso = cat_inst

                    # CASE OF CD CONTROL PLOTS
                    elif type_extrapolation == "CD":
                        if cat_inst.has_tag({"os", "iso1", njet, dm}, mode=all) and not cat_inst.has_tag("noniso2"): # cat D
                            qcd_groups[f"dm_{dm}_njet_{njet}"].os_iso = cat_inst
                        elif cat_inst.has_tag({"os", "noniso1", njet, dm}, mode=all) and not cat_inst.has_tag("noniso2"): # cat C
                            qcd_groups[f"dm_{dm}_njet_{njet}"].os_noniso = cat_inst
                    
                    # CASE OF C0D0 CONTROL PLOTS
                    elif type_extrapolation == "C0D0":
                        if cat_inst.has_tag({"os", "iso1", njet, dm}, mode=all) and not cat_inst.has_tag("noniso2"): # cat D0
                            qcd_groups[f"dm_{dm}_njet_{njet}"].os_iso = cat_inst
                        elif cat_inst.has_tag({"os", "noniso1", njet, dm}, mode=all) and not cat_inst.has_tag("noniso2"): # cat C0
                            qcd_groups[f"dm_{dm}_njet_{njet}"].os_noniso = cat_inst   

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
        path = f"{law.wlcg.WLCGFileSystem().base[0].split('root://eosuser.cern.ch')[-1]}/analysis_httcp/cf.PlotVariables1D/{config.campaign.name}/QCD"
        if not os.path.exists(path):
            os.makedirs(path)
        ratio_path = f"{path}/Ratio"
        if not os.path.exists(ratio_path):
            os.makedirs(ratio_path)
                
        for gidx, group_name in enumerate(complete_groups):
            group = qcd_groups[group_name]
            logger.info(f"Group name : {group_name}")

            #from IPython import embed; embed()
            
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
            """
            if type_extrapolation != "CD":
                # Save tne qcd histogram in a pickle file
                hname = qcd_hist.axes[2].name
                #path = "/eos/user/o/oponcet2/analysis/CP_dev/analysis_httcp/cf.PlotVariables1D/QCD"
                # Ensure the folder exists
                #if not os.path.exists(path):
                #    os.makedirs(path)
                with open(f"{path}/qcd_{hname}_{group_name}.pkl", "wb") as f:
                    pickle.dump(qcd_hist, f)

                # create a hist clone of the data_hist
                ratio_hist = data_hist.copy()

                # calultate sum_mc_hist
                mc_hists = [h for p, h in hists.items() if p.is_mc and not p.has_tag("signal")]
                mc_hist_sum = sum(mc_hists[1:], mc_hists[0].copy())

                # For inclusive region
                mc_hist_incl = mc_hist_incl + mc_hist_sum.copy()
                data_hist_incl = data_hist_incl + data_hist.copy()
                
                os_iso_mc  = hist_to_num(get_hist(mc_hist_sum, "os_iso"), "os_iso_mc")
                os_iso_data = hist_to_num(get_hist(data_hist, "os_iso"), "os_iso_data")
                
                # Calucate the DATA/MC ratio
                ratio = os_iso_data/os_iso_mc

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
                cat_axis = ratio_hist.axes["category"]
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
            
                #path = f"{path}/Ratio"
                # save the ratio in a pickle file
                with open(f"{ratio_path}/ratio_{hname}_{group_name}_{group.os_iso.id}.pkl", "wb") as f:
                    pickle.dump(ratio_hist, f)
            """        
        """
        if type_extrapolation != "CD":
            # Save the inclusive histograms
            incl_ratio = os_iso_data_incl/os_iso_mc_incl

            incl_ratio_hist_values = incl_ratio()
            incl_ratio_hist_variances = incl_ratio(sn.UP, sn.ALL, unc=True)**2

            incl_ratio_hist = mc_hist_incl.copy().reset()
            incl_ratio_hist.view().value[0, ...] = incl_ratio_hist_values
            incl_ratio_hist.view().variance[0, ...] = incl_ratio_hist_variances
        
            with open(f"{ratio_path}/RATIO_{hname}_inclusive.pkl", "wb") as f:
                pickle.dump(incl_ratio_hist, f) #data_hist_incl, f)
        """  
        return hists
    
    

    def extrapolate_fake_classifier_tautau(task, hists):
        '''
        This is a multi task function that can create extrapolate the fake process of one region to another.
        - FF x B -> A : type_extrapolation = "AB"; also create ratio plot of DATA/MC of A region for closure correction
        - FF x C -> D : type_extrapolation = "CD" it's control plots
        - FF0 x C0 -> D0 : type_extrapolation = "C0D0"; also create ratio plot of DATA/MC oof D0 region for extrapolation correction
        '''
        # Choose the type of extrapolation
        #type_extrapolation = "CD" # "AB" or "CD" or "C0D0"
        type_extrapolation = config.x.regions_to_extrapolate_fake

        
        # Check if histograms are available
        if not hists:
            print("no hists")
            return hists

        # Get the qcd proces, this will be used as the fake process
        qcd_proc = config.get_process("qcd", default=None)
        if not qcd_proc:
            print("no fake") 
            return hists

        # extract all unique category ids and verify that the axis order is exactly
        # "category -> shift -> variable" which is needed to insert values at the end
        CAT_AXIS, SHIFT_AXIS, VAR_AXIS = range(3)
        category_ids = set()
        for proc, h in hists.items():
            print(f"proc : {proc}")
            # validate axes
            assert len(h.axes) == 3
            assert h.axes[CAT_AXIS].name == "category"
            assert h.axes[SHIFT_AXIS].name == "shift"
            # get the category axis
            cat_ax = h.axes["category"]
            for cat_index in range(cat_ax.size):
                #print(cat_index, cat_ax.value(cat_index))
                #from IPython import embed; embed()
                #if cat_index >= len(cat_ax):
                #    logger.warning(f"{cat_index} not found ... skipping")
                #    continue
                category_ids.add(cat_ax.value(cat_index))

        ### Create a QCD group ofr each DM and Njet category
        qcd_groups: dict[str, dict[str, od.Category]] = defaultdict(DotDict)

        dms = ["pi_pi","pi_rho","pi_a1DM2","pi_a1DM10",
               "rho_rho","rho_a1DM2","rho_a1DM10",
               "a1DM2_a1DM2","a1DM2_a1DM10",
               "a1DM10_a1DM10"]
        xgb_nodes = ["dy_node", "fake_node", "higgs_node"]

        # Loop over all categories and create a QCD group for each DM and Njet category
        for dm in dms:
            for xgb in xgb_nodes:
                for cat_id in category_ids:
                    cat_inst = config.get_category(cat_id)
                    #if cat_inst is None:
                    #    print(f"WARNING: No category config found for ID {cat_id}")
                    #else:
                    #    print(f"{cat_id} tags: {cat_inst.tags}")

                    # CASE OF AB CLOSURE PLOTS
                    if type_extrapolation == "AB":
                        if cat_inst.has_tag({"ss", "iso1", xgb, dm}, mode=all) and not cat_inst.has_tag("noniso2"): # cat A 
                            qcd_groups[f"bdt_{xgb}_dm_{dm}"].os_iso = cat_inst
                        elif cat_inst.has_tag({"ss", "noniso1", xgb, dm}, mode=all) and not cat_inst.has_tag("noniso2"): # cat B
                            qcd_groups[f"bdt_{xgb}_dm_{dm}"].os_noniso = cat_inst
                                
                    # CASE OF CD CONTROL PLOTS
                    elif type_extrapolation == "CD":
                        if cat_inst.has_tag({"os", "iso1", xgb, dm}, mode=all) and not cat_inst.has_tag("noniso2"): # cat D
                            qcd_groups[f"bdt_{xgb}_dm_{dm}"].os_iso = cat_inst
                        elif cat_inst.has_tag({"os", "noniso1", xgb, dm}, mode=all) and not cat_inst.has_tag("noniso2"): # cat C
                            qcd_groups[f"bdt_{xgb}_dm_{dm}"].os_noniso = cat_inst
                            
                            
        # Get complete qcd groups
        complete_groups = [name for name, cats in qcd_groups.items() if len(cats) == 2]

        #from IPython import embed; embed()
        
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
        #from IPython import embed; embed()
        mc_hist_incl = mc_hist.copy().reset()
        data_hist_incl = data_hist.copy().reset()
        os_iso_mc_incl = None
        os_iso_data_incl = None
        #from IPython import embed; embed()
        law_args = law.luigi.cmdline_parser.CmdlineParser.get_instance().known_args
        path = f"{law.wlcg.WLCGFileSystem().base[0].split('root://eosuser.cern.ch')[-1]}/analysis_httcp/cf.MergeHistograms/{config.campaign.name}/qcd/{law_args.shift}/calib__{config.x.default_calibrator}/sel__{config.x.default_selector}/prod__{config.x.default_producer}/weight__{config.x.default_weight_producer}/{law_args.version}"
        #path = f"{law.wlcg.WLCGFileSystem().base[0].split('root://eosuser.cern.ch')[-1]}/analysis_httcp/cf.PlotVariables1D/{config.campaign.name}/QCD"
        if not os.path.exists(path):
            os.makedirs(path)
        #ratio_path = f"{path}/Ratio"
        #if not os.path.exists(ratio_path):
        #    os.makedirs(ratio_path)
                
        for gidx, group_name in enumerate(complete_groups):
            print("Group: ",group_name)
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

            print("Data hist (noniso) integral:", os_noniso_data().sum())
            print("MC hist (noniso) integral:", os_noniso_mc().sum())
            print("Fake hist integral:", fake_hist().sum())

            
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
                    print(f"Cat index : {cat_index} ==> {cat_axis.value(cat_index)}")
                    print("Target shape:", qcd_hist.view().value[cat_index, ...].sum())
                    print("Fake values shape:", fake_hist_values.sum())
                    qcd_hist.view().value[cat_index, ...] = fake_hist_values
                    qcd_hist.view().variance[cat_index, ...] = fake_hist_variances
                    break
            else:
                raise RuntimeError(
                    f"could not find index of bin on 'category' axis of qcd histogram {mc_hist} "
                    f"for category {group.os_iso}",
                )

        #from IPython import embed; embed()
        # Save tne qcd histogram in a pickle file
        qcd_hist_new = add_axis_with_value(qcd_hist, qcd_proc.id)
        #qcd_hist_ = hist.Hist(
        #    qcd_hist.axes[0],
        #    hist.axis.IntCategory([qcd_proc.id], growth=True, name="process"),  # New axis
        #    qcd_hist.axes[1],
        #    qcd_hist.axes[2],
        #    storage=qcd_hist._storage_type()
        #)
        hname = qcd_hist_new.axes[3].name
        with open(f"{path}/hist__{hname}.pickle", "wb") as f:
            pickle.dump(qcd_hist_new, f)
          
        return hists



    """


    def extrapolate_fake_classifier_tautau_test(task, hists):
        '''
        This is a multi task function that can create extrapolate the fake process of one region to another.
        - FF x B -> A : type_extrapolation = "AB"; also create ratio plot of DATA/MC of A region for closure correction
        - FF x C -> D : type_extrapolation = "CD" it's control plots
        - FF0 x C0 -> D0 : type_extrapolation = "C0D0"; also create ratio plot of DATA/MC oof D0 region for extrapolation correction
        '''
        # Choose the type of extrapolation
        #type_extrapolation = "CD" # "AB" or "CD" or "C0D0"
        type_extrapolation = config.x.regions_to_extrapolate_fake

        
        # Check if histograms are available
        if not hists:
            print("no hists")
            return hists

        # Get the qcd proces, this will be used as the fake process
        qcd_proc = config.get_process("qcd", default=None)
        if not qcd_proc:
            print("no fake") 
            return hists

        # extract all unique category ids and verify that the axis order is exactly
        # "category -> shift -> variable" which is needed to insert values at the end
        CAT_AXIS, SHIFT_AXIS, VAR_AXIS = range(3)
        category_ids = set()
        category_names = set()
        for proc, h in hists.items():
            print(f"proc : {proc}")
            # validate axes
            assert len(h.axes) == 3
            assert h.axes[CAT_AXIS].name == "category"
            assert h.axes[SHIFT_AXIS].name == "shift"
            # get the category axis
            cat_ax = h.axes["category"]
            category_names.update(list(cat_ax))
            for cat_index in range(cat_ax.size):
                category_ids.add(cat_ax.value(cat_index))

        ### Create a QCD group ofr each DM and Njet category
        qcd_groups: dict[str, dict[str, od.Category]] = defaultdict(DotDict)

        for cat_name in category_names:
            cat_inst = config.get_category(cat_name)
            #from IPython import embed; embed()
            if cat_inst.has_tag({"tautau", "os","noniso1","iso2"}, mode=all):
                qcd_groups[cat_inst.name].os_noniso = cat_inst
            elif cat_inst.has_tag({"tautau", "os","iso1","iso2"}, mode=all):
                qcd_groups[cat_inst.name].os_iso = cat_inst

        # Get complete qcd groups
        complete_groups = [name for name, cats in qcd_groups.items()]
        #complete_groups = [name for name, cats in qcd_groups.items() if len(cats) == 2]

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
        path = f"{law.wlcg.WLCGFileSystem().base[0].split('root://eosuser.cern.ch')[-1]}/analysis_httcp/cf.PlotVariables1D/{config.campaign.name}/QCD"
        if not os.path.exists(path):
            os.makedirs(path)

        from IPython import embed; embed()
            
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


        # Save tne qcd histogram in a pickle file
        hname = qcd_hist.axes[2].name
        with open(f"{path}/qcd_{hname}.pkl", "wb") as f:
            pickle.dump(qcd_hist, f)
            
          
        return hists

    """
    

    
    config.x.hist_hooks = {
        "fake_onthefly": fake_factors_onthefly,
        "produce_fake_factor": produce_fake_factor,
        "produce_fake_factor_DM0": produce_fake_factor_DM0,
        "extrapolate_fake": extrapolate_fake,
        "extrapolate_fake_classifier_tautau": extrapolate_fake_classifier_tautau,
    }
