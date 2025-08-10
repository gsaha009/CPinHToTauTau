# coding: utf-8

"""
Exemplary selection methods.
"""

import law

from typing import Optional
from operator import and_
from functools import reduce
from collections import defaultdict, OrderedDict

from columnflow.selection import Selector, SelectionResult, selector
from columnflow.selection.stats import increment_stats
from columnflow.selection.cms.json_filter import json_filter
from columnflow.selection.cms.met_filters import met_filters

from columnflow.production.processes import process_ids
from columnflow.production.cms.mc_weight import mc_weight
from columnflow.production.util import attach_coffea_behavior
#from columnflow.production.cms.top_pt_weight import gen_parton_top

from columnflow.util import maybe_import
from columnflow.columnar_util import optional_column as optional
from columnflow.columnar_util import EMPTY_FLOAT, Route, set_ak_column

from httcp.selection.physics_objects import *
from httcp.selection.trigger import trigger_selection
from httcp.selection.lepton_pair_etau import etau_selection
from httcp.selection.lepton_pair_mutau import mutau_selection
from httcp.selection.lepton_pair_tautau import tautau_selection
from httcp.selection.event_category import get_categories #, build_abcd_masks
#from httcp.selection.match_trigobj import match_trigobj
from httcp.selection.trigobject_matching import match_trigobj
from httcp.selection.lepton_veto import *
from httcp.selection.higgscand import higgscand, higgscandprod

#from httcp.production.main import cutflow_features
#from httcp.production.weights import scale_mc_weight
from httcp.production.stitching_NLO import process_ids_dy
from httcp.production.stitching_LO import process_ids_w
from httcp.production.extra_weights import scale_mc_weight
from httcp.production.dilepton_features import hcand_mass, mT, rel_charge #TODO: rename mutau_vars -> dilepton_vars

from httcp.util import filter_by_triggers, get_objs_p4, trigger_object_matching_deep, IF_DATASET_IS_DY, IF_DATASET_IS_W, IF_DATASET_IS_SIGNAL

from httcp.selection.debug import debug_main
from columnflow.selection.cms.jets import jet_veto_map  

logger = law.logger.get_logger(__name__)

np = maybe_import("numpy")
ak = maybe_import("awkward")
coffea = maybe_import("coffea")
maybe_import("coffea.nanoevents.methods.nanoaod")

@selector(uses={"process_id", optional("mc_weight")}) #TODO: Move it to utils
def custom_increment_stats(
    self: Selector,
    events: ak.Array,
    results: SelectionResult,
    stats: dict,
    **kwargs,
) -> ak.Array:
    """
    Unexposed selector that does not actually select objects but instead increments selection
    *stats* in-place based on all input *events* and the final selection *mask*.
    This function saves the nevents/nfiles (from campaign) per file, so that at 
    the time of merging, the total nevents produced can be used for normalization.
    instead of nevents, the sum_mc_weights will be used actually.
    ** the same nevents or sum_wt will be saved for each process id.
    e.g. In an inclusive dataset, if three different processes are there,
    all of those will have the same number. After merging, all processes
    eventually will have same numbers.
    """
    # get event masks
    event_mask = results.event

    # get a list of unique process ids present in the chunk
    unique_process_ids = np.unique(events.process_id)

    # increment plain counts
    n_evt_per_file = self.dataset_inst.aux['n_events']/self.dataset_inst.n_files # new
    sumwt_per_file = self.dataset_inst.n_events/self.dataset_inst.n_files     # new

    stats["num_events"] = int(n_evt_per_file)
    stats.setdefault(f"num_events_per_process", defaultdict(int))
    for p in unique_process_ids:
        # for splitting, each process will have the same number of events
        stats[f"num_events_per_process"][int(p)] = int(n_evt_per_file)
        
        
    if self.dataset_inst.is_mc:
        stats[f"sum_mc_weight"] = sumwt_per_file
        stats.setdefault(f"sum_mc_weight_per_process", defaultdict(float))
        for p in unique_process_ids:
            # for splitting, each process will have the same sumwt 
            stats[f"sum_mc_weight_per_process"][int(p)] = sumwt_per_file
        
    return events, results


def get_2n_pairs(etau_indices_pair,
                 mutau_indices_pair,
                 tautau_indices_pair):

    has_one_etau_pair   = ak.num(etau_indices_pair, axis=1) == 2
    has_one_mutau_pair  = ak.num(mutau_indices_pair, axis=1) == 2
    has_one_tautau_pair = ak.num(tautau_indices_pair, axis=1) == 2

    has_at_least_one_pair = (has_one_etau_pair | has_one_mutau_pair | has_one_tautau_pair)

    return has_at_least_one_pair    

    
# exposed selectors
# (those that can be invoked from the command line)
@selector(
    uses={
        "event",
        # selectors / producers called within _this_ selector
        attach_coffea_behavior,
        json_filter, 
        met_filters, 
        scale_mc_weight, 
        process_ids,
        trigger_selection,
        #IF_DATASET_IS_DY(genZ_selection),
        genZ_selection,
        muon_selection, 
        electron_selection, 
        tau_selection, 
        jet_selection,
        jet_cleaning,
        etau_selection, 
        mutau_selection, 
        tautau_selection, 
        get_categories,
        extra_lepton_veto, 
        #double_lepton_veto,
        tau_veto_from_dy,
        match_trigobj,
        increment_stats, 
        custom_increment_stats,
        higgscand,
        gentau_selection,
        higgscandprod,
        rel_charge,
        #category_ids,
        #build_abcd_masks,
        jet_veto_map,
        "PuppiMET.pt", "Jet.pt", "Jet.eta", "Jet.phi", "Jet.neEmEF", "Jet.chEmEF",
        "Flag.ecalBadCalibFilter",
    },
    produces={
        # selectors / producers whose newly created columns should be kept
        scale_mc_weight, 
        trigger_selection,
        #IF_DATASET_IS_DY(genZ_selection),
        genZ_selection,
        muon_selection, 
        electron_selection, 
        tau_selection, 
        jet_selection,
        jet_cleaning,
        etau_selection, 
        mutau_selection, 
        tautau_selection, 
        get_categories, 
        process_ids,
        extra_lepton_veto, 
        #double_lepton_veto, 
        match_trigobj,
        higgscandprod,
        gentau_selection,
        rel_charge,
        #category_ids,
        increment_stats, 
        custom_increment_stats,
        #"trigger_ids",
        "single_triggered",
        "cross_triggered",
        "single_e_triggered",
        "single_mu_triggered",
        "cross_e_triggered",
        "cross_mu_triggered",
        "cross_tau_triggered",
        "cross_tau_jet_triggered",
        #build_abcd_masks,
    },
    exposed=True,
)
def main(
    self: Selector,
    events: ak.Array,
    stats: defaultdict,
    **kwargs,
) -> tuple[ak.Array, SelectionResult]:
    
    # ensure coffea behaviors are loaded
    events = self[attach_coffea_behavior](events, **kwargs)

    # prepare the selection results that are updated at every step
    results = SelectionResult()

    if self.dataset_inst.has_tag("is_dy_m50"):
        events, tau_veto_result = self[tau_veto_from_dy](events, **kwargs)
        results += tau_veto_result
    
    # add the mc weight --> need to move to calibration main?
    if self.dataset_inst.is_mc:
        events = self[scale_mc_weight](events, **kwargs)

    ##################################### NEW ######################################
    results += SelectionResult(steps={"starts_with": np.ones(len(events), dtype=bool)})
    ##################################### NEW ######################################
        
    # filter bad data events according to golden lumi mask
    if self.dataset_inst.is_data:
        events, json_filter_results = self[json_filter](events, **kwargs)
        results += json_filter_results
    else:
        results += SelectionResult(steps={"json": np.ones(len(events), dtype=bool)})


    #if self.config_inst.campaign.x.run == 3:
    events, veto_result = self[jet_veto_map](events, **kwargs)
    results += veto_result

    # ############################ #
    #     met filter selection     #
    # ############################ #
    events, met_filter_results = self[met_filters](events, **kwargs)
    #from IPython import embed; embed()
    if self.dataset_inst.is_data:
        # BadCalibrationFilter is meant to reject events with noise which are related to bad crystals in Ecal.
        # The issue was present in some run ranges in data (end of 2022 & early 2023) only. The fraction of bad
        # data in terms of lumi is small. So, one can ignore the recipe in MC.
        # Also, DO NOT USE THIS FLAG AT ALL.
        #BadCalibrationFilter = ak.values_astype(events.Flag.ecalBadCalibFilter, bool)
        BadCalibrationFilter = events.event > 0 # all True
        met = ak.with_name(events.PuppiMET, "PtEtaPhiMLorentzVector")
        jet = ak.with_name(events.Jet, "PtEtaPhiMLorentzVector")
        BadCalibrationFilter_perjet_mask = (
            (events.PuppiMET.pt > 100)
            & (events.Jet.pt > 50)
            & (events.Jet.eta > -0.5) & (events.Jet.eta < -0.1)
            & (events.Jet.phi > -2.1) & (events.Jet.phi < -1.8)
            & ((events.Jet.neEmEF > 0.9) | (events.Jet.chEmEF > 0.9))
            & ak.all(met.metric_table(jet, metric=lambda a,b: a.delta_phi(b)) > 2.9, axis=1)
        )
        BadCalibrationFilter_mask = ak.values_astype(ak.any(BadCalibrationFilter_perjet_mask, axis=1), bool)
        # Keep the default BadCalibrationFilter if outside the run-range, else use the calibrationFilter_mask
        BadCalibrationFilter = ak.where((events.run >= 362433) & (events.run <= 367144), ~BadCalibrationFilter_mask, BadCalibrationFilter)
        
        BadCalibrationFilter_results = SelectionResult(steps={"met_filter:BadCalibration" : BadCalibrationFilter}) # negate the entire mask
        met_filter_results += BadCalibrationFilter_results
        
    results += met_filter_results    
    
    # trigger selection
    events, trigger_results = self[trigger_selection](events, **kwargs)
    results += trigger_results

    # Get genZ collection for Zpt reweighting
    #if self.dataset_inst.has_tag("is_dy"):
    if self.dataset_inst.is_mc:
        events = self[genZ_selection](events, **kwargs)
        
    # electron selection
    # e.g. ele_idx: [ [], [0,1], [], [], [1,2] ] 
    events, ele_results, good_ele_indices, veto_ele_indices, dlveto_ele_indices = self[electron_selection](events,
                                                                                                           call_force=True,
                                                                                                           **kwargs)
    results += ele_results

    
    # muon selection
    # e.g. mu_idx: [ [0,1], [], [1], [0], [] ] 
    events, muon_results, good_muon_indices, veto_muon_indices, dlveto_muon_indices = self[muon_selection](events,
                                                                                                           call_force=True,
                                                                                                           **kwargs)
    results += muon_results


    # tau selection
    # three times for three channels
    # e.g. tau_idx: [ [1], [0,1], [1,2], [], [0,1] ]
    good_tau_indices = None
    good_tau_indices_mutau = None
    good_tau_indices_etau  = None
    good_tau_indices_tautau= None

    if self.dataset_inst.is_mc:
        logger.info("MC tau_selection for mutau channel")
        _, _, good_tau_indices_mutau = self[tau_selection](events, "mutau",  call_force=True, **kwargs)
        logger.info("MC tau_selection for etau channel")
        _, _, good_tau_indices_etau       = self[tau_selection](events, "etau",   call_force=True, **kwargs)
        logger.info("MC tau_selection for tautau channel")
        events, tau_results_tau, good_tau_indices_tautau     = self[tau_selection](events, "tautau", call_force=True, **kwargs)
        results += tau_results_tau
    else:
        logger.info("DATA tau_selection for tautau channel")
        events, tau_results, good_tau_indices = self[tau_selection](events, call_force=True, **kwargs)
        results += tau_results


    # jet selection
    events, jet_results, jet_indices = self[jet_selection](events, call_force=True, **kwargs)
    

        
    # double lepton veto
    """
    events, extra_double_lepton_veto_results = self[double_lepton_veto](events,
                                                                        dlveto_ele_indices,
                                                                        dlveto_muon_indices)
    results += extra_double_lepton_veto_results
    """
        
    # e-tau pair i.e. hcand selection
    etau_results, etau_pair, etau_trig_ids, etau_trig_types = self[etau_selection](events,
                                                                                   good_ele_indices,
                                                                                   good_tau_indices_etau if self.dataset_inst.is_mc else good_tau_indices,
                                                                                   trigger_results,
                                                                                   call_force=True,
                                                                                   **kwargs)
    results += etau_results

    # mu-tau pair i.e. hcand selection
    mutau_results, mutau_pair, mutau_trig_ids, mutau_trig_types = self[mutau_selection](events,
                                                                                        good_muon_indices,
                                                                                        good_tau_indices_mutau if self.dataset_inst.is_mc else good_tau_indices,
                                                                                        trigger_results,
                                                                                        call_force=True,
                                                                                        **kwargs)
    results += mutau_results

    # tau-tau pair i.e. hcand selection
    tautau_results, tautau_pair, tautau_trig_ids, tautau_trig_types, ditaujet_jet_indices = self[tautau_selection](events,
                                                                                                                   good_tau_indices_tautau if self.dataset_inst.is_mc else good_tau_indices,
                                                                                                                   trigger_results,
                                                                                                                   jet_indices,
                                                                                                                   call_force=True,
                                                                                                                   **kwargs)
    results += tautau_results

    match_at_least_one_pair = get_2n_pairs(etau_pair.rawIdx,
                                           mutau_pair.rawIdx,
                                           tautau_pair.rawIdx)


    match_pair_result = SelectionResult(
        steps = {
            "has_at_least_1_pair"  : match_at_least_one_pair,
        },
    )
    results += match_pair_result

    
    # channel selection
    # channel_id is now in columns
    events, channel_results = self[get_categories](events,
                                                   results,
                                                   etau_pair.rawIdx,
                                                   mutau_pair.rawIdx,
                                                   tautau_pair.rawIdx)
    results += channel_results


    
    trigger_types = ak.concatenate([etau_trig_types, mutau_trig_types, tautau_trig_types], axis=1)

    # ele
    single_e_triggered = ak.any(trigger_types == 'single_e', axis=1)
    cross_e_triggered  = ak.any(trigger_types == 'cross_e_tau', axis=1)
    # mu
    single_mu_triggered = ak.any(trigger_types == 'single_mu', axis=1)
    cross_mu_triggered = ak.any(trigger_types == 'cross_mu_tau', axis=1)
    # tau
    cross_tau_triggered = ak.any(trigger_types == 'cross_tau_tau', axis=1)
    cross_tau_jet_triggered = ak.any(trigger_types == 'cross_tau_tau_jet', axis=1)

    events = set_ak_column(events, "single_e_triggered", single_e_triggered)
    events = set_ak_column(events, "single_mu_triggered", single_mu_triggered)
    events = set_ak_column(events, "cross_e_triggered", cross_e_triggered)
    events = set_ak_column(events, "cross_mu_triggered", cross_mu_triggered)
    events = set_ak_column(events, "cross_tau_triggered", cross_tau_triggered)
    events = set_ak_column(events, "cross_tau_jet_triggered", cross_tau_jet_triggered)
    
    events = set_ak_column(events, "single_triggered", (single_e_triggered | single_mu_triggered))
    events = set_ak_column(events, "cross_triggered", (cross_e_triggered | cross_mu_triggered | cross_tau_triggered | cross_tau_jet_triggered))

    
    trigger_ids = ak.concatenate([etau_trig_ids, mutau_trig_ids, tautau_trig_ids], axis=1)

    
    # save the trigger_ids column
    #events = set_ak_column(events, "trigger_ids", trigger_ids)

    # hcand pair: [ [[mu1,tau1]], [[e1,tau1],[tau1,tau2]], [[mu1,tau2]], [], [[e1,tau2]] ]
    hcand_pairs = ak.concatenate([etau_pair[:,None], mutau_pair[:,None], tautau_pair[:,None]], axis=1)

    # extra lepton veto
    # it is only applied on the events with one higgs candidate only
    events, extra_lepton_veto_results = self[extra_lepton_veto](events,
                                                                veto_ele_indices,
                                                                veto_muon_indices,
                                                                hcand_pairs)
    results += extra_lepton_veto_results

    # hcand results
    events, hcand_array, hcand_results = self[higgscand](events, hcand_pairs)
    results += hcand_results
    
    # hcand prod results
    # -------------------------------------------- #
    # Here, inside hcandprod, the replacement of raw
    # taus by calibrated taus takes place
    # -------------------------------------------- #
    events, hcandprod_results = self[higgscandprod](events, hcand_array)
    results += hcandprod_results

    # -------- Sel : b-veto -------- #
    # jet selection
    # -------------------------------------------- #
    # this is moved here, because now the jets are
    # cleaned against the tau cadidates of hacnd
    # -------------------------------------------- #

    """
    events, bjet_veto_result = self[jet_selection](events,
                                                   call_force=True, 
                                                   **kwargs)
    results += bjet_veto_result
    """
    events, jet_clean_result, ditaujet_jet_indices = self[jet_cleaning](events,
                                                                        jet_indices,
                                                                        jet_results,
                                                                        ditaujet_jet_indices,
                                                                        call_force=True, 
                                                                        **kwargs)
    results += jet_clean_result


    # gen particles info
    # ############################################ #
    # After building the higgs candidates, one can
    # switch on the production of GenTau. Those gentaus
    # will be selected which are matched to the hcand
    # hcand-gentau match = True/False (via config)
    # ############################################ #
    if self.config_inst.x.extra_tags.genmatch:
        if "is_signal" in list(self.dataset_inst.aux.keys()):
            print(" --->>> hcand-gentau matching")
            events, gentau_results = self[gentau_selection](events, True)
            results += gentau_results


    # combined event selection after all steps
    event_sel = reduce(and_, results.steps.values())
    results.event = event_sel
    
    # rel-charge
    events = self[rel_charge](events, **kwargs)

    outliers_mask_for_stitching = None
    
    if self.process_ids_dy is not None:
        events, outliers_mask_for_stitching = self[self.process_ids_dy](events, **kwargs)
    elif self.process_ids_w is not None:
        events, outliers_mask_for_stitching = self[self.process_ids_w](events, **kwargs)
    else:
        events = self[process_ids](events, **kwargs)

    if outliers_mask_for_stitching is not None:
        n_outliers = ak.sum(outliers_mask_for_stitching)
        if n_outliers > 0:
            logger.warning(f"{self.dataset_inst} has {ak.sum(outliers_mask_for_stitching)} outlier events. Removing those events only")
        outliers_result_for_stitching = SelectionResult(
            steps = {"reject_stitching_outliers": ~outliers_mask_for_stitching}
        )
        results += outliers_result_for_stitching

        
    #events, category_ids_debug_dict = self[category_ids](events, debug=True)


    events, results = self[custom_increment_stats]( 
        events,
        results,
        stats,
    )

    # On top of custom stats, the incremet_stat from columnflow is used here, just to save the
    # number of selected events per process and sum_mc_weight as well. At the time of stitching,
    # the fraction of sum_wt_selected/total_sum_wt per process can be multiplied to the total
    # sum_wt produced to get the fraction of that particular exclusive process embeded into the
    # inclusive process.
    weight_map = {
        "num_filtered_events": Ellipsis,
        "num_events_selected": event_sel,
    }
    group_map = {}
    group_combinations = []
    if self.dataset_inst.is_mc:
        weight_map["sum_filtered_mc_weight"] = events.mc_weight
        weight_map["sum_mc_weight_selected"] = (events.mc_weight, event_sel)
        # groups
        group_map = {
            **group_map,
            # per process
            "process": {
                "values": events.process_id,
                "mask_fn": (lambda v: events.process_id == v),
            },
        }
        # combinations
        #group_combinations.append(("process"))

    events, results = self[increment_stats](
        events,
        results,
        stats,
        weight_map=weight_map,
        group_map=group_map,
        #group_combinations=group_combinations,
        **kwargs,
    )
    
    # inspect cuts
    if self.config_inst.x.verbose.selection.main:
        debug_main(events,
                   results,
                   self.config_inst.x.triggers)

    return events, results


@main.init
def main_init(self: Selector) -> None:
    if getattr(self, "dataset_inst", None) is None:
        return

    self.process_ids_dy: process_ids_dy | None = None
    #from IPython import embed; embed()
    if self.dataset_inst.has_tag("is_dy_m50"):
        if self.config_inst.x.allow_dy_stitching:
            #print(f"stitching: {self.config_inst.x.dy_stitching.items()}")
            # check if this dataset is covered by any dy id producer
            for name, dy_cfg in self.config_inst.x.dy_stitching.items():
                #print(f"dataset : {name}, {dy_cfg}")
                dataset_inst = dy_cfg["inclusive_dataset"]
                # the dataset is "covered" if its process is a subprocess of that of the dy dataset
                if dataset_inst.has_process(self.dataset_inst.processes.get_first()):
                    self.process_ids_dy = process_ids_dy.derive(f"process_ids_dy_{name}", cls_dict={
                        "dy_inclusive_dataset": dataset_inst,
                        "dy_leaf_processes": dy_cfg["leaf_processes"],
                    })

                    # add it as a dependency
                    self.uses.add(self.process_ids_dy)
                    self.produces.add(self.process_ids_dy)
                    
                    # stop after the first match
                    break
            
    self.process_ids_w: process_ids_w | None = None
    if self.dataset_inst.has_tag("is_w"):
        if self.config_inst.x.allow_w_stitching:
            # check if this dataset is covered by any dy id producer
            for name, w_cfg in self.config_inst.x.w_stitching.items():
                dataset_inst = w_cfg["inclusive_dataset"]
                # the dataset is "covered" if its process is a subprocess of that of the dy dataset
                if dataset_inst.has_process(self.dataset_inst.processes.get_first()):
                    self.process_ids_w = process_ids_w.derive(f"process_ids_w_{name}", cls_dict={
                        "w_inclusive_dataset": dataset_inst,
                        "w_leaf_processes": w_cfg["leaf_processes"],
                    })

                    # add it as a dependency
                    self.uses.add(self.process_ids_w)
                    self.produces.add(self.process_ids_w)

                    # stop after the first match
                    break

    #if self.dataset_inst.has_tag("ttbar"):
    #    self.uses.add(gen_parton_top)
    #    self.produces.add(gen_parton_top)
