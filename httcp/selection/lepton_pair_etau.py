# coding: utf-8

"""
Prepare h-Candidate from SelectionResult: selected lepton indices & channel_id [trigger matched] 
"""

from typing import Optional
from columnflow.selection import Selector, SelectionResult, selector
from columnflow.selection.util import create_collections_from_masks
from columnflow.util import maybe_import
from columnflow.columnar_util import EMPTY_FLOAT, Route, set_ak_column
from columnflow.columnar_util import optional_column as optional

from httcp.util import transverse_mass
from httcp.util import IF_RUN2, IF_RUN3

np = maybe_import("numpy")
ak = maybe_import("awkward")
coffea = maybe_import("coffea")
maybe_import("coffea.nanoevents.methods.nanoaod")



from httcp.util import filter_by_triggers, get_objs_p4, trigger_matching_extra, trigger_object_matching_deep



def match_trigobjs(
        leps_pair: ak.Array,
        trigger_results: SelectionResult,
        **kwargs,
) -> tuple[ak.Array, ak.Array]:

    # extract the trigger names, types & others from trigger_results.x (aux)
    trigger_ids           = trigger_results.x.trigger_ids
    trigger_types         = trigger_results.x.trigger_types
    leg1_minpt            = trigger_results.x.leg1_minpt
    leg2_minpt            = trigger_results.x.leg2_minpt
    leg1_maxeta           = trigger_results.x.leg1_maxeta
    leg2_maxeta           = trigger_results.x.leg2_maxeta
    leg1_matched_trigobjs = trigger_results.x.leg1_matched_trigobjs
    leg2_matched_trigobjs = trigger_results.x.leg2_matched_trigobjs


    has_single_e_triggers = trigger_types == "single_e"
    has_cross_e_triggers  = trigger_types == "cross_e_tau"
    has_e_triggers = (has_single_e_triggers | has_cross_e_triggers)


    etau_pair  = filter_by_triggers(leps_pair, has_e_triggers)

    eles, taus = ak.unzip(etau_pair)

    
    # Event level masks
    # if events have etau
    has_etau_pairs = ak.fill_none(ak.num(eles, axis=1) > 0, False)

    mask_has_single_e_triggers_and_has_etau_pairs = has_single_e_triggers & has_etau_pairs
    mask_has_cross_e_triggers_and_has_etau_pairs  = has_cross_e_triggers & has_etau_pairs
    mask_has_e_triggers_and_has_etau_pairs = has_e_triggers & has_etau_pairs

    # filtering out the info based on the masks defined just above
    # for etau and mutau type of events, separate masks are created
    # for single and cross triggered events
    # for single e triggers
    single_etau_trigger_types          = trigger_types[mask_has_single_e_triggers_and_has_etau_pairs]
    single_etau_trigger_ids            = trigger_ids[mask_has_single_e_triggers_and_has_etau_pairs]
    single_etau_leg_1_minpt            = leg1_minpt[mask_has_single_e_triggers_and_has_etau_pairs]
    single_etau_leg_1_maxeta           = leg1_maxeta[mask_has_single_e_triggers_and_has_etau_pairs] 
    single_etau_leg_1_matched_trigobjs = leg1_matched_trigobjs[mask_has_single_e_triggers_and_has_etau_pairs]
    # for cross e triggers
    cross_etau_trigger_types          = trigger_types[mask_has_cross_e_triggers_and_has_etau_pairs]
    cross_etau_trigger_ids            = trigger_ids[mask_has_cross_e_triggers_and_has_etau_pairs]
    cross_etau_leg_1_minpt            = leg1_minpt[mask_has_cross_e_triggers_and_has_etau_pairs]
    cross_etau_leg_2_minpt            = leg2_minpt[mask_has_cross_e_triggers_and_has_etau_pairs] 
    cross_etau_leg_1_maxeta           = leg1_maxeta[mask_has_cross_e_triggers_and_has_etau_pairs]
    cross_etau_leg_2_maxeta           = leg2_maxeta[mask_has_cross_e_triggers_and_has_etau_pairs] 
    cross_etau_leg_1_matched_trigobjs = leg1_matched_trigobjs[mask_has_cross_e_triggers_and_has_etau_pairs]
    cross_etau_leg_2_matched_trigobjs = leg2_matched_trigobjs[mask_has_cross_e_triggers_and_has_etau_pairs]
    # concatenating single and cross info, so that it does not remain depend on the trigger hierarchy
    etau_trigger_types          = ak.concatenate([single_etau_trigger_types, cross_etau_trigger_types], axis=-1)
    etau_trigger_ids            = ak.concatenate([single_etau_trigger_ids, cross_etau_trigger_ids], axis=-1)


    # to convert the masks to event level
    # e.g. events with etau pair and pass electron triggers
    mask_has_single_e_triggers_and_has_e_evt_level = ak.any(mask_has_single_e_triggers_and_has_etau_pairs, axis=1)
    mask_has_cross_e_triggers_and_has_e_evt_level  = ak.any(mask_has_cross_e_triggers_and_has_etau_pairs, axis=1)
    mask_has_e_triggers_and_has_e_evt_level        = ak.any(mask_has_e_triggers_and_has_etau_pairs, axis=1)
    
    # dummy bool array
    trigobj_matched_mask_dummy = ak.from_regular((trigger_ids > 0)[:,:0][:,None])

    single_el_trigobj_matched_mask = trigger_object_matching_deep(eles,
                                                                  single_etau_leg_1_matched_trigobjs,
                                                                  single_etau_leg_1_minpt,
                                                                  single_etau_leg_1_maxeta,
                                                                  True)

    #print("ff")
    cross_el_trigobj_matched_mask_leg1 = trigger_object_matching_deep(eles,
                                                                      cross_etau_leg_1_matched_trigobjs,
                                                                      cross_etau_leg_1_minpt,
                                                                      cross_etau_leg_1_maxeta,
                                                                      True)
    #print("gg")
    cross_el_trigobj_matched_mask_leg2 = trigger_object_matching_deep(taus,
                                                                      cross_etau_leg_2_matched_trigobjs,
                                                                      cross_etau_leg_2_minpt,
                                                                      cross_etau_leg_2_maxeta,
                                                                      True)

    cross_el_trigobj_matched_mask = (cross_el_trigobj_matched_mask_leg1 & cross_el_trigobj_matched_mask_leg2)

    

    single_el_trigobj_matched_mask_evt_level = ak.fill_none(ak.firsts(ak.any(single_el_trigobj_matched_mask, axis=1), axis=1), False)
    cross_el_trigobj_matched_mask_evt_level = ak.fill_none(ak.firsts(ak.any(cross_el_trigobj_matched_mask, axis=1), axis=1), False) 
    

    match_single = mask_has_single_e_triggers_and_has_e_evt_level & single_el_trigobj_matched_mask_evt_level
    match_cross  = mask_has_cross_e_triggers_and_has_e_evt_level & cross_el_trigobj_matched_mask_evt_level

    
    el_trigobj_matched_mask = ak.where(match_single, single_el_trigobj_matched_mask, cross_el_trigobj_matched_mask)
    el_trigobj_matched_mask = ak.fill_none(ak.firsts(el_trigobj_matched_mask, axis=-1), False)

    new_eles = eles[el_trigobj_matched_mask]
    new_taus = taus[el_trigobj_matched_mask]

    #trigIds = ak.where(match_single, single_etau_trigger_ids, cross_etau_trigger_ids)
    #trigTypes = ak.where(match_single, single_etau_trigger_types, cross_etau_trigger_types)

    trigIds = ak.concatenate([single_etau_trigger_ids, cross_etau_trigger_ids], axis=1)
    ids = ak.values_astype(trigIds, 'int64')
    trigTypes = ak.concatenate([single_etau_trigger_types, cross_etau_trigger_types], axis=1)

    
    leps_pair = ak.zip([new_eles, new_taus])
    
    return leps_pair, ids, trigTypes



def sort_pairs(dtrpairs: ak.Array)->ak.Array:
    # Just to get the indices
    # Redundatnt as already sorted by their isolation
    sorted_idx = ak.argsort(dtrpairs["0"].pfRelIso03_all, ascending=True)
    # Sort the pairs based on pfRelIso03_all of the first object in each pair
    dtrpairs = dtrpairs[sorted_idx]

    # Check if the pfRelIso03_all values are the same for the first two objects in each pair
    where_same_iso_1 = ak.fill_none(
        ak.firsts(dtrpairs["0"].pfRelIso03_all[:,:1], axis=1) == ak.firsts(dtrpairs["0"].pfRelIso03_all[:,1:2], axis=1),
        False)

    # Sort the pairs based on pt if pfRelIso03_all is the same for the first two objects
    sorted_idx = ak.where(where_same_iso_1,
                          ak.argsort(dtrpairs["0"].pt, ascending=False),
                          sorted_idx)
    dtrpairs = dtrpairs[sorted_idx]

    # Check if the pt values are the same for the first two objects in each pair    
    where_same_pt_1 = ak.fill_none(
        ak.firsts(dtrpairs["0"].pt[:,:1], axis=1) == ak.firsts(dtrpairs["0"].pt[:,1:2], axis=1),
        False
    )
    # if so, sort the pairs with tau rawDeepTau2017v2p1VSjet
    sorted_idx = ak.where(where_same_pt_1,
                          ak.argsort(dtrpairs["1"].rawDeepTau2018v2p5VSjet, ascending=False),
                          sorted_idx)
    dtrpairs = dtrpairs[sorted_idx]
    
    # check if the first two pairs have taus with same rawDeepTau2018v2p5VSjet
    where_same_iso_2 = ak.fill_none(
        ak.firsts(dtrpairs["1"].rawDeepTau2018v2p5VSjet[:,:1], axis=1) == ak.firsts(dtrpairs["1"].rawDeepTau2018v2p5VSjet[:,1:2], axis=1),
        False
    )
    # Sort the pairs based on pt if rawDeepTau2018v2p5VSjet is the same for the first two objects
    sorted_idx = ak.where(where_same_iso_2,
                          ak.argsort(dtrpairs["1"].pt, ascending=False),
                          sorted_idx)
    # finally, the pairs are sorted
    dtrpairs = dtrpairs[sorted_idx]

    return dtrpairs



@selector(
    uses={
        # Electron
        "Electron.pt", "Electron.eta", "Electron.phi", "Electron.mass",
        "Electron.charge", "Electron.pfRelIso03_all", "Electron.rawIdx",
        # Tau
        optional("Tau.pt"),
        optional("Tau.pt_etau"),
        "Tau.eta", "Tau.phi",
        optional("Tau.mass"),
        optional("Tau.mass_etau"),
        "Tau.charge", "Tau.rawDeepTau2018v2p5VSjet", "Tau.rawIdx",
        "Tau.idDeepTau2018v2p5VSjet", "Tau.idDeepTau2018v2p5VSe", "Tau.idDeepTau2018v2p5VSmu",
        optional("Tau.genPartFlav"),
        # MET
        IF_RUN2("MET.pt", "MET.phi"),
        IF_RUN3("PuppiMET.pt", "PuppiMET.phi"),
    },
    exposed=False,
)
def etau_selection(
        self: Selector,
        events: ak.Array,
        lep1_indices: ak.Array,
        lep2_indices: ak.Array,
        trigger_results: SelectionResult,
        **kwargs,
) -> tuple[SelectionResult, ak.Array, ak.Array]:

    eles  = events.Electron[lep1_indices]
    taus  = events.Tau[lep2_indices]
    
    # Extra channel specific selections on e or tau
    tau_tagger      = self.config_inst.x.deep_tau_tagger
    tau_tagger_wps  = self.config_inst.x.deep_tau_info[tau_tagger].wp
    vs_e_wp         = self.config_inst.x.deep_tau_info[tau_tagger].vs_e["etau"]
    vs_mu_wp        = self.config_inst.x.deep_tau_info[tau_tagger].vs_m["etau"]
    vs_jet_wp       = self.config_inst.x.deep_tau_info[tau_tagger].vs_j["etau"]

    if self.dataset_inst.is_mc:
        # rename "pt_etau" and "mass_etau" to "pt" and "mass"
        taus = ak.without_field(taus, "pt")
        taus = ak.with_field(taus, taus.pt_etau, "pt")
        taus = ak.without_field(taus, "mass")
        taus = ak.with_field(taus, taus.mass_etau, "mass")
    
    is_good_tau     = (
        (taus.pt > 20.0)
        #(taus.idDeepTau2018v2p5VSjet   >= tau_tagger_wps.vs_j[vs_jet_wp])
        & (taus.idDeepTau2018v2p5VSe   >= tau_tagger_wps.vs_e[vs_e_wp])
        & (taus.idDeepTau2018v2p5VSmu  >= tau_tagger_wps.vs_m[vs_mu_wp])
    )

    taus = taus[is_good_tau]

    # puppi for Run3
    met = events.MET if self.config_inst.campaign.x.year < 2022 else events.PuppiMET

    # Sorting lep1 [Electron] by isolation [ascending]
    eles_sort_idxs = ak.argsort(eles.pfRelIso03_all, axis=-1, ascending=True)
    eles = eles[eles_sort_idxs]
    taus_sort_idx = ak.argsort(taus.rawDeepTau2018v2p5VSjet, axis=-1, ascending=False)
    taus = taus[taus_sort_idx]
    
    leps_pair  = ak.cartesian([eles, taus], axis=1)
    
    lep1, lep2         = ak.unzip(leps_pair)

        
    preselection = {
        #"etau_is_os"         : (lep1.charge * lep2.charge) < 0,
        "etau_dr_0p5"        : (1*lep1).delta_r(1*lep2) > 0.5,
        #"etau_mT_50"         : transverse_mass(lep1, met) < 50
    }

    # get preselected pairs
    good_pair_mask = lep1.rawIdx >= 0
    pair_selection_steps = {}
    pair_selection_steps["etau_starts_with"] = good_pair_mask
    for cut in preselection.keys():
        good_pair_mask = good_pair_mask & preselection[cut]
        pair_selection_steps[cut] = good_pair_mask

    good_pair_mask = ak.fill_none(good_pair_mask, False)
    leps_pair  = leps_pair[good_pair_mask]
    # check nPairs
    npair = ak.num(leps_pair["0"], axis=1)
    pair_selection_steps["etau_before_trigger_matching"] = leps_pair["0"].pt >= 0.0
    
    # sort the pairs if many
    leps_pair = ak.where(npair > 1, sort_pairs(leps_pair), leps_pair)

    #from IPython import embed; embed()
    
    # match trigger objects for all pairs
    leps_pair, trigIds, trigTypes = match_trigobjs(leps_pair, trigger_results)

    pair_selection_steps["etau_after_trigger_matching"] = leps_pair["0"].pt >= 0.0

    lep1, lep2 = ak.unzip(leps_pair)


    # take the 1st pair and 1st trigger id
    lep1 = lep1[:,:1]
    lep2 = lep2[:,:1]
    #trigId = trigIds[:,:1]
    
    # rebuild the pair with the 1st one only
    leps_pair = ak.concatenate([lep1, lep2], axis=1)

    return SelectionResult(
        aux = pair_selection_steps,
    ), leps_pair, trigIds, trigTypes
