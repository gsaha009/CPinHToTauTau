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

np = maybe_import("numpy")
ak = maybe_import("awkward")
coffea = maybe_import("coffea")
maybe_import("coffea.nanoevents.methods.nanoaod")

from httcp.util import filter_by_triggers, get_objs_p4, trigger_matching_extra, trigger_object_matching_deep, trigger_object_matching_jets_deep



def match_trigobjs(
        leps_pair: ak.Array,
        trigger_results: SelectionResult,
        jets: ak.Array,
        **kwargs,
) -> tuple[ak.Array, ak.Array]:

    # extract the trigger names, types & others from trigger_results.x (aux)
    trigger_ids           = trigger_results.x.trigger_ids
    trigger_types         = trigger_results.x.trigger_types
    leg1_minpt            = trigger_results.x.leg1_minpt
    leg2_minpt            = trigger_results.x.leg2_minpt
    leg3_minpt            = trigger_results.x.leg3_minpt # NEW
    leg1_maxeta           = trigger_results.x.leg1_maxeta
    leg2_maxeta           = trigger_results.x.leg2_maxeta
    leg3_maxeta           = trigger_results.x.leg3_maxeta # NEW
    leg1_matched_trigobjs = trigger_results.x.leg1_matched_trigobjs
    leg2_matched_trigobjs = trigger_results.x.leg2_matched_trigobjs
    leg3_matched_trigobjs = trigger_results.x.leg3_matched_trigobjs # NEW

    has_tau_triggers     = ((trigger_types == "cross_tau_tau") | (trigger_types == "cross_tau_tau_jet"))

    old_leps_pair = leps_pair

    #from IPython import embed; embed()
    # ----- WARNING : TIME CONSUMING ------ #
    leps_pair = filter_by_triggers(leps_pair, has_tau_triggers) 

    taus1, taus2 = ak.unzip(leps_pair)

    #from IPython import embed; embed()

    #jets = filter_by_triggers(jets, trigger_types == "cross_tau_tau_jet") # NEW
    #jets = events.Jet[jet_indices]
    jets = filter_by_triggers(jets, has_tau_triggers) # NEW

    # IMPORTANT : APPLY JET PT > 60 GeV BEFORE MATCHING WITH JET LEG
    jets = jets[jets.pt > 60.0]

    
    # Event level masks
    # if events have tau
    has_tau_pairs = ak.fill_none(ak.num(taus1, axis=1) > 0, False)

    # events must be fired by tau triggers and there is ta inside
    mask_has_tau_triggers_and_has_tau_pairs = has_tau_triggers & has_tau_pairs

    tautau_trigger_ids            = trigger_ids[mask_has_tau_triggers_and_has_tau_pairs]    
    tautau_trigger_types          = trigger_types[mask_has_tau_triggers_and_has_tau_pairs]
    tautau_leg_1_minpt            = leg1_minpt[mask_has_tau_triggers_and_has_tau_pairs]
    tautau_leg_2_minpt            = leg2_minpt[mask_has_tau_triggers_and_has_tau_pairs]
    tautau_leg_3_minpt            = leg3_minpt[mask_has_tau_triggers_and_has_tau_pairs] # NEW
    tautau_leg_1_maxeta           = leg1_maxeta[mask_has_tau_triggers_and_has_tau_pairs]
    tautau_leg_2_maxeta           = leg2_maxeta[mask_has_tau_triggers_and_has_tau_pairs]
    tautau_leg_3_maxeta           = leg3_maxeta[mask_has_tau_triggers_and_has_tau_pairs] # NEW
    tautau_leg_1_matched_trigobjs = leg1_matched_trigobjs[mask_has_tau_triggers_and_has_tau_pairs]
    tautau_leg_2_matched_trigobjs = leg2_matched_trigobjs[mask_has_tau_triggers_and_has_tau_pairs]
    tautau_leg_3_matched_trigobjs = leg3_matched_trigobjs[mask_has_tau_triggers_and_has_tau_pairs] # NEW

    #from IPython import embed; embed()
    mask_has_tau_triggers_and_has_tau_pairs_evt_level = ak.fill_none(ak.any(mask_has_tau_triggers_and_has_tau_pairs, axis=1), False)
    jets_dummy = jets[:,:0]
    #jets = jets[mask_has_tau_triggers_and_has_tau_pairs]
    jets = ak.where(mask_has_tau_triggers_and_has_tau_pairs_evt_level, jets, jets_dummy)
    

    # tau1            : [ [    t11,        t12      ], [  t11  ] ]
    # tau2            : [ [    t21,        t22      ], [  t21  ] ]

    # trigger_ids     : [ [    1120,       1121     ], [  1120 ] ]
    # leg1/2_minpt    : [ [     40 ,        40      ], [   40  ] ]
    # leg1/2_trigobjs : [ [[o1,o2,o3], [o1,o2,o3,o4]], [[o1,o2]] ]

    # In [67]: taus1.pt[1044]
    # Out[67]: <Array [111, 111, 58.7] type='3 * float32'>

    # In [72]: ak.to_list(tautau_leg_1_matched_trigobjs.pt[1044])
    # Out[72]: 
    # [[110.203125, 53.0, 61.3828125, 210.625],
    #  [110.203125, 53.0, 61.3828125, 210.625, 35.0703125]]
    
    # In [70]: ak.to_list(ak.any(dr[1044] < 0.5, axis=-1))
    # Out[70]: [[True, True], [True, True], [True, True]]

    # same for tau1-legs2 and for tau2
    
    # final mask may look like this:
    # [[True, True], [True, True], [False, True]]

    # is tau1 matched to leg1?
    pass_tau1_leg1_triglevel = trigger_object_matching_deep(taus1,
                                                            tautau_leg_1_matched_trigobjs,
                                                            tautau_leg_1_minpt,
                                                            tautau_leg_1_maxeta,
                                                            True)
    # is tau1 matched to leg2?
    pass_tau1_leg2_triglevel = trigger_object_matching_deep(taus1,
                                                            tautau_leg_2_matched_trigobjs,
                                                            tautau_leg_2_minpt,
                                                            tautau_leg_2_maxeta,
                                                            True)

    # is tau2 matched to leg1?
    pass_tau2_leg1_triglevel = trigger_object_matching_deep(taus2,
                                                            tautau_leg_1_matched_trigobjs,
                                                            tautau_leg_1_minpt,
                                                            tautau_leg_1_maxeta,
                                                            True)
    # is tau2 matched to leg2?
    pass_tau2_leg2_triglevel = trigger_object_matching_deep(taus2,
                                                            tautau_leg_2_matched_trigobjs,
                                                            tautau_leg_2_minpt,
                                                            tautau_leg_2_maxeta,
                                                            True)

    # is any jet matched to leg3?
    pass_jets_leg3_triglevel, jets = trigger_object_matching_jets_deep(jets,
                                                                       tautau_leg_3_matched_trigobjs,
                                                                       tautau_leg_3_minpt,
                                                                       tautau_leg_3_maxeta,
                                                                       True)
    
    #from IPython import embed; embed()





    

    # tau1 to leg1 & tau2 to leg2 or, tau1 to leg2 & tau2 to leg1
    pass_taus_legs = (pass_tau1_leg1_triglevel & pass_tau2_leg2_triglevel) | (pass_tau1_leg2_triglevel & pass_tau2_leg1_triglevel)

    #pass_jet_leg_jet_level = ak.any(pass_jets_leg3_triglevel, axis=-1)
    #pass_jet_leg_jet_level_1 = pass_jet_leg_jet_level[:,:,0:1]
    #pass_jet_leg_jet_level_2 = pass_jet_leg_jet_level[:,:,1:2]
    #pass_jet_leg_any_1 = ak.any(pass_jet_leg_jet_level_1, axis=1)[:,None]
    #pass_jet_leg_any_2 = ak.any(pass_jet_leg_jet_level_2, axis=1)[:,None]
    #pass_jet_leg_tau_pair_level = ak.concatenate([pass_jet_leg_any_1, pass_jet_leg_any_2], axis=1)

    pass_jet_leg_jet_level = ak.any(pass_jets_leg3_triglevel, axis=-1)


    #dummy_true = ak.values_astype(ak.ones_like(tautau_trigger_ids), np.bool)
    #ak.where(tautau_trigger_ids == 15152, pass_jet_leg_trg_level, dummy_true)


    

    #pass_jet_leg_2 = ak.any(pass_jets_leg3_triglevel, axis=1)
    #temp = ak.values_astype(ak.ones_like(tautau_trigger_ids), np.bool)
    #dummy = temp[:,:0]
    #num_mask = ak.num(tautau_trigger_ids, axis=1) > 0
    #temp2 = ak.where(num_mask, temp[:,None], dummy)
    #pass_jet_leg_2 = ak.where(num_mask, pass_jet_leg[:,None], dummy)

    #from IPython import embed; embed()

    #pass_jet_leg2 = ak.where(num_mask, pass_jet_leg_2[:,None], dummy)
    #dummy_2 = tautau_trigger_ids[:,:0]
    #ids_2 = ak.where(num_mask, tautau_trigger_ids[:,None], dummy_2)
    
    
    #ditaujet_mask = (tautau_trigger_ids == 15152)
    #ak.where(ditaujet_mask, pass_jet_leg, temp2)
    
    
    #mask_has_tau_triggers_and_has_tau_pairs_evt_level = ak.fill_none(ak.any(mask_has_tau_triggers_and_has_tau_pairs, axis=1), False)
    trigobj_matched_mask_dummy = ak.from_regular((trigger_ids > 0)[:,:0][:,None])
    
    pass_taus_legs = ak.where(mask_has_tau_triggers_and_has_tau_pairs_evt_level, pass_taus_legs, trigobj_matched_mask_dummy)
    pass_taus_legs = ak.enforce_type(ak.values_astype(pass_taus_legs, "bool"), "var * var * bool") # 100000 * var * var * bool
    
    
    pass_taus = ak.fill_none(ak.any(pass_taus_legs, axis=-1), False)

    #tautau_trigger_ids_dummy = tautau_trigger_ids[:,:0]
    #tautau_trigger_ids = ak.where(mask_has_tau_triggers_and_has_tau_pairs_evt_level, tautau_trigger_ids, tautau_trigger_ids_dummy)
    
    
    tautau_trigger_ids_brdcst, _ = ak.broadcast_arrays(tautau_trigger_ids[:,None], pass_taus)
    #tautau_trigger_ids_brdcst = ak.where(mask_has_tau_triggers_and_has_tau_pairs_evt_level,
    #                                     tautau_trigger_ids_brdcst, 
    #                                     trigobj_matched_mask_dummy)

    tautau_trigger_ids_brdcst = tautau_trigger_ids_brdcst[pass_taus]
    # -------------> Order of triggers is important. Mind the hierarchy
    # this ids are the trig-obj match ids

    tautau_trigger_ids_brdcst_dummy = ak.from_regular(tautau_trigger_ids_brdcst[:,:0][:,None])
    #_tautau_trigger_ids_brdcst = ak.fill_none(ak.firsts(tautau_trigger_ids_brdcst, axis=1), 0)
    tautau_trigger_ids_brdcst = ak.where(ak.num(tautau_trigger_ids_brdcst) > 0, tautau_trigger_ids_brdcst, tautau_trigger_ids_brdcst_dummy)
    
    #ids = ak.fill_none(ak.firsts(tautau_trigger_ids_brdcst, axis=-1), -1)

    ids = ak.Array(ak.to_list(ak.firsts(tautau_trigger_ids_brdcst, axis=1))) # BAD Practice !!!
    
    #ids_dummy = ak.from_regular((trigger_ids > 0)[:,:0])
    #ids = ak.where(mask_has_tau_triggers_and_has_tau_pairs_evt_level, ids, ids_dummy)
    
    ids = ids[ak.fill_none(ak.firsts(pass_taus_legs, axis=1), False)]
    ids = ak.values_astype(ids, 'int64')

    new_taus1 = taus1[pass_taus]
    new_taus2 = taus2[pass_taus]
    #ids = ids[pass_taus]
    
    leps_pair = ak.zip([new_taus1, new_taus2])

    tautau_trigger_types_brdcst, _ = ak.broadcast_arrays(tautau_trigger_types[:,None], pass_taus)
    #tautau_trigger_ids_brdcst = ak.where(mask_has_tau_triggers_and_has_tau_pairs_evt_level,
    #                                     tautau_trigger_ids_brdcst, 
    #                                     trigobj_matched_mask_dummy)

    tautau_trigger_types_brdcst = tautau_trigger_types_brdcst[pass_taus]
    # -------------> Order of triggers is important. Mind the hierarchy
    # this ids are the trig-obj match ids


    #from IPython import embed; embed()

    tautau_trigger_types_brdcst_dummy = ak.from_regular(tautau_trigger_types_brdcst[:,:0][:,None])
    tautau_trigger_types_brdcst = ak.where(ak.num(tautau_trigger_types_brdcst) > 0, tautau_trigger_types_brdcst, tautau_trigger_types_brdcst_dummy)
    
    
    types = ak.Array(ak.to_list(ak.firsts(tautau_trigger_types_brdcst, axis=1))) # BAD Practice !!!
    types = types[ak.fill_none(ak.firsts(pass_taus_legs, axis=1), False)]
    
    #types = ak.fill_none(ak.firsts(tautau_trigger_types_brdcst, axis=-1), "")

    #from IPython import embed; embed()


    #pass_jet_leg_2 = pass_jet_leg[:,None]
    #dummy = pass_jet_leg_jet_level[:,:0]
    #pass_jet_leg_2 = ak.where(mask_has_tau_triggers_and_has_tau_pairs_evt_level, pass_jet_leg_2, dummy)

    #ids_mask_ditaujet = ak.sort(ids == 15152, ascending=False)
    #mask_ditaujet = ak.fill_none(ak.firsts(ids_mask_ditaujet, axis=1), False)
    #ids_mask_ditau = ids == 15151
    #mask_ditau = ak.fill_none(ak.firsts(ids_mask_ditaujet, axis=1), False)
    #mask_ditaujet_with_jetmatch = mask_ditaujet & pass_jet_leg

    #mask_all_matched = mask_ditau | mask_ditaujet_with_jetmatch

    #ids_ditau = 

    #from IPython import embed; embed()


    
    #-----ids = ak.enforce_type(ids, "var * int64")
    #-----types = ak.enforce_type(types, "var * string")
    """
    jets_idx = jets[pass_jet_leg_jet_level].rawIdx[:,:1]
    """
    jets_idx = jets.rawIdx[:,:1]
    #jets_idx = ak.local_index(jets.pt)[pass_jet_leg_jet_level][:,:1]
    
    return leps_pair, ids, types, jets_idx




def sort_pairs(dtrpairs: ak.Array)->ak.Array:
    sorted_idx = ak.argsort(dtrpairs["0"].rawDeepTau2018v2p5VSjet, ascending=False)
    dtrpairs = dtrpairs[sorted_idx]

    # if the deep tau val of tau-0 is the same for the first two pair
    where_same_iso_1 = ak.fill_none(
        ak.firsts(dtrpairs["0"].rawDeepTau2018v2p5VSjet[:,:1], axis=1) == ak.firsts(dtrpairs["0"].rawDeepTau2018v2p5VSjet[:,1:2], axis=1),
        False
    )

    # if so, sort the pairs according to the deep tau of the 2nd tau
    sorted_idx = ak.where(where_same_iso_1,
                          ak.argsort(dtrpairs["1"].rawDeepTau2018v2p5VSjet, ascending=False),
                          sorted_idx)
    dtrpairs = dtrpairs[sorted_idx]

    """
    # if the deep tau val of tau-1 is the same for the first two pair 
    where_same_iso_2 = ak.fill_none(
        ak.firsts(dtrpairs["1"].rawDeepTau2018v2p5VSjet[:,:1], axis=1) == ak.firsts(dtrpairs["1"].rawDeepTau2018v2p5VSjet[:,1:2], axis=1),
        False
    )
    where_same_iso_2 = ak.fill_none(where_same_iso_2, False)

    # sort them with the pt of the 1st tau
    sorted_idx = ak.where(where_same_iso_2,
                          ak.argsort(dtrpairs["0"].pt, ascending=False),
                          sorted_idx)
    dtrpairs = dtrpairs[sorted_idx]

    
    # check if the first two pairs have the second tau with same rawDeepTau2017v2p1VSjet
    where_same_pt_1 = ak.fill_none(
        ak.firsts(dtrpairs["0"].pt[:,:1], axis=1) == ak.firsts(dtrpairs["0"].pt[:,1:2], axis=1),
        False
    )
    
    # if so, sort the taus with their pt
    sorted_idx = ak.where(where_same_pt_1,
                          ak.argsort(dtrpairs["1"].pt, ascending=False),
                          sorted_idx)

    # finally, the pairs are sorted
    dtrpairs = dtrpairs[sorted_idx]

    #lep1 = ak.singletons(ak.firsts(dtrpairs["0"], axis=1))
    #lep2 = ak.singletons(ak.firsts(dtrpairs["1"], axis=1))

    #dtrpair    = ak.concatenate([lep1, lep2], axis=1) 
    """

    
    return dtrpairs



@selector(
    uses={
        optional("Tau.pt"),
        optional("Tau.pt_tautau"),
        optional("Tau.mass"),
        optional("Tau.mass_tautau"),
        "Tau.eta", "Tau.phi",
        "Tau.rawIdx", optional("Tau.genPartFlav"),
        "Tau.charge", "Tau.rawDeepTau2018v2p5VSjet",
        "Tau.idDeepTau2018v2p5VSjet", "Tau.idDeepTau2018v2p5VSe", "Tau.idDeepTau2018v2p5VSmu",
        "Jet.pt","Jet.eta","Jet.phi","Jet.mass",
    },
    exposed=False,
)
def tautau_selection(
        self: Selector,
        events: ak.Array,
        lep_indices: ak.Array,
        trigger_results: SelectionResult,
        jet_indices : ak.Array,
        **kwargs,
) -> tuple[SelectionResult, ak.Array, ak.Array]:

    taus            = events.Tau[lep_indices]
    # Extra channel specific selections on tau
    tau_tagger      = self.config_inst.x.deep_tau_tagger
    tau_tagger_wps  = self.config_inst.x.deep_tau_info[tau_tagger].wp
    vs_e_wp         = self.config_inst.x.deep_tau_info[tau_tagger].vs_e["tautau"]
    vs_mu_wp        = self.config_inst.x.deep_tau_info[tau_tagger].vs_m["tautau"]
    vs_jet_wp       = self.config_inst.x.deep_tau_info[tau_tagger].vs_j["tautau"]

    # rename the {channel}_pt/mass to pt/mass
    if self.dataset_inst.is_mc:
        taus = ak.without_field(taus, "pt")
        taus = ak.with_field(taus, taus.pt_tautau, "pt")
        taus = ak.without_field(taus, "mass")
        taus = ak.with_field(taus, taus.mass_tautau, "mass")
    
    is_good_tau     = (
        (taus.pt > 20.0)
        #(taus.idDeepTau2018v2p5VSjet   >= tau_tagger_wps.vs_j[vs_jet_wp])
        & (taus.idDeepTau2018v2p5VSe   >= tau_tagger_wps.vs_e[vs_e_wp])
        & (taus.idDeepTau2018v2p5VSmu  >= tau_tagger_wps.vs_m[vs_mu_wp])
    )

    taus = taus[is_good_tau]

    
    # Sorting leps [Tau] by deeptau [descending]
    taus_sort_idx = ak.argsort(taus.rawDeepTau2018v2p5VSjet, axis=-1, ascending=False)
    taus = taus[taus_sort_idx]

    leps_pair  = ak.combinations(taus, 2, axis=1)    
    lep1, lep2 = ak.unzip(leps_pair)

    #from IPython import embed; embed()
    
    preselection = {
        #"tautau_tau1_iso"      : (lep1.idDeepTau2018v2p5VSjet >= tau_tagger_wps.vs_j[vs_jet_wp]),
        "tautau_is_pt_35"      : (lep1.pt > 35.0) & (lep2.pt > 35.0), # just changed 40.0 to 35.0 (19.12.2024)
        "tautau_is_eta_2p1"    : (np.abs(lep1.eta) < 2.1) & (np.abs(lep2.eta) < 2.1),
        #"tautau_is_os"         : (lep1.charge * lep2.charge) < 0,
        "tautau_dr_0p5"        : (1*lep1).delta_r(1*lep2) > 0.5,  #deltaR(lep1, lep2) > 0.5,
        "tautau_invmass_40"    : (1*lep1 + 1*lep2).mass > 40.0, # invariant_mass(lep1, lep2) > 40
    }

    
    good_pair_mask = lep1.rawIdx >= 0
    pair_selection_steps = {}
    category_selections = {}
    pair_selection_steps["tautau_starts_with"] = good_pair_mask
    for cut in preselection.keys():
        good_pair_mask = good_pair_mask & preselection[cut]
        pair_selection_steps[cut] = good_pair_mask

    good_pair_mask = ak.fill_none(good_pair_mask, False)
    leps_pair = leps_pair[good_pair_mask]

    # check nPairs
    npair = ak.num(leps_pair["0"], axis=1)
    pair_selection_steps["tautau_before_trigger_matching"] = leps_pair["0"].pt >= 0.0
    
    # sort the pairs if many
    leps_pair = ak.where(npair > 1, sort_pairs(leps_pair), leps_pair)

    
    # match trigger objects for all pairs
    leps_pair, trigIds, trigTypes, jet_idx = match_trigobjs(leps_pair, trigger_results, events.Jet[jet_indices])

    
    pair_selection_steps["tautau_after_trigger_matching"] = leps_pair["0"].pt >= 0.0
    
    lep1, lep2 = ak.unzip(leps_pair)

    # take the 1st pair and 1st trigger id
    lep1 = lep1[:,:1]
    lep2 = lep2[:,:1]
    #trigId = trigIds[:,:1]
    #trigTypes = trigTypes[:,:1]


    #trigIds_ = ak.sort(trigIds, ascending=False)[:,:1]
    #is_ditaujet = trigIds_ == 15152 #         [ [], [True], [], [False] ]
                                    # jetidx  [ [], [1],    [], [0]     ]

    #has_jet_ = ak.num(jet_idx, axis=1) > 0
    ditau_ids = trigIds[trigIds == 15151]
    ditaujet_ids =	trigIds[trigIds == 15152]

    has_jet = jet_idx >= 0
    a = ak.concatenate([(ditaujet_ids == 15152), has_jet], axis=1)
    b = ak.sum(a, axis=1) == 2
    ditaujet_ids_matched = ak.where(b, ditaujet_ids, ditaujet_ids[:,:0])

    trig_ids_matched = ak.concatenate([ditau_ids,ditaujet_ids_matched], axis=1)
    trig_types_matched = trigTypes[trig_ids_matched > 0]
    
    
    # rebuild the pair with the 1st one only
    leps_pair = ak.concatenate([lep1, lep2], axis=1)

    sort_idx = ak.argsort(leps_pair.pt, ascending=False)
    leps_pair = leps_pair[sort_idx]

    leps_pair_dummy = leps_pair[:,:0]
    leps_pair_matched = ak.where(ak.num(trig_ids_matched) > 0, leps_pair, leps_pair_dummy)

    #from IPython import embed; embed()
    

    #return SelectionResult(
    #    aux = pair_selection_steps,
    #), leps_pair, trigIds, trigTypes

    return SelectionResult(
        aux = pair_selection_steps,
    ), leps_pair_matched, trig_ids_matched, trig_types_matched, jet_idx
