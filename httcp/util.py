# coding: utf-8

"""
Collection of helpers
"""

from __future__ import annotations


import law
import order as od
from typing import Any, Optional
from functools import wraps
from collections import defaultdict, OrderedDict

from columnflow.util import maybe_import
from columnflow.columnar_util import ArrayFunction, deferred_column
from columnflow.selection import Selector, SelectionResult, selector
from columnflow.columnar_util import optional_column as optional

np = maybe_import("numpy")
ak = maybe_import("awkward")
coffea = maybe_import("coffea")
maybe_import("coffea.nanoevents.methods.nanoaod")

logger = law.logger.get_logger(__name__)


def get_objs_p4(collobj):
    return ak.zip(
        {
            "pt"  : collobj.pt,
            "eta" : collobj.eta,
            "phi" : collobj.phi,
            "mass": collobj.mass if "mass" in collobj.fields else ak.zeros_like(collobj.pt),
        },
        with_name="PtEtaPhiMLorentzVector",
        behavior=coffea.nanoevents.methods.vector.behavior, 
    )


def filter_by_triggers(lep_pair, mask):
    dummy = lep_pair[:,:0]
    any_mask = ak.any(mask, axis=1)
    out_pair = ak.where(any_mask, lep_pair, dummy)
    return out_pair

@deferred_column
def IF_RUN2(self, func: ArrayFunction)  -> Any | set[Any]:
    return self.get() if func.config_inst.campaign.x.year < 2022 else None

@deferred_column
def IF_RUN3(self, func: ArrayFunction)  -> Any | set[Any]:
    return self.get() if func.config_inst.campaign.x.year >= 2022 else None

@deferred_column
def IF_NANO_V9(self, func: ArrayFunction) -> Any | set[Any]:
    return self.get() if func.config_inst.campaign.x.version == 9 else None

@deferred_column
def IF_NANO_V11(self, func: ArrayFunction) -> Any | set[Any]:
    return self.get() if func.config_inst.campaign.x.version >= 10 else None

@deferred_column
def IF_DATASET_HAS_LHE_WEIGHTS(
    self: ArrayFunction.DeferredColumn,
    func: ArrayFunction,
) -> Any | set[Any]:
    if getattr(func, "dataset_inst", None) is None:
        return self.get()

    return None if func.dataset_inst.has_tag("no_lhe_weights") else self.get()

@deferred_column
def IF_DATASET_IS_TT(
    self: ArrayFunction.DeferredColumn,
    func: ArrayFunction,
) -> Any | set[Any]:
    if getattr(func, "dataset_inst", None) is None:
        return self.get()
    return None if not func.dataset_inst.has_tag("is_tt") else self.get()

@deferred_column
def IF_DATASET_IS_DY(
    self: ArrayFunction.DeferredColumn,
    func: ArrayFunction,
) -> Any | set[Any]:
    if getattr(func, "dataset_inst", None) is None:
        return self.get()
    return None if not func.dataset_inst.has_tag("is_dy") else self.get()

@deferred_column
def IF_DATASET_IS_DY_M50(
    self: ArrayFunction.DeferredColumn,
    func: ArrayFunction,
) -> Any | set[Any]:
    if getattr(func, "dataset_inst", None) is None:
        return self.get()
    return None if not func.dataset_inst.has_tag("is_dy_m50") else self.get()


@deferred_column
def IF_DATASET_IS_DY_TAUTAU(
    self: ArrayFunction.DeferredColumn,
    func: ArrayFunction,
) -> Any | set[Any]:
    if getattr(func, "dataset_inst", None) is None:
        return self.get()
    return None if not func.dataset_inst.has_tag("is_dy_tautau") else self.get()


@deferred_column
def IF_DATASET_IS_W(
    self: ArrayFunction.DeferredColumn,
    func: ArrayFunction,
) -> Any | set[Any]:
    if getattr(func, "dataset_inst", None) is None:
        return self.get()
    return None if not func.dataset_inst.has_tag("is_w") else self.get()

@deferred_column
def IF_DATASET_IS_SIGNAL(
    self: ArrayFunction.DeferredColumn,
    func: ArrayFunction,
) -> Any | set[Any]:
    if getattr(func, "dataset_inst", None) is None:
        return self.get()
    return None if not (func.dataset_inst.has_tag("is_ggf_signal")
                        | func.dataset_inst.has_tag("is_vh_signal")
                        | func.dataset_inst.has_tag("is_vbf_signal")) else self.get()

@deferred_column
def IF_ALLOW_STITCHING(
        self: ArrayFunction.DeferredColumn,
        func: ArrayFunction,
) -> Any | set[Any]:
    if getattr(func, "dataset_inst", None) is None:
        return self.get()
    allow_dy = func.dataset_inst.has_tag("is_dy_m50") & func.config_inst.x.allow_dy_stitching
    allow_w  = func.dataset_inst.has_tag("is_w") & func.config_inst.x.allow_w_stitching
    #return None if not (func.dataset_inst.has_tag("is_w") | func.dataset_inst.has_tag("is_dy")) else self.get()
    return None if not (allow_dy | allow_w) else self.get()

@deferred_column
def IF_GENMATCH_ON_FOR_SIGNAL(
    self: ArrayFunction.DeferredColumn,
    func: ArrayFunction,
) -> Any | set[Any]:
    if getattr(func, "dataset_inst", None) is None:
        return self.get()
    dogenmatch = False
    if "is_signal" in list(func.dataset_inst.aux.keys()):
        if func.config_inst.x.extra_tags.genmatch == True:
            dogenmatch = True
    return self.get() if dogenmatch else None



def transverse_mass(lepton: ak.Array, met: ak.Array) -> ak.Array:
    dphi_lep_met = lepton.delta_phi(met)
    mt = np.sqrt(2 * lepton.pt * met.pt * (1 - np.cos(dphi_lep_met)))
    return mt


#def trigger_object_matching(
#    vectors1: ak.Array,
#    vectors2: ak.Array,
#    threshold: float = 0.5,
#    axis: int = 2,
#) -> ak.Array:
#    """
#    Helper to check per object in *vectors1* if there is at least one object in *vectors2* that
#    leads to a delta R metric below *threshold*. The final reduction is applied over *axis* of the
#    resulting metric table containing the full combinatorics. When *return_all_matches* is *True*,
#    the matrix with all matching decisions is returned as well.
#    """
#    # delta_r for all combinations
#    dr = vectors1.metric_table(vectors2)
#    # check per element in vectors1 if there is at least one matching element in vectors2
#    any_match = ak.any(dr < threshold, axis=axis)##
#
#    return any_match



def trigger_object_matching(
    vectors1: ak.Array,
    vectors2: ak.Array,
    threshold: float = 0.5,
    axis: int = -1,
) -> ak.Array:
    """
    Helper to check per object in *vectors1* if there is at least one object in *vectors2* that
    leads to a delta R metric below *threshold*. The final reduction is applied over *axis* of the
    resulting metric table containing the full combinatorics. When *return_all_matches* is *True*,
    the matrix with all matching decisions is returned as well.
    """
    # delta_r for all combinations
    dr = vectors1.metric_table(vectors2)
    # check per element in vectors1 if there is at least one matching element in vectors2
    any_match = ak.fill_none(ak.any(dr < threshold, axis=axis), False)
    any_match = ak.enforce_type(ak.values_astype(any_match, "bool"), "var * var * bool")

    return any_match




def trigger_matching_extra(
        objs: ak.Array,
        legminpt: ak.Array,
        legmaxeta: ak.Array,
):
    # check pt and eta match
    pt  = objs.pt
    eta = np.abs(objs.eta)

    pt_brdcst, minpt_brdcst   = ak.broadcast_arrays(pt,  legminpt[:,None], depth_limit=-1)
    eta_brdcst, maxeta_brdcst = ak.broadcast_arrays(eta, legmaxeta[:,None], depth_limit=-1)        
    #match_pt = (obj_pt_brdcst - minpt_brdcst) >= 0.0 # 100000 * var * var * ?bool
    match_pt  = pt_brdcst >= minpt_brdcst   # 100000 * var * var * ?bool
    match_eta = eta_brdcst <= maxeta_brdcst # 100000 * var * var * ?bool

    any_match = match_pt & match_eta
    any_match = ak.enforce_type(ak.values_astype(any_match, "bool"), "var * var * bool")
    
    return any_match




def trigger_object_matching_deep(
        vectors1: ak.Array,
        vectors2: ak.Array,
        legminpt: ak.Array,
        legmaxeta: ak.Array,
        checkdr: bool = True,
        threshold: float = 0.5,
        axis: int = 1,
) -> ak.Array:
    """
    trigger object matching with the reconstructed objects
    This function is named as deep because it makes sure of both dr and pt matching
    vectors1 would always be the p4s of e/mu/tau
    vectors2 would always be the p4s of trigger objects
    dr is calculated using metric_table which brings all combinatorics.
    Finally, if any of the trigger objects matches the e/mu/tau, it will be True
    Next, the trigger leg minpt is provided as another argument.
    the reco object pt should be greater than the above
    Finally, both of these conditions must be satisfied for the selection of the reco object
    """

    # check leg min pt and leg max eta
    obj_pt = vectors1.pt
    obj_eta = np.abs(vectors1.eta)

    obj_pt_brdcst, minpt_brdcst = ak.broadcast_arrays(obj_pt, legminpt[:,None], depth_limit=-1)
    obj_eta_brdcst, maxeta_brdcst = ak.broadcast_arrays(obj_eta, legmaxeta[:,None], depth_limit=-1)        

    match_pt  = obj_pt_brdcst >= minpt_brdcst   # 100000 * var * var * ?bool
    match_eta = obj_eta_brdcst <= maxeta_brdcst # 100000 * var * var * ?bool
    
    match = match_pt & match_eta

    if checkdr:
        # delta_r for all combinations
        dr  = vectors1.metric_table(vectors2)  # 100000 * var * var * option[var * float32]
        match_dr = dr < threshold              # 100000 * var * var * option[var * bool]
        match_dr = ak.any(match_dr, axis=-1)   # 100000 * var * var * ?bool

        match = match & match_dr
    
    match = ak.fill_none(match, False) # 100000 * var * var * bool

    return match




def trigger_object_matching_jets_deep(
        vectors1: ak.Array,
        vectors2: ak.Array,
        legminpt: ak.Array,
        legmaxeta: ak.Array,
        checkdr: bool = True,
        threshold: float = 0.5,
        axis: int = 1,
) -> ak.Array:
    """
    trigger object matching with the reconstructed objects
    This function is named as deep because it makes sure of both dr and pt matching
    vectors1 would always be the p4s of e/mu/tau
    vectors2 would always be the p4s of trigger objects
    dr is calculated using metric_table which brings all combinatorics.
    Finally, if any of the trigger objects matches the e/mu/tau, it will be True
    Next, the trigger leg minpt is provided as another argument.
    the reco object pt should be greater than the above
    Finally, both of these conditions must be satisfied for the selection of the reco object
    """

    # check leg min pt and leg max eta
    obj_pt = vectors1.pt
    obj_eta = np.abs(vectors1.eta)

    obj_pt_brdcst, minpt_brdcst = ak.broadcast_arrays(obj_pt, legminpt[:,None], depth_limit=-1)
    obj_eta_brdcst, maxeta_brdcst = ak.broadcast_arrays(obj_eta, legmaxeta[:,None], depth_limit=-1)        

    match_pt  = obj_pt_brdcst >= minpt_brdcst   # 100000 * var * var * ?bool
    match_eta = obj_eta_brdcst <= maxeta_brdcst # 100000 * var * var * ?bool
    
    match = match_pt & match_eta

    match = ak.fill_none(ak.any(match, axis=-1), False)
    
    if checkdr:
        vectors2 = ak.drop_none(vectors2)
        # delta_r for all combinations
        dr  = vectors1.metric_table(vectors2)  # 100000 * var * var * option[var * float32]
        mindr = ak.fill_none(ak.min(ak.min(dr, axis=-1), axis=-1), 10.0)

        match_dr = mindr < threshold
        match = match & match_dr

        mindr = mindr[match]
        vectors1 = vectors1[match]

        mindr_sort_idx = ak.argsort(mindr, ascending=True)        
        vectors1 = vectors1[mindr_sort_idx]
        
    match = ak.fill_none(match, False) # 100000 * var * var * bool

    return match, vectors1





def get_dataset_lfns(
        dataset_inst: od.Dataset,
        shift_inst: od.Shift,
        dataset_key: str,
) -> list[str]:
    # destructure dataset_key into parts and create the lfn base directory
    lfn_base = law.wlcg.WLCGDirectoryTarget(
        dataset_key,
        fs=f"local",
    )
    # loop though files and interpret paths as lfns
    paths = [lfn_base.child(basename, type="f").path for basename in lfn_base.listdir(pattern="*.root")]

    return paths


def getGenTauDecayMode(prod: ak.Array):
    pids = prod.pdgId

    is_ele  = np.abs(pids) == 11
    is_muon = np.abs(pids) == 13
    #is_charged = ((np.abs(pids) == 211) | (np.abs(pids) == 321)) # | (np.abs(pids) == 323) | (np.abs(pids) == 325) | (np.abs(pids) == 327) | (np.abs(pids) == 329) )
    #is_neutral = ((pids == 111) | (pids == 311) | (pids == 130) | (pids == 310)) # | (pids == 313) | (pids == 315) | (pids == 317) | (pids == 319))
    is_charged = ((np.abs(pids) == 211) | (np.abs(pids) == 321) | (np.abs(pids) == 323) | (np.abs(pids) == 10321) | (np.abs(pids) == 10211))
    is_neutral = ((pids == 111) | (pids == 311) | (pids == 130) | (pids == 310))

    edecay = ak.sum(is_ele,  axis=-1) > 0
    mdecay = ak.sum(is_muon, axis=-1) > 0
    #hdecay = (ak.sum(is_charged, axis=-1) > 0) | (ak.sum(is_neutral, axis=-1) >= 0)
    hdecay = (
        (ak.sum(is_charged, axis=-1) > 0)
        & (ak.sum(is_neutral, axis=-1) >= 0)
    ) # either, one prong i.e. only one charged hadron or more than one prong where at least one charged hadron with zero or more neutral hadrons

    Nc = ak.sum(is_charged, axis=-1)
    Np = ak.sum(is_neutral, axis=-1)

    dm = ak.where(edecay, 
                  -1, 
                  ak.where(mdecay, 
                           -2, 
                           ak.where(hdecay, 
                                    (5 * (Nc - 1) + Np),
                                    -9)
                       )
              )
    
    dm = ak.fill_none(dm, -3) # just to make ?int64 to int64
    return dm



def enforce_hcand_type(hcand_pair_concat, field_type_dict):
    temp = {}
    for field, typename in field_type_dict.items():
        # 2022PreEE tt_dl --branch=8 has one single nan value for the tau in hcand. So, applying ak.nan_to_num for safety !!!
        # But, why nan?? Babushcha knows
        # event : 6784092 hcand_IPx: [[-0.000813, nan]]
        if field not in hcand_pair_concat.fields:
            logger.warning(f"{field} not found in enforce_hcand_type: 1 is used instead")
            field_in = ak.ones_like(hcand_pair_concat["pt"])
        else:
            field_in = hcand_pair_concat[field]

        temp[field] = ak.enforce_type(ak.values_astype(ak.nan_to_num(field_in, 0.0),typename),
                                      f"var * var * {typename}")
        #from IPython import embed; embed()
        
    hcand_array = ak.zip(temp)
    return hcand_array
    

@selector(
    uses={
        "process_id", optional("mc_weight")
    },
)
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
    """
    # get event masks
    event_mask = results.event

    # get a list of unique process ids present in the chunk
    unique_process_ids = np.unique(events.process_id)
    # increment plain counts
    n_evt_per_file = self.dataset_inst.n_events/self.dataset_inst.n_files
    stats["num_events"] = n_evt_per_file
    stats["num_events_selected"] += ak.sum(event_mask, axis=0)
    if self.dataset_inst.is_mc:
        stats[f"sum_mc_weight"] = n_evt_per_file
        stats.setdefault(f"sum_mc_weight_per_process", defaultdict(float))
        for p in unique_process_ids:
            stats[f"sum_mc_weight_per_process"][int(p)] = n_evt_per_file
        
    # create a map of entry names to (weight, mask) pairs that will be written to stats
    weight_map = OrderedDict()
    if self.dataset_inst.is_mc:
        # mc weight for selected events
        weight_map["mc_weight_selected"] = (events.mc_weight, event_mask)

    # get and store the sum of weights in the stats dictionary
    for name, (weights, mask) in weight_map.items():
        joinable_mask = True if mask is Ellipsis else mask

        # sum of different weights in weight_map for all processes
        stats[f"sum_{name}"] += ak.sum(weights[mask])
        # sums per process id
        stats.setdefault(f"sum_{name}_per_process", defaultdict(float))
        for p in unique_process_ids:
            stats[f"sum_{name}_per_process"][int(p)] += ak.sum(
                weights[(events.process_id == p) & joinable_mask],
            )

    return events, results


def get_trigger_id_map(triggers):
    tmap = {}
    for trigger in triggers:
        tmap[trigger.name] = trigger.id
    return tmap


def call_once_on_config(include_hash=False):
    """
    Parametrized decorator to ensure that function *func* is only called once for the config *config*
    """
    def outer(func):
        @wraps(func)
        def inner(config, *args, **kwargs):
            tag = f"{func.__name__}_called"
            if include_hash:
                tag += f"_{func.__hash__()}"

            if config.has_tag(tag):
                return

            config.add_tag(tag)
            return func(config, *args, **kwargs)
        return inner
    return outer


def setp4(name="LorentzVector", *args, verbose: Optional[bool] = False):
    if len(args) < 4:
        raise RuntimeError ("Need at least four components")

    if name == "PtEtaPhiMLorentzVector":
        if verbose:
            print(f" --- pT  : {args[0]}")
            print(f" --- eta : {args[1]}")
            print(f" --- phi : {args[2]}")
            print(f" --- mass: {args[3]}")
        return ak.zip(
            {
                "pt": args[0],
                "eta": args[1],
                "phi": args[2],
                "mass": args[3],
            },
            with_name=name,
            behavior=coffea.nanoevents.methods.vector.behavior
        )
    else:
        if verbose:
            print(f" --- px     : {args[0]}")
            print(f" --- py     : {args[1]}")
            print(f" --- pz     : {args[2]}")
            print(f" --- energy : {args[3]}")
        return ak.zip(
            {
                "x": args[0],
                "y": args[1],
                "z": args[2],
                "t": args[3],
            },
            with_name=name,
            behavior=coffea.nanoevents.methods.vector.behavior
        )
