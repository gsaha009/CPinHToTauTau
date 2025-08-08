# coding: utf-8

"""
Exemplary calibration methods.
"""
import law
import time
import functools
import itertools

from columnflow.calibration import Calibrator, calibrator
from columnflow.calibration.util import propagate_met
from columnflow.production.cms.seeds import deterministic_seeds
from columnflow.util import maybe_import, InsertableDict
from columnflow.columnar_util import set_ak_column, flat_np_view


np = maybe_import("numpy")
ak = maybe_import("awkward")

set_ak_column_f32 = functools.partial(set_ak_column, value_type=np.float32)

logger = law.logger.get_logger(__name__)

@calibrator(
    uses={"nTau"} | {f"Tau.{var}" for var in ["pt","eta","phi","mass","charge","decayMode", "decayModePNet", "genPartFlav"]} | {f"PuppiMET.{var}" for var in ["pt", "phi"]},
    produces={
        "Tau.pt_no_tes", "Tau.mass_no_tes",
        "Tau.pt_etau", "Tau.mass_etau",
        "Tau.pt_mutau", "Tau.mass_mutau",
        "Tau.pt_tautau", "Tau.mass_tautau",
        "PuppiMET.pt_no_tes", "PuppiMET.phi_no_tes",
        "PuppiMET.pt", "PuppiMET.phi",
    },
    mc_only=True,
)
def tau_energy_scale(self: Calibrator, events: ak.Array, **kwargs) -> ak.Array:
    start = time.time()
    # fail when running on data
    if self.dataset_inst.is_data:
        raise ValueError("attempt to apply tau energy corrections in data")
    
    # the correction tool only supports flat arrays, so convert inputs to flat np view first
    pt = flat_np_view(events.Tau.pt, axis=1)
    mass = flat_np_view(events.Tau.mass, axis=1)
    abseta = flat_np_view(abs(events.Tau.eta), axis=1)
    dm = flat_np_view(events.Tau.decayModePNet, axis=1)
    match = flat_np_view(events.Tau.genPartFlav, axis=1)

    tausum_before = events.Tau.sum(axis=1)
    
    #from IPython import embed; embed(); exit()
    
    syst = "nom" # TODO define this systematics inside config file
    #Get working points of the DeepTau tagger
    deep_tau_tagger    = self.config_inst.x.deep_tau_tagger
   
    #Create get energy scale correction for each tau
    tes_nom = np.ones_like(pt, dtype=np.float32)
    #Calculate tau ID scale factors for genuine taus
    # pt, eta, dm, genmatch, deep_tau_id, jet_wp, e_wp, syst

    arr_shape = ak.num(events.Tau.pt, axis=1)
    
    mask2prong = ((dm >= 0) & (dm != 5) & (dm != 6))

    #logger.critical("TES_Nom = 0 will be replaced by TES_Nom = 0.8 for safety == Bug in 2023PostBPix DM0 both WPs Tight for 20.0 < pT < 30.0 GeV")

    tes_args = lambda events, mask, deep_tau_tagger, syst: (pt[mask],
                                                            abseta[mask],
                                                            dm[mask],
                                                            match[mask],
                                                            deep_tau_tagger,
                                                            self.config_inst.x.deep_tau_info[deep_tau_tagger].vs_j["mutau"],
                                                            self.config_inst.x.deep_tau_info[deep_tau_tagger].vs_e["mutau"],
                                                            syst)
    tes_nom[mask2prong] = self.tes_corrector.evaluate(*tes_args(events, mask2prong, deep_tau_tagger, syst))
    tes_nom     = np.asarray(tes_nom)
    #tau_pt      = np.asarray(ak.flatten(events.Tau.pt))
    #tau_mass    = np.asarray(ak.flatten(events.Tau.mass))
    mask_mutau_tes_0 = (pt > 20.0) & (tes_nom == 0)
    n = np.sum(mask_mutau_tes_0)
    if n > 0:
        logger.warning(f"Calibration with mu-tau WPs : {n} taus with nom tes = 0 for tau pt > 20 GeV")
    tes_nom = np.where(tes_nom == 0.0, 0.8, tes_nom)
    
    events = set_ak_column_f32(events, "Tau.pt_no_tes", ak.unflatten(pt, arr_shape))
    events = set_ak_column_f32(events, "Tau.mass_no_tes", ak.unflatten(mass, arr_shape))

    events = set_ak_column_f32(events, "Tau.pt_mutau", ak.unflatten(pt * tes_nom, arr_shape))
    events = set_ak_column_f32(events, "Tau.mass_mutau", ak.unflatten(mass * tes_nom, arr_shape))
    # CAREFUL
    #events = set_ak_column_f32(events, "Tau.pt_mutau", ak.unflatten(pt, arr_shape))
    #events = set_ak_column_f32(events, "Tau.mass_mutau", ak.unflatten(mass, arr_shape))

    tes_nom = np.ones_like(pt, dtype=np.float32)
    tes_args = lambda events, mask, deep_tau_tagger, syst: (pt[mask],
                                                            abseta[mask],
                                                            dm[mask],
                                                            match[mask],
                                                            deep_tau_tagger,
                                                            self.config_inst.x.deep_tau_info[deep_tau_tagger].vs_j["etau"],
                                                            self.config_inst.x.deep_tau_info[deep_tau_tagger].vs_e["etau"],
                                                            syst)
    tes_nom[mask2prong] = self.tes_corrector.evaluate(*tes_args(events, mask2prong, deep_tau_tagger, syst))
    tes_nom     = np.asarray(tes_nom)
    #tau_pt      = np.asarray(ak.flatten(events.Tau.pt))
    #tau_mass    = np.asarray(ak.flatten(events.Tau.mass))
    mask_etau_tes_0 = (pt > 20.0) & (tes_nom == 0)
    n = np.sum(mask_etau_tes_0)
    if n > 0:
        logger.warning(f"Calibration with e-tau WPs : {n} taus with nom tes = 0 for tau pt > 20 GeV")
    tes_nom = np.where(tes_nom == 0.0, 0.8, tes_nom)

    events = set_ak_column_f32(events, "Tau.pt_etau", ak.unflatten(pt * tes_nom, arr_shape))
    events = set_ak_column_f32(events, "Tau.mass_etau", ak.unflatten(mass * tes_nom, arr_shape))
    # CAREFUL
    #events = set_ak_column_f32(events, "Tau.pt_etau", ak.unflatten(pt, arr_shape))
    #events = set_ak_column_f32(events, "Tau.mass_etau", ak.unflatten(mass, arr_shape))

    tes_nom = np.ones_like(pt, dtype=np.float32)
    tes_args = lambda events, mask, deep_tau_tagger, syst: (pt[mask],
                                                            abseta[mask],
                                                            dm[mask],
                                                            match[mask],
                                                            deep_tau_tagger,
                                                            self.config_inst.x.deep_tau_info[deep_tau_tagger].vs_j["tautau"],
                                                            self.config_inst.x.deep_tau_info[deep_tau_tagger].vs_e["tautau"],
                                                            syst)
    tes_nom[mask2prong] = self.tes_corrector.evaluate(*tes_args(events, mask2prong, deep_tau_tagger, syst))
    tes_nom     = np.asarray(tes_nom)
    #tau_pt      = np.asarray(ak.flatten(events.Tau.pt))
    #tau_mass    = np.asarray(ak.flatten(events.Tau.mass))
    mask_tautau_tes_0 = (pt > 20.0) & (tes_nom == 0)
    n = np.sum(mask_tautau_tes_0)
    if n > 0:
        logger.warning(f"Calibration with tau-tau WPs : {n} taus with nom tes = 0 for tau pt > 20 GeV")
    tes_nom = np.where(tes_nom == 0.0, 0.8, tes_nom)

    pt_tautau = ak.unflatten(pt * tes_nom, arr_shape)
    events = set_ak_column_f32(events, "Tau.pt_tautau", pt_tautau)
    events = set_ak_column_f32(events, "Tau.mass_tautau", ak.unflatten(mass * tes_nom, arr_shape))
    # CAREFUL
    #events = set_ak_column_f32(events, "Tau.pt_tautau", ak.unflatten(pt, arr_shape))
    #events = set_ak_column_f32(events, "Tau.mass_tautau", ak.unflatten(mass, arr_shape))

    # propagate to met
    events = set_ak_column_f32(events, "PuppiMET.pt_no_tes", events.PuppiMET.pt)
    events = set_ak_column_f32(events, "PuppiMET.phi_no_tes", events.PuppiMET.phi)
    
    #from IPython import embed; embed()
    # propagate changes to MET
    met_pt, met_phi = propagate_met(
        tausum_before.pt,
        tausum_before.phi,
        events.Tau.pt_tautau,
        events.Tau.phi,
        events.PuppiMET.pt,
        events.PuppiMET.phi,
    )
    events = set_ak_column_f32(events, "PuppiMET.pt", met_pt)
    events = set_ak_column_f32(events, "PuppiMET.phi", met_phi)

    
    stop = time.time()
    if self.config_inst.x.verbose.calibration.tau:
        print(f"tau energy correction takes : {round((stop - start)/60.0, 3)} min")
    
    return events


@tau_energy_scale.requires
def tau_energy_scale_requires(self: Calibrator, reqs: dict) -> None:
    if "external_files" in reqs:
        return
    
    from columnflow.tasks.external import BundleExternalFiles
    reqs["external_files"] = BundleExternalFiles.req(self.task)
    
@tau_energy_scale.setup
def tau_energy_scale_setup(
    self: Calibrator,
    reqs: dict,
    inputs: dict,
    reader_targets: InsertableDict,
) -> None:
    bundle = reqs["external_files"]
    import correctionlib
    correctionlib.highlevel.Correction.__call__ = correctionlib.highlevel.Correction.evaluate
    
    correction_set = correctionlib.CorrectionSet.from_string(
        #bundle.files.tau_sf.load(formatter="gzip").decode("utf-8"),
        bundle.files.tes_sf.load(formatter="gzip").decode("utf-8"),
    )
    #tagger_name = self.config_inst.x.deep_tau.tagger
    tagger_name = self.config_inst.x.deep_tau_tagger
    #self.tes_corrector = correction_set["tau_energy_scale"]
    self.tes_corrector = correction_set[f"tau_es_dm_{tagger_name}_{self.config_inst.campaign.x.year}_{self.config_inst.campaign.x.postfix}"]


@calibrator(
    uses={
        "channel_id", "hcand",
        "Tau.pt_no_tes", "Tau.mass_no_tes",
        "Tau.pt", "Tau.mass",
    },
    produces={"Tau.pt", "Tau.mass"},
    mc_only=True,
)
def insert_calibrated_taus(
        self: Calibrator,
        events: ak.Array,
        **kwargs
) -> ak.Array:
    """
      Thanks Jim: https://github.com/scikit-hep/awkward/discussions/3250
      A difficult task to do. 
      Only needed if one wants to run all channles simultaneously
    """
    start = time.time()

    etau_id   = self.config_inst.get_channel("etau").id
    mutau_id  = self.config_inst.get_channel("mutau").id
    tautau_id = self.config_inst.get_channel("tautau").id

    # set Muon, Electron and Tau from hcands
    dummy_idx = events.hcand.rawIdx[:,:0]
    dummy_pt = events.Tau.pt[:,:0]
    dummy_mass = events.Tau.mass[:,:0]

    raw_tau_pt = events.Tau.pt
    # check the channel_id and keep only the selected taus from hcand
    # here keeping the tau pt from hcand
    hcand_tau_pt = ak.where(events.channel_id == etau_id,
                            events.hcand.pt[:,1:2],
                            ak.where(events.channel_id == mutau_id,
                                     events.hcand.pt[:,1:2],
                                     ak.where(events.channel_id == tautau_id,
                                              events.hcand.pt,
                                              dummy_pt)
                                     )
                            )

    raw_tau_mass = events.Tau.mass
    # here keeping the tau mass from hcand
    hcand_tau_mass = ak.where(events.channel_id == etau_id,
                              events.hcand.mass[:,1:2],
                              ak.where(events.channel_id == mutau_id,
                                       events.hcand.mass[:,1:2],
                                       ak.where(events.channel_id == tautau_id,
                                                events.hcand.mass,
                                                dummy_mass)
                                       )
                              )

    # getting the raw tau idx
    # e.g. say 3 events
    # [ [0,1,2],  [0,1,2,3],  [0,1] ]
    raw_tau_idx = events.Tau.rawIdx
    # extract the tau indices from hcand using rawIdx
    # e.g. for the same 3 events
    # [ [0,2],  [3],  [] ]
    hcand_tau_idx = ak.where(events.channel_id == etau_id,
                             events.hcand.rawIdx[:,1:2],
                             ak.where(events.channel_id == mutau_id,
                                      events.hcand.rawIdx[:,1:2],
                                      ak.where(events.channel_id == tautau_id,
                                               events.hcand.rawIdx,
                                               dummy_idx)
                                      )
                             )

    # creating cartesian pairs
    # e.g.
    # [
    #   [ [(0,0),(0,2)], [(1,0),(1,2)], [(2,0),(2,2)] ],
    #   [ [(0,3)], [(1,3)], [(2,3)], [(3,3)] ],
    #   [ [], [] ]
    # ]
    comb_idx_raw_tau_hcand_tau = ak.cartesian([raw_tau_idx, hcand_tau_idx], nested=True, axis=1)
    # unzipping
    # 0 : [ [[0,0], [1,1], [2,2]],  [[0], [1], [2] ,[3]],  [[], []] ]
    # 1 : [ [[0,2], [0,2], [0,2]],  [[3], [3], [3], [3]],  [[], []] ]
    raw_idx_unzip, hcand_idx_unzip  = ak.unzip(comb_idx_raw_tau_hcand_tau)

    # get the mask to know the same positions
    # "0" == "1" :
    # [ [[T,F], [F,F], [F,T]],  [[F], [F], [F], [T]],  [[], []] ]
    equal_mask = raw_idx_unzip == hcand_idx_unzip
    # get this with ak.any
    # [ [T,F,T], [F,F,F,T], [F,F] ]
    # mind the shape of this boolean mask is the same as the shape of raw tau idx
    replace_mask = ak.any(equal_mask, axis=-1)

    # To reverse the pointers, we need to know: at what position in the axis=2 lists
    # do we have a match? We could use ak.argmax to find the position where x == y is
    # true (1), rather than false (0).
    equal_mask_to_int = ak.where(equal_mask, 1, 0)
    to_inplace_hcand_idx = ak.argmax(equal_mask_to_int, axis=-1)
    # place the hcand taus in the samne place of raw tau
    inplace_hcand_idx = to_inplace_hcand_idx.mask[replace_mask]

    pt_calib = ak.fill_none(hcand_tau_pt[inplace_hcand_idx], -1.0)
    mass_calib = ak.fill_none(hcand_tau_mass[inplace_hcand_idx], -1.0)

    tau_pt   = ak.where(replace_mask, pt_calib, raw_tau_pt)
    tau_mass = ak.where(replace_mask, mass_calib, raw_tau_mass)


    #from IPython import embed; embed()

    events = set_ak_column_f32(events, "Tau.pt", tau_pt)
    events = set_ak_column_f32(events, "Tau.mass", tau_mass)

    stop = time.time()
    if self.config_inst.x.verbose.calibration.tau:
        print(f"replacing raw taus by calibrated taus takes : {round((stop - start)/60.0, 3)} min")
    
    return events
