import law
import re
import functools
from columnflow.production import Producer, producer
from columnflow.util import maybe_import, safe_div, InsertableDict, dev_sandbox
from columnflow.columnar_util import set_ak_column, has_ak_column, EMPTY_FLOAT, Route, flat_np_view, optional_column as optional
from columnflow.production.util import attach_coffea_behavior

from httcp.util import get_trigger_id_map

ak     = maybe_import("awkward")
np     = maybe_import("numpy")
pd     = maybe_import("pandas")
coffea = maybe_import("coffea")
warn   = maybe_import("warnings")

# helper
set_ak_column_f32 = functools.partial(set_ak_column, value_type=np.float32)

logger = law.logger.get_logger(__name__)

# ------------------------------------------------- #
# Assign MC weight [gen Weight / LHEWeight]
# ------------------------------------------------- #
@producer(
    uses={"genWeight", optional("LHEWeight.originalXWGTUP")},
    produces={"mc_weight"},
    mc_only=True,
)
def scale_mc_weight(self: Producer, events: ak.Array, **kwargs) -> ak.Array:
    """
    Reads the genWeight and LHEWeight columns and makes a decision about which one to save. This
    should have been configured centrally [1] and stored in genWeight, but there are some samples
    where this failed.

    Strategy:

      1. Use LHEWeight.originalXWGTUP when it exists and genWeight is always 1.
      2. In all other cases, use genWeight.

    [1] https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookNanoAOD?rev=99#Weigths
    """
    # determine the mc_weight
    mc_weight = np.sign(events.genWeight)
    if has_ak_column(events, "LHEWeight.originalXWGTUP") and ak.all(events.genWeight == 1.0):
        mc_weight = np.sign(events.LHEWeight.originalXWGTUP)

    # store the column
    events = set_ak_column(events, "mc_weight", mc_weight, value_type=np.float32)

    return events


# ------------------------------------------------- #
# Calculate Zpt weights
# ------------------------------------------------- #

@producer(
    uses={
        "GenZ.pt", "GenZ.mass",
    },
    produces={
        "zpt_reweight", "zpt_reweight_up", "zpt_reweight_down",
    },
    mc_only=True,
)
def zpt_reweight(
        self: Producer,
        events: ak.Array,
        **kwargs,
) :
    # if no gen particle found, all fields of GenZ will be zero
    # Those events will have nominal zpt rewt = 1.0
    # the same will be applied also for the evnts outside the range

    is_outside_range = (
        ((events.GenZ.pt == 0.0) & (events.GenZ.mass == 0.0))
        | ((events.GenZ.pt >= 600.0) | (events.GenZ.mass >= 1000.0))
    )

    # for safety
    zm  = ak.where(events.GenZ.mass > 1000.0, 999.99, events.GenZ.mass)
    zpt = ak.where(events.GenZ.pt > 600.0, 599.99, events.GenZ.pt)

    sf_nom = ak.where(is_outside_range,
                      1.0,
                      self.zpt_corrector.evaluate(zm,zpt))

    events = set_ak_column(events, "zpt_reweight", sf_nom, value_type=np.float32)

    sf_up   = ak.where(is_outside_range, 1.0, 1.67 * sf_nom)
    events = set_ak_column(events, "zpt_reweight_up", sf_up, value_type=np.float32)
    sf_down = ak.where(is_outside_range, 1.0, 0.33 * sf_nom)
    events = set_ak_column(events, "zpt_reweight_down", sf_down, value_type=np.float32)    
    
    return events


@zpt_reweight.requires
def zpt_reweight_requires(self: Producer, reqs: dict) -> None:
    if "external_files" in reqs:
        return
    
    from columnflow.tasks.external import BundleExternalFiles
    reqs["external_files"] = BundleExternalFiles.req(self.task)

@zpt_reweight.setup
def zpt_reweight_setup(
    self: Producer,
    reqs: dict,
    inputs: dict,
    reader_targets: InsertableDict,
) -> None:
    bundle = reqs["external_files"]
    import correctionlib
    correctionlib.highlevel.Correction.__call__ = correctionlib.highlevel.Correction.evaluate
    
    correction_set = correctionlib.CorrectionSet.from_string(
        bundle.files.zpt_rewt_v1_sf.load(formatter="gzip").decode("utf-8"),
    )
    self.zpt_corrector    = correction_set["zptreweight"]



# ------------------------------------------------- #
# Calculate Zpt weights (Dennis)
# ------------------------------------------------- #

@producer(
    uses={
        "GenZ.pt", "GenZ.mass",
    },
    produces={
        "zpt_reweight", "zpt_reweight_up", "zpt_reweight_down",
    },
    mc_only=True,
)
def zpt_reweight_v2(
        self: Producer,
        events: ak.Array,
        **kwargs,
) :
    # if no gen particle found, all fields of GenZ will be zero
    # Those events will have nominal zpt rewt = 1.0
    # the same will be applied also for the evnts outside the range

    is_outside_range = (
        (events.GenZ.pt == 0.0) & (events.GenZ.mass == 0.0)
    )

    zpt = events.GenZ.pt

    dataset_type = "LO" if "madgraph" in self.dataset_inst.name else "NLO"
    #era = f"{self.config_inst.campaign.x.year}{self.config_inst.campaign.x.postfix}_{dataset_type}"
    era = f"{self.config_inst.campaign.x.year}{self.config_inst.campaign.x.postfix}"
    
    for syst in ["nom", "up1", "down1"]:
        tag = re.match(r'([a-zA-Z]+)\d*', syst).group(1)
        name = "zpt_reweight" if syst=="nom" else f"zpt_reweight_{tag}"
        #from IPython import embed; embed()
        sf = ak.where(is_outside_range,
                      1.0,
                      self.zpt_corrector.evaluate(era,
                                                  dataset_type,
                                                  zpt,
                                                  syst)
                      )
        
        events = set_ak_column(events, name, sf, value_type=np.float32)
    
    return events


@zpt_reweight_v2.requires
def zpt_reweight_v2_requires(self: Producer, reqs: dict) -> None:
    if "external_files" in reqs:
        return
    
    from columnflow.tasks.external import BundleExternalFiles
    reqs["external_files"] = BundleExternalFiles.req(self.task)

@zpt_reweight_v2.setup
def zpt_reweight_v2_setup(
    self: Producer,
    reqs: dict,
    inputs: dict,
    reader_targets: InsertableDict,
) -> None:
    bundle = reqs["external_files"]
    import correctionlib
    correctionlib.highlevel.Correction.__call__ = correctionlib.highlevel.Correction.evaluate
    
    correction_set = correctionlib.CorrectionSet.from_string(
        bundle.files.zpt_rewt_v2_sf.load(formatter="gzip").decode("utf-8"),
    )
    self.zpt_corrector    = correction_set["DY_pTll_reweighting"]






# ------------------------------------------------- #
# Calculate FF weights (Dummy)
# ------------------------------------------------- #

@producer(
    uses={
        "channel_id",
        "category_ids",
        "hcand.pt", "hcand.decayMode", "n_jet",
        "met_var_qcd_h1",
    },
    produces={
        "ff_weight",
        #"ff_cls_corr_weight",
    },
    mc_only=False,
)
def ff_weight(
        self: Producer,
        events: ak.Array,
        **kwargs,
) :
    # Leading candidate
    hcand1 = events.hcand[:,0] 

    is_A = (
        (events.channel_id == 4)
        & ~events.is_os
        & events.is_real_1
        & events.is_iso_2
        & events.is_iso_1
    )
    is_B = (
        (events.channel_id == 4)
        & ~events.is_os
        & events.is_real_1
        & events.is_iso_2
        & ~events.is_iso_1
    )
    is_C0 = (
        (events.channel_id == 4)
        & events.is_os
        & events.is_real_1
        & ~events.is_iso_2
        & ~events.is_iso_1
    )
    is_C = (
        (events.channel_id == 4)
        & events.is_os
        & events.is_real_1
        & events.is_iso_2
        & ~events.is_iso_1
    )

    is_A_category = ak.to_numpy(is_A)
    is_B_category = ak.to_numpy(is_B)
    is_C0_category = ak.to_numpy(is_C0)
    is_C_category = ak.to_numpy(is_C)

    
    # Get the pt of the leading candidate
    pt1 = flat_np_view(hcand1.pt[:,None])
    metvarqcdh1 = flat_np_view(events.met_var_qcd_h1[:,None])
    
    njet = ak.where(events.n_jet > 2, 2, events.n_jet)
    njet = ak.to_numpy(njet)
    dm = ak.where(hcand1.decayMode < 0, 0, hcand1.decayMode)
    dm = ak.where(dm == 11, 10, dm)
    dm = ak.to_numpy(dm)
        
    fake_factors_nom = self.ff_corrector.evaluate(
        pt1,
        dm,
        njet,
        "nom",
    )
    ##from IPython import embed; embed()
    #closure_corr_nom = self.cls_corrector.evaluate(
    #    metvarqcdh1,
    #    #dm,
    #    njet,
    #    "nom",
    #)
    
    # Apply the fake factor only for the C category
    ff_nom  = np.where((is_C_category | is_B_category), fake_factors_nom, 1.0)
    #cls_nom = np.ones_like(ff_nom)
    #cls_nom = np.where((is_C_category | is_B_category), closure_corr_nom, cls_nom)
    ##ff_nom = np.where(is_C0_category, fake_0_factors_nom, ff_nom)

    ##ext_corr_nom = np.where((is_C_category | is_C0_category), ext_corr_nom, 1.0)
    
    ##closure_nom = np.where(is_C_category, closure_nom, 1.0)

    # Add the column to the events
    events = set_ak_column(events, "ff_weight", ff_nom, value_type=np.float32)
    #events = set_ak_column(events, "ff_cls_corr_weight", cls_nom, value_type=np.float32)
    ##events = set_ak_column(events, "ff_ext_corr_weight", ext_corr_nom, value_type=np.float32)
    
    return events


@ff_weight.requires
def ff_weight_requires(self: Producer, reqs: dict) -> None:
    if "external_files" in reqs:
        return
    
    from columnflow.tasks.external import BundleExternalFiles
    reqs["external_files"] = BundleExternalFiles.req(self.task)

@ff_weight.setup
def ff_weight_setup(
    self: Producer,
    reqs: dict,
    inputs: dict,
    reader_targets: InsertableDict,
) -> None:
    bundle = reqs["external_files"]
    import correctionlib
    correctionlib.highlevel.Correction.__call__ = correctionlib.highlevel.Correction.evaluate

    correction_set = correctionlib.CorrectionSet.from_file(
        bundle.files.tautau_ff.path,
    ) 
    self.ff_corrector = correction_set["fake_factors_fit"]

    closure_correction_set = correctionlib.CorrectionSet.from_file(
        bundle.files.tautau_ff_closure.path,
    ) 
    self.cls_corrector = closure_correction_set["closure_corrections_fit"]
    
    ## extrapolation correction on FF
    #ext_correction_set = correctionlib.CorrectionSet.from_file(
    #    bundle.files.tautau_ext_corr.path,
    #)
    #self.ext_corrector = ext_correction_set["extrapolation_correction"]

    
    #correction_set_0 = correctionlib.CorrectionSet.from_file(
    #    bundle.files.tautau_ff0.path,
    #)
    #self.ff0_corrector = correction_set["fake_factors_fit"]

    ## extrapolation correction on FF
    #ext_correction_set = correctionlib.CorrectionSet.from_file(
    #    bundle.files.tautau_ext_corr.path,
    #)
    #self.ext_corrector = ext_correction_set["extrapolation_correction"]
    


# ------------------------------------------------- #
# Calculate MET Recoil Corrections
# ------------------------------------------------- #

@producer(
    uses={
        "metRecoilJet.pt",
        "GenZ.pt", "GenZ.phi",
        "GenZvis.pt", "GenZvis.phi",
        "PuppiMET.pt", "PuppiMET.phi",
    },
    produces={
        "PuppiMET.pt", "PuppiMET.phi",
    },
    mc_only=True,
)
def met_recoil_corr(
        self: Producer,
        events: ak.Array,
        **kwargs,
) :
    """
      Ref:
        - https://gitlab.cern.ch/dwinterb/HiggsDNA/-/blob/master/higgs_dna/tools/ditau/recoil_corrector.py#L62
    """
    n_jets = ak.to_numpy(ak.num(events.metRecoilJet.pt, axis=1)).astype(float)
    
    met_pt = events.PuppiMET.pt
    met_phi = events.PuppiMET.phi
    gen_boson_pt = events.GenZ.pt
    gen_boson_phi = events.GenZ.phi
    gen_boson_visible_pt = events.GenZvis.pt
    gen_boson_visible_phi = events.GenZvis.phi

    met_x = met_pt * np.cos(met_phi)
    met_y = met_pt * np.sin(met_phi)
    gen_boson_x = gen_boson_pt * np.cos(gen_boson_phi)
    gen_boson_y = gen_boson_pt * np.sin(gen_boson_phi)
    gen_boson_visible_x = gen_boson_visible_pt * np.cos(gen_boson_visible_phi)
    gen_boson_visible_y = gen_boson_visible_pt * np.sin(gen_boson_visible_phi)

    U_x = met_x + gen_boson_visible_x - gen_boson_x
    U_y = met_y + gen_boson_visible_y - gen_boson_y

    U_pt = np.sqrt(U_x**2 + U_y**2)
    U_phi = np.arctan2(U_y, U_x)

    Upara = U_pt * np.cos(U_phi - gen_boson_phi)
    Uperp = U_pt * np.sin(U_phi - gen_boson_phi)


    dataset_type = "LO" if "madgraph" in self.dataset_inst.name else "NLO"
    era = f"{self.config_inst.campaign.x.year}{self.config_inst.campaign.x.postfix}"


    Upara_corr = self.met_recoil_corrector.evaluate(era, dataset_type, n_jets, gen_boson_pt, "Upara", Upara)
    Uperp_corr = self.met_recoil_corrector.evaluate(era, dataset_type, n_jets, gen_boson_pt, "Uperp", Uperp)

    U_pt_corr = np.sqrt(Upara_corr**2 + Uperp_corr**2)
    U_phi_corr = np.arctan2(Uperp_corr, Upara_corr) + gen_boson_phi
    U_x_corr = U_pt_corr * np.cos(U_phi_corr)
    U_y_corr = U_pt_corr * np.sin(U_phi_corr)
    met_x_corr = U_x_corr - gen_boson_visible_x + gen_boson_x
    met_y_corr = U_y_corr - gen_boson_visible_y + gen_boson_y

    met_pt_corr = np.sqrt(met_x_corr**2 + met_y_corr**2)
    met_phi_corr = np.arctan2(met_y_corr, met_x_corr)    
    
    # Add the column to the events
    #events = set_ak_column(events, "PuppiMET.pt_pre_recoil", met_pt)
    #events = set_ak_column(events, "PuppiMET.phi_pre_recoil", met_phi)

    events = set_ak_column(events, "PuppiMET.pt", met_pt_corr)
    events = set_ak_column(events, "PuppiMET.phi", met_phi_corr)

    #from IPython import embed; embed()
    
    return events


@met_recoil_corr.requires
def met_recoil_corr_requires(self: Producer, reqs: dict) -> None:
    if "external_files" in reqs:
        return
    
    from columnflow.tasks.external import BundleExternalFiles
    reqs["external_files"] = BundleExternalFiles.req(self.task)

@met_recoil_corr.setup
def met_recoil_corr_setup(
    self: Producer,
    reqs: dict,
    inputs: dict,
    reader_targets: InsertableDict,
) -> None:
    bundle = reqs["external_files"]
    import correctionlib
    correctionlib.highlevel.Correction.__call__ = correctionlib.highlevel.Correction.evaluate

    #from IPython import embed; embed()
    correction_set = correctionlib.CorrectionSet.from_file(
        bundle.files.met_recoil.path,
    ) 
    self.met_recoil_corrector = correction_set["Recoil_correction_QuantileMapHist"]




@producer(
    uses={
        "GenTop_pt",
    },
    produces={
        "top_pt_weight", "top_pt_weight_up", "top_pt_weight_down",
    },
    get_top_pt_config=(lambda self: self.config_inst.x.top_pt_reweighting_params),
    mc_only=True,
)
def top_pt_weight(self: Producer, events: ak.Array, **kwargs) -> ak.Array:
    """
    Compute SF to be used for top pt reweighting.

    The *GenPartonTop.pt* column can be produced with the :py:class:`gen_parton_top` Producer.

    The SF should *only be applied in ttbar MC* as an event weight and is computed
    based on the gen-level top quark transverse momenta.

    The function is skipped when the dataset is data or when it does not have the tag *is_ttbar*.

    The top pt reweighting parameters should be given as an auxiliary entry in the config:

    .. code-block:: python

        cfg.x.top_pt_reweighting_params = {
            "a": 0.0615,
            "a_up": 0.0615 * 1.5,
            "a_down": 0.0615 * 0.5,
            "b": -0.0005,
            "b_up": -0.0005 * 1.5,
            "b_down": -0.0005 * 0.5,
        }

    *get_top_pt_config* can be adapted in a subclass in case it is stored differently in the config.

    :param events: awkward array containing events to process
    """

    # get SF function parameters from config
    params = self.get_top_pt_config()

    # check the number of gen tops
    if ak.any(ak.num(events.GenTop_pt, axis=1) != 2):
        logger.warning("There are events with != 2 GenPartonTops. This producer should only run for ttbar")

    # clamp top pT < 500 GeV
    pt_clamped = ak.where(events.GenTop_pt > 500.0, 500.0, events.GenTop_pt)
    for variation in ("", "_up", "_down"):
        # evaluate SF function
        #sf = np.exp(params[f"a{variation}"] + params[f"b{variation}"] * pt_clamped)

        sf_13p0 = params[f"a{variation}"] * np.exp(params[f"b{variation}"] * pt_clamped) + params[f"c{variation}"] * pt_clamped + params[f"d{variation}"]
        sf_ext_13p6 = params[f"e{variation}"] + params[f"f{variation}"] * pt_clamped

        sf = sf_13p0 * sf_ext_13p6
        
        # compute weight from SF product for top and anti-top
        weight = np.sqrt(np.prod(sf, axis=1))

        # write out weights
        events = set_ak_column(events, f"top_pt_weight{variation}", ak.fill_none(weight, 1.0))

    return events


# ------------------------------------------------- #
# Calculate Classifier Score
# ------------------------------------------------- #

@producer(
    uses={
        "channel_id",
        "hcand.*", "PuppiMET.*",
        "dphi_met_h1", "dphi_met_h2",
        "hcand_dr", "hcand_invm", "hcand_dphi",
        "pt_tt", "pt_vis",
        "Jet.*",
        "bJet.pt",
        "hcand_invm_fastMTT",
    },
    produces={
        "classifier_score",
    },
    #sandbox=dev_sandbox("bash::$HTTCP_BASE/sandboxes/venv_columnar_xgb.sh"),
    mc_only=False,
)
def classify_events(
        self: Producer,
        events: ak.Array,
        **kwargs,
) :
    is_even = events.event % 2 == 0
    is_odd  = ~is_even

    is_etau = events.channel_id == 1
    is_mutau = events.channel_id == 2
    is_tautau = events.channel_id == 4

    is_etau_even = ak.to_numpy(is_etau & is_even)
    is_etau_odd  = ak.to_numpy(is_etau & is_odd)
    is_mutau_even = ak.to_numpy(is_mutau & is_even)
    is_mutau_odd  = ak.to_numpy(is_mutau & is_odd)
    is_tautau_even = ak.to_numpy(is_tautau & is_even)
    is_tautau_odd  = ak.to_numpy(is_tautau & is_odd)


    # add more Jet features
    default_val = -9999.0
    njets = ak.num(events.Jet.pt, axis=1)
    nbjets = ak.num(events.bJet.pt, axis=1)
    
    jets = ak.with_name(events.Jet, "PtEtaPhiMLorentzVector")
    jet1 = jets[:,0:1]
    jet2 = jets[:,1:2]
    # dummyp4
    dummyp4 = jet1[:,:0]
    # make two jets array the same structure
    jet1 = ak.where(njets >= 2, jet1, dummyp4)
    jet2 = ak.where(njets >= 2, jet2, dummyp4)

    dijet = jet1 + jet2
    
    jpt_1 = ak.from_regular(ak.fill_none(ak.firsts(jet1.pt, axis=1), default_val)[:,None])
    jpt_2 = ak.from_regular(ak.fill_none(ak.firsts(jet1.pt, axis=1), default_val)[:,None])
    jeta_1 = ak.from_regular(ak.fill_none(ak.firsts(jet1.eta, axis=1), default_val)[:,None])
    jeta_2 = ak.from_regular(ak.fill_none(ak.firsts(jet1.eta, axis=1), default_val)[:,None])
    # mjj
    mjj = (jet1 + jet2).mass
    mjj = ak.fill_none(ak.firsts(mjj, axis=1), default_val)
    mjj = ak.from_regular(ak.where(njets >= 2, mjj, default_val)[:,None])
    # jdeta
    jdeta = np.abs(jet1.eta - jet2.eta)
    jdeta = ak.fill_none(ak.firsts(jdeta, axis=1), default_val)
    jdeta = ak.from_regular(ak.where(njets >= 2, jdeta, default_val)[:,None])

    # dijetpt
    dijetpt = ak.fill_none(ak.firsts(dijet.pt, axis=1), default_val)
    dijetpt = ak.from_regular(ak.where(njets >= 2, dijetpt, default_val)[:,None])

    
    # Create fartures
    features = {
        "pt_1"         : events.hcand.pt[:,0:1],
        "pt_2"         : events.hcand.pt[:,1:2],
        "met_pt"       : events.PuppiMET.pt[:,None],
        "met_dphi_1"   : events.dphi_met_h1,
        "met_dphi_2"   : events.dphi_met_h2,
        "dR"           : events.hcand_dr,
        "dphi"         : events.hcand_dphi,
        "pt_tt"        : events.pt_tt,
        "m_vis"        : events.hcand_invm,
        "pt_vis"       : events.pt_vis,
        "FastMTT_mass" : events.hcand_invm_fastMTT,
        "mt_1"         : events.mt_1,
        "mt_2"         : events.mt_2,
        "mt_lep"       : events.mt_lep,
        "mt_tot"       : events.mt_tot,
        "jpt_1"        : jpt_1,
        "jpt_2"        : jpt_2,
        "jeta_1"       : jeta_1,
        "jeta_2"       : jeta_2,
        "mjj"          : mjj,
        "jdeta"        : jdeta,
        "dijetpt"      : dijetpt,
        "n_jets"       : njets[:,None],
        #"n_bjets"      : nbjets[:,None],
    }
    # create arrays of feats
    x_flat_feats = {f: ak.to_numpy(a).flatten() for f,a in features.items()}
    # create pd df
    df = pd.DataFrame(x_flat_feats)

    # categorize events
    df_etau_even = df[is_etau_even].copy()
    df_etau_odd  = df[is_etau_odd].copy()
    df_mutau_even = df[is_mutau_even].copy()
    df_mutau_odd  = df[is_mutau_odd].copy()
    df_tautau_even = df[is_tautau_even].copy()
    df_tautau_odd  = df[is_tautau_odd].copy()
    
    #from IPython import embed; embed()

    #score_etau = -99.9 * np.ones((np.array(events.channel_id).shape[0], 3), dtype=float)
    #score_mutau = -99.9 * np.ones((np.array(events.channel_id).shape[0], 3), dtype=float)
    #score_tautau = -99.9 * np.ones((np.array(events.channel_id).shape[0], 3), dtype=float)

    score = -99.9 * np.ones((np.array(events.channel_id).shape[0], 3), dtype=float)


    h1_abs_eta = ak.to_numpy(np.abs(events.hcand.eta[:,0]))
    # Evaluation
    # ===>> etau : Even
    if df_etau_even.shape[0] > 0:
        # adding h1 eta to df
        df_etau_even.insert(2, 'abs_eta_1', h1_abs_eta[is_etau_even])
        score_etau_even = self.model_et_even.predict_proba(df_etau_even)
        score[is_etau_even] = score_etau_even
    else:
        logger.warning(f"0 events in df_etau_even")
    # Odd
    if df_etau_odd.shape[0] > 0:
        df_etau_odd.insert(2, 'abs_eta_1', h1_abs_eta[is_etau_odd])
        score_etau_odd  = self.model_et_odd.predict_proba(df_etau_odd)
        score[is_etau_odd]  = score_etau_odd
    else:
        logger.warning(f"0 events in df_etau_odd")

    # ===>> mutau : Even
    if df_mutau_even.shape[0] > 0:
        df_mutau_even.insert(2, 'abs_eta_1', h1_abs_eta[is_mutau_even])
        score_mutau_even = self.model_mt_even.predict_proba(df_mutau_even)
        score[is_mutau_even] = score_mutau_even
    else:
        logger.warning(f"0 events in df_etau_even")
    # Odd
    if df_mutau_odd.shape[0] > 0:
        df_mutau_odd.insert(2, 'abs_eta_1', h1_abs_eta[is_mutau_odd])        
        score_mutau_odd  = self.model_mt_odd.predict_proba(df_mutau_odd)
        score[is_mutau_odd]  = score_mutau_odd
    else:
        logger.warning(f"0 events in df_etau_odd")

    # ===>> tautau : Even 
    if df_tautau_even.shape[0] > 0:
        score_tautau_even = self.model_tt_even.predict_proba(df_tautau_even)
        score[is_tautau_even] = score_tautau_even
    else:
        logger.warning(f"0 events in df_tautau_even")
    # ===>> Odd 
    if df_tautau_odd.shape[0] > 0:
        score_tautau_odd  = self.model_tt_odd.predict_proba(df_tautau_odd)
        score[is_tautau_odd]  = score_tautau_odd
    else:
        logger.warning(f"0 events in df_tautau_odd")

    #from IPython import embed; embed()
        
    #score_tautau[is_tautau_even] = score_tautau_even
    #score_tautau[is_tautau_odd]  = score_tautau_odd

    """
    # Score tautau
    # 0: tau, 1: higgs, 2: fake
    score_tautau = ak.from_regular(ak.Array(score_tautau))
    score_dummy = score_tautau[:,:0]

    # DO the same for mutau and etau

    score = ak.where(is_tautau, score_tautau, score_dummy)
    events = set_ak_column(events, "classifier_score", score)
    """

    score = ak.from_regular(ak.Array(score)) 
    events = set_ak_column(events, "classifier_score", score) 
    
    return events


@classify_events.requires
def classify_events_requires(self: Producer, reqs: dict) -> None:
    if "external_files" in reqs:
        return
    
    from columnflow.tasks.external import BundleExternalFiles
    reqs["external_files"] = BundleExternalFiles.req(self.task)

@classify_events.setup
def classify_events_setup(
    self: Producer,
    reqs: dict,
    inputs: dict,
    reader_targets: InsertableDict,
) -> None:
    bundle = reqs["external_files"]

    model_etau_even_json = bundle.files.model_et_EVEN.path
    model_etau_odd_json  = bundle.files.model_et_ODD.path
    model_mutau_even_json = bundle.files.model_mt_EVEN.path
    model_mutau_odd_json  = bundle.files.model_mt_ODD.path
    model_tautau_even_json = bundle.files.model_tt_EVEN.path
    model_tautau_odd_json  = bundle.files.model_tt_ODD.path

    import xgboost as xgb

    self.model_et_even = xgb.XGBClassifier()
    self.model_et_even.load_model(model_etau_even_json)
    self.model_et_odd = xgb.XGBClassifier()
    self.model_et_odd.load_model(model_etau_odd_json)

    self.model_mt_even = xgb.XGBClassifier()
    self.model_mt_even.load_model(model_mutau_even_json)
    self.model_mt_odd = xgb.XGBClassifier()
    self.model_mt_odd.load_model(model_mutau_odd_json)    
    
    self.model_tt_even = xgb.XGBClassifier()
    self.model_tt_even.load_model(model_tautau_even_json)
    self.model_tt_odd = xgb.XGBClassifier()
    self.model_tt_odd.load_model(model_tautau_odd_json)
