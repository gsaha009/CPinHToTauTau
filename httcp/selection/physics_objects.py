# coding: utf-8

"""
Exemplary selection methods.
"""

from collections import defaultdict

from typing import Optional
from columnflow.selection import Selector, SelectionResult, selector
from columnflow.selection.cms.jets import jet_veto_map
from columnflow.selection.util import sorted_indices_from_mask
from columnflow.util import maybe_import, DotDict
from columnflow.columnar_util import EMPTY_FLOAT, Route, set_ak_column
from columnflow.columnar_util import optional_column as optional
from columnflow.production.util import attach_coffea_behavior

from httcp.util import IF_NANO_V9, IF_NANO_V11, IF_RUN2, IF_RUN3, getGenTauDecayMode

np = maybe_import("numpy")
ak = maybe_import("awkward")
coffea = maybe_import("coffea")

# ------------------------------------------------------------------------------------------------------- #
# Muon Selection
# Reference:
#   https://cms.cern.ch/iCMS/analysisadmin/cadilines?id=2325&ancode=HIG-20-006&tp=an&line=HIG-20-006
#   http://cms.cern.ch/iCMS/jsp/openfile.jsp?tp=draft&files=AN2019_192_v15.pdf
# ------------------------------------------------------------------------------------------------------- #
@selector(
    uses={
        f"Muon.{var}" for var in [
            "pt", "eta", "phi", "dxy", "dz", "mediumId", 
            "pfRelIso04_all", "isGlobal", "isPFcand", 
            "IPx", "IPy", "IPz", "sip3d",
            "isTracker",
        ]
    } | {optional("genPartFlav")},
    produces={
        f"Muon.{var}" for var in [
            "rawIdx", "decayMode", "IPsig", "isolation",
            "decayModeHPS", "SVx", "SVy", "SVz",
        ]
    },
    exposed=False,
)
def muon_selection(
        self: Selector,
        events: ak.Array,
        **kwargs
) -> tuple[ak.Array, SelectionResult, ak.Array, ak.Array, ak.Array]:
    """
    Muon selection returning two sets of indidces for default and veto muons.
    
    References:
      - Isolation working point: https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2?rev=59
      - ID und ISO : https://twiki.cern.ch/twiki/bin/view/CMS/MuonUL2017?rev=15
    """
    # Adding new columns to the muon collection for convenience
    events = set_ak_column(events, "Muon.rawIdx",    ak.local_index(events.Muon))
    events = set_ak_column(events, "Muon.decayMode", -2)
    events = set_ak_column(events, "Muon.decayModeHPS", -2)
    # make sure that there is no nan sip3d
    ipsig_dummy = 0.0 #-999.9
    events = set_ak_column(events, "Muon.IPsig", ak.nan_to_num(events.Muon.sip3d, nan=ipsig_dummy))
    events = set_ak_column(events, "Muon.isolation", events.Muon.pfRelIso04_all, value_type=np.float32)

    events = set_ak_column(events, "Muon.SVx", 0.0)
    events = set_ak_column(events, "Muon.SVy", 0.0)
    events = set_ak_column(events, "Muon.SVz", 0.0)
    
    # pt sorted indices for converting masks to indices
    sorted_indices = ak.argsort(events.Muon.pt, axis=-1, ascending=False)
    muons = events.Muon[sorted_indices]

    good_selections = {
        "muon_pt_20"          : muons.pt > 20,
        "muon_eta_2p4"        : abs(muons.eta) < 2.4,
        "muon_mediumID"       : muons.mediumId == 1,
        "muon_dxy_0p045"      : abs(muons.dxy) < 0.045,
        "muon_dz_0p2"         : abs(muons.dz) < 0.2,
        #"muon_iso_0p15"       : muons.pfRelIso04_all < 0.15, # For SR, moved to categorization
        "muon_iso_0p5"        : muons.pfRelIso04_all < 0.5, # new categorization following DESY setup 
        # not before applying correction from IPsig calibration
        #"muon_ipsig_1p0"      : np.abs(muons.IPsig) > 1.0,
    }
    single_veto_selections = {
        "muon_pt_10"          : muons.pt > 10,
        "muon_eta_2p4"        : abs(muons.eta) < 2.4,
        "muon_mediumID"       : muons.mediumId == 1,
        "muon_dxy_0p045"      : abs(muons.dxy) < 0.045,
        "muon_dz_0p2"         : abs(muons.dz) < 0.2,
        "muon_iso_0p3"        : muons.pfRelIso04_all < 0.3,
    }
    double_veto_selections = {
        "muon_pt_15"          : muons.pt > 15,
        "muon_eta_2p4"        : abs(muons.eta) < 2.4,
        "muon_isGlobal"       : muons.isGlobal == True,
        "muon_isPF"           : muons.isPFcand == True,
        "muon_isTracker"      : muons.isTracker ==True,
        "muon_dxy_0p045"      : abs(muons.dxy) < 0.045,
        "muon_dz_0p2"         : abs(muons.dz) < 0.2,
        "muon_iso_0p3"        : muons.pfRelIso04_all < 0.3,
    }

    muon_mask      = ak.local_index(muons) >= 0
    
    good_muon_mask = muon_mask
    single_veto_muon_mask = muon_mask
    double_veto_muon_mask = muon_mask

    good_selection_steps = {}
    single_veto_selection_steps = {}
    double_veto_selection_steps = {}
    
    good_selection_steps        = {"muon_starts_with": good_muon_mask}
    single_veto_selection_steps = {"muon_starts_with": good_muon_mask}
    double_veto_selection_steps = {"muon_starts_with": good_muon_mask}
            
    for cut in good_selections.keys():
        good_muon_mask = good_muon_mask & ak.fill_none(good_selections[cut], False)
        good_selection_steps[cut] = good_muon_mask

    for cut in single_veto_selections.keys():
        single_veto_muon_mask = single_veto_muon_mask & ak.fill_none(single_veto_selections[cut], False)
        single_veto_selection_steps[cut] = single_veto_muon_mask

    for cut in double_veto_selections.keys():
        double_veto_muon_mask = double_veto_muon_mask & ak.fill_none(double_veto_selections[cut], False)
        double_veto_selection_steps[cut] = double_veto_muon_mask


    # filtered muons
    good_muons = muons[good_muon_mask]
    single_veto_muons = muons[single_veto_muon_mask]
    double_veto_muons = muons[double_veto_muon_mask]

    # indices
    raw_muon_indices = events.Muon.rawIdx
    sorted_muon_indices = muons.rawIdx
    good_muon_indices = good_muons.rawIdx
    single_veto_muon_indices = single_veto_muons.rawIdx
    double_veto_muon_indices = double_veto_muons.rawIdx
    
    return events, SelectionResult(
        objects={
            "Muon": {
                "RawMuon": raw_muon_indices,
                "SortedMuon": sorted_muon_indices,
                "GoodMuon": good_muon_indices,
                "VetoMuon": single_veto_muon_indices,
                "DoubleVetoMuon": double_veto_muon_indices,
            },
        },
        aux={
            "muon_good_selection": good_selection_steps,
            "muon_single_veto_selection": single_veto_selection_steps,
            "muon_double_veto_selection": double_veto_selection_steps,
        }
    ), good_muon_indices, single_veto_muon_indices, double_veto_muon_indices



# ------------------------------------------------------------------------------------------------------- #
# Electron Selection
# Reference:
#   https://cms.cern.ch/iCMS/analysisadmin/cadilines?id=2325&ancode=HIG-20-006&tp=an&line=HIG-20-006
#   http://cms.cern.ch/iCMS/jsp/openfile.jsp?tp=draft&files=AN2019_192_v15.pdf
# ------------------------------------------------------------------------------------------------------- #
@selector(
    uses={
        "Electron.pt", "Electron.eta", "Electron.phi", "Electron.mass", "Electron.dxy", "Electron.dz",
        "Electron.pfRelIso03_all", "Electron.convVeto", #"lostHits",
        IF_NANO_V9("Electron.mvaFall17V2Iso_WP80", "Electron.mvaFall17V2Iso_WP90", "Electron.mvaFall17V2noIso_WP90"),
        IF_NANO_V11("Electron.mvaIso_WP80", "Electron.mvaIso_WP90", "Electron.mvaNoIso_WP90"),
        "Electron.cutBased",
        "Electron.IPx", "Electron.IPy", "Electron.IPz", "Electron.sip3d",
        optional("Electron.genPartFlav"),
    },
    produces={
        f"Electron.{var}" for var in [
            "rawIdx", "decayMode", "IPsig", "isolation",
            "decayModeHPS", "SVx", "SVy", "SVz",
        ]
    },
    exposed=False,
)
def electron_selection(
        self: Selector,
        events: ak.Array,
        **kwargs
) -> tuple[ak.Array, SelectionResult, ak.Array, ak.Array, ak.Array]:
    """
    Electron selection returning two sets of indidces for default and veto muons.
    
    References:
      - https://twiki.cern.ch/twiki/bin/view/CMS/EgammaNanoAOD?rev=4
    """
    # Adding new columns to the ele collection for convenience
    events = set_ak_column(events, "Electron.rawIdx",    ak.local_index(events.Electron))
    events = set_ak_column(events, "Electron.decayMode", -1)
    events = set_ak_column(events, "Electron.decayModeHPS", -1)
    # make sure that there is no nan sip3d
    ipsig_dummy = 0.0 #-999.9
    events = set_ak_column(events, "Electron.IPsig", ak.nan_to_num(events.Electron.sip3d, nan=ipsig_dummy))
    #events = set_ak_column(events, "Electron.isolation", -1.0)
    events = set_ak_column(events, "Electron.isolation", events.Electron.pfRelIso03_all, value_type=np.float32)

    events = set_ak_column(events, "Electron.SVx", 0.0)
    events = set_ak_column(events, "Electron.SVy", 0.0)
    events = set_ak_column(events, "Electron.SVz", 0.0)

    # pt sorted indices for converting masks to indices
    sorted_indices = ak.argsort(events.Electron.pt, axis=-1, ascending=False)
    electrons = events.Electron[sorted_indices]

    # >= nano v10
    mva_iso_wp80 = electrons.mvaIso_WP80
    mva_iso_wp90 = electrons.mvaIso_WP90
    mva_noniso_wp90 = electrons.mvaNoIso_WP90

    good_selections = {
        "electron_pt_25"          : electrons.pt > 25,
        "electron_eta_2p5"        : abs(electrons.eta) < 2.5,
        "electron_dxy_0p045"      : abs(electrons.dxy) < 0.045,
        "electron_dz_0p2"         : abs(electrons.dz) < 0.2,
        "electron_mva_iso_wp80"   : mva_iso_wp80 == 1,
        "electron_iso_0p15"       : electrons.pfRelIso03_all < 0.5,
        # not before applying correction from IPsig calibration        
        #"electron_ipsig_1p0"      : np.abs(electrons.IPsig) > 1.0,
    }
    single_veto_selections = {
        "electron_pt_10"          : electrons.pt > 10,
        "electron_eta_2p5"        : abs(electrons.eta) < 2.5,
        "electron_dxy_0p045"      : abs(electrons.dxy) < 0.045,
        "electron_dz_0p2"         : abs(electrons.dz) < 0.2,
        "electron_mva_noniso_wp90": mva_noniso_wp90 == 1,
        "electron_convVeto"       : electrons.convVeto == 1,
        #"electron_lostHits"       : electrons.lostHits <= 1,
        "electron_iso_0p3"        : electrons.pfRelIso03_all < 0.3,
    }
    double_veto_selections = {
        "electron_pt_15"          : electrons.pt > 15,
        "electron_eta_2p5"        : abs(electrons.eta) < 2.5,
        "electron_dxy_0p045"      : abs(electrons.dxy) < 0.045,
        "electron_dz_0p2"         : abs(electrons.dz) < 0.2,
        "electron_cutBased"       : electrons.cutBased >= 1, # probably == 1 is wrong !!!
        "electron_iso_0p3"        : electrons.pfRelIso03_all < 0.3,
    }

    electron_mask  = ak.local_index(electrons) >= 0

    good_electron_mask        = electron_mask
    single_veto_electron_mask = electron_mask
    double_veto_electron_mask = electron_mask

    good_selection_steps = {}
    single_veto_selection_steps = {}
    double_veto_selection_steps = {}

    good_selection_steps        = {"electron_starts_with": good_electron_mask}
    single_veto_selection_steps = {"electron_starts_with": good_electron_mask}
    double_veto_selection_steps = {"electron_starts_with": good_electron_mask}
    
    for cut in good_selections.keys():
        good_electron_mask = good_electron_mask & ak.fill_none(good_selections[cut], False)
        good_selection_steps[cut] = good_electron_mask

    for cut in single_veto_selections.keys():
        single_veto_electron_mask = single_veto_electron_mask & ak.fill_none(single_veto_selections[cut], False)
        single_veto_selection_steps[cut] = single_veto_electron_mask

    for cut in double_veto_selections.keys():
        double_veto_electron_mask = double_veto_electron_mask & ak.fill_none(double_veto_selections[cut], False)
        double_veto_selection_steps[cut] = double_veto_electron_mask


    # filtered electrons
    good_electrons = electrons[good_electron_mask]
    single_veto_electrons = electrons[single_veto_electron_mask]
    double_veto_electrons = electrons[double_veto_electron_mask]

    # indices
    raw_electron_indices = events.Electron.rawIdx
    sorted_electron_indices = electrons.rawIdx
    good_electron_indices = good_electrons.rawIdx
    single_veto_electron_indices = single_veto_electrons.rawIdx
    double_veto_electron_indices = double_veto_electrons.rawIdx

    return events, SelectionResult(
        objects={
            "Electron": {
                "RawElectron": raw_electron_indices,
                "SortedElectron": sorted_electron_indices,
                "GoodElectron": good_electron_indices,
                "VetoElectron": single_veto_electron_indices,
                "DoubleVetoElectron": double_veto_electron_indices,
            },
        },
        aux={
            "electron_good_selection": good_selection_steps,
            "electron_single_veto_selection": single_veto_selection_steps,
            "electron_double_veto_selection": double_veto_selection_steps,
        }
    ), good_electron_indices, single_veto_electron_indices, double_veto_electron_indices


# ------------------------------------------------------------------------------------------------------- #
# Tau Selection
# Reference:
#   https://cms.cern.ch/iCMS/analysisadmin/cadilines?id=2325&ancode=HIG-20-006&tp=an&line=HIG-20-006
#   http://cms.cern.ch/iCMS/jsp/openfile.jsp?tp=draft&files=AN2019_192_v15.pdf
# ------------------------------------------------------------------------------------------------------- #
@selector(
    uses={
        f"Tau.{var}" for var in [
            "eta", "phi", "dz",
            "idDeepTau2018v2p5VSe", "idDeepTau2018v2p5VSmu", "idDeepTau2018v2p5VSjet",
            "decayMode", "decayModePNet",
            "IPx","IPy","IPz", "ipLengthSig","hasRefitSV",
        ]
    } | {optional("Tau.pt"), optional("Tau.pt_etau"), optional("Tau.pt_mutau"), optional("Tau.pt_tautau"), optional("genPartFlav")},
    produces={
        f"Tau.{var}" for var in [
            "rawIdx", "decayMode", "decayModeHPS", "IPsig", "isolation",
            "SVx","SVy","SVz",
        ]
    },
    exposed=False,
)
def tau_selection(
        self: Selector,
        events: ak.Array,
        channel_: Optional[str] = "",
        electron_indices: Optional[ak.Array]=None,
        muon_indices    : Optional[ak.Array]=None,
        **kwargs
) -> tuple[ak.Array, SelectionResult, ak.Array]:
    """
    Tau selection returning two sets of indidces for default and veto muons.    
    References:
      - 
    """
    tau_local_indices = ak.local_index(events.Tau)

    events = set_ak_column(events, "Tau.rawIdx", tau_local_indices)
    # to get rid of any nan values
    #from IPython import embed; embed()
    #ipsig_dummy = ak.max(events.Tau.ipLengthSig) + np.abs(ak.min(events.Tau.ipLengthSig))
    ipsig_dummy = 0.0
    events = set_ak_column(events, "Tau.IPsig",  ak.nan_to_num(events.Tau.ipLengthSig, nan=ipsig_dummy))
    events = set_ak_column(events, "Tau.isolation", events.Tau.idDeepTau2018v2p5VSjet, value_type=np.float32)

    events = set_ak_column(events, "Tau.SVx", ak.nan_to_num(events.Tau.refitSVx, 0.0))
    events = set_ak_column(events, "Tau.SVy", ak.nan_to_num(events.Tau.refitSVy, 0.0))
    events = set_ak_column(events, "Tau.SVz", ak.nan_to_num(events.Tau.refitSVz, 0.0))
    
    # https://cms-nanoaod-integration.web.cern.ch/integration/cms-swmaster/data106Xul17v2_v10_doc.html#Tau
    tau_vs_e = DotDict(vvloose=2, vloose=3)
    tau_vs_mu = DotDict(vloose=1, tight=4)
    tau_vs_jet = DotDict(vvloose=2, loose=4, medium=5)
    
    tau_tagger         = self.config_inst.x.deep_tau_tagger
    tau_tagger_wps     = self.config_inst.x.deep_tau_info[tau_tagger].wp

    if "decayModeHPS" not in events.Tau.fields:
        events = set_ak_column(events, "Tau.decayModeHPS", events.Tau.decayMode)      # explicitly renaming decayMode to decayModeHPS
        events = set_ak_column(events, "Tau.decayMode",    events.Tau.decayModePNet)  # set decayModePNet as decayMode
    
    tau_pt = None
    if channel_ == "mutau":
        tau_pt = events.Tau.pt_mutau
    elif channel_ == "etau":
        tau_pt = events.Tau.pt_etau
    elif channel_ == "tautau":
        tau_pt = events.Tau.pt_tautau
    else:
        tau_pt = events.Tau.pt
        
    # [20.6] [19.4] [19.7]

    sorted_indices = ak.argsort(tau_pt, axis=-1, ascending=False)
    taus = events.Tau[sorted_indices]

    good_selections = {
        "tau_pt_15"     : taus.pt > 15.0, # 20 GeV is in pair selection
        "tau_eta_2p5"   : abs(taus.eta) < 2.5, # 2.3
        "tau_dz_0p2"    : abs(taus.dz) < 0.2,
        # have to make them channel-specific later
        #                  e-tau  mu-tau  tau-tau     SafeHere
        #   DeepTauVSjet : Tight  Medium  Medium  --> Medium  
        #   DeepTauVSe   : Tight  VVLoose VVLoose --> VVLoose 
        #   DeepTauVSmu  : Loose  Tight   VLoose  --> VLoose  
        "tau_DeepTauVSjet"  : taus.idDeepTau2018v2p5VSjet >= tau_tagger_wps.vs_j.VVVLoose, # for tautau fake region
        "tau_DeepTauVSe"    : taus.idDeepTau2018v2p5VSe   >= tau_tagger_wps.vs_e.VVLoose,
        "tau_DeepTauVSmu"   : taus.idDeepTau2018v2p5VSmu  >= tau_tagger_wps.vs_m.VLoose,
        "tau_HPSDMveto_5or6" : ((taus.decayModeHPS != 5) & (taus.decayModeHPS != 6)),
        "tau_no_undefinedPNetDM": (taus.decayMode != -1),
        "tau_DecayMode"  : (
            (taus.decayMode ==  0)
            | ((taus.decayMode ==  1) & (taus.decayModeHPS == 1))
            | ((taus.decayMode ==  2) & (taus.decayModeHPS == 1))
            | ((taus.decayMode ==  10) & taus.hasRefitSV)
            #| ((taus.decayMode ==  11) & taus.hasRefitSV)
        ),
        "tau_DM0_IPsig_1p25" : ak.where(taus.decayMode == 0, np.abs(taus.IPsig) >= 1.25, True),
    }
    
    tau_mask = ak.local_index(taus) >= 0
    
    good_tau_mask = tau_mask
    selection_steps = {}

    selection_steps = {"tau_starts_with": good_tau_mask}
    for cut in good_selections.keys():
        good_tau_mask = good_tau_mask & ak.fill_none(good_selections[cut], False)
        selection_steps[cut] = good_tau_mask
        
    if electron_indices is not None:
        good_tau_mask = good_tau_mask & ak.all(taus.metric_table(events.Electron[electron_indices]) > 0.2, axis=2)
        selection_steps["tau_clean_against_electrons"] = good_tau_mask 

    if muon_indices is not None:
        good_tau_mask = good_tau_mask & ak.all(taus.metric_table(events.Muon[muon_indices]) > 0.2, axis=2)
        selection_steps["tau_clean_against_muons"] = good_tau_mask

    good_taus = taus[good_tau_mask]

    # indices
    raw_tau_indices = events.Tau.rawIdx
    sorted_tau_indices = taus.rawIdx
    good_tau_indices = good_taus.rawIdx
    
    
    return events, SelectionResult(
        objects={
            "Tau": {
                "RawTau": raw_tau_indices,
                "SortedTau": sorted_tau_indices,
                "GoodTau": good_tau_indices,
            },
        },
        aux=selection_steps,
    ), good_tau_indices


@tau_selection.init
def tau_selection_init(self: Selector) -> None:
    # register tec shifts
    self.shifts |= {
        shift_inst.name
        for shift_inst in self.config_inst.shifts
        if shift_inst.has_tag("tec")
    }


    
# ------------------------------------------------------------------------------------------------------- #
# Jet Selection
# Reference:
#   https://cms.cern.ch/iCMS/analysisadmin/cadilines?id=2325&ancode=HIG-20-006&tp=an&line=HIG-20-006
#   http://cms.cern.ch/iCMS/jsp/openfile.jsp?tp=draft&files=AN2019_192_v15.pdf
# ------------------------------------------------------------------------------------------------------- #
@selector(
    uses={ f"Jet.{var}" for var in 
        [
            "pt", "eta", "phi", "mass",
            "jetId", "btagDeepFlavB",
            "neHEF", "neEmEF", "chMultiplicity", "neMultiplicity",
            "chHEF", "chMultiplicity", "muEF", "chEmEF",
        ]}
    | {optional("Jet.puId")}
    | {optional("Jet.genJetIdx")}
    | {IF_RUN3(jet_veto_map)}
    | {optional("GenJet.*")}
    | {attach_coffea_behavior},
    produces={"Jet.rawIdx"},
    exposed=False,
)
def jet_selection(
        self: Selector,
        events: ak.Array,
        **kwargs
) -> tuple[ak.Array, SelectionResult, ak.Array]:
    """
    This function vetoes b-jets with sufficiently high pt and incide eta region of interest
    """
    is_run2 = self.config_inst.campaign.x.run == 2
    is_run3 = self.config_inst.campaign.x.run == 3

    jet_mask  = ak.local_index(events.Jet.pt) >= 0 #Create a mask filled with ones

    # Redefination of JetID because of the bug in NanoAOD v12-v15
    # https://gitlab.cern.ch/cms-jetmet/coordination/coordination/-/issues/117
    passJetIdTight = ak.where(np.abs(events.Jet.eta) <= 2.6,
                              ((events.Jet.neHEF < 0.99)
                               & (events.Jet.neEmEF < 0.9)
                               & (events.Jet.chMultiplicity + events.Jet.neMultiplicity > 1)
                               & (events.Jet.chHEF > 0.01)
                               & (events.Jet.chMultiplicity > 0)),  # Tight criteria for abs_eta <= 2.6
                              ak.where((np.abs(events.Jet.eta) > 2.6) & (np.abs(events.Jet.eta) <= 2.7),
                                       ((events.Jet.neHEF < 0.9)
                                        & (events.Jet.neEmEF < 0.99)),  # Tight criteria for 2.6 < abs_eta <= 2.7
                                       ak.where((np.abs(events.Jet.eta) > 2.7) & (np.abs(events.Jet.eta) <= 3.0),
                                                events.Jet.neHEF < 0.99,  # Tight criteria for 2.7 < abs_eta <= 3.0
                                                ((events.Jet.neMultiplicity >= 2) & (events.Jet.neEmEF < 0.4))  # Tight criteria for abs_eta > 3.0
                                                )
                                       )
                              )
    
    # Default tight lepton veto
    passJetIdTightLepVeto = ak.where(
        np.abs(events.Jet.eta) <= 2.7,
        (passJetIdTight & (events.Jet.muEF < 0.8) & (events.Jet.chEmEF < 0.8)),  # add lepton veto for abs_eta <= 2.7
        passJetIdTight  # No lepton veto for 2.7 < abs_eta
    )
    
        
    # nominal selection
    good_selections = {
        "jet_pt_20"               : events.Jet.pt > 20.0,
        "jet_eta_4p7"             : abs(events.Jet.eta) <= 4.7,
        "jet_special_for_PU"      : ak.where(((abs(events.Jet.eta) >= 2.5) & (abs(events.Jet.eta) < 3.0)), events.Jet.pt > 50.0, jet_mask),
        "jet_forward"             : ak.where((abs(events.Jet.eta) >= 3.0), events.Jet.pt > 30.0, jet_mask),
        # use the newly defined TightLepVeto ID
        #"jet_id"                  : events.Jet.jetId >= 2,  # Jet ID flag: bit2 is tight, bit3 is tightLepVeto            
                                                            # So, 0000010 : 2**1 = 2 : pass tight, fail lep-veto          
                                                            #     0000110 : 2**1 + 2**2 = 6 : pass both tight and lep-veto
        "jet_id"                  : passJetIdTightLepVeto,
    }
        
    events = set_ak_column(events, "Jet.rawIdx", ak.local_index(events.Jet.pt))    

    selection_steps = {}
    selection_steps = {"jet_starts_with": jet_mask}
    for cut in good_selections.keys():
        jet_mask = jet_mask & ak.fill_none(good_selections[cut], False)
        selection_steps[cut] = jet_mask
    
    sorted_indices = ak.argsort(events.Jet.pt, axis=-1, ascending=False)
    good_jet_indices = sorted_indices[jet_mask[sorted_indices]]
    good_jet_indices = ak.values_astype(good_jet_indices, np.int32)
    
    results = SelectionResult(
        # bveto selection will be required later for ABCD
        objects = {
            "Jet": {
                "RawJet": events.Jet.rawIdx,
                "SortedJet": sorted_indices,
            },
        },
        aux = selection_steps,
    )
        
    return events, results, good_jet_indices


@jet_selection.init
def jet_selection_init(self: Selector) -> None:
    # register shifts
    self.shifts |= {
        shift_inst.name
        for shift_inst in self.config_inst.shifts
        if shift_inst.has_tag(("jec", "jer"))
    }




@selector(
    uses={ f"Jet.{var}" for var in 
           [
            "pt", "eta", "phi", "mass",
            "jetId", "btagDeepFlavB"
           ]} | {optional("hcand.pt"), optional("hcand.eta"), optional("hcand.phi"),
                 optional("hcand.mass"), optional("hcand.decayMode")} | {"cross_tau_jet_triggered", "cross_tau_triggered"} | {"channel_id"},
    #produces={"Jet.pass_ditaujet"},
    exposed=False,
)
def jet_cleaning(
        self: Selector,
        events: ak.Array,
        good_jet_indices: ak.Array,
        jet_selection_results: SelectionResult,
        ditaujet_jet_indices: Optional[ak.Array]=None,
        **kwargs
) -> tuple[ak.Array, SelectionResult]:
    """
    This function vetoes b-jets with sufficiently high pt and incide eta region of interest
    """

    selection_steps = jet_selection_results.aux

    jet_mask  = good_jet_indices >= 0

    if "hcand" in events.fields:
        good_jets = ak.with_name(events.Jet[good_jet_indices], "PtEtaPhiMLorentzVector")
        # should we clean the jets wrt the tauh only or with e/mu/tauh?
        #hcand = ak.with_name(events.hcand[events.hcand.decayMode >= 0], "PtEtaPhiMLorentzVector") # to make sure the presence of tauh only
        hcand = ak.with_name(events.hcand, "PtEtaPhiMLorentzVector")
        dr_jets_hcand = good_jets.metric_table(hcand)
        jet_is_closed_to_hcand = ak.any(dr_jets_hcand < 0.5, axis=-1)
        selection_steps["jet_isclean"] = ~jet_is_closed_to_hcand
        good_clean_jets = good_jets[~jet_is_closed_to_hcand]
        good_jet_indices = good_clean_jets.rawIdx

    # b-tagged jets, tight working point
    btag_wp = self.config_inst.x.btag_working_points.deepjet.medium
    b_jet_mask = (np.abs(events.Jet[good_jet_indices].eta) < 2.5) & (events.Jet[good_jet_indices].btagDeepFlavB >= btag_wp)
    b_jet_indices = good_jet_indices[b_jet_mask]
    selection_steps["jet_isbtag"] = ak.fill_none(b_jet_mask, False)

    # bjet veto
    bjet_veto = ak.sum(b_jet_mask, axis=1) == 0

    temp = ak.cartesian([good_jet_indices, ditaujet_jet_indices], axis=1)
    good_jet_indices_, ditaujet_jet_indices_ = ak.unzip(temp)
    temp2 = ak.any((good_jet_indices_ - ditaujet_jet_indices_ == 0), axis=1)

    trigger_mask = (events.channel_id == 4) & (ak.num(ditaujet_jet_indices_) > 0)  & events.cross_tau_jet_triggered & ~events.cross_tau_triggered
    mask = ak.where(trigger_mask, temp2, True)

    met_recoil_mask = ak.where((np.abs(events.Jet[good_jet_indices].eta) < 2.5), (events.Jet[good_jet_indices].pt > 30), True)
    met_recoil_mask = ak.where((np.abs(events.Jet[good_jet_indices].eta) > 4.7), (events.Jet[good_jet_indices].pt > 50), met_recoil_mask) 
    met_recoil_jet_indices = good_jet_indices[met_recoil_mask]

    
    results = SelectionResult(
        # check if ditaujet matched jet is in cleaned jets
        steps={
            "ditaujet_is_in_selected_jets": mask,
        },
        objects = {
            "Jet": {
                "Jet": good_jet_indices,
                "bJet": b_jet_indices,
                "trigJet": ditaujet_jet_indices,
                "metRecoilJet": met_recoil_jet_indices,
            },
        },
        aux = selection_steps,
    )
        
    return events, results, ditaujet_jet_indices #, bjet_veto


    
    
# ------------------------------------------------------------------------------------------------------- #
# GenTau Selection
# Reference:
#   GenPart =  ['eta', 'mass', 'phi', 'pt', 'genPartIdxMother', 'pdgId', 'status', 'statusFlags', 
#               'genPartIdxMotherG', 'distinctParentIdxG', 'childrenIdxG', 'distinctChildrenIdxG', 
#               'distinctChildrenDeepIdxG']
# ------------------------------------------------------------------------------------------------------- #
@selector(
    uses={
        "GenPart.*",
        "hcand.pt", "hcand.eta", "hcand.phi", "hcand.mass",
        "PV.x", "PV.y", "PV.z",
    },
    produces={
        'GenPart.rawIdx',
        'GenTau.rawIdx', 'GenTau.eta', 'GenTau.mass', 'GenTau.phi', 'GenTau.pt', 'GenTau.pdgId', 'GenTau.decayMode', 
        'GenTau.charge', 'GenTau.IPx', 'GenTau.IPy', 'GenTau.IPz', 'GenTau.IPmag',
        'GenTauProd.rawIdx', 'GenTauProd.eta', 'GenTauProd.mass', 'GenTauProd.phi', 'GenTauProd.pt', 'GenTauProd.pdgId',
        'GenTauProd.charge',
    },
    mc_only=True,
    exposed=False,
)
def gentau_selection(
        self: Selector,
        events: ak.Array,
        match: Optional[bool]=True,
        **kwargs
) -> tuple[ak.Array, SelectionResult]:
    """
    Selecting the generator level taus only, no martching here
    select the gen tau decay products as well
    
    References:
      - 
    """
    genpart_indices = ak.local_index(events.GenPart.pt)
    events = set_ak_column(events, "GenPart.rawIdx", genpart_indices)

    # masks to select gen tau+ and tau-    
    good_selections = {
        "genpart_pdgId"           : np.abs(events.GenPart.pdgId) == 15,
        "genpart_status"          : events.GenPart.status == 2,
        "genpart_status_flags"    : events.GenPart.hasFlags(["isPrompt", "isFirstCopy"]),
        "genpart_pt_10"           : events.GenPart.pt > 5.0, # CHANGE 10.0
        "genpart_eta_2p5"         : np.abs(events.GenPart.eta) < 3.0, # CHANGE 2.5
        "genpart_momid_25"        : events.GenPart[events.GenPart.distinctParent.genPartIdxMother].pdgId == 25,
        "genpart_mom_status_22"   : events.GenPart[events.GenPart.distinctParent.genPartIdxMother].status == 22,
    }
    
    gen_mask  = genpart_indices >= 0
    good_gen_mask = gen_mask

    selection_steps = {"genpart_starts_with": good_gen_mask}
    for cut in good_selections.keys():
        good_gen_mask = good_gen_mask & ak.fill_none(good_selections[cut], False)
        selection_steps[cut] = good_gen_mask

    gentau_indices = genpart_indices[good_gen_mask]
    
    gentaus = ak.with_name(events.GenPart[gentau_indices], "PtEtaPhiMLorentzVector")
    hcands  = ak.with_name(events.hcand, "PtEtaPhiMLorentzVector")

    # check the matching
    matched_gentaus = hcands.nearest(gentaus, threshold=0.8) if match else gentaus # CHANGE 0.5

    # nearest method can include None if a particle is not matched
    # so, taking care of the none values before adding it as a new column
    matched_gentaus_dummy = matched_gentaus[:,:0]
    """
    is_none = ak.sum(ak.is_none(matched_gentaus, axis=1), axis=1) > 0
    #matched_gentaus = ak.where(is_none, gentaus[:,:0], matched_gentaus)
    matched_gentaus = ak.where(is_none, matched_gentaus_dummy, matched_gentaus) # new
    """
    matched_gentaus = ak.drop_none(matched_gentaus) # CHANGE
    
    """
    has_full_match           = ~is_none
    """
    #has_two_matched_gentaus  = has_full_match & ak.fill_none(ak.num(matched_gentaus.rawIdx, axis=1) == 2, False)
    has_two_matched_gentaus = ak.fill_none(ak.num(matched_gentaus.rawIdx, axis=1) == 2, False) # CHANGE
    gentaus_of_opposite_sign = ak.fill_none(ak.sum(matched_gentaus.pdgId, axis=1) == 0, False)

    # new
    # filter, again
    """
    matched_gentaus = ak.where((has_full_match & has_two_matched_gentaus & gentaus_of_opposite_sign),
                               matched_gentaus,
                               matched_gentaus_dummy)
    """
    matched_gentaus = ak.where((has_two_matched_gentaus & gentaus_of_opposite_sign),
                               matched_gentaus,
                               matched_gentaus_dummy) 


    
    # Get gentau decay products 
    # hack: _apply_global_index [todo: https://github.com/columnflow/columnflow/discussions/430]
    # get decay modes for the GenTaus
    decay_gentau_indices = matched_gentaus.distinctChildrenIdxG
    decay_gentaus = events.GenPart._apply_global_index(decay_gentau_indices)
    gentaus_dm = getGenTauDecayMode(decay_gentaus)

    mask_genmatchedtaus_1 = ak.fill_none(ak.firsts(((gentaus_dm[:,:1]   == -2)
                                                    |(gentaus_dm[:,:1]  == -1)
                                                    |(gentaus_dm[:,:1]  ==  0) 
                                                    |(gentaus_dm[:,:1]  ==  1)
                                                    |(gentaus_dm[:,:1]  ==  2)
                                                    |(gentaus_dm[:,:1]  == 10)
                                                    |(gentaus_dm[:,:1]  == 11)), axis=1), False) # ele/mu/had
    mask_genmatchedtaus_2 = ak.fill_none(ak.firsts(((gentaus_dm[:,1:2]  ==  0) 
                                                    |(gentaus_dm[:,1:2] ==  1)
                                                    |(gentaus_dm[:,1:2] ==  2)
                                                    |(gentaus_dm[:,1:2] == 10)
                                                    |(gentaus_dm[:,1:2] == 11)), axis=1), False) # had only

    mask_genmatchedtaus   = mask_genmatchedtaus_1 & mask_genmatchedtaus_2

    # check decaymodes
    # make sure that the decay mode of hcand is the same as decay mode of GenTau
    dm_match_evt_mask = ak.num(hcands.decayMode, axis=1) == 2
    if match:
        has_2         = (ak.num(hcands.decayMode, axis=1) == 2) & (ak.num(gentaus_dm, axis=1) == 2)
        _gentaus_dm   = ak.where(has_2, gentaus_dm, gentaus_dm[:,:0])
        _hcands_dm    = ak.where(has_2, hcands.decayMode, hcands.decayMode[:,:0])
        dm_match_mask = _hcands_dm == _gentaus_dm
        dm_match_evt_mask = ak.sum(dm_match_mask, axis=1) == 2
    
    # -- N E W
    has_finite_decay_prods = ak.prod(ak.num(decay_gentaus.pdgId, axis=-1), axis=1) > 1
        
    matched_gentaus =  ak.where((mask_genmatchedtaus & dm_match_evt_mask & has_finite_decay_prods),
                                matched_gentaus, matched_gentaus_dummy)
    decay_gentaus_dummy = decay_gentaus[:,:0]
    decay_gentaus = ak.where((mask_genmatchedtaus & dm_match_evt_mask & has_finite_decay_prods),
                             decay_gentaus, decay_gentaus_dummy)
    gentaus_dm_dummy = gentaus_dm[:,:0]
    gentaus_dm = ak.where((mask_genmatchedtaus & dm_match_evt_mask & has_finite_decay_prods),
                          gentaus_dm, gentaus_dm_dummy)
    # -- N E W
    
    
    # creating a proper array to save it as a new column
    dummy_decay_gentaus = decay_gentaus[:,:0][:,None]
    decay_1             = decay_gentaus[:,:1]
    decay_2             = decay_gentaus[:,1:2]
    decay_gentaus       = ak.concatenate([decay_1, decay_2], axis=1)
    
    # WARNING: Not a smart way to convert ak.Array -> List -> ak.Array
    # Must use ak.enforce_type: NOT WORKING here, but it was good for hcand selection
    events = set_ak_column(events, "GenTau",           ak.Array(ak.to_list(matched_gentaus)))
    events = set_ak_column(events, "GenTau.decayMode", gentaus_dm)
    events = set_ak_column(events, "GenTau.mass",      ak.ones_like(events.GenTau.mass) * 1.777)
    events = set_ak_column(events, "GenTau.charge",    ak.where(events.GenTau.pdgId > 0,
                                                                -1, 
                                                                ak.where(events.GenTau.pdgId < 0, 
                                                                         1, 
                                                                         0)))
    events = set_ak_column(events, "GenTauProd",       ak.Array(ak.to_list(decay_gentaus)))
    events = set_ak_column(events, "GenTauProd.charge",ak.where(events.GenTauProd.pdgId > 0,
                                                                -1,
                                                                ak.where(events.GenTauProd.pdgId < 0,
                                                                         1, 0)))

    is_mu = np.abs(events.GenTauProd.pdgId) == 13
    is_e  = np.abs(events.GenTauProd.pdgId) == 11

    # new implementation for GenTau : IPx,IPy,IPz
    # Ref: https://indico.cern.ch/event/1451226/contributions/6253060/attachments/2976341/5239440/Pi_CP_28_11_24.pdf

    is_pi = (np.abs(events.GenTauProd.pdgId) == 211) | (np.abs(events.GenTauProd.pdgId) == 321) | (np.abs(events.GenTauProd.pdgId) == 323) | (np.abs(events.GenTauProd.pdgId) == 10321) | (np.abs(events.GenTauProd.pdgId) == 10211)
    is_pi_oneprong = (events.GenTau.decayMode == 0) & is_pi # 1-prong
    prod_e = events.GenTauProd[is_e]
    prod_mu = events.GenTauProd[is_mu]
    prod_pi = events.GenTauProd[is_pi_oneprong]
    prod = events.GenTauProd[is_e | is_mu | is_pi_oneprong]
    
    ### IMPACT PARAMETER ###
    # from IPython import embed; embed()
    # Displacement vector L components
    Lx = prod.vx - events.PV.x
    Ly = prod.vy - events.PV.y
    Lz = prod.vz - events.PV.z
    # Direction vector d components (normalized momentum of GenTau)
    # Calculate px, py, pz from pt, eta, phi
    prod_px = prod.pt * np.cos(prod.phi)
    prod_py = prod.pt * np.sin(prod.phi)
    prod_pz = prod.pt * np.sinh(prod.eta)
    d_norm = np.sqrt(prod_px**2 + prod_py**2 + prod_pz**2)
    dx = prod_px / d_norm
    dy = prod_py / d_norm
    dz = prod_pz / d_norm
    # Projection of L onto d (L_par components)
    L_dot_d = Lx * dx + Ly * dy + Lz * dz
    Lpar_x = L_dot_d * dx
    Lpar_y = L_dot_d * dy
    Lpar_z = L_dot_d * dz
    # Impact parameter vector IP = L - L_par (L_perp components)
    IPx = Lx - Lpar_x
    IPy = Ly - Lpar_y
    IPz = Lz - Lpar_z
    IPx = ak.firsts(IPx, axis=-1)
    IPy = ak.firsts(IPy, axis=-1)
    IPz = ak.firsts(IPz, axis=-1)
    # Replace None by -999 
    IPx = ak.fill_none(IPx, -999)
    IPy = ak.fill_none(IPy, -999)
    IPz = ak.fill_none(IPz, -999)
    # Set new columns for impact parameter components
    events = set_ak_column(events, "GenTau.IPx", IPx)
    events = set_ak_column(events, "GenTau.IPy", IPy)
    events = set_ak_column(events, "GenTau.IPz", IPz)
    # Optional: magnitude of the impact parameter vector
    IP_magnitude = np.sqrt(IPx**2 + IPy**2 + IPz**2)
    events = set_ak_column(events, "GenTau.IPmag", IP_magnitude)

    return events, SelectionResult(
        steps = {
            "has_two_matched_gentaus"  : has_two_matched_gentaus,
            "gentaus_of_opposite sign" : gentaus_of_opposite_sign,
            "valid_decay_products"     : mask_genmatchedtaus,
            "has_finite_decay_products": has_finite_decay_prods,
            "gen_DMs_same_as_hcands"   : dm_match_evt_mask, # probably not necessary, will test later
        },
        aux = selection_steps,
    )



# ------------------------------------------------------------------------------------------------------- #
# GenZ selection
# will be used for Zpt reweighting
# https://github.com/danielwinterbottom/ICHiggsTauTau/blob/UL_ditau/Analysis/HiggsTauTauRun2/src/HTTWeights.cc#L2079-L2114
# ------------------------------------------------------------------------------------------------------- #
def get_gen_p4_array(gen_part):
    # form LV
    gen_part = ak.Array(gen_part, behavior=coffea.nanoevents.methods.nanoaod.behavior)
    p4_gen_part = ak.with_name(gen_part, "PtEtaPhiMLorentzVector")    
    p4_gen_part = ak.zip(
        {
            "x" : p4_gen_part.px,
            "y" : p4_gen_part.py,
            "z" : p4_gen_part.pz,
            "t": p4_gen_part.energy,
        },
        with_name="LorentzVector",
        behavior=coffea.nanoevents.methods.vector.behavior,
    )
    #sum_p4_gen_part = ak.sum(p4_gen_part, axis=1)
    sum_p4_gen_part = ak.zip(
        {
            "x" : ak.sum(p4_gen_part.x, axis=1),
            "y" : ak.sum(p4_gen_part.y, axis=1),
            "z" : ak.sum(p4_gen_part.z, axis=1),
            "t" : ak.sum(p4_gen_part.energy, axis=1),
        },
        with_name="LorentzVector",
        behavior=coffea.nanoevents.methods.vector.behavior,
    )

    p4_gen_part_array = ak.zip(
        {
            "pt"  : ak.nan_to_num(sum_p4_gen_part.pt, 0.0),  # nan values are for empty genpart
            "eta" : ak.nan_to_num(sum_p4_gen_part.eta, 0.0), 
            "phi" : ak.nan_to_num(sum_p4_gen_part.phi, 0.0),
            "mass": ak.nan_to_num(sum_p4_gen_part.mass, 0.0),
        }
    )
    return p4_gen_part_array



@selector(
    uses={
        "GenPart.*",
    },
    produces={
        "GenZ.pt", "GenZ.eta", "GenZ.phi", "GenZ.mass",
        "GenZvis.pt", "GenZvis.eta", "GenZvis.phi", "GenZvis.mass",
        "GenTop_pt",
    },
    mc_only=True,
    exposed=False,
)
def genZ_selection(
        self: Selector,
        events: ak.Array,
        **kwargs
) -> ak.Array:
    """
    
    References:
      - https://gitlab.cern.ch/dwinterb/HiggsDNA/-/blob/master/higgs_dna/tools/ditau/add_gen_info.py#L33
    """
    genpart_indices = ak.local_index(events.GenPart.pt)

    """
    sel_gen_part = (
        (((np.abs(events.GenPart.pdgId) >= 11) & (np.abs(events.GenPart.pdgId) <= 16))
         & (events.GenPart.hasFlags(["isHardProcess"]))
         & (events.GenPart.status == 1)) | (events.GenPart.hasFlags(["isDirectHardProcessTauDecayProduct"]))
    ) 
    sel_gen_nus = (
        (((np.abs(events.GenPart.pdgId) == 12) | (np.abs(events.GenPart.pdgId) == 14) | (np.abs(events.GenPart.pdgId) == 16))
         & (events.GenPart.hasFlags(["isHardProcess"]))
         & (events.GenPart.status == 1)) | (events.GenPart.hasFlags(["isDirectHardProcessTauDecayProduct"]))
    )
    sel_gen_vis_part = sel_gen_part & ~sel_gen_nus

    sel_top = 
    """
    genParts = events.GenPart
    # Define cuts to select generator-level particles:
    # get all decay products of leptonic W, Z, or H decays
    sel_gen_part = ((np.abs(genParts.pdgId) >= 11) & (np.abs(genParts.pdgId) <= 16) & (genParts.statusFlags & 256 != 0) & (genParts.status == 1)) | (genParts.statusFlags & 1024 != 0)

    # gen neutrinos from tau decays
    sel_gen_nus = ((np.abs(genParts.pdgId) == 12) | (np.abs(genParts.pdgId) == 14) | (np.abs(genParts.pdgId) == 16)) & (genParts.status == 1) & (genParts.statusFlags & 1024 != 0)

    # get visible decay products (pdgId 11, 13, 15)
    sel_gen_vis_part = sel_gen_part & ~sel_gen_nus

    # top quarks for top pt reweighting
    sel_gen_top = ((np.abs(genParts.pdgId) == 6) & (genParts.statusFlags & 256 != 0) & (genParts.statusFlags & 8192 != 0))

    sel_gen_ids = genpart_indices[sel_gen_part]
    gen_part = events.GenPart[sel_gen_ids]

    sel_gen_vis_ids = genpart_indices[sel_gen_vis_part]
    gen_vis_part = events.GenPart[sel_gen_vis_ids]

    sel_gen_top_ids = genpart_indices[sel_gen_top]
    gen_top = events.GenPart[sel_gen_top_ids]

    #from IPython import embed; embed()
    
    p4_gen_part_array = get_gen_p4_array(gen_part)
    p4_gen_vis_part_array = get_gen_p4_array(gen_vis_part)

    events = set_ak_column(events, "GenZ", p4_gen_part_array)
    events = set_ak_column(events, "GenZvis", p4_gen_vis_part_array)
    events = set_ak_column(events, "GenTop_pt", gen_top.pt)
    
    return events    
