# coding: utf-8

"""
main calibration script
"""

import functools

import law

from columnflow.calibration import Calibrator, calibrator
from columnflow.calibration.cms.jets import jets, jec, jec_nominal, jer
from columnflow.production.cms.seeds import deterministic_seeds
from columnflow.util import maybe_import
from columnflow.columnar_util import set_ak_column
from columnflow.columnar_util import optional_column as optional
from columnflow.production.util import attach_coffea_behavior

from httcp.calibration.electron import electron_smearing_scaling
from httcp.calibration.tau import tau_energy_scale
from httcp.util import IF_RUN2, IF_RUN3

np = maybe_import("numpy")
ak = maybe_import("awkward")

set_ak_column_f32 = functools.partial(set_ak_column, value_type=np.float32)

# derive calibrators to add settings
#jec_full = jec.derive("jec_full", cls_dict={"mc_only": True, "nominal_only": True})
logger = law.logger.get_logger(__name__)


@calibrator(
    uses={
        deterministic_seeds,
        "Jet.*",
        "RawPuppiMET.*",
        "PuppiMET.*",
        jets,
        electron_smearing_scaling, # comment for 2023
        tau_energy_scale,
    } | {optional("GenJet.*")} | {attach_coffea_behavior},
    produces={
        deterministic_seeds,
        "Jet.pt_no_corr", "Jet.phi_no_corr", "Jet.eta_no_corr", "Jet.mass_no_corr",
        "Jet.isgenmatched",
        "PuppiMET.pt_no_corr", "PuppiMET.phi_no_corr",
        "Electron.pt_no_ss", # comment for 2023
        jets,
        electron_smearing_scaling, # comment for 2023
        tau_energy_scale,
    },
)
def main(self: Calibrator, events: ak.Array, **kwargs) -> ak.Array:
    events = self[deterministic_seeds](events, **kwargs)

    events = set_ak_column_f32(events, "Jet.pt_no_corr", events.Jet.pt)
    events = set_ak_column_f32(events, "Jet.phi_no_corr", events.Jet.phi)
    events = set_ak_column_f32(events, "Jet.eta_no_corr", events.Jet.eta)
    events = set_ak_column_f32(events, "Jet.mass_no_corr", events.Jet.mass)
    
    #PuppiMET variables before applying energy corrections
    events = set_ak_column_f32(events, "PuppiMET.pt_no_corr", events.PuppiMET.pt)
    events = set_ak_column_f32(events, "PuppiMET.phi_no_corr", events.PuppiMET.phi)

    genjet_matched = ak.values_astype(ak.ones_like(events.Jet.pt), np.int32)
    if self.dataset_inst.is_mc:
        events = self[attach_coffea_behavior](events, collections=["Jet", "GenJet"], **kwargs)
        jet_is_matched_to_genjet = events.Jet.matched_gen
        temp_ = ak.fill_none(jet_is_matched_to_genjet.pt, -1.0)
        genjet_matched = ak.values_astype(ak.where(temp_ < 0.0, 0, 1), np.int32)

    events = set_ak_column(events, "Jet.isgenmatched", genjet_matched)
    
    if self.config_inst.campaign.x.run == 3:
        events = ak.with_field(events, events.RawPuppiMET, "RawMET")
        events = ak.with_field(events, events.PuppiMET, "MET")
        events = ak.without_field(events, "RawPuppiMET")
        events = ak.without_field(events, "PuppiMET")

    logger.warning("JER is not going to be applied to jets in MC with 2.5 < eta < 3.0 (NEW) and has no genjet match : L716-L722 (columnflow.calibration.jets.py)")
    events = self[jets](events, **kwargs)

    if self.config_inst.campaign.x.run == 3:
        events = ak.with_field(events, events.RawMET, "RawPuppiMET")
        events = ak.with_field(events, events.MET, "PuppiMET")
        events = ak.without_field(events, "RawMET")
        events = ak.without_field(events, "MET")

    # BE CAREFUL!!!
    # TES CAN BE SWITCHED OFF FOR SYNC STUDY
    if self.dataset_inst.is_mc: 
        events = self[tau_energy_scale](events, **kwargs)

    # electron scale and smearing correction
    events = set_ak_column_f32(events, "Electron.pt_no_ss", events.Electron.pt) # comment for 2023
    events = self[electron_smearing_scaling](events, **kwargs) # comment for 2023
        
    return events
