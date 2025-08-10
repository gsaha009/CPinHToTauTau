# coding: utf-8

"""
Definition of triggers
"""

import order as od

from httcp.config.trigger_util import Trigger, TriggerLeg



# ----------------------------------------------- #
#                   Run2 2016 UL                  #
# ----------------------------------------------- #
def add_triggers_UL2016(config: od.Config, postfix: str) -> None:
    """
    Adds all triggers to a *config*. For the conversion from filter names to trigger bits, see
    https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/NanoAOD/python/triggerObjects_cff.py.
    """
    config.x.triggers = od.UniqueObjectIndex(Trigger, [
        # ===>>> single electron
        Trigger(
            name="HLT_Ele25_eta2p1_WPTight_Gsf",
            id=111,
            legs=[
                TriggerLeg(
                    pdg_id=11,
                    min_pt=26.0,
                    # filter names:
                    # hltEle25erWPTightGsfTrackIsoFilter
                    trigger_bits=2,
                ),
            ],
            tags={"single_trigger", "single_e", "channel_e_tau"},
        ),
        # ===>>> single muon
        Trigger(
            name="HLT_IsoMu22",
            id=131,
            legs=[
                TriggerLeg(
                    pdg_id=13,
                    min_pt=23.0,
                    # filter names:
                    trigger_bits=2,
                ),
            ],
            tags={"single_trigger", "single_mu", "channel_mu_tau"},
        ),
        Trigger(
            name="HLT_IsoMu22_eta2p1",
            id=132,
            legs=[
                TriggerLeg(
                    pdg_id=13,
                    min_pt=23.0,
                    # filter names:
                    # hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07
                    trigger_bits=2,
                ),
            ],
            tags={"single_trigger", "single_mu", "channel_mu_tau"},
        ),
        Trigger(
            name="HLT_IsoTkMu22",
            id=133,
            legs=[
                TriggerLeg(
                    pdg_id=13,
                    min_pt=23.0,
                    # filter names:
                    # hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07
                    trigger_bits=2,
                ),
            ],
            tags={"single_trigger", "single_mu", "channel_mu_tau"},
        ),
        Trigger(
            name="HLT_IsoTkMu22_eta2p1",
            id=134,
            legs=[
                TriggerLeg(
                    pdg_id=13,
                    min_pt=23.0,
                    # filter names:
                    # hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07
                    trigger_bits=2,
                ),
            ],
            tags={"single_trigger", "single_mu", "channel_mu_tau"},
        ),
        # ===>>> e-tauh

        # ===>>> mu-tauh
        Trigger(
            name="HLT_IsoMu19_eta2p1_LooseIsoPFTau20",
            id=13151,
            legs=[
                TriggerLeg(
                    pdg_id=13,
                    min_pt=20.0,
                    # filter names:
                    # hltL3crIsoL1sMu18erTau24erIorMu20erTau24erL1f0L2f10QL3f20QL3trkIsoFiltered0p07
                    # hltOverlapFilterIsoMu20LooseChargedIsoPFTau27L1Seeded
                    trigger_bits=2 + 64,
                ),
                TriggerLeg(
                    pdg_id=15,
                    min_pt=21.0,
                    # filter names:
                    # hltSelectedPFTau27LooseChargedIsolationAgainstMuonL1HLTMatched or
                    # hltOverlapFilterIsoMu20LooseChargedIsoPFTau27L1Seeded
                    trigger_bits=1024 + 512,
                ),
            ],
            tags={"cross_trigger", "cross_mu_tau", "channel_mu_tau"},
        ),
        Trigger(
            name="HLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1",
            id=13152,
            legs=[
                TriggerLeg(
                    pdg_id=13,
                    min_pt=20.0,
                    # filter names:
                    # hltL3crIsoL1sMu18erTau24erIorMu20erTau24erL1f0L2f10QL3f20QL3trkIsoFiltered0p07
                    # hltOverlapFilterIsoMu20LooseChargedIsoPFTau27L1Seeded
                    trigger_bits=2 + 64,
                ),
                TriggerLeg(
                    pdg_id=15,
                    min_pt=21.0,
                    # filter names:
                    # hltSelectedPFTau27LooseChargedIsolationAgainstMuonL1HLTMatched or
                    # hltOverlapFilterIsoMu20LooseChargedIsoPFTau27L1Seeded
                    trigger_bits=1024 + 512,
                ),
            ],
            tags={"cross_trigger", "cross_mu_tau", "channel_mu_tau"},
        ),
        # ===>>> tauh-tauh
        Trigger(
            name="HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg",
            id=15151,
            legs=[
                TriggerLeg(
                    pdg_id=15,
                    min_pt=36.0,
                    # filter names:
                    # hltDoublePFTau35TrackPt1MediumChargedIsolationAndTightOOSCPhotonsDz02Reg
                    trigger_bits=64,
                ),
                TriggerLeg(
                    pdg_id=15,
                    min_pt=36.0,
                    # filter names:
                    # hltDoublePFTau35TrackPt1MediumChargedIsolationAndTightOOSCPhotonsDz02Reg
                    trigger_bits=64,
                ),
            ],
            applies_to_dataset=(lambda dataset_inst: dataset_inst.is_mc or dataset_inst.is_mc) if postfix == "preVFP" else (lambda dataset_inst: dataset_inst.x.era < "H"),
            tags={"cross_trigger", "cross_tau_tau", "channel_tau_tau"},
        ),
        Trigger(
            name="HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg",
            id=15152,
            legs=[
                TriggerLeg(
                    pdg_id=15,
                    min_pt=36.0,
                    # filter names:
                    # hltDoublePFTau35TrackPt1TightChargedIsolationAndTightOOSCPhotonsDz02Reg
                    trigger_bits=64,
                ),
                TriggerLeg(
                    pdg_id=15,
                    min_pt=36.0,
                    # filter names:
                    # hltDoublePFTau35TrackPt1TightChargedIsolationAndTightOOSCPhotonsDz02Reg
                    trigger_bits=64,
                ),
            ],
            applies_to_dataset=(lambda dataset_inst: dataset_inst.is_mc) if postfix == "preVFP" else (lambda dataset_inst: dataset_inst.is_mc or dataset_inst.x.era == "H"),
            tags={"cross_trigger", "cross_tau_tau", "channel_tau_tau"},
        ),
    ])    



# ----------------------------------------------- #
#                   Run2 2017 UL                  #
# ----------------------------------------------- #
def add_triggers_UL2017(config: od.Config) -> None:
    """
    Adds all triggers to a *config*. For the conversion from filter names to trigger bits, see
    https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/NanoAOD/python/triggerObjects_cff.py.
    """
    config.x.triggers = od.UniqueObjectIndex(Trigger, [
        # ===>>> single electron
        Trigger(
            name="HLT_Ele27_WPTight_Gsf",
            id=111,
            legs=[
                TriggerLeg(
                    pdg_id=11,
                    min_pt=28.0,
                    # filter names:
                    # hltEle32L1DoubleEGWPTightGsfTrackIsoFilter
                    # hltEGL1SingleEGOrFilter
                    trigger_bits=2,
                ),
            ],
            tags={"single_trigger", "single_e", "channel_e_tau"},
        ),
        Trigger(
            name="HLT_Ele32_WPTight_Gsf",
            id=112,
            legs=[
                TriggerLeg(
                    pdg_id=11,
                    min_pt=28.0,
                    # filter names:
                    # hltEle32WPTightGsfTrackIsoFilter
                    trigger_bits=2 + 1024,
                ),
            ],
            applies_to_dataset=(lambda dataset_inst: dataset_inst.is_mc or dataset_inst.x.era >= "D"),
            tags={"single_trigger", "single_e", "channel_e_tau"},
        ),
        # ===>>> single muon
        Trigger(
            name="HLT_IsoMu24",
            id=131,
            legs=[
                TriggerLeg(
                    pdg_id=13,
                    min_pt=25.0,
                    # filter names:
                    # hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p07
                    trigger_bits=2,
                ),
            ],
            tags={"single_trigger", "single_mu", "channel_mu_tau"},
        ),
        Trigger(
            name="HLT_IsoMu27",
            id=132,
            legs=[
                TriggerLeg(
                    pdg_id=13,
                    min_pt=25.0,
                    # filter names:
                    # hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07
                    trigger_bits=2,
                ),
            ],
            tags={"single_trigger", "single_mu", "channel_mu_tau"},
        ),
        # ===>>> e-tauh
        Trigger(
            name="HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1",
            id=11151,
            legs=[
                TriggerLeg(
                    pdg_id=11,
                    min_pt=25.0,
                    # filter names:
                    # hltEle24erWPTightGsfTrackIsoFilterForTau
                    # hltOverlapFilterIsoEle24WPTightGsfLooseIsoPFTau30
                    trigger_bits=2 + 64,
                ),
                TriggerLeg(
                    pdg_id=15,
                    min_pt=32.0,
                    # filter names:
                    # hltSelectedPFTau30LooseChargedIsolationL1HLTMatched
                    # hltOverlapFilterIsoEle24WPTightGsfLooseIsoPFTau30
                    trigger_bits=1024 + 256,
                ),
            ],
            tags={"cross_trigger", "cross_e_tau", "channel_e_tau"},
        ),
        # ===>>> mu-tauh
        Trigger(
            name="HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1",
            id=13151,
            legs=[
                TriggerLeg(
                    pdg_id=13,
                    min_pt=21.0,
                    # filter names:
                    # hltL3crIsoL1sMu18erTau24erIorMu20erTau24erL1f0L2f10QL3f20QL3trkIsoFiltered0p07
                    # hltOverlapFilterIsoMu20LooseChargedIsoPFTau27L1Seeded
                    trigger_bits=2 + 64,
                ),
                TriggerLeg(
                    pdg_id=15,
                    min_pt=32.0,
                    # filter names:
                    # hltSelectedPFTau27LooseChargedIsolationAgainstMuonL1HLTMatched or
                    # hltOverlapFilterIsoMu20LooseChargedIsoPFTau27L1Seeded
                    trigger_bits=1024 + 512,
                ),
            ],
            tags={"cross_trigger", "cross_mu_tau", "channel_mu_tau"},
        ),
        # ===>>> tauh-tauh
        Trigger(
            name="HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg",
            id=15151,
            legs=[
                TriggerLeg(
                    pdg_id=15,
                    min_pt=42.0,
                    # filter names:
                    # hltDoublePFTau35TrackPt1MediumChargedIsolationAndTightOOSCPhotonsDz02Reg
                    trigger_bits=64,
                ),
                TriggerLeg(
                    pdg_id=15,
                    min_pt=42.0,
                    # filter names:
                    # hltDoublePFTau35TrackPt1MediumChargedIsolationAndTightOOSCPhotonsDz02Reg
                    trigger_bits=64,
                ),
            ],
            tags={"cross_trigger", "cross_tau_tau", "channel_tau_tau"},
        ),
        Trigger(
            name="HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg",
            id=15152,
            legs=[
                TriggerLeg(
                    pdg_id=15,
                    min_pt=42.0,
                    # filter names:
                    # hltDoublePFTau35TrackPt1TightChargedIsolationAndTightOOSCPhotonsDz02Reg
                    trigger_bits=64,
                ),
                TriggerLeg(
                    pdg_id=15,
                    min_pt=42.0,
                    # filter names:
                    # hltDoublePFTau35TrackPt1TightChargedIsolationAndTightOOSCPhotonsDz02Reg
                    trigger_bits=64,
                ),
            ],
            #applies_to_dataset=(lambda dataset_inst: dataset_inst.is_data),
            tags={"cross_trigger", "cross_tau_tau", "channel_tau_tau"},
        ),
        Trigger(
            name="HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg",
            id=15153,
            legs=[
                TriggerLeg(
                    pdg_id=15,
                    min_pt=37.0,
                    # filter names:
                    # hltDoublePFTau40TrackPt1MediumChargedIsolationAndTightOOSCPhotonsDz02Reg
                    trigger_bits=64,
                ),
                TriggerLeg(
                    pdg_id=15,
                    min_pt=37.0,
                    # filter names:
                    # hltDoublePFTau40TrackPt1MediumChargedIsolationAndTightOOSCPhotonsDz02Reg
                    trigger_bits=64,
                ),
            ],
            #applies_to_dataset=(lambda dataset_inst: dataset_inst.is_data),
            tags={"cross_trigger", "cross_tau_tau", "channel_tau_tau"},
        ),
    ])    



# ----------------------------------------------- #
#                   Run2 2018 UL                  #
# ----------------------------------------------- #
def add_triggers_UL2018(config: od.Config) -> None:
    """
    Adds all triggers to a *config*. For the conversion from filter names to trigger bits, see
    https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/NanoAOD/python/triggerObjects_cff.py.
    """
    config.x.triggers = od.UniqueObjectIndex(Trigger, [
        # ===>>> single muon
        Trigger(
            name="HLT_IsoMu24",
            id=131000,
            legs=[
                TriggerLeg(
                    pdg_id=13,
                    min_pt=25.0,
                    max_abseta=None,
                    # filter names:
                    # hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p07
                    trigger_bits=2,
                ),
            ],
            tags={"single_trigger", "single_mu", "channel_mu_tau"},
        ),
        Trigger(
            name="HLT_IsoMu27",
            id=132000,
            legs=[
                TriggerLeg(
                    pdg_id=13,
                    min_pt=25.0,
                    max_abseta=None,
                    # filter names:
                    # hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07
                    trigger_bits=2,
                ),
            ],
            tags={"single_trigger", "single_mu", "channel_mu_tau"},
        ),
        # ===>>> mu-tauh
        Trigger(
            name="HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1",
            id=13151,
            run_range=(None,317509), #315974),
            legs=[
                TriggerLeg(
                    pdg_id=13,
                    min_pt=22.0, # 21.0
                    max_abseta=2.1, # None
                    # filter names:
                    # hltL3crIsoL1sMu18erTau24erIorMu20erTau24erL1f0L2f10QL3f20QL3trkIsoFiltered0p07
                    # hltOverlapFilterIsoMu20LooseChargedIsoPFTau27L1Seeded
                    trigger_bits=64, #2 + 64,
                ),
                TriggerLeg(
                    pdg_id=15, 
                    min_pt=32.0,
                    max_abseta=2.1,  #None,
                    # filter names:
                    # hltSelectedPFTau27LooseChargedIsolationAgainstMuonL1HLTMatched or
                    # hltOverlapFilterIsoMu20LooseChargedIsoPFTau27L1Seeded
                    trigger_bits=64, #1024 + 512,
                ),
            ],
            applies_to_dataset=(lambda dataset_inst: dataset_inst.is_data and dataset_inst.x.era <= "B"),
            tags={"cross_trigger", "cross_mu_tau", "channel_mu_tau"},
        ),
        #Trigger(
        #    name="HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1",
        #    id=13152,
        #    run_range=(315974,317509),
        #    legs=[
        #        TriggerLeg(
        #            pdg_id=13,
        #            min_pt=21.0,
        #            max_abseta=None,
        #            # filter names:
        #            # hltL3crIsoL1sMu18erTau24erIorMu20erTau24erL1f0L2f10QL3f20QL3trkIsoFiltered0p07
        #            # hltOverlapFilterIsoMu20LooseChargedIsoPFTau27L1Seeded
        #            trigger_bits=2 + 64,
        #        ),
        #        TriggerLeg(
        #            pdg_id=15,
        #            min_pt=32.0,
        #            max_abseta=None,
        #            # filter names:
        #            # hltSelectedPFTau27LooseChargedIsolationAgainstMuonL1HLTMatched or
        #            # hltOverlapFilterIsoMu20LooseChargedIsoPFTau27L1Seeded
        #            trigger_bits=1024 + 512,
        #        ),
        #    ],
        #    tags={"cross_trigger", "cross_mu_tau", "channel_mu_tau"},
        #),
        Trigger(
            name="HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1",
            id=13153,
            run_range=(317508,None),
            legs=[
                TriggerLeg(
                    pdg_id=13,
                    min_pt=22.0,
                    max_abseta=2.1,
                    # filter names:
                    # hltL3crIsoL1sMu18erTau24erIorMu20erTau24erL1f0L2f10QL3f20QL3trkIsoFiltered0p07
                    # hltOverlapFilterIsoMu20LooseChargedIsoPFTau27L1Seeded
                    trigger_bits=64, #2 + 64,
                ),
                TriggerLeg(
                    pdg_id=15,
                    min_pt=32.0,
                    max_abseta=2.1,
                    # filter names:
                    # hltSelectedPFTau27LooseChargedIsolationAgainstMuonL1HLTMatched or
                    # hltOverlapFilterIsoMu20LooseChargedIsoPFTau27L1Seeded
                    trigger_bits=512, #1024 + 512,
                ),
            ],
            applies_to_dataset=(lambda dataset_inst: dataset_inst.is_mc or dataset_inst.x.era >= "B"),
            tags={"cross_trigger", "cross_mu_tau", "channel_mu_tau"},
        ),        
        # ===>>> single electron
        Trigger(
            name="HLT_Ele32_WPTight_Gsf",
            id=111000,
            legs=[
                TriggerLeg(
                    pdg_id=11,
                    min_pt=33.0,
                    max_abseta=None,
                    # filter names:
                    # hltEle32WPTightGsfTrackIsoFilter
                    trigger_bits=2,
                ),
            ],
            tags={"single_trigger", "single_e", "channel_e_tau"},
        ),
        Trigger(
            name="HLT_Ele35_WPTight_Gsf",
            id=112000,
            legs=[
                TriggerLeg(
                    pdg_id=11,
                    min_pt=33.0,
                    max_abseta=None,
                    # filter names:
                    # hltEle35noerWPTightGsfTrackIsoFilter
                    trigger_bits=2,
                ),
            ],
            tags={"single_trigger", "single_e", "channel_e_tau"},
        ),
        # ===>>> e-tauh
        Trigger(
            name="HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1",
            id=11151,
            run_range=(None,317509),
            legs=[
                TriggerLeg(
                    pdg_id=11,
                    min_pt=26.0, #25.0,
                    max_abseta=2.1, #None,
                    # filter names:
                    # hltEle24erWPTightGsfTrackIsoFilterForTau
                    # hltOverlapFilterIsoEle24WPTightGsfLooseIsoPFTau30
                    trigger_bits=64, #2 + 64,
                ),
                TriggerLeg(
                    pdg_id=15,
                    min_pt=35.0,
                    max_abseta=2.1, #None,
                    # filter names:
                    # hltSelectedPFTau30LooseChargedIsolationL1HLTMatched
                    # hltOverlapFilterIsoEle24WPTightGsfLooseIsoPFTau30
                    trigger_bits=256, #1024 + 256,
                ),
            ],
            applies_to_dataset=(lambda dataset_inst: dataset_inst.is_data and dataset_inst.x.era <= "B"),
            tags={"cross_trigger", "cross_e_tau", "channel_e_tau"},
        ),
        Trigger(
            name="HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_CrossL1",
            id=11152,
            run_range=(317508,None),
            legs=[
                TriggerLeg(
                    pdg_id=11,
                    min_pt=26.0, #25.0,
                    max_abseta=2.1, #None,
                    # filter names:
                    # hltEle24erWPTightGsfTrackIsoFilterForTau
                    # hltOverlapFilterIsoEle24WPTightGsfLooseIsoPFTau30
                    trigger_bits=2 + 64,
                ),
                TriggerLeg(
                    pdg_id=15,
                    min_pt=35.0,
                    max_abseta=2.1, #None,
                    # filter names:
                    # hltSelectedPFTau30LooseChargedIsolationL1HLTMatched
                    # hltOverlapFilterIsoEle24WPTightGsfLooseIsoPFTau30
                    trigger_bits=1024 + 256,
                ),
            ],
            applies_to_dataset=(lambda dataset_inst: dataset_inst.is_mc or dataset_inst.x.era >= "B"),
            tags={"cross_trigger", "cross_e_tau", "channel_e_tau"},
        ),        
        # ===>>> tauh-tauh 
        Trigger(
            name="HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg",
            id=15151,
            run_range=(None,317509), # cover up to 317509
            legs=[
                TriggerLeg(
                    pdg_id=15,
                    min_pt=40.0,
                    max_abseta=2.1,
                    # filter names:
                    # hltDoublePFTau35TrackPt1TightChargedIsolationAndTightOOSCPhotonsDz02Reg
                    trigger_bits=64,
                ),
                TriggerLeg(
                    pdg_id=15,
                    min_pt=40.0,
                    max_abseta=2.1,
                    # filter names:
                    # hltDoublePFTau35TrackPt1TightChargedIsolationAndTightOOSCPhotonsDz02Reg
                    trigger_bits=64,
                ),
            ],
            applies_to_dataset=(lambda dataset_inst: dataset_inst.is_data and dataset_inst.x.era <= "B"),
            tags={"cross_trigger", "cross_tau_tau", "channel_tau_tau"},
        ),
        Trigger(
            name="HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg",
            id=15152,
            run_range=(None,317509), # cover up to 317509
            legs=[
                TriggerLeg(
                    pdg_id=15,
                    min_pt=40.0,
                    max_abseta=2.1,
                    # filter names:
                    # hltDoublePFTau35TrackPt1MediumChargedIsolationAndTightOOSCPhotonsDz02Reg
                    trigger_bits=64,
                ),
                TriggerLeg(
                    pdg_id=15,
                    min_pt=45.0,
                    max_abseta=None,
                    # filter names:
                    # hltDoublePFTau35TrackPt1MediumChargedIsolationAndTightOOSCPhotonsDz02Reg
                    trigger_bits=64,
                ),
            ],
            applies_to_dataset=(lambda dataset_inst: dataset_inst.is_data and dataset_inst.x.era <= "B"),
            tags={"cross_trigger", "cross_tau_tau", "channel_tau_tau"},
        ),
        Trigger(
            name="HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg",
            id=15153,
            run_range=(None,317509), # cover up to 317509
            legs=[
                TriggerLeg(
                    pdg_id=15,
                    min_pt=40.0,
                    max_abseta=2.1,
                    # filter names:
                    # hltDoublePFTau35TrackPt1MediumChargedIsolationAndTightOOSCPhotonsDz02Reg
                    trigger_bits=64,
                ),
                TriggerLeg(
                    pdg_id=15,
                    min_pt=40.0,
                    max_abseta=2.1,
                    # filter names:
                    # hltDoublePFTau35TrackPt1MediumChargedIsolationAndTightOOSCPhotonsDz02Reg
                    trigger_bits=64,
                ),
            ],
            applies_to_dataset=(lambda dataset_inst: dataset_inst.is_data and dataset_inst.x.era <= "B"),
            tags={"cross_trigger", "cross_tau_tau", "channel_tau_tau"},
        ),
        #Trigger(
        #    name="HLT_DoubleTightChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg",
        #    id=15152,
        #    run_range=(None,317510), # cover up to 317509
        #    legs=[
        #        TriggerLeg(
        #            pdg_id=15,
        #            min_pt=45.0,
        #            max_abseta=None,
        #            # filter names:
        #            # hltDoublePFTau35TrackPt1MediumChargedIsolationAndTightOOSCPhotonsDz02Reg
        #            trigger_bits=64,
        #        ),
        #        TriggerLeg(
        #            pdg_id=15,
        #            min_pt=45.0,
        #            max_abseta=None,
        #            # filter names:
        #            # hltDoublePFTau35TrackPt1MediumChargedIsolationAndTightOOSCPhotonsDz02Reg
        #            trigger_bits=64,
        #        ),
        #    ],
        #    applies_to_dataset=(lambda dataset_inst: dataset_inst.is_data and dataset_inst.x.era <= "B"),
        #    tags={"cross_trigger", "cross_tau_tau", "channel_tau_tau"},
        #),
        Trigger(
            name="HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg",
            id=15154,
            run_range=(317508,None), # cover 317509
            legs=[
                TriggerLeg(
                    pdg_id=15,
                    min_pt=40.0,
                    max_abseta=2.1,
                    # filter names:
                    # hltDoublePFTau40TrackPt1MediumChargedIsolationAndTightOOSCPhotonsDz02Reg
                    trigger_bits=64,
                ),
                TriggerLeg(
                    pdg_id=15,
                    min_pt=40.0,
                    max_abseta=2.1,
                    # filter names:
                    # hltDoublePFTau40TrackPt1MediumChargedIsolationAndTightOOSCPhotonsDz02Reg
                    trigger_bits=64,
                ),
            ],
            tags={"cross_trigger", "cross_tau_tau", "channel_tau_tau"},
        ),
    ])



# ----------------------------------------------- #
#                   Run3 2022                     #
# ----------------------------------------------- #
def add_triggers_run3_2022(config: od.Config, postfix: str) -> None:
    """
    # https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauTrigger#Tau_Triggers_in_NanoAOD
    Adds all triggers to a *config*. For the conversion from filter names to trigger bits, see
    https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/NanoAOD/python/triggerObjects_cff.py.
    ** Tau Trigger: https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauTrigger#Trigger_Table_for_2022
    ** Electron Trigger: https://twiki.cern.ch/twiki/bin/view/CMS/EgHLTRunIIISummary
    ** Muon Trigger: https://twiki.cern.ch/twiki/bin/view/CMS/MuonHLT2022
    ** HiggsDNA : https://gitlab.cern.ch/dwinterb/HiggsDNA/-/blob/NewProduction_Run3-2022/higgs_dna/metaconditions/Era2022_ditau.json?ref_type=heads
    PreEE:
      /afs/cern.ch/work/g/gsaha/public/IPHC/Work/ColumnFlowAnalyses/CPinHToTauTau/yamls/HLTlog_2022PreEE.log
    PostEE:
      /afs/cern.ch/work/g/gsaha/public/IPHC/Work/ColumnFlowAnalyses/CPinHToTauTau/yamls/HLTlog_2022PostEE.log
    """
    config.x.triggers = od.UniqueObjectIndex(Trigger,[
        ## ===>>> single electron
        Trigger(
            name="HLT_Ele30_WPTight_Gsf",
            id=111000,
            legs=[
                TriggerLeg(
                    pdg_id=11,
                    min_pt=31.0,
                    min_pt_online=30.0,
                    max_abseta=2.1,
                    max_abseta_online=2.3,
                    # filter names:
                    # hltEle30WPTightGsfTrackIsoFilter
                    trigger_bits=2**1,  # 1e (WPTight) (bit 1)
                ),
            ],
            tags=["single_trigger", "single_e", "channel_e_tau"],
        ),
        ## ===>>> e-tauh
        Trigger(
            name="HLT_Ele24_eta2p1_WPTight_Gsf_LooseDeepTauPFTauHPS30_eta2p1_CrossL1",
            id=11151,
            legs=[
                TriggerLeg(
                    pdg_id=11,
                    min_pt=25.0, # 26 -> Imperial
                    min_pt_online=24.0,
                    max_abseta=2.1,
                    max_abseta_online=2.1,
                    # filter names:
                    # hltOverlapFilterIsoEle24IsoTau30WPTightGsfCaloJet5
                    # hltHpsOverlapFilterIsoEle24WPTightGsfLooseETauWPDeepTauPFTau30
                    trigger_bits=2**7,
                ),
                TriggerLeg(
                    pdg_id=15,
                    min_pt=35.0,
                    min_pt_online=30.0,
                    max_abseta=2.1,
                    max_abseta_online=2.1,
                    # filter names:
                    # hltOverlapFilterIsoEle24IsoTau30WPTightGsfCaloJet5
                    # hltHpsOverlapFilterIsoEle24WPTightGsfLooseETauWPDeepTauPFTau30
                    trigger_bits=2**12,
                ),
            ],
            tags={"cross_trigger", "cross_e_tau", "channel_e_tau"},
        ),
        ## ===>>> single muon
        # https://cmshltinfo.app.cern.ch/path/HLT_IsoMu24_v#state=53332ee4-249a-4fd3-8152-5ddbcec178c6&session_state=4d32aa9a-e3a3-4947-8c1a-6d96d0a49833&code=90d18d3f-5685-4f64-899f-02aa42d33b2c.4d32aa9a-e3a3-4947-8c1a-6d96d0a49833.1363e04b-e180-4d83-92b3-3aca653d1d8d
        Trigger(
            name="HLT_IsoMu24",
            id=131000,
            legs=[
                TriggerLeg(
                    pdg_id=13,
                    min_pt=26.0,
                    min_pt_online=24.0,
                    max_abseta=2.4,
                    max_abseta_online=2.4,
                    # filter names:
                    # hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p08
                    trigger_bits=2**1 + 2**3,  # Iso (bit 1) + 1mu (bit 3)
                ),
            ],
            tags={"single_trigger", "single_mu", "channel_mu_tau"},
        ),
        ## ===>>> mu-tauh
        # https://cmshltinfo.app.cern.ch/path/HLT_IsoMu20_eta2p1_LooseDeepTauPFTauHPS27_eta2p1_CrossL1_v
        # mu : https://github.com/cms-sw/cmssw/blob/CMSSW_13_0_X/PhysicsTools/NanoAOD/python/triggerObjects_cff.py#L118C66-L118C74
        # tau : https://github.com/cms-sw/cmssw/blob/CMSSW_13_0_X/PhysicsTools/NanoAOD/python/triggerObjects_cff.py#L143
        # https://cms-nanoaod-integration.web.cern.ch/autoDoc/NanoAODv14/2024Prompt/doc_EGamma1_Run2024D-PromptReco-v1.html#TrigObj
        Trigger(
            name="HLT_IsoMu20_eta2p1_LooseDeepTauPFTauHPS27_eta2p1_CrossL1",
            id=13151,
            legs=[
                TriggerLeg(
                    pdg_id=13,
                    min_pt=21.0,
                    min_pt_online=20.0,
                    max_abseta=2.1,
                    max_abseta_online=2.1,
                    # filter names:
                    # 
                    trigger_bits=2**6,
                ),
                TriggerLeg(
                    pdg_id=15,
                    min_pt=32.0,
                    min_pt_online=27.0,
                    max_abseta=2.1,
                    max_abseta_online=2.1,
                    # filter names:
                    # hltHpsOverlapFilterIsoMu20LooseMuTauWPDeepTauPFTau27L1Seeded
                    trigger_bits=2**13,
                ),
            ],
            tags={"cross_trigger", "cross_mu_tau", "channel_mu_tau"},
        ),
        ## ===>>> tauh-tauh
        # https://cmshltinfo.app.cern.ch/path/HLT_DoubleMediumDeepTauPFTauHPS35_L2NN_eta2p1_v
        # https://github.com/cms-sw/cmssw/blob/CMSSW_13_0_X/PhysicsTools/NanoAOD/python/triggerObjects_cff.py#L141C73-L141C87
        # https://cms-nanoaod-integration.web.cern.ch/autoDoc/NanoAODv14/2024Prompt/doc_EGamma1_Run2024D-PromptReco-v1.html#TrigObj
        Trigger(
            name="HLT_DoubleMediumDeepTauPFTauHPS35_L2NN_eta2p1",
            id=15151,
            #run_range=[355862,362760],
            legs=[
                TriggerLeg(
                    pdg_id=15,
                    min_pt=40.0,
                    min_pt_online=35.0,
                    max_abseta=2.1,
                    max_abseta_online=2.1,
                    # filter names:
                    # hltHpsSelectedPFTausMediumDitauWPDeepTau
                    # hltHpsDoublePFTau35MediumDitauWPDeepTauL1HLTMatched
                    # 3 => DeepTau no spec WP, 11 => di-tau
                    trigger_bits=2**3 + 2**11,
               ),
                TriggerLeg(
                    pdg_id=15,
                    min_pt=40.0,
                    min_pt_online=35.0,
                    max_abseta=2.1,
                    max_abseta_online=2.1,
                    # filter names:
                    # hltHpsDoublePFTau35MediumDitauWPDeepTauDz02 (Deeptau + HPS)
                    # 3 => DeepTau no spec WP, 11 => di-tau 
                    trigger_bits=2**3 + 2**11,
                ),
            ],
            tags={"cross_trigger", "cross_tau_tau", "channel_tau_tau"},
        ),
        # ==> diTau + Jet
        # https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauTrigger#Trigger_Table_for_2022
        # https://github.com/cms-sw/cmssw/blob/CMSSW_13_0_X/PhysicsTools/NanoAOD/python/triggerObjects_cff.py#L148
        # https://cms-nanoaod-integration.web.cern.ch/autoDoc/NanoAODv14/2024Prompt/doc_EGamma1_Run2024D-PromptReco-v1.html#TrigObj
        # same as https://gitlab.cern.ch/dwinterb/HiggsDNA/-/blob/NewProduction_Run3-2022/higgs_dna/metaconditions/Era2022_ditau.json?ref_type=heads#L146
        Trigger(
            name="HLT_DoubleMediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet60",
            id=15152,
            #run_range=[355862,362760],
            legs=[
                TriggerLeg(
                    pdg_id=15,
                    min_pt=35.0,
                    min_pt_online=30.0,
                    max_abseta=2.1,
                    max_abseta_online=2.1,
                    # filter names:
                    # hltHpsOverlapFilterDeepTauDoublePFTau30PFJet60
                    # 14 => di-tau + PFJet
                    trigger_bits=2**3 + 2**14,
                ),
                TriggerLeg(
                    pdg_id=15,
                    min_pt=35.0,
                    min_pt_online=30.0,
                    max_abseta=2.1,
                    max_abseta_online=2.1,
                    # filter names:
                    # hltHpsOverlapFilterDeepTauDoublePFTau30PFJet60
                    # 14 => di-tau + PFJet
                    trigger_bits=2**3 + 2**14,
                ),
                TriggerLeg(
                    pdg_id=1,
                    min_pt=60.0,
                    min_pt_online=60.0,
                    max_abseta=4.7,
                    max_abseta_online=4.9,
                    # filter names:
                    # hltHpsOverlapFilterDeepTauDoublePFTau30PFJet60
                    # 14 => di-tau + PFJet
                    trigger_bits=2**17,
                ),                
            ],
            tags={"cross_trigger", "cross_tau_tau_jet", "channel_tau_tau"},
        ),
    ])


# ----------------------------------------------- #
#                   Run3 2023                     #
# ----------------------------------------------- #
def add_triggers_run3_2023(config: od.Config, postfix: str) -> None:
    """
    Adds all triggers to a *config*. For the conversion from filter names to trigger bits, see
    https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/NanoAOD/python/triggerObjects_cff.py.
    PreBPix:
    """
    config.x.triggers = od.UniqueObjectIndex(Trigger,[
        # ===>>> single electron
        Trigger(
            name="HLT_Ele30_WPTight_Gsf",
            id=111000,
            legs=[
                TriggerLeg(
                    pdg_id=11,
                    min_pt=31.0,
                    min_pt_online=30.0,
                    max_abseta=2.1,
                    max_abseta_online=2.3,
                    # filter names:
                    # 
                    trigger_bits=2**1,
                ),
            ],
            #applies_to_dataset=(lambda dataset_inst: dataset_inst.is_mc or dataset_inst.x.era >= "D"),
            tags={"single_trigger", "single_e", "channel_e_tau"},
        ),
        # ===>>> e-tauh
        Trigger(
            name="HLT_Ele24_eta2p1_WPTight_Gsf_LooseDeepTauPFTauHPS30_eta2p1_CrossL1",
            id=11151,
            legs=[
                TriggerLeg(
                    pdg_id=11,
                    min_pt=25.0,
                    min_pt_online=24.0,
                    max_abseta=2.1,
                    max_abseta_online=2.1,
                    # filter names:
                    # hltEle24erWPTightGsfTrackIsoFilterForTau
                    # hltOverlapFilterIsoEle24WPTightGsfLooseIsoPFTau30
                    trigger_bits=2**7,
                ),
                TriggerLeg(
                    pdg_id=15,
                    min_pt=35.0,
                    min_pt_online=30.0,
                    max_abseta=2.1,
                    max_abseta_online=2.1,
                    # filter names:
                    # hltSelectedPFTau30LooseChargedIsolationL1HLTMatched
                    # hltOverlapFilterIsoEle24WPTightGsfLooseIsoPFTau30
                    trigger_bits=2**12,
                ),
            ],
            tags={"cross_trigger", "cross_e_tau", "channel_e_tau"},
        ),
        # ===>>> single muon
        Trigger(
            name="HLT_IsoMu24",
            id=131000,
            legs=[
                TriggerLeg(
                    pdg_id=13,
                    min_pt=26.0,
                    min_pt_online=24.0,
                    max_abseta=2.4,
                    max_abseta_online=2.4,
                    # filter
                    # hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3trkIsoFiltered
                    # https://cms-nanoaod-integration.web.cern.ch/autoDoc/NanoAODv14/2024Prompt/doc_TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8_RunIII2024Summer24NanoAOD-140X_mcRun3_2024_realistic_v26-v2.html#TrigObj
                    # https://github.com/cms-sw/cmssw/blob/CMSSW_13_0_X/PhysicsTools/NanoAOD/python/triggerObjects_cff.py#L115
                    trigger_bits=2**1 + 2**3, # Iso (bit 1) + 1mu (bit 3)
                ),
            ],
            tags={"single_trigger", "single_mu", "channel_mu_tau"},
        ),
        # ===>>> mu-tauh
        Trigger(
            name="HLT_IsoMu20_eta2p1_LooseDeepTauPFTauHPS27_eta2p1_CrossL1",
            id=13151,
            legs=[
                TriggerLeg(
                    pdg_id=13,
                    min_pt=21.0,
                    min_pt_online=20.0,
                    max_abseta=2.1,
                    max_abseta_online=2.1,
                    # filter names:
                    # hltL3crIsoL1sMu18erTau24erIorMu20erTau24erL1f0L2f10QL3f20QL3trkIsoFiltered0p07
                    # hltOverlapFilterIsoMu20LooseChargedIsoPFTau27L1Seeded
                    trigger_bits=2**6,
                ),
                TriggerLeg(
                    pdg_id=15,
                    min_pt=32.0,
                    min_pt_online=27.0,
                    max_abseta=2.1,
                    max_abseta_online=2.1,
                    # filter names:
                    # hltSelectedPFTau27LooseChargedIsolationAgainstMuonL1HLTMatched or
                    # hltOverlapFilterIsoMu20LooseChargedIsoPFTau27L1Seeded
                    trigger_bits=2**13,
                ),
            ],
            tags={"cross_trigger", "cross_mu_tau", "channel_mu_tau"},
        ),
        # ===>>> tauh-tauh
        Trigger(
            name="HLT_DoubleMediumDeepTauPFTauHPS35_L2NN_eta2p1",
            id=15151,
            #run_range=[355862,362760],
            legs=[
                TriggerLeg(
                    pdg_id=15,
                    min_pt=40.0,
                    min_pt_online=35.0,
                    max_abseta=2.1,
                    max_abseta_online=2.1,
                    # filter names:
                    # hltHpsSelectedPFTausMediumDitauWPDeepTau
                    # hltHpsDoublePFTau35MediumDitauWPDeepTauL1HLTMatched
                    # 3 => DeepTau no spec WP, 11 => di-tau
                    trigger_bits=2**3 + 2**11,
               ),
                TriggerLeg(
                    pdg_id=15,
                    min_pt=40.0,
                    min_pt_online=35.0,
                    max_abseta=2.1,
                    max_abseta_online=2.1,
                    # filter names:
                    # hltHpsDoublePFTau35MediumDitauWPDeepTauDz02 (Deeptau + HPS)
                    # 3 => DeepTau no spec WP, 11 => di-tau 
                    trigger_bits=2**3 + 2**11,
                ),
            ],
            tags={"cross_trigger", "cross_tau_tau", "channel_tau_tau"},
        ),
        # ==> diTau + Jet
        # https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauTrigger#Trigger_Table_for_2022
        # https://github.com/cms-sw/cmssw/blob/CMSSW_13_0_X/PhysicsTools/NanoAOD/python/triggerObjects_cff.py#L148
        # https://cms-nanoaod-integration.web.cern.ch/autoDoc/NanoAODv14/2024Prompt/doc_EGamma1_Run2024D-PromptReco-v1.html#TrigObj
        # same as https://gitlab.cern.ch/dwinterb/HiggsDNA/-/blob/NewProduction_Run3-2022/higgs_dna/metaconditions/Era2022_ditau.json?ref_type=heads#L146
        Trigger(
            name="HLT_DoubleMediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet60",
            id=15152,
            #run_range=[355862,362760],
            legs=[
                TriggerLeg(
                    pdg_id=15,
                    min_pt=35.0,
                    min_pt_online=30.0,
                    max_abseta=2.1,
                    max_abseta_online=2.1,
                    # filter names:
                    # hltHpsOverlapFilterDeepTauDoublePFTau30PFJet60
                    # 14 => di-tau + PFJet
                    trigger_bits=2**3 + 2**14,
                ),
                TriggerLeg(
                    pdg_id=15,
                    min_pt=35.0,
                    min_pt_online=30.0,
                    max_abseta=2.1,
                    max_abseta_online=2.1,
                    # filter names:
                    # hltHpsOverlapFilterDeepTauDoublePFTau30PFJet60
                    # 14 => di-tau + PFJet
                    trigger_bits=2**3 + 2**14,
                ),
                TriggerLeg(
                    pdg_id=1,
                    min_pt=60.0,
                    min_pt_online=60.0,
                    max_abseta=4.7,
                    max_abseta_online=4.9,
                    # filter names:
                    # hltHpsOverlapFilterDeepTauDoublePFTau30PFJet60
                    # 14 => di-tau + PFJet
                    trigger_bits=2**17,
                ),
            ],
            tags={"cross_trigger", "cross_tau_tau_jet", "channel_tau_tau"},
        ),
    ])





    

# ----------------------------------------------- #
#                      EXTRAS                     #
# ----------------------------------------------- #
def add_triggers_run3_2022_tau_tau_postEE(config: od.Config) -> None:
    """
    Adds all triggers to a *config*. For the conversion from filter names to trigger bits, see
    https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/NanoAOD/python/triggerObjects_cff.py.

    HLT_DoubleMediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet60
    HLT_DoubleMediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet75
    HLT_DoubleMediumDeepTauPFTauHPS35_L2NN_eta2p1
    """
    config.x.triggers = od.UniqueObjectIndex(Trigger,[
        #
        # cross-tau
        #
        Trigger(
            name="HLT_DoubleMediumChargedIsoPFTauHPS40_Trk1_eta2p1",
            id=151,
            legs=[
                TriggerLeg(
                    pdg_id=15,
                    min_pt=40.0,
                    # filter names:
                    # hltDoublePFTau35TrackPt1TightChargedIsolationAndTightOOSCPhotonsDz02Reg
                    trigger_bits=64,
                ),
                TriggerLeg(
                    pdg_id=15,
                    min_pt=40.0,
                    # filter names:
                    # hltDoublePFTau35TrackPt1TightChargedIsolationAndTightOOSCPhotonsDz02Reg
                    trigger_bits=64,
                ),
            ],
            tags={"cross_trigger", "cross_tau_tau", "channel_tau_tau"},
        ),

        Trigger(
            name="HLT_DoubleTightChargedIsoPFTauHPS35_Trk1_eta2p1",
            id=152,
            legs=[
                TriggerLeg(
                    pdg_id=15,
                    min_pt=40.0,
                    # filter names:
                    # hltDoublePFTau40TrackPt1MediumChargedIsolationAndTightOOSCPhotonsDz02Reg
                    trigger_bits=64,
                ),
                TriggerLeg(
                    pdg_id=15,
                    min_pt=40.0,
                    # filter names:
                    # hltDoublePFTau40TrackPt1MediumChargedIsolationAndTightOOSCPhotonsDz02Reg
                    trigger_bits=64,
                ),
            ],
            tags={"cross_trigger", "cross_tau_tau", "channel_tau_tau"},
        ),
    ])


def add_triggers_run2_UL2017_mu_tau(config: od.Config) -> None:
    """
    Adds all triggers to a *config*. For the conversion from filter names to trigger bits, see
    https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/NanoAOD/python/triggerObjects_cff.py.
    """
    config.x.triggers = od.UniqueObjectIndex(Trigger,[
        
        #
        # single muon
        #
        Trigger(
            name="HLT_IsoMu24",
            id=131, #13 is for muon pdg_id, 1 because it's first muon trigger
            legs=[
                TriggerLeg(
                    pdg_id=13,
                    min_pt=25.0,
                    trigger_bits=2,
                ),
            ],
            tags={"single_trigger", "single_mu", "channel_mu_tau"},
        ),
    ])
