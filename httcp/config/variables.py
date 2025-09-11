### This config is used for listing the variables used in the analysis ###

from columnflow.config_util import add_category

import order as od

from columnflow.columnar_util import EMPTY_FLOAT,EMPTY_INT
from columnflow.util import DotDict, maybe_import
from columnflow.columnar_util import ColumnCollection

np = maybe_import("numpy")
ak = maybe_import("awkward")

def add_common_features(cfg: od.config) -> None:
    """
    Adds common features
    """
    cfg.add_variable(
        name="event",
        expression="event",
        binning=(1, 0.0, 1.0e6),
        x_title="Event number",
        discrete_x=True,
    )
    cfg.add_variable(
        name="run",
        expression="run",
        binning=(1, 100000.0, 500000.0),
        x_title="Run number",
        discrete_x=True,
    )
    cfg.add_variable(
        name="lumi",
        expression="luminosityBlock",
        binning=(1, 0.0, 5000.0),
        x_title="Luminosity block",
        discrete_x=True,
    )


def add_lepton_features(cfg: od.Config) -> None:
    """
    Adds lepton features only , ex electron_1_pt
    """
    for obj in ["Electron", "Muon", "Tau"]:
        for i in range(2):
            cfg.add_variable(
                name=f"{obj.lower()}_{i+1}_pt",
                expression=f"{obj}.pt[:,{i}]",
                null_value=EMPTY_FLOAT,
                binning=(40, 0., 200.),
                unit="GeV",
                x_title=obj + r" $p_{T}$",
            )
            cfg.add_variable(
                name=f"{obj.lower()}_{i+1}_phi",
                expression=f"{obj}.phi[:,{i}]",
                null_value=EMPTY_FLOAT,
                binning=(32, -3.2, 3.2),
                x_title=obj + r" $\phi$",
            )
            cfg.add_variable(
                name=f"{obj.lower()}_{i+1}_eta",
                expression=f"{obj}.eta[:,{i}]",
                null_value=EMPTY_FLOAT,
                binning=(25, -2.5, 2.5),
                x_title=obj + r" $\eta$",
            )
            cfg.add_variable(
                name=f"{obj.lower()}_{i+1}_IPsig",
                expression=f"{obj}.IPsig[:,{i}]",
                null_value=EMPTY_FLOAT,
                binning=(80, -10.0, 10.0) if obj == "Tau" else (40, 0.0, 10.0),
                x_title=obj + r" $IP Significance$",
            )
        cfg.add_variable(
            name=f"{obj.lower()}_mT",
            expression=f"{obj}.mT",
            null_value=EMPTY_FLOAT,
            binning=(40, 0.0, 200.0),
            unit="GeV",
            x_title=obj + r"$m_{T}$",
    )

def add_gen_features(cfg: od.Config) -> None:
    for i in range(2):
        cfg.add_variable(
            name=f"gentau_{i+1}_IPx",
            expression=f"GenTau.IPx[:,{i}]",
            null_value=EMPTY_FLOAT,
            binning=(60, -60.0, 60.0),
            unit="GeV",
            x_title=f"GenTau[{i+1}]" + r" $IPx$",
        )
        cfg.add_variable(
            name=f"gentau_{i+1}_IPy",
            expression=f"GenTau.IPy[:,{i}]",
            null_value=EMPTY_FLOAT,
            binning=(60, -60.0, 60.0),
            x_title=f"GenTau[{i+1}]" + r" $IPy$",
        )
        cfg.add_variable(
            name=f"gentau_{i+1}_IPz",
            expression=f"GenTau.IPz[:,{i}]",
            null_value=EMPTY_FLOAT,
            binning=(60, -60.0, 60.0),
            x_title=f"GenTau[{i+1}]" + r" $IPz$",
        )
        

def add_jet_features(cfg: od.Config) -> None:
    """
    Adds jet features only
    """
    cfg.add_variable(
        name="n_jet",
        expression="n_jet",
        #binning=(11, -0.5, 10.5),
        binning=(10, 0, 10),
        x_title="Number of jets",
        discrete_x=True,
    )
    cfg.add_variable(
        name="jets_pt",
        expression="Jet.pt",
        binning=(40, 0.0, 400.0),
        unit="GeV",
        x_title=r"$p_{T} of all jets$",
    )
    for i in range(2):
        cfg.add_variable(
            name=f"jet_{i+1}_pt",
            expression=f"Jet.pt[:,{i}]",
            null_value=EMPTY_FLOAT,
            binning=(50, 0.0, 200.0),
            unit="GeV",
            x_title=r"Jet $p_{T}$",
        )
        cfg.add_variable(
            name=f"jet_{i+1}_eta",
            expression=f"Jet.eta[:,{i}]",
            null_value=EMPTY_FLOAT,
            binning=(50, -5.0, 5.0),
            x_title=r"Jet $\eta$",
        )
        cfg.add_variable(
            name=f"jet_{i+1}_phi",
            expression=f"Jet.phi[:,{i}]",
            null_value=EMPTY_FLOAT,
            binning=(32, -3.2, 3.2),
            x_title=r"Jet $\phi$",
        )
    cfg.add_variable(
        name="hT",
        expression="hT",
        binning=(60, 20.0, 320.0),
        #null_value=EMPTY_FLOAT,
        unit="GeV",
        x_title="HT",
    )
    cfg.add_variable(
        name="hT_binvar",
        expression="hT",
        binning=[20, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 1000, 1500, 2500],
        #null_value=EMPTY_FLOAT,
        unit="GeV",
        x_title="HT",
    )
    cfg.add_variable(
        name="jet_raw_DeepJetFlavB",
        expression="Jet.btagDeepFlavB",
        null_value=EMPTY_FLOAT,
        binning=(30, 0,1),
        x_title=r"raw DeepJetFlawB",
    )


def add_highlevel_features(cfg: od.Config) -> None:    
    """
    Adds MET and other high-level features
    """
    cfg.add_variable(
        name="met",
        expression="MET.pt",
        null_value=EMPTY_FLOAT,
        binning=(20, 0.0, 100.0),
        x_title=r"MET",
    )
    cfg.add_variable(
        name="puppi_met_pt",
        expression="PuppiMET.pt",
        null_value=EMPTY_FLOAT,
        binning=(50, 0,100),
        unit="GeV",
        x_title=r"PuppiMET p$_{T}$",
    )
    cfg.add_variable(
        name="puppi_met_phi",
        expression="PuppiMET.phi",
        null_value=EMPTY_FLOAT,
        binning=(32, -3.2, 3.2),
        x_title=r"PuppiMET $\phi$",
    )
    cfg.add_variable(
        name="hcand_mass",
        expression="hcand_obj.mass",
        null_value=EMPTY_FLOAT,
        binning=(40, 0.0, 200.0),
        unit="GeV",
        x_title=r"$m_{vis}$",
    )
    # classifier scores
    cfg.add_variable(
        name="bdt_score_dy",
        expression="classifier_score[:,0]",
        null_value=EMPTY_FLOAT,
        binning=(25, 0.0, 1.0),
        unit="",
        x_title="BDT Score (DYTau prob)",
    )
    cfg.add_variable(
        name="bdt_score_higgs",
        expression="classifier_score[:,1]",
        null_value=EMPTY_FLOAT,
        binning=(25, 0.0, 1.0),
        unit="",
        x_title="BDT Score (Higgs prob)",
    )
    cfg.add_variable(
        name="bdt_score_higgs_binvar",
        expression="classifier_score[:,1]",
        null_value=EMPTY_FLOAT,
        binning=[0., 0.55, 0.65, 0.80, 0.90, 1.0],
        unit="",
        x_title="BDT Score (Higgs prob)",
    )
    cfg.add_variable(
        name="bdt_score_fake",
        expression="classifier_score[:,2]",
        null_value=EMPTY_FLOAT,
        binning=(25, 0.0, 1.0),
        unit="",
        x_title="BDT Score (Fake prob)",
    )
    

def add_weight_features(cfg: od.Config) -> None:
    """
    Adds weights
    """
    weights = ["mc", "pu",
               "muon_id", "muon_iso", "muon_IsoMu24_trigger", "muon_xtrig",
               "electron_idiso", "electron_Ele30_WPTight_trigger", "electron_xtrig",
               "tau", "tau_trigger"]

    for wt in weights:
        cfg.add_variable(
            name=f"{wt}_weight",
            expression=f"{wt}_weight",
            null_value=EMPTY_FLOAT,
            binning=(60, -3, 3),
            x_title=f"{wt} weight",
            tags={"mc_only"},
        )
        cfg.add_variable(
            name=f"{wt}_weight_log",
            expression=f"{wt}_weight",
            null_value=EMPTY_FLOAT,
            binning=list(np.logspace(-3, 3, 100)),
            log_x=True,
            x_title=f"{wt} weight",
            tags={"mc_only"},
        )

    

def add_cutflow_features(cfg: od.Config) -> None:
    """
    Adds cf features
    """
    cfg.add_variable(
        name="cf_jet1_pt",
        expression="cutflow.jet1_pt",
        binning=(40, 0.0, 400.0),
        unit="GeV",
        x_title=r"Jet 1 $p_{T}$",
    )


def add_hcand_features(cfg: od.Config) -> None:
    """
    Adds h lepton features only
    """
    for i in range(2):
        cfg.add_variable(
            name=f"hcand_{i+1}_pt",
            expression=f"hcand.pt[:,{i}]",
            null_value=EMPTY_FLOAT,
            binning=(40, 0., 200.),
            unit="GeV",
            #x_title=f"hcand[{i+1}]" + r" $p_{T}$",
            x_title=r"Leading tau p$_{T}$" if i == 0 else r"Subleading tau p$_{T}$",
        )
        cfg.add_variable(
            name=f"hcand_{i+1}_pt_MediumWP_binvar",
            expression=f"hcand.pt[:,{i}]",
            null_value=EMPTY_FLOAT,
            binning=[35,40,45,50,55,60,65,70,80,90,100,120,140,200],
            unit="GeV",
            x_title=f"hcand[{i+1}]" + r" $p_{T}$",
        )        
        cfg.add_variable(
            name=f"hcand_{i+1}_pt_VTightWP_binvar",
            expression=f"hcand.pt[:,{i}]",
            null_value=EMPTY_FLOAT,
            binning=[35,40,45,50,55,60,65,70,80,120,200],
            unit="GeV",
            x_title=f"hcand[{i+1}]" + r" $p_{T}$",
        )        
        cfg.add_variable(
            name=f"hcand_{i+1}_pt_VTightWP_binvar_v2",
            expression=f"hcand.pt[:,{i}]",
            null_value=EMPTY_FLOAT,
            binning=[35,40,45,50,55,60,65,70,80,100,140,200],
            unit="GeV",
            x_title=f"hcand[{i+1}]" + r" $p_{T}$",
        )
        cfg.add_variable(
            name=f"hcand_{i+1}_pt_fastMTT",
            expression=f"hcand.pt_fastMTT[:,{i}]",
            null_value=EMPTY_FLOAT,
            binning=(40, 0., 200.),
            unit="GeV",
            x_title=r"Leading tau p$_{T}$ (fastMTT)" if i == 0 else r"Subleading tau p$_{T}$ (fastMTT)",  
            #x_title=f"hcand[{i+1}]" + r" $p_{T}$ (fastMTT)",
        )
        cfg.add_variable(
            name=f"hcand_{i+1}_phi",
            expression=f"hcand.phi[:,{i}]",
            null_value=EMPTY_FLOAT,
            binning=(32, -3.2, 3.2),
            x_title=r"Leading tau $\phi$" if i == 0 else r"Subleading tau $\phi$",
            #x_title=f"hcand[{i+1}]" + r" $\phi$",
        )
        cfg.add_variable(
            name=f"hcand_{i+1}_phi_fastMTT",
            expression=f"hcand.phi_fastMTT[:,{i}]",
            null_value=EMPTY_FLOAT,
            binning=(32, -3.2, 3.2),
            x_title=f"hcand[{i+1}]" + r" $\phi$ (fastMTT)",
        )
        cfg.add_variable(
            name=f"hcand_{i+1}_eta",
            expression=f"hcand.eta[:,{i}]",
            null_value=EMPTY_FLOAT,
            binning=(25, -2.5, 2.5),
            x_title=r"Leading tau $\eta$" if i == 0 else r"Subleading tau $\eta$",
            #x_title=f"hcand[{i+1}]" + r" $\eta$",
        )
        cfg.add_variable(
            name=f"hcand_{i+1}_eta_fastMTT",
            expression=f"hcand.eta_fastMTT[:,{i}]",
            null_value=EMPTY_FLOAT,
            binning=(25, -2.5, 2.5),
            x_title=f"hcand[{i+1}]" + r" $\eta$ (fastMTT)",
        )
        cfg.add_variable(
            name=f"hcand_{i+1}_mass",
            expression=f"hcand.mass[:,{i}]",
            null_value=EMPTY_FLOAT,
            binning=(30, 0., 3.0),
            unit="GeV",
            x_title=r"Leading tau mass" if i == 0 else r"Subleading tau mass",
            #x_title=f"hcand[{i+1}]" + " mass",
        )
        cfg.add_variable(
            name=f"hcand_{i+1}_mass_fastMTT",
            expression=f"hcand.mass_fastMTT[:,{i}]",
            null_value=EMPTY_FLOAT,
            binning=(30, 0., 3.0),
            unit="GeV",
            x_title=f"hcand[{i+1}]" + " mass (fastMTT)",
        )
        cfg.add_variable(
            name=f"hcand_{i+1}_decayMode",
            expression=f"hcand.decayMode[:,{i}]",
            #null_value=EMPTY_INT,
            binning=(12, -0.5, 11.5),
            x_title=r"Leading tau DM" if i == 0 else r"Subleading tau DM",            
            #x_title=f"hcand[{i+1}]" + r" $DM (PNet)$",
        )
        cfg.add_variable(
            name=f"hcand_{i+1}_IPx",
            expression=f"hcand.IPx[:,{i}]",
            null_value=EMPTY_FLOAT,
            binning=(30, -0.015, 0.015),
            x_title=f"hcand[{i+1}]" + r" $IP_{x}$",
        )
        cfg.add_variable(
            name=f"hcand_{i+1}_IPy",
            expression=f"hcand.IPy[:,{i}]",
            null_value=EMPTY_FLOAT,
            binning=(30, -0.015, 0.015),
            x_title=f"hcand[{i+1}]" + r" $IP_{y}$",
        )
        cfg.add_variable(
            name=f"hcand_{i+1}_IPz",
            expression=f"hcand.IPz[:,{i}]",
            null_value=EMPTY_FLOAT,
            binning=(30, -0.015, 0.015),
            x_title=f"hcand[{i+1}]" + r" $IP_{z}$",
        )
        cfg.add_variable(
            name=f"hcand_{i+1}_IPsig",
            expression=f"hcand.IPsig[:,{i}]",
            null_value=EMPTY_FLOAT,
            binning=(40, 0.0, 10),
            x_title=f"hcand[{i+1}]" + r" $IP Significance$",
        )
        cfg.add_variable(
            name=f"dphi_met_h{i+1}",
            expression=f"dphi_met_h{i+1}",
            null_value=EMPTY_FLOAT,
            binning=(32, 0, 3.2),
            x_title=f"hcand[{i+1}]" + r" $-MET \Delta_{phi}$",
        )
        cfg.add_variable(
            name=f"met_var_qcd_h{i+1}",
            expression=f"met_var_qcd_h{i+1}",
            null_value=EMPTY_FLOAT,
            binning=(30, -1.5, 1.5),
            x_title=r"$MET var QCD$",
        )    

    cfg.add_variable(
        name="hcand_invm",
        expression="hcand_invm",
        null_value=EMPTY_FLOAT,
        binning=(50, 0.0, 400.0),
        #binning=(50, 0.0, 200.0),
        unit="GeV",
        x_title=r"$\tau_{h}\tau_{h}$ visible mass",
    )
    cfg.add_variable(
        name="hcand_invm_1bin",
        expression="hcand_invm",
        null_value=EMPTY_FLOAT,
        binning=(1, 0.0, 10000.0),
        unit="GeV",
        x_title=r"$visible mass$",
    )
    cfg.add_variable(
        name="hcand_invm_10GeV",
        expression="hcand_invm",
        null_value=EMPTY_FLOAT,
        binning=(40, 0.0, 400.0),
        unit="GeV",
        x_title=r"$visible mass$",
    )
    cfg.add_variable(
        name="hcand_invm_fastMTT",
        expression="hcand_invm_fastMTT",
        null_value=EMPTY_FLOAT,
        binning=(50, 0.0, 400.0),
        unit="GeV",
        x_title=r"$\tau\tau$ invariant mass (fastMTT)",
    )
    cfg.add_variable(
        name="hcand_invm_fastMTT_limit",
        expression="hcand_invm_fastMTT",
        null_value=EMPTY_FLOAT,
        binning=(4000, 110.0, 150.0),
        unit="GeV",
        x_title=r"$invariant mass (fastMTT)$",
    )
    cfg.add_variable(
        name="hcand_dr",
        expression="hcand_dr",
        null_value=EMPTY_FLOAT,
        binning=(40, 0.0, 5.0),
        x_title=r"$\Delta R (\tau_{h}^{1},\tau_{h}^{2})$",
    )
    cfg.add_variable(
        name="hcand_dphi",
        expression="hcand_dphi",
        null_value=EMPTY_FLOAT,
        binning=(32, 0.0, 3.2),
        x_title=r"$\Delta \phi(h1,h2)$",
    )
    cfg.add_variable(
        name="pt_vis",
        expression="pt_vis",
        null_value=EMPTY_FLOAT,
        binning=(50, 0.0, 200.0),
        unit="GeV",
        x_title=r"$p_{T}$ (visible)",
    )
    cfg.add_variable(
        name="pt_tt",
        expression="pt_tt",
        null_value=EMPTY_FLOAT,
        binning=(50, 0.0, 200.0),
        unit="GeV",
        x_title=r"$p_{T}$ (transverse)",
    )
    cfg.add_variable(
        name="mt_1",
        expression="mt_1",
        null_value=EMPTY_FLOAT,
        binning=(50, 0.0, 200.0),
        unit="GeV",
        x_title=r"$m_{T} ~(\tau_{1}, E_{T})$",
    )
    cfg.add_variable(
        name="mt_2",
        expression="mt_2",
        null_value=EMPTY_FLOAT,
        binning=(50, 0.0, 200.0),
        unit="GeV",
        x_title=r"$m_{T} ~(\tau_{2}, E_{T})$",
    )
    cfg.add_variable(
        name="mt_lep",
        expression="mt_lep",
        null_value=EMPTY_FLOAT,
        binning=(50, 0.0, 200.0),
        unit="GeV",
        x_title=r"$m_{T} ~(\tau_{1}, \tau_{2})$",
    )
    cfg.add_variable(
        name="mt_tot",
        expression="mt_tot",
        null_value=EMPTY_FLOAT,
        binning=(50, 0.0, 200.0),
        unit="GeV",
        x_title=r"$m_{T} ~(\tau_{1}, \tau_{2}, E_{T})$",
    )
    
    # PhiCP - Det
    cfg.add_variable(
        name="PhiCP_IPIP",
        expression="PhiCP_IPIP",
        null_value=EMPTY_FLOAT,
        binning=(int(2*np.pi/(0.1*np.pi)), 0, 2*np.pi),
        x_title=r"$\Phi_{CP}^{IP-IP}$ (rad)",
    )
    cfg.add_variable(
        name="PhiCP_DPDP",
        expression="PhiCP_DPDP",
        null_value=EMPTY_FLOAT,
        binning=(int(2*np.pi/(0.1*np.pi)), 0, 2*np.pi),
        x_title=r"$\Phi_{CP}^{DP-DP}$ (rad)",
    )
    cfg.add_variable(
        name="PhiCP_PVPV",
        expression="PhiCP_PVPV",
        null_value=EMPTY_FLOAT,
        binning=(int(2*np.pi/(0.1*np.pi)), 0, 2*np.pi),
        x_title=r"$\Phi_{CP}^{PV-PV}$ (rad)",
    )
    cfg.add_variable(
        name="PhiCP_IPDP",
        expression="PhiCP_IPDP",
        null_value=EMPTY_FLOAT,
        binning=(int(2*np.pi/(0.1*np.pi)), 0, 2*np.pi),
        x_title=r"$\Phi_{CP}^{IP-DP}$ (rad)",
    )
    cfg.add_variable(
        name="PhiCP_IPPV",
        expression="PhiCP_IPPV",
        null_value=EMPTY_FLOAT,
        binning=(int(2*np.pi/(0.1*np.pi)), 0, 2*np.pi),
        x_title=r"$\Phi_{CP}^{IP-PV}$ (rad)",
    )
    # PhiCP - Gen
    cfg.add_variable(
        name="PhiCPGen_IPIP",
        expression="PhiCPGen_IPIP",
        null_value=EMPTY_FLOAT,
        binning=(int(2*np.pi/(0.1*np.pi)), 0, 2*np.pi),
        x_title=r"$\Phi_{CP}^{IP-IP}$ (rad) [Gen level]",
    )
    cfg.add_variable(
        name="PhiCPGen_DPDP",
        expression="PhiCPGen_DPDP",
        null_value=EMPTY_FLOAT,
        binning=(int(2*np.pi/(0.1*np.pi)), 0, 2*np.pi),
        x_title=r"$\Phi_{CP}^{DP-DP}$ (rad) [Gen level]",
    )
    cfg.add_variable(
        name="PhiCPGen_PVPV",
        expression="PhiCPGen_PVPV",
        null_value=EMPTY_FLOAT,
        binning=(int(2*np.pi/(0.1*np.pi)), 0, 2*np.pi),
        x_title=r"$\Phi_{CP}^{PV-PV}$ (rad) [Gen level]",
    )
    cfg.add_variable(
        name="PhiCPGen_IPDP",
        expression="PhiCPGen_IPDP",
        null_value=EMPTY_FLOAT,
        binning=(int(2*np.pi/(0.1*np.pi)), 0, 2*np.pi),
        x_title=r"$\Phi_{CP}^{IP-DP}$ (rad) [Gen level]",
    )
    cfg.add_variable(
        name="PhiCPGen_IPPV",
        expression="PhiCPGen_IPPV",
        null_value=EMPTY_FLOAT,
        binning=(int(2*np.pi/(0.1*np.pi)), 0, 2*np.pi),
        x_title=r"$\Phi_{CP}^{IP-PV}$ (rad) [Gen level]",
    )
    # alpha minus Det
    cfg.add_variable(
        name="Alpha",
        expression="alpha",
        null_value=EMPTY_FLOAT,
        binning=(20, 0, np.pi/2.0),
        x_title=r"$\alpha_{-}$",
    )
    # alpha minus Gen
    cfg.add_variable(
        name="AlphaGen",
        expression="alphaGen",
        null_value=EMPTY_FLOAT,
        binning=(20, 0, np.pi/2.0),
        x_title=r"$\alpha_{-} [GenLevel]$",
    )
    # conditional PhiCP
    # TODO:
    cfg.add_variable(
        name="PhiCP_IPIP_alpha_lt_piby4",
        expression="PhiCP_IPIP_alpha_lt_piby4",
        null_value=EMPTY_FLOAT,
        binning=(int(2*np.pi/(0.1*np.pi)), 0, 2*np.pi),
        x_title=r"$\Phi_{CP}^{IP-IP}$ (rad) [a- < pi/4]",
    )
    cfg.add_variable(
        name="PhiCP_IPIP_alpha_gt_piby4",
        expression="PhiCP_IPIP_alpha_gt_piby4",
        null_value=EMPTY_FLOAT,
        binning=(int(2*np.pi/(0.1*np.pi)), 0, 2*np.pi),
        x_title=r"$\Phi_{CP}^{IP-IP}$ (rad) [a- > pi/4]",
    )
    cfg.add_variable(
        name="PhiCP_IPDP_alpha_lt_piby4",
        expression="PhiCP_IPDP_alpha_lt_piby4",
        null_value=EMPTY_FLOAT,
        binning=(int(2*np.pi/(0.1*np.pi)), 0, 2*np.pi),
        x_title=r"$\Phi_{CP}^{IP-NP}$ (rad) [a- < pi/4]",
    )
    cfg.add_variable(
        name="PhiCP_IPDP_alpha_gt_piby4",
        expression="PhiCP_IPDP_alpha_gt_piby4",
        null_value=EMPTY_FLOAT,
        binning=(int(2*np.pi/(0.1*np.pi)), 0, 2*np.pi),
        x_title=r"$\Phi_{CP}^{IP-NP}$ (rad) [a- > pi/4]",
    )
    cfg.add_variable(
        name="PhiCPGen_IPIP_alpha_lt_piby4",
        expression="PhiCPGen_IPIP_alpha_lt_piby4",
        null_value=EMPTY_FLOAT,
        binning=(int(2*np.pi/(0.1*np.pi)), 0, 2*np.pi),
        x_title=r"$\Phi_{CP}^{IP-IP}$ (rad) [a- < pi/4] [Gen Level]",
    )
    cfg.add_variable(
        name="PhiCPGen_IPIP_alpha_gt_piby4",
        expression="PhiCPGen_IPIP_alpha_gt_piby4",
        null_value=EMPTY_FLOAT,
        binning=(int(2*np.pi/(0.1*np.pi)), 0, 2*np.pi),
        x_title=r"$\Phi_{CP}^{IP-IP}$ (rad) [a- > pi/4] [Gen Level]",
    )
    cfg.add_variable(
        name="PhiCPGen_IPDP_alpha_lt_piby4",
        expression="PhiCPGen_IPDP_alpha_lt_piby4",
        null_value=EMPTY_FLOAT,
        binning=(int(2*np.pi/(0.1*np.pi)), 0, 2*np.pi),
        x_title=r"$\Phi_{CP}^{IP-NP}$ (rad) [a- < pi/4] [Gen Level]",
    )
    cfg.add_variable(
        name="PhiCPGen_IPDP_alpha_gt_piby4",
        expression="PhiCPGen_IPDP_alpha_gt_piby4",
        null_value=EMPTY_FLOAT,
        binning=(int(2*np.pi/(0.1*np.pi)), 0, 2*np.pi),
        x_title=r"$\Phi_{CP}^{IP-NP}$ (rad) [a- > pi/4] [Gen Level]",
    )

    

def add_test_variables(cfg: od.Config) -> None:
        cfg.add_variable(
            name="tau_pt",
            expression="Tau.pt",
            null_value=EMPTY_FLOAT,
            binning=(20, 0, 100),
            unit="GeV",
            x_title=r"tau $p_{T}$ (TES)",
        )
        cfg.add_variable(
            name="tau_mass",
            expression="Tau.mass",
            null_value=EMPTY_FLOAT,
            binning=(30, 25, 85),
            unit="GeV",
            x_title=r"tau $M$ (TES)",
        )

        cfg.add_variable(
            name="tau_pt_no_tes",
            expression="Tau.pt_no_tes",
            null_value=EMPTY_FLOAT,
            binning=(20, 0, 100),
            unit="GeV",
            x_title=r"tau $p_{T}$ (no TES)",
        )
        cfg.add_variable(
            name="tau_mass_no_tes",
            expression="Tau.mass_no_tes",
            null_value=EMPTY_FLOAT,
            binning=(30, 25, 85),
            unit="GeV",
            x_title=r"tau $M$ (TES)",
        )

        cfg.add_variable(
            name="mutau_mass_no_tes",
            expression="mutau_mass_no_tes",
            null_value=EMPTY_FLOAT,
            binning=(40, 0.0, 200.0),
            unit="GeV",
            x_title=r"$m_{vis}$(no TES)",
        )
    
         #single bin variables for transfer factor calculation
        cfg.add_variable(
            name="muon_eta_1bin",
            expression="Muon.eta",
            null_value=EMPTY_FLOAT,
            binning=(1, -3.0, 3.0),
            x_title=r"muon $\eta$",
        )
        cfg.add_variable(
            name="muon_pt_1bin",
            expression="Muon.pt",
            null_value=EMPTY_FLOAT,
            binning=(1, 20.0, 80.0),
            unit="GeV",
            x_title=r"muon $p_{T}$",
        )
        cfg.add_variable(
            name="muon_phi_1bin",
            expression="Muon.phi",
            null_value=EMPTY_FLOAT,
            binning=(1, -3.14159, 3.14159),
            x_title=r"muon $\varphi$",
        )
        cfg.add_variable(
            name="mutau_mass_1bin",
            expression="mutau_mass",
            null_value=EMPTY_FLOAT,
            binning=(1, 0.0, 200.0),
            unit="GeV",
            x_title=r"$m_{vis}$",
        )
        

# ############################ #
#  main add_variables function #
# ############################ #
def add_variables(cfg: od.Config) -> None:
    """
    Adds all variables to a *config*.
    """
    add_common_features(cfg)
    add_lepton_features(cfg)
    add_jet_features(cfg)
    add_highlevel_features(cfg)
    add_hcand_features(cfg)
    add_weight_features(cfg)
    add_cutflow_features(cfg)
    add_test_variables(cfg)
    add_gen_features(cfg)
