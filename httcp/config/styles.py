# coding: utf-8

"""
Style definitions.
"""

import order as od

from columnflow.util import DotDict


def stylize_processes(config: od.Config) -> None:
    """
    Adds process colors and adjust labels.
    """
    cfg = config

    # recommended cms colors
    cfg.x.colors = DotDict(
        col_h_ggf_htt="#690301",
        col_h_ggf_htt_cpeven="#D72638",
        col_h_ggf_htt_cpmix="#1E56A0",
        col_h_ggf_htt_cpodd="#3CAEA3",
        col_vh_htt="#1b07f0",
        col_h_vbf_htt="#014d06",
        col_tt="#998ec3",
        col_st="#5ab4ac",
        col_vv="#5a3a1a",
        col_vvv="#b3987d",
        col_w="#b06eb5",#"#d78a7e",
        col_dy="#fec44f",
        col_dy_lep="#3690c0",
        col_dy_tau="#fec44f",
        col_dy_jet="#db9af5",
        col_dy_lm="#357591",
        col_h="#f768a1",
        light_blue="#bbd3f3",
        bright_blue="#0aa4f6",
        dark_blue="#10376c",
        light_purple="#cbabfa",
        bright_purple="#bf5ffa",
        purple="#761fc3",
        aubergine="#431f67",
        yellow="#ebe341",
        dark_yellow="#b1b111",
        bright_orange="#e96a0c",
        dark_orange="#8c4106",
        red="#FF0000",
        teal="#008080",
        grey="#808080",
        maroon="#800000",
        brown="#8f3007",
        black="#000000",
        fake="#fbb4ae",
    )

    if (p := config.get_process("h_ggf_htt", default=None)):
        p.color1 = cfg.x.colors.col_h_ggf_htt
        p.label = r"$H_{ggf} \rightarrow \tau\tau$"

    if (p := config.get_process("h_ggf_htt_cpeven", default=None)):
        p.color1 = cfg.x.colors.col_h_ggf_htt_cpeven
        p.label = r"$H_{ggf} \rightarrow \tau\tau$ (CP-even)"

    if (p := config.get_process("h_ggf_htt_cpmix", default=None)):
        p.color1 = cfg.x.colors.col_h_ggf_htt_cpmix
        p.label = r"$H_{ggf} \rightarrow \tau\tau$ (CP-mix)"

    if (p := config.get_process("h_ggf_htt_cpodd", default=None)):
        p.color1 = cfg.x.colors.col_h_ggf_htt_cpodd
        p.label = r"$H_{ggf} \rightarrow \tau\tau$ (CP-odd)"        
        
    if (p := config.get_process("h_vbf_htt", default=None)):
        p.color1 = cfg.x.colors.col_h_vbf_htt
        p.label = r"$H_{vbf} \rightarrow \tau\tau$"
        
    if (p := config.get_process("vh_htt", default=None)):
        p.color1 = cfg.x.colors.col_vh_htt
        p.label = r"$(V)H \rightarrow \tau\tau$"

    if (p := config.get_process("h", default=None)):
        p.color1 = cfg.x.colors.col_h
        p.label = r"$Higgs$"

    if (p := config.get_process("tt", default=None)):
        p.color1 = cfg.x.colors.col_tt
        p.label = r"$t\bar{t}$"

    if (p := config.get_process("st", default=None)):
        p.color1 = cfg.x.colors.col_st
        p.label = r"$Single ~t(W)$"

    if (p := config.get_process("top", default=None)):
        p.color1 = cfg.x.colors.col_st
        p.label = "Top"
        
        
    #if (p := config.get_process("dy", default=None)):
    #    p.color1 = cfg.x.colors.brown
    if (p := config.get_process("dy", default=None)):
        p.color1 = cfg.x.colors.col_dy
        p.label = r"$Z\to \ell^+ \ell^- (\ell \equiv e/\mu/\tau)$"

    if (p := config.get_process("dy_m50toinf", default=None)):
        p.color1 = cfg.x.colors.col_dy
        p.label = r"$Z\to \ell^+ \ell^- (\ell \equiv e/\mu/\tau)$"
        
    if (p := config.get_process("dy_m50toinf_lep", default=None)):
        p.color1 = cfg.x.colors.col_dy_lep
        p.label = r"$Z \to \ell^+ \ell^- (\ell \equiv e/\mu)$"
        
    if (p := config.get_process("dy_m50toinf_tau", default=None)):
        p.color1 = cfg.x.colors.col_dy_tau
        p.label = r"$Z \to \ell^+ \ell^- (\ell \equiv \tau)$"

    if (p := config.get_process("dy_m50toinf_jet", default=None)):
        p.color1 = cfg.x.colors.col_dy_jet
        p.label = r"$Z \to \ell^+ \ell^- (\ell \equiv j)$"
        
    if (p := config.get_process("dy_m10to50", default=None)):
        p.color1 = cfg.x.colors.col_dy_lm
        p.label = r"$Z\to \ell^+ \ell^- (m < 50)$"

    if (p := config.get_process("vv", default=None)):
        p.color1 = cfg.x.colors.col_vv

    if (p := config.get_process("vvv", default=None)):
        p.color1 = cfg.x.colors.col_vvv

    if (p := config.get_process("multiboson", default=None)):
        p.color1 = cfg.x.colors.col_vv

    if (p := config.get_process("w", default=None)):
        p.color1 = cfg.x.colors.col_w
        p.label = "W + jets"
        
    if (p := config.get_process("w_lnu", default=None)):
        p.color1 = cfg.x.colors.col_w
        p.label = "W + jets"

    if (p := config.get_process("ewk", default=None)):
        p.color1 = cfg.x.colors.brown

    if (p := config.get_process("ttv", default=None)):
        p.color1 = cfg.x.colors.grey
        p.label = r"$t\bar{t} + V$"

    if (p := config.get_process("ttvv", default=None)):
        p.color1 = cfg.x.colors.grey
        p.label = r"$t\bar{t} + VV$"

    if (p := config.get_process("tt_multiboson", default=None)):
        p.color1 = cfg.x.colors.grey

    if (p := config.get_process("qcd", default=None)):
        p.color1 = cfg.x.colors.fake
        p.label = r"$j \to \tau_h fakes$"
