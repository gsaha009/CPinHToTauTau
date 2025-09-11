# coding: utf-8

"""
Main categories file for the Higgs CP analysis
"""

from columnflow.categorization import Categorizer, categorizer
from columnflow.util import maybe_import

ak = maybe_import("awkward")


# ---------------------------------------------------------- #
#                          Channels                          #
# ---------------------------------------------------------- #

# inclusive
@categorizer(uses={"event"})
def cat_incl(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, ak.ones_like(events.event) == 1

# etau
@categorizer(uses={"channel_id"})
def cat_etau(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    ch = self.config_inst.get_channel("etau")
    return events, events["channel_id"] == ch.id

# mutau
@categorizer(uses={"channel_id"})
def cat_mutau(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    ch = self.config_inst.get_channel("mutau")
    return events, events["channel_id"] == ch.id

# tautau
@categorizer(uses={"channel_id"})
def cat_tautau(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    ch = self.config_inst.get_channel("tautau")
    return events, events["channel_id"] == ch.id


# ---------------------------------------------------------- #
#                          For nJets                         #
# ---------------------------------------------------------- #
@categorizer(uses={"has_0jet"})
def cat_0j(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.has_0jet

@categorizer(uses={"has_1jet"})
def cat_1j(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.has_1jet

@categorizer(uses={"has_2jet"})
def cat_2j(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.has_2jet

# ---------------------------------------------------------- #
#                          For PhiCP                         #
# ---------------------------------------------------------- #

# for tau-tau
# tau -> pi
@categorizer(uses={"is_pi_1"})
def cat_pi_1(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.is_pi_1
# tau -> rho
@categorizer(uses={"is_rho_1"})
def cat_rho_1(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.is_rho_1
#@categorizer(uses={"is_pi_1", "is_rho_1"})
#def cat_pi_rho_1(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
#    return events, (events.is_pi_1 | events.is_rho_1)
# tau -> a1 (DM2)
@categorizer(uses={"is_a1_1pr_2pi0_1"})
def cat_a1dm2_1(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.is_a1_1pr_2pi0_1
# tau -> a1 (DM10)
@categorizer(uses={"is_a1_3pr_0pi0_1"})
def cat_a1dm10_1(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.is_a1_3pr_0pi0_1
# tau -> a1 (DM11)
@categorizer(uses={"is_a1_3pr_1pi0_1"})
def cat_a1dm11_1(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.is_a1_3pr_1pi0_1


@categorizer(uses={"channel_id", "is_real_1", "has_1jet", "is_rho_1", "is_os", "is_iso_1", "is_iso_2", "is_b_veto"})
def cat_tautau_test11(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    ch = self.config_inst.get_channel("tautau")
    return events, (events["channel_id"] == ch.id) & events.is_real_1 & events.has_1jet & events.is_rho_1 & events.is_os & events.is_iso_1 & events.is_iso_2 & events.is_b_veto



# ---------- >>> for e/mu-tauh
# tau -> pi
@categorizer(uses={"is_pi_2"})
def cat_pi_2(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.is_pi_2
# tau -> rho
@categorizer(uses={"is_rho_2"})
def cat_rho_2(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.is_rho_2
# tau -> a1 (DM2)
@categorizer(uses={"is_a1_1pr_2pi0_2"})
def cat_a1dm2_2(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.is_a1_1pr_2pi0_2
# tau -> a1 (DM10)
@categorizer(uses={"is_a1_3pr_0pi0_2"})
def cat_a1dm10_2(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.is_a1_3pr_0pi0_2
# tau -> a1 (DM11)
@categorizer(uses={"is_a1_3pr_1pi0_2"})
def cat_a1dm11_2(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.is_a1_3pr_1pi0_2

# IPSig
@categorizer(uses={"is_ipsig_0to1_1"})
def cat_ipsig_0to1_1(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.is_ipsig_0to1_1
@categorizer(uses={"is_ipsig_0to1_1"})
def cat_ipsig_1toany_1(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, ~events.is_ipsig_0to1_1


# ----- >>>> to make less combinatorics for tautau
@categorizer(uses={"is_pi_1", "is_pi_2"})
def cat_pi_pi(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.is_pi_1 & events.is_pi_2
# tau -> pi, tau -> rho, vice versa
@categorizer(uses={"is_pi_1", "is_rho_2",
                   "is_pi_2", "is_rho_1"})
def cat_pi_rho(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, ((events.is_pi_1 & events.is_rho_2) | (events.is_pi_2 & events.is_rho_1))
# tau -> pi, tau -> a1 DM2, vice versa
@categorizer(uses={"is_pi_1", "is_a1_1pr_2pi0_2",
                   "is_pi_2", "is_a1_1pr_2pi0_1"})
def cat_pi_a1dm2(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, ((events.is_pi_1 & events.is_a1_1pr_2pi0_2) | (events.is_pi_2 & events.is_a1_1pr_2pi0_1))
# tau -> pi, tau -> a1 DM10, vice versa
@categorizer(uses={"is_pi_1", "is_a1_3pr_0pi0_2",
                   "is_pi_2", "is_a1_3pr_0pi0_1"})
def cat_pi_a1dm10(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, ((events.is_pi_1 & events.is_a1_3pr_0pi0_2) | (events.is_pi_2 & events.is_a1_3pr_0pi0_1))
# tau -> pi, tau -> a1 DM11, vice versa
@categorizer(uses={"is_pi_1", "is_a1_3pr_1pi0_2",
                   "is_pi_2", "is_a1_3pr_1pi0_1"})
def cat_pi_a1dm11(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, ((events.is_pi_1 & events.is_a1_3pr_1pi0_2) | (events.is_pi_2 & events.is_a1_3pr_1pi0_1))
# tau -> rho, tau -> rho
@categorizer(uses={"is_rho_1", "is_rho_2"})
def cat_rho_rho(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.is_rho_1 & events.is_rho_2
# tau -> rho, tau -> a1 DM2, vice versa
@categorizer(uses={"is_rho_1", "is_a1_1pr_2pi0_2",
                   "is_rho_2", "is_a1_1pr_2pi0_1"})
def cat_rho_a1dm2(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, ((events.is_rho_1 & events.is_a1_1pr_2pi0_2) | (events.is_rho_2 & events.is_a1_1pr_2pi0_1))
# tau -> rho, tau -> a1 DM10, vice versa
@categorizer(uses={"is_rho_1", "is_a1_3pr_0pi0_2",
                   "is_rho_2", "is_a1_3pr_0pi0_1"})
def cat_rho_a1dm10(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, ((events.is_rho_1 & events.is_a1_3pr_0pi0_2) | (events.is_rho_2 & events.is_a1_3pr_0pi0_1))
# tau -> rho, tau -> a1 DM11, vice versa
@categorizer(uses={"is_rho_1", "is_a1_3pr_1pi0_2",
                   "is_rho_2", "is_a1_3pr_1pi0_1"})
def cat_rho_a1dm11(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, ((events.is_rho_1 & events.is_a1_3pr_1pi0_2) | (events.is_rho_2 & events.is_a1_3pr_1pi0_1))
# tau -> a1 DM2, tau -> a1 DM2
@categorizer(uses={"is_a1_1pr_2pi0_1",
                   "is_a1_1pr_2pi0_2"})
def cat_a1dm2_a1dm2(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.is_a1_1pr_2pi0_1 & events.is_a1_1pr_2pi0_2
# tau -> a1 DM2, tau -> a1 DM10, vice versa
@categorizer(uses={"is_a1_1pr_2pi0_1", "is_a1_3pr_0pi0_2",
                   "is_a1_1pr_2pi0_2", "is_a1_3pr_0pi0_1"})
def cat_a1dm2_a1dm10(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, ((events.is_a1_1pr_2pi0_1 & events.is_a1_3pr_0pi0_2) | (events.is_a1_1pr_2pi0_2 & events.is_a1_3pr_0pi0_1))
# tau -> a1 DM2, tau -> a1 DM11, vice versa
@categorizer(uses={"is_a1_1pr_2pi0_1", "is_a1_3pr_1pi0_2",
                   "is_a1_1pr_2pi0_2", "is_a1_3pr_1pi0_1"})
def cat_a1dm2_a1dm11(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, ((events.is_a1_1pr_2pi0_1 & events.is_a1_3pr_1pi0_2) | (events.is_a1_1pr_2pi0_2 & events.is_a1_3pr_1pi0_1))
# tau -> a1 DM10, tau -> a1 DM10
@categorizer(uses={"is_a1_3pr_0pi0_1",
                   "is_a1_3pr_0pi0_2"})
def cat_a1dm10_a1dm10(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.is_a1_3pr_0pi0_1 & events.is_a1_3pr_0pi0_2
# tau -> a1 DM10, tau -> a1 DM11, vice versa
@categorizer(uses={"is_a1_3pr_0pi0_1", "is_a1_3pr_1pi0_2",
                   "is_a1_3pr_0pi0_2", "is_a1_3pr_1pi0_1"})
def cat_a1dm10_a1dm11(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, ((events.is_a1_3pr_0pi0_1 & events.is_a1_3pr_1pi0_2) | (events.is_a1_3pr_0pi0_2 & events.is_a1_3pr_1pi0_1))
# tau -> a1 DM11, tau -> a1 DM11
@categorizer(uses={"is_a1_3pr_1pi0_1", 
                   "is_a1_3pr_1pi0_2"})
def cat_a1dm11_a1dm11(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.is_a1_3pr_1pi0_1 & events.is_a1_3pr_1pi0_2


# to know true or fake tau
@categorizer(uses={"is_real_1"})
def cat_real_1(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.is_real_1
@categorizer(uses={"is_fake_1"})
def cat_fake_1(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.is_fake_1
@categorizer(uses={"is_real_2"})
def cat_real_2(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.is_real_2
@categorizer(uses={"is_fake_2"})
def cat_fake_2(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.is_fake_2

# TODO : need to fix this later
@categorizer(uses={"channel_id","is_real_1"})
def cat_tautau_real_1(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    ch = self.config_inst.get_channel("tautau")
    ch_mask = events["channel_id"] == ch.id
    return events, ch_mask & events.is_real_1
@categorizer(uses={"channel_id","is_fake_1"})
def cat_tautau_fake_1(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    ch = self.config_inst.get_channel("tautau")
    ch_mask = events["channel_id"] == ch.id
    return events, ch_mask & events.is_fake_1

@categorizer(uses={"channel_id","is_real_2"})
def cat_etau_real_2(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    ch = self.config_inst.get_channel("etau")
    ch_mask = events["channel_id"] == ch.id
    return events, ch_mask & events.is_real_2
@categorizer(uses={"channel_id","is_fake_2"})
def cat_etau_fake_2(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    ch = self.config_inst.get_channel("etau")
    ch_mask = events["channel_id"] == ch.id
    return events, ch_mask & events.is_fake_2

@categorizer(uses={"channel_id","is_real_2"})
def cat_mutau_real_2(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    ch = self.config_inst.get_channel("mutau")
    ch_mask = events["channel_id"] == ch.id
    return events, ch_mask & events.is_real_2
@categorizer(uses={"channel_id","is_fake_2"})
def cat_mutau_fake_2(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    ch = self.config_inst.get_channel("mutau")
    ch_mask = events["channel_id"] == ch.id
    return events, ch_mask & events.is_fake_2

# ---------------------------------------------------------- #
#                   CLASSIFIER [tau-tau]                     #
# ---------------------------------------------------------- #

# Tau or DY (Z to TauTau)
@categorizer(uses={"is_tautau_dy"})
def cat_tautau_node_dy(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.is_tautau_dy

# Fake
@categorizer(uses={"is_tautau_fake"})
def cat_tautau_node_fake(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.is_tautau_fake

# Higgs
@categorizer(uses={"is_tautau_higgs"})
def cat_tautau_node_higgs(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.is_tautau_higgs

@categorizer(uses={"is_tautau_higgs_bin_1"})
def cat_tautau_node_higgs_bin_1(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.is_tautau_higgs_bin_1

@categorizer(uses={"is_tautau_higgs_bin_2"})
def cat_tautau_node_higgs_bin_2(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.is_tautau_higgs_bin_2

@categorizer(uses={"is_tautau_higgs_bin_3"})
def cat_tautau_node_higgs_bin_3(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.is_tautau_higgs_bin_3

@categorizer(uses={"is_tautau_higgs_bin_4"})
def cat_tautau_node_higgs_bin_4(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.is_tautau_higgs_bin_4

@categorizer(uses={"is_tautau_higgs_bin_5"})
def cat_tautau_node_higgs_bin_5(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.is_tautau_higgs_bin_5

# ---------------------------------------------------------- #
#                          For ABCD                          #
# ---------------------------------------------------------- #


## --- tautau --->>>

#      tau2-iso          tau2-noniso            tau2-iso      
# << ------------ | --------------------- | -------------- >> 
# |------------------------------------------------------- | ^
# |               |           |           |                | ^
# |               |           |           |                | | tau-1
# |       A       |    A0     |     D0    |       D        | |  iso
# |               |           |           |                | |
# |               |           |           |                | |
# |------------------------------------------------------- | -
# |               |           |           |                | |
# |               |           |           |                | |
# |       B       |    B0     |     C0    |       C        | |  tau-1
# |               |           |           |                | | antiIso
# |               |           |           |                | |
# |------------------------------------------------------- | v
# << ------------------------ | ------------------------- >> v
#               SS                           OS                 
#
# A   : tautau [ss__iso1__iso2__nobjet]
# B   : tautau [ss__noniso1__iso2__nobjet]
# A0  : tautau [ss__iso1__noniso2__nobjet]
# B0  : tautau [ss__noniso1__noniso2__nobjet]
# D0  : tautau [os__iso1__noniso2__nobjet]
# C0  : tautau [os__noniso1__noniso2__nobjet]
# D   : tautau [os__iso1__iso2__nobjet]
# C   : tautau [os__noniso1__iso2__nobjet]

# ##################################################################### #
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ #
# W A R N I N G
# REMOVING b-veto to be on the same page of IMPERIAL
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ #
# ##################################################################### #

# A
@categorizer(uses={"is_os", "is_iso_1", "is_iso_2", "is_b_veto"})
def cat_ss_iso1_iso2_bveto(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, ~events.is_os & events.is_iso_1 & events.is_iso_2 #& events.is_b_veto
# B
@categorizer(uses={"is_os", "is_iso_1", "is_iso_2", "is_b_veto"})
def cat_ss_noniso1_iso2_bveto(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, ~events.is_os & ~events.is_iso_1 & events.is_iso_2 #& events.is_b_veto
# A0
@categorizer(uses={"is_os", "is_iso_1", "is_iso_2", "is_b_veto"})
def cat_ss_iso1_noniso2_bveto(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, ~events.is_os & events.is_iso_1 & ~events.is_iso_2 #& events.is_b_veto
# B0
@categorizer(uses={"is_os", "is_iso_1", "is_iso_2", "is_b_veto"})
def cat_ss_noniso1_noniso2_bveto(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, ~events.is_os & ~events.is_iso_1 & ~events.is_iso_2 #& events.is_b_veto
# D0
@categorizer(uses={"is_os", "is_iso_1", "is_iso_2", "is_b_veto"})
def cat_os_iso1_noniso2_bveto(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.is_os & events.is_iso_1 & ~events.is_iso_2 #& events.is_b_veto
# C0
@categorizer(uses={"is_os", "is_iso_1", "is_iso_2", "is_b_veto"})
def cat_os_noniso1_noniso2_bveto(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.is_os & ~events.is_iso_1 & ~events.is_iso_2 #& events.is_b_veto
# D
@categorizer(uses={"is_os", "is_iso_1", "is_iso_2", "is_b_veto"})
def cat_os_iso1_iso2_bveto(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.is_os & events.is_iso_1 & events.is_iso_2 #& events.is_b_veto
# C
@categorizer(uses={"is_os", "is_iso_1", "is_iso_2", "is_b_veto"})
def cat_os_noniso1_iso2_bveto(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.is_os & ~events.is_iso_1 & events.is_iso_2 #& events.is_b_veto


## --- e/mutau --->>>

#        bveto      no-bveto                 bveto                 
# << ------------ | --------- | -------------------------------- >>
# |-------------------------------------------------------------- | ^
# |               |           |           |                       | ^
# |               |           |           |                       | | tau-2
# |      QCD      |     T     |     SR    |          W            | |  iso
# |       A       |    A0     |     D     |          A1           | |
# |               |           |           |                       | |
# |-------------------------------------------------------------- | -
# |               |           |           |                       | |
# |               |           |           |                       | |
# |      QCD      |     T     |     SR    |          W            | |  tau-2
# |       B       |    B0     |     C     |          B1           | | antiIso
# |               |           |           |                       | |
# |-------------------------------------------------------------- | v
# << ------------ | -------------------------------------------- >> v
#        SS                           OS                           
# << ------------------------------------ | -------------------- >>
#                  mT < 50                          mT > 50        
#
# A   : e/mutau [ss__iso2__nobjet__lowmt]
# B   : e/mutau [ss__noniso2__nobjet__lowmt]
# A0  : e/mutau [os__iso2__bjet__lowmt]
# B0  : e/mutau [os__noniso2__bjet__lowmt]
# A1  : e/mutau [os__iso2__nobjet__highmt]
# B1  : e/mutau [os__noniso2__nobjet__highmt]
# D   : e/mutau [os__iso2__nobjet__lowmt]
# C   : e/mutau [os__noniso2__nobjet__lowmt]

# A
@categorizer(uses={"is_os", "is_iso_2", "is_b_veto", "is_low_mt"})
def cat_ss_iso2_bveto_lowmt(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, ~events.is_os & events.is_iso_2 & events.is_b_veto & events.is_low_mt 
# B
@categorizer(uses={"is_os", "is_iso_2", "is_b_veto", "is_low_mt"})
def cat_ss_noniso2_bveto_lowmt(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, ~events.is_os & ~events.is_iso_2 & events.is_b_veto & events.is_low_mt 
# A0
@categorizer(uses={"is_os", "is_iso_2", "is_b_veto", "is_low_mt"})
def cat_os_iso2_nobveto_lowmt(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.is_os & events.is_iso_2 & ~events.is_b_veto & events.is_low_mt 
# B0
@categorizer(uses={"is_os", "is_iso_2", "is_b_veto", "is_low_mt"})
def cat_os_noniso2_nobveto_lowmt(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.is_os & ~events.is_iso_2 & ~events.is_b_veto & events.is_low_mt 
# A1
@categorizer(uses={"is_os", "is_iso_2", "is_b_veto", "is_low_mt"})
def cat_os_iso2_bveto_highmt(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.is_os & events.is_iso_2 & events.is_b_veto & ~events.is_low_mt 
# B1
@categorizer(uses={"is_os", "is_iso_2", "is_b_veto", "is_low_mt"})
def cat_os_noniso2_bveto_highmt(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.is_os & ~events.is_iso_2 & events.is_b_veto & ~events.is_low_mt 
# D
@categorizer(uses={"is_os", "is_iso_2", "is_b_veto", "is_low_mt"})
def cat_os_iso2_bveto_lowmt(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.is_os & events.is_iso_2 & events.is_b_veto & events.is_low_mt 
# C
@categorizer(uses={"is_os", "is_iso_2", "is_b_veto", "is_low_mt"})
def cat_os_noniso2_bveto_lowmt(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.is_os & ~events.is_iso_2 & events.is_b_veto & events.is_low_mt


# DESY (Stepan)
# SR
@categorizer(uses={"is_os", "is_iso_1", "is_iso_2", "is_low_mt"})
def cat_os_iso1_iso2_lowmt(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.is_os & events.is_iso_1 & events.is_iso_2 & events.is_low_mt 
# AR
@categorizer(uses={"is_os", "is_iso_1", "is_iso_2", "is_low_mt"})
def cat_ss_iso1_iso2_lowmt(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, ~events.is_os & events.is_iso_1 & events.is_iso_2 & events.is_low_mt 
# DR_Num
@categorizer(uses={"is_os", "is_iso_1", "is_iso_2", "is_low_mt"})
def cat_os_noniso1_iso2_lowmt(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.is_os & ~events.is_iso_1 & events.is_iso_2 & events.is_low_mt 
# DR_Den
@categorizer(uses={"is_os", "is_iso_1", "is_iso_2", "is_low_mt"})
def cat_ss_noniso1_iso2_lowmt(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, ~events.is_os & ~events.is_iso_1 & events.is_iso_2 & events.is_low_mt





# test
@categorizer(uses={"channel_id", "has_0jet"})
def cat_tautau_0j(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, (events.channel_id == 4) & events.has_0jet
@categorizer(uses={"channel_id", "has_1jet"})
def cat_tautau_1j(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, (events.channel_id == 4) & events.has_1jet
@categorizer(uses={"channel_id", "has_2jet"})
def cat_tautau_2j(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, (events.channel_id == 4) & events.has_2jet
