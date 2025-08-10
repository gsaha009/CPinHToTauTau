# coding: utf-8

# To handle many categories:
# https://github.com/columnflow/columnflow/commit/3a104d633fa47a8efc789f7aba054ed967017347#diff-01da7ecfbc4b8bb83460821147201da604f9825e32cbd51e365bdfcbc7cb0912R61
# https://github.com/columnflow/columnflow/issues/547
# https://github.com/columnflow/columnflow/issues/559

import law
import order as od

from columnflow.config_util import add_category, create_category_combinations
from httcp.util import call_once_on_config

logger = law.logger.get_logger(__name__)


# ############################################################### #
# To create combinations of categories                            #
# the entire root -> leaf categories will be in the category_ids  #
# other combinatorics will be produced in createHistograms task,  #
# if mentioned                                                    #
# ############################################################### #
def name_fn(root_categories):
    catlist = [cat.name for cat in root_categories.values() if cat]
    catname = "__".join(cat.name for cat in root_categories.values() if cat)
    return catname

def kwargs_fn(root_categories):
    return {
        "id": sum([c.id for c in root_categories.values()]),
        #"label": ",".join([c.label for c in root_categories.values()]),
        "label": "+".join([c.label for c in root_categories.values()]),
        "tags": set.union(*[cat.tags for cat in root_categories.values() if cat]),
    }

@call_once_on_config()
def add_common_categories(config: od.Config) -> None:
    """
    add_category(config,
                 name="incl",
                 id=9000000,
                 selection="cat_incl",
                 label=r"$Inclusive$",
                 tags={"incl"})
    """
    add_category(config,
                 name="etau",
                 id=10000000,
                 selection="cat_etau",
                 label=r"$e\tau_{h}$",
                 tags={"etau"})
    add_category(config,
                 name="mutau",
                 id=20000000,
                 selection="cat_mutau",
                 label=r"$\mu\tau_{h}$",
                 tags={"mutau"})
    add_category(config,
                 name="tautau",
                 id=40000000,
                 selection="cat_tautau",
                 label=r"$\tau_{h} \tau_{h}$",
                 tags={"tautau"})

@call_once_on_config()
def add_RealOrFake_categories(config: od.Config) -> None:
    #add_category(config, name="real_1", id=1000000, selection="cat_real_1", label=r"$t_1$", tags={"tau1isRealMC"})
    #add_category(config, name="fake_1", id=2000000, selection="cat_fake_1", label=r"$f_1$", tags={"tau1isFakeMC"})
    #add_category(config, name="real_2", id=3000000, selection="cat_real_2", label=r"$t_2$", tags={"tau2isRealMC"})
    #add_category(config, name="fake_2", id=4000000, selection="cat_fake_2", label=r"$f_2$", tags={"tau2isFakeMC"})
    add_category(config, name="real_1", id=1000000, selection="cat_real_1", label="prompt",     tags={"tau1isRealMC"})
    add_category(config, name="fake_1", id=2000000, selection="cat_fake_1", label="non-prompt", tags={"tau1isFakeMC"})
    add_category(config, name="real_2", id=3000000, selection="cat_real_2", label="prompt",     tags={"tau2isRealMC"})
    add_category(config, name="fake_2", id=4000000, selection="cat_fake_2", label="non-prompt", tags={"tau2isFakeMC"})


    #add_category(config, name="ipsig_0to1_1",   id=57, selection="cat_ipsig_0to1_1",   label=r"0 < IP_significance < 1", tags={"tau1_ipsig_0to1"})
    #add_category(config, name="ipsig_1toany_1", id=59, selection="cat_ipsig_1toany_1", label=r"IP_significance > 1",     tags={"tau1_ipsig_1toany"})

    
@call_once_on_config()
def add_njet_categories(config: od.Config) -> None:
    add_category(config, name="has_0j", id=100000, selection="cat_0j", label=r"$0$ jet",           tags={"has0j"})
    add_category(config, name="has_1j", id=200000, selection="cat_1j", label=r"$1$ jet",           tags={"has1j"})
    add_category(config, name="has_2j", id=300000, selection="cat_2j", label=r"$\geq{2}$ jets",    tags={"has2j"})
    

@call_once_on_config()
def add_ABCD_categories(config: od.Config) -> None:
    # hardonic
    add_category(config,name="hadA",  id=1000,  selection="cat_ss_iso1_iso2_bveto",       label="A",   tags={"ss","iso1",   "iso2",   "bveto"})
    add_category(config,name="hadB",  id=3000,  selection="cat_ss_noniso1_iso2_bveto",    label="B",   tags={"ss","noniso1","iso2",   "bveto"})
    add_category(config,name="hadA0", id=5000,  selection="cat_ss_iso1_noniso2_bveto",    label="A0",  tags={"ss","iso1",   "noniso2","bveto"})
    add_category(config,name="hadB0", id=7000,  selection="cat_ss_noniso1_noniso2_bveto", label="B0",  tags={"ss","noniso1","noniso2","bveto"})
    add_category(config,name="hadD0", id=9000,  selection="cat_os_iso1_noniso2_bveto",    label="D0",  tags={"os","iso1",   "noniso2","bveto"})
    add_category(config,name="hadC0", id=11000, selection="cat_os_noniso1_noniso2_bveto", label="C0",  tags={"os","noniso1","noniso2","bveto"})
    add_category(config,name="hadD",  id=13000, selection="cat_os_iso1_iso2_bveto",       label="Signal Region",  tags={"os","iso1",   "iso2",   "bveto"})
    add_category(config,name="hadC",  id=15000, selection="cat_os_noniso1_iso2_bveto",    label="Application Region",  tags={"os","noniso1","iso2",   "bveto"})
    # leptonic
    #add_category(config,name="lepA",  id=2000,  selection="cat_ss_iso2_bveto_lowmt",      label="A",   tags={"ss","iso2",   "bveto",  "lowmt" })
    #add_category(config,name="lepB",  id=4000,  selection="cat_ss_noniso2_bveto_lowmt",   label="B",   tags={"ss","noniso2","bveto",  "lowmt" })
    #add_category(config,name="lepA0", id=6000,  selection="cat_os_iso2_nobveto_lowmt",    label="A0",  tags={"os","iso2",   "nobveto","lowmt" })
    #add_category(config,name="lepB0", id=8000,  selection="cat_os_noniso2_nobveto_lowmt", label="B0",  tags={"os","noniso2","nobveto","lowmt" })
    #add_category(config,name="lepA1", id=10000, selection="cat_os_iso2_bveto_highmt",     label="A1",  tags={"os","iso2",   "bveto",  "highmt"})
    #add_category(config,name="lepB1", id=12000, selection="cat_os_noniso2_bveto_highmt",  label="B1",  tags={"os","noniso2","bveto",  "highmt"}) 
    #add_category(config,name="lepD",  id=14000, selection="cat_os_iso2_bveto_lowmt",      label="Signal Region",  tags={"os","iso2",   "bveto",  "lowmt" })
    #add_category(config,name="lepC",  id=16000, selection="cat_os_noniso2_bveto_lowmt",   label="Application Region",  tags={"os","noniso2","bveto",  "lowmt" })

    # DESY
    add_category(config,name="DRnum",  id=2000,  selection="cat_os_noniso1_iso2_lowmt",      label="dr_num",   tags={"os","noniso1", "iso2", "lowmt" })
    add_category(config,name="DRden",  id=4000,  selection="cat_ss_noniso1_iso2_lowmt",      label="dr_den",   tags={"ss","noniso1", "iso2", "lowmt" })
    add_category(config,name="AR",     id=6000,  selection="cat_ss_iso1_iso2_lowmt",         label="ar",       tags={"ss","iso1",    "iso2", "lowmt" })
    add_category(config,name="SR",     id=8000,  selection="cat_os_iso1_iso2_lowmt",         label="sr",       tags={"os","iso1",    "iso2", "lowmt" })
    

@call_once_on_config()
def add_classifier_categories(config: od.Config) -> None:
    # hardonic
    add_category(config,name="nodeDY_tautau",    id=100,  selection="cat_tautau_node_dy",      label="DY_node",    tags={"dy_node"})
    add_category(config,name="nodeFake_tautau",  id=300,  selection="cat_tautau_node_fake",    label="Fake_node",  tags={"fake_node"})
    add_category(config,name="nodeHiggs_tautau", id=500,  selection="cat_tautau_node_higgs",   label="Higgs_node", tags={"higgs_node"})    


    
@call_once_on_config()
def add_DM_categories(config: od.Config) -> None:
    # hadronic
    add_category(config, name="pi_1",          id=1,  selection="cat_pi_1",            label=r"$\tau_{h}^{1}\to\pi$",                         tags={"tau1pi"         })  # h1 -> pi
    add_category(config, name="rho_1",         id=3,  selection="cat_rho_1",           label=r"$\tau_{h}^{1}\to\rho$",                        tags={"tau1rho"        })  # h1 -> rho
    add_category(config, name="a1dm2_1",       id=5,  selection="cat_a1dm2_1",         label=r"$\tau_{h}^{1}\to a_{1}(1\pi-2\pi^{0})$",       tags={"tau1a1DM2"      })  # h1 -> a1
    add_category(config, name="a1dm10_1",      id=7,  selection="cat_a1dm10_1",        label=r"$\tau_{h}^{1}\to a_{1}(3\pi-0\pi^{0})$",       tags={"tau1a1DM10"     })  # h1 -> a1
    add_category(config, name="a1dm11_1",      id=9,  selection="cat_a1dm11_1",        label=r"$\tau_{h}^{1}\to a_{1}(3\pi-1\pi^{0})$",       tags={"tau1a1DM11"     })  # h1 -> a1
    #add_category(config, name="pi_rho_1",      id=11, selection="cat_pi_rho_1",        label=r"$\pi/\rho$",                            tags={"tau1piorrho"    })  # h1 -> a1
    add_category(config, name="pi_pi",         id=11, selection="cat_pi_pi",           label=r"$\tau_{h}\to\pi-\tau_{h}\to\pi$",                                       tags={"pi_pi"        })  # 
    add_category(config, name="pi_rho",        id=13, selection="cat_pi_rho",          label=r"$\tau_{h}\to\pi-\tau_{h}\to\rho$",                                      tags={"pi_rho"       })  # 
    add_category(config, name="pi_a1dm2",      id=15, selection="cat_pi_a1dm2",        label=r"$\tau_{h}\to\pi-\tau_{h}\to a_{1}(1\pi-2\pi^{0})$",                     tags={"pi_a1DM2"     })  # 
    add_category(config, name="pi_a1dm10",     id=17, selection="cat_pi_a1dm10",       label=r"$\tau_{h}\to\pi-\tau_{h}\to a_{1}(3\pi-0\pi^{0})$",                     tags={"pi_a1DM10"    })  # 
    add_category(config, name="pi_a1dm11",     id=19, selection="cat_pi_a1dm11",       label=r"$\tau_{h}\to\pi-\tau_{h}\to a_{1}(3\pi-1\pi^{0})$",                     tags={"pi_a1DM11"    })  # 
    add_category(config, name="rho_rho",       id=21, selection="cat_rho_rho",         label=r"$\tau_{h}\to\rho-\tau_{h}\to\rho$",                                     tags={"rho_rho"      })  # 
    add_category(config, name="rho_a1dm2",     id=23, selection="cat_rho_a1dm2",       label=r"$\tau_{h}\to\rho-\tau_{h}\to a_{1}(1\pi-2\pi^{0})$",                    tags={"rho_a1DM2"    })  # 
    add_category(config, name="rho_a1dm10",    id=25, selection="cat_rho_a1dm10",      label=r"$\tau_{h}\to\rho-\tau_{h}\to a_{1}(3\pi-0\pi^{0})$",                    tags={"rho_a1DM10"   })  # 
    add_category(config, name="rho_a1dm11",    id=27, selection="cat_rho_a1dm11",      label=r"$\tau_{h}\to\rho-\tau_{h}\to a_{1}(3\pi-1\pi^{0})$",                    tags={"rho_a1DM11"   })  # 
    add_category(config, name="a1dm2_a1dm2",   id=29, selection="cat_a1dm2_a1dm2",     label=r"$\tau_{h}\to a_{1}(1\pi-2\pi^{0})-\tau_{h}\to a_{1}(1\pi-2\pi^{0})$",   tags={"a1DM2_a1DM2"  })  # 
    add_category(config, name="a1dm2_a1dm10",  id=31, selection="cat_a1dm2_a1dm10",    label=r"$\tau_{h}\to a_{1}(1\pi-2\pi^{0})-\tau_{h}\to a_{1}(3\pi-0\pi^{0})$",   tags={"a1DM2_a1DM10" })  # 
    add_category(config, name="a1dm2_a1dm11",  id=33, selection="cat_a1dm2_a1dm11",    label=r"$\tau_{h}\to a_{1}(1\pi-2\pi^{0})-\tau_{h}\to a_{1}(3\pi-1\pi^{0})$",   tags={"a1DM2_a1DM11" })  # 
    add_category(config, name="a1dm10_a1dm10", id=35, selection="cat_a1dm10_a1dm10",   label=r"$\tau_{h}\to a_{1}(3\pi-0\pi^{0})-\tau_{h}\to a_{1}(3\pi-0\pi^{0})$",   tags={"a1DM10_a1DM10"})  # 
    add_category(config, name="a1dm10_a1dm11", id=37, selection="cat_a1dm10_a1dm11",   label=r"$\tau_{h}\to a_{1}(3\pi-0\pi^{0})-\tau_{h}\to a_{1}(3\pi-1\pi^{0})$",   tags={"a1DM10_a1DM11"})  # 
    add_category(config, name="a1dm11_a1dm11", id=39, selection="cat_a1dm11_a1dm11",   label=r"$\tau_{h}\to a_{1}(3\pi-1\pi^{0})-\tau_{h}\to a_{1}(3\pi-1\pi^{0})$",   tags={"a1DM11_a1DM11"})  # 
    # leptonic
    add_category(config, name="pi_2",          id=2,  selection="cat_pi_2",            label=r"$\tau_{h}\to\pi$",                                    tags={"tau2pi"    }) # h2 -> pi
    add_category(config, name="rho_2",         id=4,  selection="cat_rho_2",           label=r"$\tau_{h}\to\rho$",                                   tags={"tau2rho"   }) # h2 -> rho
    add_category(config, name="a1dm2_2",       id=6,  selection="cat_a1dm2_2",         label=r"$\tau_{h}\to a_{1}(1\pi-2\pi^{0})$",                  tags={"tau2a1DM2" }) # h2 -> a1
    add_category(config, name="a1dm10_2",      id=8,  selection="cat_a1dm10_2",        label=r"$\tau_{h}\to a_{1}(3\pi-0\pi^{0})$",                  tags={"tau2a1DM10"}) # h2 -> a1
    add_category(config, name="a1dm11_2",      id=10, selection="cat_a1dm11_2",        label=r"$\tau_{h}\to a_{1}(3\pi-1\pi^{0})$",                  tags={"tau2a1DM11"}) # h2 -> a1

    
@call_once_on_config()
def add_mutau_categories(config: od.Config) -> None:
    """
    categories = {
        "channel": [config.get_category("mutau")],
        "RorF"   : [config.get_category("real_2")],
        "abcd"   : [config.get_category("lepA"),  config.get_category("lepB"),
                    config.get_category("lepA0"), config.get_category("lepB0"),
                    config.get_category("lepA1"), config.get_category("lepB1"),
                    config.get_category("lepC"),  config.get_category("lepD")],
        "cp"     : [config.get_category("pi_2"),
                    config.get_category("rho_2"),
                    config.get_category("a1dm2_2"),
                    config.get_category("a1dm10_2"),                    
                    config.get_category("a1dm11_2")],
    }
    """
    categories = {
        "channel": [config.get_category("mutau")],
        "RorF"   : [config.get_category("real_2")],
        "abcd"   : [
            config.get_category("DRnum"),
            config.get_category("DRden"),
            config.get_category("AR"),
            config.get_category("SR")],
        "cp"     : [config.get_category("pi_2"),
                    config.get_category("rho_2"),
                    config.get_category("a1dm2_2"),
                    config.get_category("a1dm10_2"),                    
                    config.get_category("a1dm11_2")],
    }
    logger.info("mutau_categories")
    n = create_category_combinations(config,
                                     categories,
                                     name_fn=name_fn,
                                     kwargs_fn=kwargs_fn,
                                     skip_existing=False)
    logger.info(f"{n} categories have been created")



@call_once_on_config()
def add_etau_categories(config: od.Config) -> None:
    """
    categories = {
        "channel": [config.get_category("etau")],
        "RorF"   : [config.get_category("real_2")],
        "abcd"   : [config.get_category("lepA"),  config.get_category("lepB"),
                    config.get_category("lepA0"), config.get_category("lepB0"),
                    config.get_category("lepA1"), config.get_category("lepB1"),
                    config.get_category("lepC"),  config.get_category("lepD")],
        "cp"     : [config.get_category("pi_2"),
                    config.get_category("rho_2"),
                    config.get_category("a1dm2_2"),
                    config.get_category("a1dm10_2"),                    
                    config.get_category("a1dm11_2")],
    }
    """
    categories = {
        "channel": [config.get_category("etau")],
        "RorF"   : [config.get_category("real_2")],
        "abcd"   : [
            config.get_category("DRnum"),
            config.get_category("DRden"),
            config.get_category("AR"),
            config.get_category("SR")],
        "cp"     : [config.get_category("pi_2"),
                    config.get_category("rho_2"),
                    config.get_category("a1dm2_2"),
                    config.get_category("a1dm10_2"),                    
                    config.get_category("a1dm11_2")],
    }
    logger.info("etau_categories")
    n = create_category_combinations(config,
                                     categories,
                                     name_fn=name_fn,
                                     kwargs_fn=kwargs_fn,
                                     skip_existing=False)
    logger.info(f"{n} categories have been created")


    
"""
@call_once_on_config()
def add_tautau_categories(config: od.Config) -> None:
    categories = {
        "channel": [config.get_category("tautau")],
        "RorF"   : [config.get_category("real_1")],
        "abcd"   : [
            config.get_category("hadA"),  config.get_category("hadB"),
            config.get_category("hadA0"), config.get_category("hadB0"),
            config.get_category("hadC0"), config.get_category("hadD0"),
            config.get_category("hadC"),  config.get_category("hadD"),
        ],
        "nj"     : [
            config.get_category("has_0j"),
            config.get_category("has_1j"),
            config.get_category("has_2j"),
        ],        
        "cp"     : [
            config.get_category("pi_pi"),
            config.get_category("pi_rho"),
            config.get_category("pi_a1dm2"),
            config.get_category("pi_a1dm10"),
            config.get_category("pi_a1dm11"),
            config.get_category("rho_rho"),
            config.get_category("rho_a1dm2"),
            config.get_category("rho_a1dm10"),
            config.get_category("rho_a1dm11"),
            config.get_category("a1dm2_a1dm10"),
            config.get_category("a1dm2_a1dm11"),
            config.get_category("a1dm10_a1dm10"),
            config.get_category("a1dm10_a1dm11"),
            config.get_category("a1dm11_a1dm11"),
        ],
    }
    
    logger.info("tautau_categories")
    n = create_category_combinations(config,
                                     categories,
                                     name_fn=name_fn,
                                     kwargs_fn=kwargs_fn,
                                     skip_existing=False)
    logger.info(f"{n} categories have been created")
"""

    
@call_once_on_config()
def add_tautau_categories(config: od.Config) -> None:
    categories = {
        "channel": [config.get_category("tautau")],
        "RorF"   : [config.get_category("real_1")],
        "abcd"   : [
            config.get_category("hadA"),  config.get_category("hadB"),
            config.get_category("hadA0"), config.get_category("hadB0"),
            config.get_category("hadC0"), config.get_category("hadD0"),
            config.get_category("hadC"),  config.get_category("hadD"),
        ],
        "nj"     : [
            config.get_category("has_0j"),
            config.get_category("has_1j"),
            config.get_category("has_2j"),
        ],
        #"xgb"    : [
        #    config.get_category("nodeDY_tautau"),
        #    config.get_category("nodeFake_tautau"),
        #    config.get_category("nodeHiggs_tautau"),
        #],
        "dm"     : [
            config.get_category("pi_1"),
            config.get_category("rho_1"),
            config.get_category("a1dm2_1"),
            config.get_category("a1dm10_1"),
            #config.get_category("a1dm11_1"),
        ],
    }
    
    logger.info("tautau_categories")
    n = create_category_combinations(config,
                                     categories,
                                     name_fn=name_fn,
                                     kwargs_fn=kwargs_fn,
                                     skip_existing=False)
    logger.info(f"{n} categories have been created")


@call_once_on_config()
def add_tautau_phiCP_categories(config: od.Config) -> None:
    categories = {
        "channel": [
            config.get_category("tautau"),
        ],
        "RorF"   : [
            config.get_category("real_1"),
        ],
        "abcd"   : [
            config.get_category("hadC"),
            config.get_category("hadD"),
        ],
        "xgb"    : [
            config.get_category("nodeDY_tautau"),
            config.get_category("nodeFake_tautau"),
            config.get_category("nodeHiggs_tautau"),
        ],
        "cp"     : [
            config.get_category("pi_pi"),
            config.get_category("pi_rho"),
            config.get_category("pi_a1dm2"),
            config.get_category("pi_a1dm10"),
            #config.get_category("pi_a1dm11"),
            config.get_category("rho_rho"),
            config.get_category("rho_a1dm2"),
            config.get_category("rho_a1dm10"),
            #config.get_category("rho_a1dm11"),
            config.get_category("a1dm2_a1dm2"),            
            config.get_category("a1dm2_a1dm10"),
            #config.get_category("a1dm2_a1dm11"),
            config.get_category("a1dm10_a1dm10"),
            #config.get_category("a1dm10_a1dm11"),
            #config.get_category("a1dm11_a1dm11"),
        ],
    }
    
    logger.info("tautau_real_categories")
    n = create_category_combinations(config,
                                     categories,
                                     name_fn=name_fn,
                                     kwargs_fn=kwargs_fn,
                                     skip_existing=False)
    logger.info(f"{n} categories have been created")


# ################### #
# test categorization #
# ################### #
"""
@call_once_on_config()
def add_test_categories(config: od.config) -> None:
    categories = {
        "channel": [config.get_category("tautau")],
        "TorF"   : [config.get_category("real_1")],
        "abcd"   : [
            config.get_category("hadA"),   config.get_category("hadB"),
            config.get_category("hadA0"),  config.get_category("hadB0"),
            config.get_category("hadC0"),  config.get_category("hadD0"),
            config.get_category("hadC"),   config.get_category("hadD"),
        ],
        "IP"     : [config.get_category("ipsig_0to1_1"), config.get_category("ipsig_1toany_1")],
        "cp"     : [config.get_category("pi_rho_1")],
    }

    logger.info("test_categories")
    n = create_category_combinations(config,
                                     categories,
                                     name_fn=name_fn,
                                     kwargs_fn=kwargs_fn,
                                     skip_existing=True)
    logger.info(f"{n} categories have been created")
"""
    
# ################### #
# main categorization #
# ################### #

@call_once_on_config()    
def add_categories(config: od.Config) -> None:
    """
    Adds all categories to a *config*.
    """
    add_common_categories(config)
    
    add_njet_categories(config)
    add_RealOrFake_categories(config)
    
    add_ABCD_categories(config)
    add_DM_categories(config)

    #add_classifier_categories(config)
    
    add_mutau_categories(config)
    add_tautau_categories(config)
    ##add_tautau_real_categories(config)
    #add_tautau_phiCP_categories(config)
    ###add_test_categories(config)
    
    #all_cats = [cat.name for cat, _, _ in config.walk_categories()]
    #print(f"List of all cateogries in config: \n{all_cats}")
    """
    add_category(config, name="test0", id=9997, selection="cat_tautau_0j", label=r"test0", tags={"test0"})
    add_category(config, name="test1", id=9998, selection="cat_tautau_1j", label=r"test1", tags={"test1"})
    add_category(config, name="test2", id=9999, selection="cat_tautau_2j", label=r"test2", tags={"test2"})
    """
    #add_category(config, name="tautau_test11", id=9999, selection="cat_tautau_test11", label=r"test1", tags={"test1"})
