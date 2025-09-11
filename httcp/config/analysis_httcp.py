# coding: utf-8

"""
Configuration of the CPinHToTauTau analysis.
"""
import law
import order as od

# ------------------------ #
# The main analysis object #
# ------------------------ #
analysis_httcp = ana = od.Analysis(
    name="analysis_httcp",
    id=1,
)

# analysis-global versions
# (see cfg.x.versions below for more info)
ana.x.versions = {}

# files of bash sandboxes that might be required by remote tasks
# (used in cf.HTCondorWorkflow)
ana.x.bash_sandboxes = ["$CF_BASE/sandboxes/cf.sh"]
default_sandbox = law.Sandbox.new(law.config.get("analysis", "default_columnar_sandbox"))
if default_sandbox.sandbox_type == "bash" and default_sandbox.name not in ana.x.bash_sandboxes:
    ana.x.bash_sandboxes.append(default_sandbox.name)

# files of cmssw sandboxes that might be required by remote tasks
# (used in cf.HTCondorWorkflow)
ana.x.cmssw_sandboxes = [
    "$CF_BASE/sandboxes/cmssw_default.sh",
]

# config groups for conveniently looping over certain configs
# (used in wrapper_factory)
ana.x.config_groups = {}

logger = law.logger.get_logger(__name__) 

# ------------- #
# setup configs #
# ------------- #

# ------------------------------------------------------------- #
#                               Run3                            #
# ------------------------------------------------------------- #

from httcp.config.config_run3 import add_config as add_config_run3

from cmsdb.campaigns.run3_2022_preEE_nano_cp_tau_v14_nanoprod_2024_v2 import campaign_run3_2022_preEE_nano_cp_tau_v14_nanoprod_2024_v2
from cmsdb.campaigns.run3_2022_postEE_nano_cp_tau_v14_nanoprod_2024_v2 import campaign_run3_2022_postEE_nano_cp_tau_v14_nanoprod_2024_v2
from cmsdb.campaigns.run3_2023_preBPix_nano_cp_tau_v14_nanoprod_2024_v2 import campaign_run3_2023_preBPix_nano_cp_tau_v14_nanoprod_2024_v2
from cmsdb.campaigns.run3_2023_postBPix_nano_cp_tau_v14_nanoprod_2024_v2 import campaign_run3_2023_postBPix_nano_cp_tau_v14_nanoprod_2024_v2

# ################################################################## #
# TO-DO --->>>                                                       #
# Modify this dictionary to choose campaigns                         #
# comment out or uncomment keys according to the requirements        #
# islimited to load the limited config, i.e. with one ROOT file      #
# and isfull for the full real config                                #
# ################################################################## #
campaign_dict = {
    "2022preEE"    : {"campaign" : campaign_run3_2022_preEE_nano_cp_tau_v14_nanoprod_2024_v2,    "islimited": False, "isfull": True},
    "2022postEE"   : {"campaign" : campaign_run3_2022_postEE_nano_cp_tau_v14_nanoprod_2024_v2,   "islimited": False, "isfull": True},
    "2023preBPix"  : {"campaign" : campaign_run3_2023_preBPix_nano_cp_tau_v14_nanoprod_2024_v2,  "islimited": False, "isfull": True},
    "2023postBPix" : {"campaign" : campaign_run3_2023_postBPix_nano_cp_tau_v14_nanoprod_2024_v2, "islimited": False, "isfull": True},
}

for key,val in campaign_dict.items():
    _campaign  = val["campaign"]
    _islimited = val["islimited"]
    _isfull    = val["isfull"]
    _name = _campaign.name
    _id   = _campaign.id
    if _isfull:
        logger.warning(f"Full Campaign for {key} : <{_name}> - <{_id}>")
        add_config_run3(
            analysis_httcp,
            _campaign.copy(),
            config_name=_name,
            config_id=int(_id))
    if _islimited:
        logger.warning(f"Limited Campaign for {key} : <{_name}> - <{_id+10}> - % only 1 root file per dataset will be considered %")
        add_config_run3(
            analysis_httcp,
            _campaign.copy(),
            config_name=f"{_name}_limited",
            config_id=int(_id)+10,
            limit_dataset_files=1)
        
