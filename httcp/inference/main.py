# coding: utf-8

"""
Example inference model.
"""

from columnflow.inference import inference_model, ParameterType, ParameterTransformation


@inference_model
def main(self):

    #
    # categories
    #
    self.add_category(
        f"cat__tautau__SR__nodeHiggs__pi_pi",
        config_category=f"tautau__real_1__hadD__nodeHiggs_tautau__pi_pi",
        config_variable="PhiCP_IPIP",
        config_data_datasets=["data_tau_C"],
        mc_stats=True,
    )
    self.add_category(
        f"cat__tautau__SR__nodeHiggs__pi_rho",
        config_category=f"tautau__real_1__hadD__nodeHiggs_tautau__pi_rho",
        config_variable="PhiCP_IPDP",
        config_data_datasets=["data_tau_C"],
        mc_stats=True,
    )
    self.add_category(
        f"cat__tautau__SR__nodeHiggs__pi_a1dm2",
        config_category=f"tautau__real_1__hadD__nodeHiggs_tautau__pi_a1dm2",
        config_variable="PhiCP_IPDP",
        config_data_datasets=["data_tau_C"],
        mc_stats=True,
    )
    self.add_category(
        f"cat__tautau__SR__nodeHiggs__pi_a1dm10",
        config_category=f"tautau__real_1__hadD__nodeHiggs_tautau__pi_a1dm10",
        config_variable="PhiCP_IPDP",
        config_data_datasets=["data_tau_C"],
        mc_stats=True,
    )
    self.add_category(
        f"cat__tautau__SR__nodeHiggs__rho_rho",
        config_category=f"tautau__real_1__hadD__nodeHiggs_tautau__rho_rho",
        config_variable="PhiCP_DPDP",
        config_data_datasets=["data_tau_C"],
        mc_stats=True,
    )
    self.add_category(
        f"cat__tautau__SR__nodeHiggs__rho_a1dm2",
        config_category=f"tautau__real_1__hadD__nodeHiggs_tautau__rho_a1dm2",
        config_variable="PhiCP_DPDP",
        config_data_datasets=["data_tau_C"],
        mc_stats=True,
    )
    self.add_category(
        f"cat__tautau__SR__nodeHiggs__rho_a1dm10",
        config_category=f"tautau__real_1__hadD__nodeHiggs_tautau__rho_a1dm10",
        config_variable="PhiCP_DPDP",
        config_data_datasets=["data_tau_C"],
        mc_stats=True,
    )
    self.add_category(
        f"cat__tautau__SR__nodeHiggs__a1dm2_a1dm2",
        config_category=f"tautau__real_1__hadD__nodeHiggs_tautau__a1dm2_a1dm2",
        config_variable="PhiCP_DPDP",
        config_data_datasets=["data_tau_C"],
        mc_stats=True,
    )
    self.add_category(
        f"cat__tautau__SR__nodeHiggs__a1dm2_a1dm10",
        config_category=f"tautau__real_1__hadD__nodeHiggs_tautau__a1dm2_a1dm10",
        config_variable="PhiCP_DPDP",
        config_data_datasets=["data_tau_C"],
        mc_stats=True,
    )
    self.add_category(
        f"cat__tautau__SR__nodeHiggs__a1dm10_a1dm10",
        config_category=f"tautau__real_1__hadD__nodeHiggs_tautau__a1dm10_a1dm10",
        config_variable="PhiCP_DPDP",
        config_data_datasets=["data_tau_C"],
        mc_stats=True,
    )

    #
    # processes
    #

    # backgrounds
    self.add_process(
        "DY",
        config_process="dy",
        config_mc_datasets=["dy_lep_m50*_amcatnlo",
                            "dy_2tau_m50_0j_amcatnlo",
                            "dy_2tau_m50_1j_amcatnlo",
                            "dy_2tau_m50_2j_amcatnlo"]
    )
    self.add_process(
        "WJ",
        config_process="w_lnu",
        config_mc_datasets=["wj_*_madgraph"]
    )
    self.add_process(
        "Top",
        config_process="top",
        config_mc_datasets=["tt_*","st_*"],
    )
    self.add_process(
        "VV",
        config_process="multiboson",
        config_mc_datasets=["ww","wz","zz","www","wwz","wzz","zzz"],
    )
    self.add_process(
        "QCD",
        config_process="qcd",
        config_mc_datasets=["qcd"],
    )
    # signals
    # -- ggf -- #
    self.add_process(
        "h_ggf",
        is_signal=True,
        config_process="h_ggf_htt",
        config_mc_datasets=["h_ggf_tautau_*"],
    )
    self.add_process(
        "h_ggf__tauspinnerUp",
        is_signal=True,
        config_process="h_ggf_htt",
        config_mc_datasets=["h_ggf_tautau_*"],
    )
    self.add_process(
        "h_ggf__tauspinnerDown",
        is_signal=True,
        config_process="h_ggf_htt",
        config_mc_datasets=["h_ggf_tautau_*"],
    )


    #
    # parameters
    #

    # groups
    #self.add_parameter_group("experiment")
    #self.add_parameter_group("theory")

    # lumi
    lumi = self.config_inst.x.luminosity
    for unc_name in lumi.uncertainties:
        self.add_parameter(
            unc_name,
            type=ParameterType.rate_gauss,
            effect=lumi.get(names=unc_name, direction=("down", "up"), factor=True),
            transformations=[ParameterTransformation.symmetrize],
        )

    # tune uncertainty
    #for proc in ["h_ggf", "h_vbf", "vh"]:
    self.add_parameter(
        "tauspinner",
        process="h_ggf",
        type=ParameterType.shape,
        config_shift_source="tauspinner",
    )
    #self.add_parameter(
    #    "minbias_xs",
    #    process=["*"],
    #    type=ParameterType.shape,
    #    config_shift_source="minbias_xs",
    #)
    #self.add_parameter(
    #    "tau",
    #    process=["*"],
    #    type=ParameterType.shape,
    #    config_shift_source="tau",
    #)
    #self.add_parameter(
    #    "tau_trig",
    #    process=["*"],
    #    type=ParameterType.shape,
    #    config_shift_source="tau_trig",
    #)

"""
@inference_model
def main_no_shapes(self):
    # same initialization as "example" above
    main.init_func.__get__(self, self.__class__)()

    #
    # remove all shape parameters
    #

    for category_name, process_name, parameter in self.iter_parameters():
        if parameter.type.is_shape or any(trafo.from_shape for trafo in parameter.transformations):
            self.remove_parameter(parameter.name, process=process_name, category=category_name)
"""
