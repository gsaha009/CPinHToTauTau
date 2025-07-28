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
        "cat__tautau__SR__nodeHiggs__pi_pi",
        config_category="tautau__real_1__hadD__nodeHiggs_tautau__pi_pi",
        config_variable="PhiCP_IPIP",
        config_data_datasets=["data_tau_C"],
        mc_stats=True,
    )
    self.add_category(
        "cat__tautau__SR__nodeHiggs__pi_rho",
        config_category="tautau__real_1__hadD__nodeHiggs_tautau__pi_rho",
        config_variable="PhiCP_IPDP",
        config_data_datasets=["data_tau_C"],
        mc_stats=True,
    )
    self.add_category(
        "cat__tautau__SR__nodeHiggs__pi_a1dm2",
        config_category="tautau__real_1__hadD__nodeHiggs_tautau__pi_a1dm2",
        config_variable="PhiCP_IPDP",
        config_data_datasets=["data_tau_C"],
        mc_stats=True,
    )
    self.add_category(
        "cat__tautau__SR__nodeHiggs__pi_a1dm10",
        config_category="tautau__real_1__hadD__nodeHiggs_tautau__pi_a1dm10",
        config_variable="PhiCP_IPDP",
        config_data_datasets=["data_tau_C"],
        mc_stats=True,
    )
    self.add_category(
        "cat__tautau__SR__nodeHiggs__rho_rho",
        config_category="tautau__real_1__hadD__nodeHiggs_tautau__rho_rho",
        config_variable="PhiCP_DPDP",
        config_data_datasets=["data_tau_C"],
        mc_stats=True,
    )
    self.add_category(
        "cat__tautau__SR__nodeHiggs__rho_a1dm2",
        config_category="tautau__real_1__hadD__nodeHiggs_tautau__rho_a1dm2",
        config_variable="PhiCP_DPDP",
        config_data_datasets=["data_tau_C"],
        mc_stats=True,
    )
    self.add_category(
        "cat__tautau__SR__nodeHiggs__rho_a1dm10",
        config_category="tautau__real_1__hadD__nodeHiggs_tautau__rho_a1dm10",
        config_variable="PhiCP_DPDP",
        config_data_datasets=["data_tau_C"],
        mc_stats=True,
    )
    self.add_category(
        "cat__tautau__SR__nodeHiggs__a1dm2_a1dm2",
        config_category="tautau__real_1__hadD__nodeHiggs_tautau__a1dm2_a1dm2",
        config_variable="PhiCP_DPDP",
        config_data_datasets=["data_tau_C"],
        mc_stats=True,
    )
    self.add_category(
        "cat__tautau__SR__nodeHiggs__a1dm2_a1dm10",
        config_category="tautau__real_1__hadD__nodeHiggs_tautau__a1dm2_a1dm10",
        config_variable="PhiCP_DPDP",
        config_data_datasets=["data_tau_C"],
        mc_stats=True,
    )
    self.add_category(
        "cat__tautau__SR__nodeHiggs__a1dm10_a1dm10",
        config_category="tautau__real_1__hadD__nodeHiggs_tautau__a1dm10_a1dm10",
        config_variable="PhiCP_DPDP",
        config_data_datasets=["data_tau_C"],
        mc_stats=True,
    )

    #
    # processes
    #

    self.add_process(
        "DY",
        config_process="dy",
        config_mc_datasets=["dy_lep_m*"],
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
    self.add_parameter(
        "tauspinner",
        process="h_ggf",
        type=ParameterType.shape,
        config_shift_source="tauspinner",
    )

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
