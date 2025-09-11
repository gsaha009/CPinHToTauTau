# coding: utf-8

"""
Default event weight definitions.
"""
import law
from columnflow.weight import WeightProducer, weight_producer
from columnflow.columnar_util import Route
from columnflow.util import maybe_import, pattern_matcher

ak = maybe_import("awkward")
np = maybe_import("numpy")


logger = law.logger.get_logger(__name__)

@weight_producer(
    #uses={"channel_id", "is_os", "is_low_mt", "is_b_veto", "is_real_2", "is_fake_2", "is_iso_2", "hcand.*"},
    # both produced columns and dependent shifts are defined in init below
    # only run on mc
    mc_only=False,
    # options to keep or drop specific weights
    keep_weights=None,
    drop_weights=None,
)
def main(self: WeightProducer, events: ak.Array, **kwargs) -> ak.Array:
    # build the full event weight
    weight = ak.Array(np.ones(len(events), dtype=np.float32))
    logger.warning("The following weights will be applied for this dataset")
    for column in self.weight_columns:
        logger.info(column)
        weight = weight * Route(column).apply(events)

    logger.warning("And, these are the shifts")
    logger.info(self.shifts)

    return events, weight


@main.init
def main_init(self: WeightProducer) -> None:
    # use the config's auxiliary event_weights, drop some of them based on drop_weights, and on this
    # weight producer instance, store weight_columns, used columns, and shifts
    self.weight_columns = []

    # helpers to match to kept or dropped weights
    do_keep = pattern_matcher(self.keep_weights) if self.keep_weights else (lambda _: True)
    do_drop = pattern_matcher(self.drop_weights) if self.drop_weights else (lambda _: False)

    for weight_name in self.config_inst.x.event_weights:
        if not do_keep(weight_name) or do_drop(weight_name):
            continue

        # manually skip weights for samples that do not have lhe info
        if getattr(self, "dataset_inst", None) is not None:

            #if (weight_name != "ff_weight" or weight_name != "closure_weight") and self.dataset_inst.is_data:    
            #    continue
            if self.dataset_inst.is_data:
                if not weight_name in ["ff_weight","ff_cls_corr_weight","ff_ext_corr_weight"]:
                    continue

            # skip pdf weights for samples that dont have lhe weight
            is_lhe_weight = any(
                shift_inst.has_tag("pdf_weight")
                for shift_inst in self.config_inst.x.event_weights[weight_name]
            )
            if self.dataset_inst.has_tag("no_lhe_weights"):
                if weight_name in ["pdf_weight"] or is_lhe_weight:
                    continue

            # zpt weight only for DY samples
            is_zpt_reweight = any(
                shift_inst.has_tag("zpt_reweight")
                for shift_inst in self.config_inst.x.event_weights[weight_name]
            )
            if not self.dataset_inst.has_tag("is_dy"):
                if weight_name in ["zpt_reweight"] or is_zpt_reweight:
                    continue

            # top pt weight only for ttbar samples
            is_top_pt_reweight = any(
                shift_inst.has_tag("top_pt_weight")
                for shift_inst in self.config_inst.x.event_weights[weight_name]
            )
            if not self.dataset_inst.has_tag("is_tt"):
                if weight_name in ["top_pt_weight"]:
                    continue

            # tau-spinner weights are only for signal samples
            is_tauspinner_weight = any(
                shift_inst.has_tag("tauspinner_weight")
                for shift_inst in self.config_inst.x.event_weights[weight_name] 
            )
            if not (self.dataset_inst.has_tag("is_ggf_signal") or self.dataset_inst.has_tag("is_vh_signal") or self.dataset_inst.has_tag("is_vbf_signal")):
                if weight_name in ["tauspinner_weight"] or is_tauspinner_weight:
                    continue

        self.weight_columns.append(weight_name)
        self.uses.add(weight_name)
        self.shifts |= {
            shift_inst.name
            for shift_inst in self.config_inst.x.event_weights[weight_name]
        }


normalization_only = main.derive(
    "normalization_only",
    cls_dict={"keep_weights": "normalization_weight"},
)
