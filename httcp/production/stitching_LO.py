# coding: utf-8

"""
Process ID producer relevant for the stitching of the DY samples.
"""

import functools

import law

from columnflow.production import Producer, producer
from columnflow.util import maybe_import, InsertableDict
from columnflow.columnar_util import set_ak_column

from httcp.util import transverse_mass, IF_RUN2, IF_RUN3, IF_DATASET_IS_DY_M50,IF_DATASET_IS_W

np = maybe_import("numpy")
ak = maybe_import("awkward")
sp = maybe_import("scipy")
maybe_import("scipy.sparse")


logger = law.logger.get_logger(__name__)

NJetsRange = tuple[int, int]
PtRange = tuple[float, float]

set_ak_column_i64 = functools.partial(set_ak_column, value_type=np.int64)



# ################################## #
#     Stitching WJetBinned samples   #
# ################################## #
@producer(
    uses={
        IF_DATASET_IS_W("LHE.NpLO"), IF_DATASET_IS_W("LHE.Njets"),
    },
    produces={
        IF_DATASET_IS_W("process_id"),
    },
)
def process_ids_w(self: Producer, events: ak.Array, **kwargs) -> ak.Array:
    """
    Assigns each dy event a single process id, based on the number of jets and the di-lepton pt of
    the LHE record. This is used for the stitching of the DY samples.
    """
    # as always, we assume that each dataset has exactly one process associated to it
    if len(self.dataset_inst.processes) != 1:
        raise NotImplementedError(
            f"dataset {self.dataset_inst.name} has {len(self.dataset_inst.processes)} processes "
            "assigned, which is not yet implemented",
        )
    process_inst = self.dataset_inst.processes.get_first()

    # get the number of nlo jets and the di-lepton pt
    #njets = events.LHE.NpNLO
    njets = events.LHE.Njets

    outliers_mask = events.event < 0 # all False
    # raise a warning if a datasets was already created for a specific "bin" (leaf process),
    # but actually does not fit
    njets_range = process_inst.x("njets", None)
    if njets_range is not None:
        outliers = (njets < njets_range[0]) | (njets >= njets_range[1])
        outliers_mask = (outliers_mask | outliers)
        if ak.any(outliers):
            logger.warning(
                f"dataset {self.dataset_inst.name} is meant to contain njet values in the range "
                f"[{njets_range[0]}, {njets_range[0]}), but found {ak.sum(outliers)} events "
                "outside this range",
            )

    # lookup the id and check for invalid values
    process_ids = np.squeeze(np.asarray(self.id_table[self.key_func(njets)].todense()))
    invalid_mask = process_ids == 0
    if ak.any(invalid_mask):
        raise ValueError(
            f"found {sum(invalid_mask)} w events that could not be assigned to a process",
        )

    # store them
    events = set_ak_column_i64(events, "process_id", process_ids)

    return events, outliers_mask


@process_ids_w.setup
def process_ids_w_setup(
    self: Producer,
    reqs: dict,
    inputs: dict,
    reader_targets: InsertableDict,
) -> None:
    # define stitching ranges for the DY datasets covered by this producer's dy_inclusive_dataset
    stitching_ranges: list[NJetsRange] = []
    for proc in self.w_leaf_processes:
        logger.warning(f"proc : {proc}")
        njets = proc.x.njets
        stitching_ranges.append(njets)

    # sort by the first element of the ptll range
    sorted_stitching_ranges: list[NJetsRange] = [
        nj_range for nj_range in sorted(stitching_ranges, key=lambda nj_range: nj_range[0])
    ]

    # define a key function that maps njets and pt to a unique key for use in a lookup table
    def key_func(njets):
        # potentially convert single values into arrays
        single = False
        if isinstance(njets, int):
            njets = np.array([njets], dtype=np.int32)
            single = True

        # map into bins (index 0 means no binning)
        nj_bins = np.zeros(len(njets), dtype=np.int32)
        for nj_bin, nj_range in enumerate(sorted_stitching_ranges, 1):
            # nj_bin
            nj_mask = (nj_range[0] <= njets) & (njets < nj_range[1])
            nj_bins[nj_mask] = nj_bin

        return nj_bins[0] if single else nj_bins

    self.key_func = key_func

    # define the lookup table
    max_nj_bin = len(sorted_stitching_ranges)

    self.id_table = sp.sparse.lil_matrix((max_nj_bin + 1,1), dtype=np.int64)

    # fill it
    for proc in self.w_leaf_processes:
        key = key_func(proc.x.njets[0])
        self.id_table[key] = proc.id
    logger.info(f"id_table : \n{self.id_table}")


# ################################## #
#    Stitching DYJetBinned samples   #
# ################################## #
@producer(
    uses={
        IF_DATASET_IS_DY_M50("LHE.NpLO"), IF_DATASET_IS_DY_M50("LHE.Njets"),
    },
    produces={
        IF_DATASET_IS_DY_M50("process_id"),
    },
)
def process_ids_dy(self: Producer, events: ak.Array, **kwargs) -> ak.Array:
    """
    Assigns each dy event a single process id, based on the number of jets and the di-lepton pt of
    the LHE record. This is used for the stitching of the DY samples.
    """
    # as always, we assume that each dataset has exactly one process associated to it
    if len(self.dataset_inst.processes) != 1:
        raise NotImplementedError(
            f"dataset {self.dataset_inst.name} has {len(self.dataset_inst.processes)} processes "
            "assigned, which is not yet implemented",
        )
    process_inst = self.dataset_inst.processes.get_first()

    # get the number of nlo jets and the di-lepton pt
    njets = events.LHE.Njets

    outliers_mask = events.event < 0 # all False
    # raise a warning if a datasets was already created for a specific "bin" (leaf process),
    # but actually does not fit
    njets_range = process_inst.x("njets", None)
    if njets_range is not None:
        outliers = (njets < njets_range[0]) | (njets >= njets_range[1])
        outliers_mask = (outliers_mask | outliers)
        if ak.any(outliers):
            logger.warning(
                f"dataset {self.dataset_inst.name} is meant to contain njet values in the range "
                f"[{njets_range[0]}, {njets_range[0]}), but found {ak.sum(outliers)} events "
                "outside this range",
            )

    # lookup the id and check for invalid values
    process_ids = np.squeeze(np.asarray(self.id_table[self.key_func(njets)].todense()))
    invalid_mask = process_ids == 0
    if ak.any(invalid_mask):
        raise ValueError(
            f"found {sum(invalid_mask)} dy events that could not be assigned to a process",
        )

    # store them
    events = set_ak_column_i64(events, "process_id", process_ids)

    return events, outliers_mask


@process_ids_dy.setup
def process_ids_dy_setup(
    self: Producer,
    reqs: dict,
    inputs: dict,
    reader_targets: InsertableDict,
) -> None:
    # define stitching ranges for the DY datasets covered by this producer's dy_inclusive_dataset
    stitching_ranges: list[NJetsRange] = []
    for proc in self.dy_leaf_processes:
        logger.warning(f"proc : {proc}")
        njets = proc.x.njets
        stitching_ranges.append(njets)

    # sort by the first element of the ptll range
    sorted_stitching_ranges: list[NJetsRange] = [
        nj_range for nj_range in sorted(stitching_ranges, key=lambda nj_range: nj_range[0])
    ]

    # define a key function that maps njets and pt to a unique key for use in a lookup table
    def key_func(njets):
        # potentially convert single values into arrays
        single = False
        if isinstance(njets, int):
            njets = np.array([njets], dtype=np.int32)
            single = True

        # map into bins (index 0 means no binning)
        nj_bins = np.zeros(len(njets), dtype=np.int32)
        for nj_bin, nj_range in enumerate(sorted_stitching_ranges, 1):
            # nj_bin
            nj_mask = (nj_range[0] <= njets) & (njets < nj_range[1])
            nj_bins[nj_mask] = nj_bin

        return nj_bins[0] if single else nj_bins

    self.key_func = key_func

    # define the lookup table
    max_nj_bin = len(sorted_stitching_ranges)

    self.id_table = sp.sparse.lil_matrix((max_nj_bin + 1,1), dtype=np.int64)

    # fill it
    for proc in self.dy_leaf_processes:
        key = key_func(proc.x.njets[0])
        self.id_table[key] = proc.id
    #logger.info(f"id_table : \n{self.id_table}")
