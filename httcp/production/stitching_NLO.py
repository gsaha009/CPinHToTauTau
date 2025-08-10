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
#    Stitching DYJetBinned samples   #
# ################################## #
#@producer(
#    uses={
#        IF_DATASET_IS_DY_M50("LHE.NpNLO"),
#    },
#    produces={
#        IF_DATASET_IS_DY_M50("process_id"),
#    },
#)
#def process_ids_dy(self: Producer, events: ak.Array, **kwargs) -> ak.Array:
#    """
#    Assigns each dy event a single process id, based on the number of jets and the di-lepton pt of
#    the LHE record. This is used for the stitching of the DY samples.
#    """
#    # as always, we assume that each dataset has exactly one process associated to it
#    if len(self.dataset_inst.processes) != 1:
#        raise NotImplementedError(
#            f"dataset {self.dataset_inst.name} has {len(self.dataset_inst.processes)} processes "
#            "assigned, which is not yet implemented",
#        )
#    process_inst = self.dataset_inst.processes.get_first()

#    # get the number of nlo jets and the di-lepton pt
#    njets = events.LHE.NpNLO

#    outliers_mask = events.event < 0 # all False
#    # raise a warning if a datasets was already created for a specific "bin" (leaf process),
#    # but actually does not fit
#    njets_range = process_inst.x("njets", None)
#    if njets_range is not None:
#        outliers = (njets < njets_range[0]) | (njets >= njets_range[1])
#        outliers_mask = (outliers_mask | outliers)
#        if ak.any(outliers):
#            logger.warning(
#                f"dataset {self.dataset_inst.name} is meant to contain njet values in the range "
#                f"[{njets_range[0]}, {njets_range[0]}), but found {ak.sum(outliers)} events "
#                "outside this range",
#            )

#    # lookup the id and check for invalid values
#    process_ids = np.squeeze(np.asarray(self.id_table[self.key_func(njets)].todense()))
#    invalid_mask = process_ids == 0
#    if ak.any(invalid_mask):
#        raise ValueError(
#            f"found {sum(invalid_mask)} dy events that could not be assigned to a process",
#        )

#    # store them
#    events = set_ak_column_i64(events, "process_id", process_ids)

#    return events, outliers_mask


#@process_ids_dy.setup
#def process_ids_dy_setup(
#    self: Producer,
#    reqs: dict,
#    inputs: dict,
#    reader_targets: InsertableDict,
#) -> None:
#    # define stitching ranges for the DY datasets covered by this producer's dy_inclusive_dataset
#    stitching_ranges: list[NJetsRange] = []
#    for proc in self.dy_leaf_processes:
#        logger.warning(f"proc : {proc}")
#        njets = proc.x.njets
#        stitching_ranges.append(njets)

#    # sort by the first element of the ptll range
#    sorted_stitching_ranges: list[NJetsRange] = [
#        nj_range for nj_range in sorted(stitching_ranges, key=lambda nj_range: nj_range[0])
#    ]

#    # define a key function that maps njets and pt to a unique key for use in a lookup table
#    def key_func(njets):
#        # potentially convert single values into arrays
#        single = False
#        if isinstance(njets, int):
#            njets = np.array([njets], dtype=np.int32)
#            single = True

#        # map into bins (index 0 means no binning)
#        nj_bins = np.zeros(len(njets), dtype=np.int32)
#        for nj_bin, nj_range in enumerate(sorted_stitching_ranges, 1):
#            # nj_bin
#            nj_mask = (nj_range[0] <= njets) & (njets < nj_range[1])
#            nj_bins[nj_mask] = nj_bin

#        return nj_bins[0] if single else nj_bins

#    self.key_func = key_func

#    # define the lookup table
#    max_nj_bin = len(sorted_stitching_ranges)

#    self.id_table = sp.sparse.lil_matrix((max_nj_bin + 1,1), dtype=np.int64)

#    # fill it
#    for proc in self.dy_leaf_processes:
#        key = key_func(proc.x.njets[0])
#        self.id_table[key] = proc.id
#    #logger.info(f"id_table : \n{self.id_table}")



# ################################## #
#      Stitching DYJets samples      #
# ################################## #

@producer(
    uses={
        IF_DATASET_IS_DY_M50("LHE.NpNLO", "LHE.Vpt"),
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
    #logger.warning(f"dataste : {self.dataset_inst}")
    if len(self.dataset_inst.processes) != 1:
        raise NotImplementedError(
            f"dataset {self.dataset_inst.name} has {len(self.dataset_inst.processes)} processes "
            "assigned, which is not yet implemented",
        )
    process_inst = self.dataset_inst.processes.get_first()

    #from IPython import embed; embed()
    # get the number of nlo jets and the di-lepton pt
    njets = events.LHE.NpNLO
    pt = events.LHE.Vpt

    
    outliers_mask = events.event < 0 # all False
    # raise a warning if a datasets was already created for a specific "bin" (leaf process)
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
    pt_range = process_inst.x("ptll", None)
    if pt_range is not None:
        outliers = (pt < pt_range[0]) | (pt >= pt_range[1])
        outliers_mask =	(outliers_mask | outliers)
        if ak.any(outliers):
            logger.warning(
                f"dataset {self.dataset_inst.name} is meant to contain ptll values in the range "
                f"[{pt_range[0]}, {pt_range[1]}), but found {ak.sum(outliers)} events outside this "
                "range",
            )
            
    # lookup the id and check for invalid values
    process_ids = np.squeeze(np.asarray(self.id_table[self.key_func(njets, pt)].todense()))

    # modifying the process id for outlier events
    if process_inst in self.dy_leaf_processes:
        if ak.sum(outliers_mask) > 0:
            # outlier : dy_lep_m50_1j_pt40to100_amcatnlo : branch 5, chunk 3, idx 5392 : pt = 100.0 GeV for 1 event only
            logger.warning(f"{self.dataset_inst.name} has {ak.sum(outliers_mask)} outlier events, and their process_ids are intentionally set to {process_inst.id} \nMake sure that the oulier events are vetoed in the selection task")
            process_ids = process_inst.id * ak.ones_like(process_ids)
            
    
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
    stitching_ranges: dict[NJetsRange, list[PtRange]] = {}
    for proc in self.dy_leaf_processes:
        njets = proc.x.njets
        stitching_ranges.setdefault(njets, [])
        if proc.has_aux("ptll"):
            stitching_ranges[njets].append(proc.x.ptll)

    # sort by the first element of the ptll range
    sorted_stitching_ranges: list[tuple[NJetsRange, list[PtRange]]] = [
        (nj_range, sorted(stitching_ranges[nj_range], key=lambda ptll_range: ptll_range[0]))
        for nj_range in sorted(stitching_ranges.keys(), key=lambda nj_range: nj_range[0])
    ]

    # define a key function that maps njets and pt to a unique key for use in a lookup table
    def key_func(njets, pt):
        # potentially convert single values into arrays
        single = False
        if isinstance(njets, int):
            assert isinstance(pt, (int, float))
            njets = np.array([njets], dtype=np.int32)
            pt = np.array([pt], dtype=np.float32)
            single = True

        # map into bins (index 0 means no binning)
        nj_bins = np.zeros(len(njets), dtype=np.int32)
        pt_bins = np.zeros(len(pt), dtype=np.int32)

        #from IPython import embed; embed()
        #logger.warning(f"sorted_stitching_ranges : {sorted_stitching_ranges}")
        for nj_bin, (nj_range, pt_ranges) in enumerate(sorted_stitching_ranges, 1):
            # nj_bin
            #logger.info(f"nj_bin : {nj_bin}")
            #logger.warning(f"nj_range : {nj_range}, pt_ranges : {pt_ranges}")
            nj_mask = (nj_range[0] <= njets) & (njets < nj_range[1])
            nj_bins[nj_mask] = nj_bin
            # pt_bin
            for pt_bin, (pt_min, pt_max) in enumerate(pt_ranges, 1):
                pt_mask = (pt_min <= pt) & (pt < pt_max)
                pt_bins[nj_mask & pt_mask] = pt_bin

        #from IPython import embed; embed()        
        return (nj_bins[0], pt_bins[0]) if single else (nj_bins, pt_bins)

    self.key_func = key_func

    # define the lookup table
    max_nj_bin = len(sorted_stitching_ranges)
    max_pt_bin = max(map(len, stitching_ranges.values()))
    self.id_table = sp.sparse.lil_matrix((max_nj_bin + 1, max_pt_bin + 1), dtype=np.int64)

    # fill it
    for proc in self.dy_leaf_processes:
        key = key_func(proc.x.njets[0], proc.x("ptll", [-1])[0])
        self.id_table[key] = proc.id
    #logger.info(f"id_table : \n{self.id_table}")
    
        


# ################################## #
#     Stitching WJetBinned samples   #
# ################################## #
@producer(
    uses={
        IF_DATASET_IS_W("LHE.NpNLO"),
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
    njets = events.LHE.NpNLO
    #njets = events.LHE.Njets

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
