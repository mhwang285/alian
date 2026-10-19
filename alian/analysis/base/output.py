from collections import UserDict
from functools import singledispatchmethod
from pathlib import Path

import numpy as np
import ROOT
from RooUnfold import RooUnfoldResponse

from .logs import set_up_logger
from .utils import linbins, logbins, read_yaml

ROOT.TH1.AddDirectory(False)
ROOT.TH1.SetDefaultSumw2()
ROOT.TH2.SetDefaultSumw2()

registry = {
    ROOT.TH1F: [1, '1', '1F', '1f'],
    ROOT.TH2F: [2, '2', '2F', '2f'],
    ROOT.TH3F: [3, '3', '3F', '3f'],
    ROOT.TH1D: ['1D', '1d'],
    ROOT.TH2D: ['2D', '2d'],
    ROOT.TH3D: ['3D', '3d'],
    ROOT.TH1I: ['1I', '1i'],
    ROOT.TH2I: ['2I', '2i'],
    ROOT.TH3I: ['3I', '3i'],
}

mapping_double = {
    ROOT.TH1F: ROOT.TH2F,
    ROOT.TH1D: ROOT.TH2D,
    ROOT.TH1I: ROOT.TH2I
}

# reverse registry mapping to map key to specific ROOT class
HISTOGRAM_REGISTRY = {
    key: value
    for value, keys in registry.items()
    for key in keys
}

class Collection(UserDict):
    def __init__(self, *args, **kwargs):
        nsplit = kwargs.pop("nsplit")
        super().__init__(*args, **kwargs)
        if nsplit > 0:
            # give a sensible initial split so one can query histogram properties
            # prior to analysis
            self.split = "1_"
        else:
            self.split = ""

    def __getitem__(self, key):
        return super().__getitem__(f"{self.split}{key}")
    def set_split(self, split):
        self.split = split
    def get_direct(self, key):
        return self.data[key]
    def save(self, file):
        for obj in self.data.values():
            file.WriteTObject(obj)


class Output:
    def __init__(self, bins, histograms, branches, trees, responses, split = None):
        self.logger = set_up_logger(__name__)
        self._bins_init_done = False
        self.names = []
        self.logger.info(f"Output splitting: {split}")
        if split is None:
            self.split = split
            self.nsplit = 0
            self.split_ids = [""]
        elif isinstance(split, int):
            if split > 1:
                # split the dataset in a 1 : split - 1 ratio
                # for example, split = 5 will split the data into 1/5 and 4/5
                self.split = ["1_"] + ["2_"] * (split - 1)
                self.nsplit = split
                self.split_ids = ["1_", "2_"]
            else:
                raise ValueError(f"Output splitting parameter must be greater than 1 if int (was given: {split})")
        elif isinstance(split, list):
            if all(i > 0 for i in split):
                # split the dataset into the given chunks
                # for example, split = [1, 3, 3] will split the data into 1/7, 3/7 and 3/7
                self.split = []
                self.split_ids = []
                for idx, num in enumerate(split, start = 1):
                    self.split += [f"{idx}_"] * num
                    self.split_ids = [f"{idx}_"]
                self.nsplit = len(self.split)
            else:
                raise ValueError(f"Output splitting parameters must all greater than 0 if list (was given: {split})")
        else:
            raise TypeError(f"Could not parse splitting parameter of type {type(split)}, must be None, int or list")

        self.logger.info("Initializing histograms...", stacklevel = 2)
        self._init_bins(bins)
        self._init_histograms(histograms)
        self.logger.info("Histograms initialized.", stacklevel = 2)
        self.logger.info("Initializing trees...", stacklevel = 2)
        self._init_branches(branches)
        self._init_trees(trees)
        self.logger.info("Trees initialized.", stacklevel = 2)
        self.logger.info("Initializing response matrices...", stacklevel = 2)
        if not self._bins_init_done:
            self._init_bins(bins)
        self._init_responses(responses)
        self.logger.info("Response matrices initialized.", stacklevel = 2)

    @singledispatchmethod
    @classmethod
    def load(cls, file):
        raise NotImplementedError(f'Cannot configure output from type {type(file)}.')
    @load.register(dict)
    @classmethod
    def _load(cls, cfg):
        if "output" not in cfg:
            raise KeyError("Output must be configured in a 'output' block in the YAML configuration!")
        if cfg["output"] is None:
            raise KeyError("The 'output' block in the YAML configuration cannot be empty, remove or fill it!")

        cfg_output = {}
        fields = ["bins", "histograms", "branches", "trees", "responses"]
        for field in fields:
            if field not in cfg["output"] or cfg["output"][field] is None:
                cfg_output[field] = {}
            else:
                cfg_output[field] = cfg["output"][field]
        if "split" in cfg["output"]:
            cfg_output["split"] = cfg["output"]["split"]
        # if split not in config, let it be unset and use the default value

        return cls(**cfg_output)
    @load.register(str)
    @load.register(Path)
    @classmethod
    def _load(cls, file):
        cfg = read_yaml(file)
        return cls.load(cfg)

    def _init_bins(self, bins_cfg):
        self._bins = {}
        self._nbins = {}
        for name, arglist in bins_cfg.items():
            self._bins[name] = self._calculate_bins(*arglist)
            self._nbins[name] = len(self._bins[name]) - 1
        self._bins_init_done = True

    def _calculate_bins(self, binning, *bin_info):
        match binning:
            case 'logarithmic' | 'logarithm' | 'log' | 'lg':
                binmin, binmax, nbins = bin_info
                return logbins(binmin, binmax, nbins)
            case 'linear' | 'line' | 'lin' | 'ln':
                binmin, binmax, nbins = bin_info
                return linbins(binmin, binmax, nbins)
            case 'custom' | 'cst' | 'c':
                return np.array(bin_info[0], dtype = float)
            case _:
                raise ValueError(f"Binning type {binning} invalid.")

    def _init_histograms(self, hists_cfg):
        self.hists = Collection(nsplit = self.nsplit)
        # for now, doesn't support nested structures, but could in the future
        for split in self.split_ids:
            for tag, cfg in hists_cfg.items():
                htype, name, title, *bin_names = cfg
                self._validate_bin_names(bin_names)
                tag = split + tag
                name = split + name
                # interlace number of bins and bin arrays for each axis
                binnings = [val for name in bin_names for val in (self._nbins[name], self._bins[name])]
                root_hist_cls = HISTOGRAM_REGISTRY[htype]
                self.hists[tag] = root_hist_cls(name, title, *binnings)
                self.names.append(name)
        self._check_for_duplicates()

    def _validate_bin_names(self, bin_names):
        """Check bin names from Histogram configs are defined."""
        invalid_names = [name for name in bin_names if name not in self._bins]
        if invalid_names:
            raise KeyError(f"The following binnings are undefined: {', '.join(invalid_names)}")

    def _init_branches(self, branches):
        # TODO: implement TTrees as output objects
        pass

    def _init_trees(self, trees):
        # TODO: implement TTrees as output objects
        self.trees = {}

    def _init_responses(self, responses_cfg):
        self.responses = Collection(nsplit = self.nsplit)
        for split in self.split_ids:
            for tag, cfg in responses_cfg.items():
                rtype, name, title, *bin_names = cfg
                self._validate_bin_names(bin_names)
                tag = split + tag
                name = split + name
                title = split + title
                # interlace number of bins and bin arrays for each axis
                if len(bin_names) % 2 != 0:
                    raise ValueError("Number of axes defined in this response matrix is not even!")
                ndims = len(bin_names) // 2
                root_hist_cls = HISTOGRAM_REGISTRY[rtype]

                helper_det_name = f"{name}_helper_det"
                binnings_det = [val for name in bin_names[:ndims] for val in (self._nbins[name], self._bins[name])]
                helper_hist_det = root_hist_cls(helper_det_name, helper_det_name, *binnings_det)
                helper_gen_name = f"{name}_helper_gen"
                binnings_gen = [val for name in bin_names[ndims:] for val in (self._nbins[name], self._bins[name])]
                helper_hist_gen = root_hist_cls(helper_gen_name, helper_gen_name, *binnings_gen)

                if ndims == 1:
                    # workaround for 3.0.5 where not specifying the 2D matrix will force bin edges to 0 -> 1
                    # is fixed in 3.1.0 but then 3.1.0 has its own (worse) issues...
                    root_hist_cls_double = mapping_double[root_hist_cls]
                    helper_resp_name = f"{name}_helper_resp"
                    helper_hist_resp = root_hist_cls_double(helper_resp_name, helper_resp_name, *binnings_det, *binnings_gen)
                    self.responses[tag] = RooUnfoldResponse(helper_hist_det, helper_hist_gen, helper_hist_resp, name, title)
                else:
                    self.responses[tag] = RooUnfoldResponse(helper_hist_det, helper_hist_gen, name, title)
                self.names.append(name)
        self._check_for_duplicates()

    def _check_for_duplicates(self):
        duplicates = {name: self.names.count(name) for name in self.names if self.names.count(name) > 1}
        if duplicates:
            raise ValueError(f"Multiple ROOT objects defined with the same name, they will overwrite each other: {list(duplicates.keys())}")

    def save(self, tfile):
        """Write output to an *already-opened* TFile."""
        self.hists.save(tfile)
        self.responses.save(tfile)
        for tree in self.trees.values():
            tfile.WriteTObject(tree)

    def set_split(self, isplit):
        isplit %= self.nsplit
        split_id = self.split[isplit]
        self.hists.set_split(split_id)
        self.responses.set_split(split_id)