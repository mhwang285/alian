from collections import defaultdict
from time import perf_counter

from alian.io import data_io
from ROOT import TFile

import heppyy

from .event import (
    Event,
    EventMC,
    get_particles,
    get_selected_clusters,
    get_selected_tracks,
    get_selected_tracks_mc,
)
from .jet_finder import JetFinder
from .logs import set_up_logger
from .output import Output
from .selector import AnalysisSelector
from .utils import is_slurm, read_yaml, delta_R

alian = heppyy.load_cppyy("alian")

class BaseAnalysis(heppyy.GenericObject):
    _defaults = {}

    def __init__(self, **kwargs):
        super(BaseAnalysis, self).__init__(**kwargs)
        for k, val in self.__class__._defaults.items():
            if not hasattr(self, k) or getattr(self, k) is None:
                setattr(self, k, val)

class AnalysisBase:
    """A base class for analysis tasks.

    This base class describes the basic event analysis structure.
    At minimum, the method `analyze_event()` must be overridden. Depending
    on the selectors defined in the configuration, clusters may or may not be
    loaded into the event.
    """
    _defaults = {}

    def __init__(self, input_file, output_file, cfg, tree_struct, nev = -1, lhc_run = 3):
        self.input_file = input_file
        self.output_file = output_file
        self.cfg_file = cfg
        self.tree_struct = tree_struct
        self.nev = nev
        self.lhc_run = lhc_run

    def run(self):
        self.start_time = perf_counter()
        # set up logger
        self.logger = set_up_logger(__name__)
        self.logger.info("Starting analysis!")

        self.logger.info("Configuring analysis pipeline...")
        self.init()
        self.logger.info("Analysis pipeline configured.")
        self.analyze_events()
        self.finalize()
        self.save()

        self.note_time("Analysis complete")

    def init(self):
        self.cfg = read_yaml(self.cfg_file)
        # validate YAML config
        self._validate_cfg(self.cfg)
        self.logger.info("Config schema validated.")
        # initialize data input
        self.data_source = data_io.DataInput(self.input_file,
                                             lhc_run = self.lhc_run,
                                             yaml_file = self.tree_struct,
                                             n_events = self.nev)
        # initialize selectors for events, tracks, clusters
        self.logger.info("Configuring selectors...")
        self.selector = AnalysisSelector.load(self.cfg)
        self.selector.dump()
        self.logger.info("Selectors configured.")
        if self.selector.is_active("cluster"):
            self.logger.info("Cluster selector found, clusters will be loaded into events.")
            self.load_clusters = True
        else:
            self.logger.info("Cluster selector not found, clusters will NOT be loaded into events.")
            self.load_clusters = False
        # initialize jet finder
        if 'jet_finder' in self.cfg:
            self.init_jet_finder()
        else:
            self.logger.info("No configuration for jet finder, no jets will be loaded.")
            self.load_jets = False
        self.logger.info("Configuring output...")
        self.init_output()
        self.logger.info("Output configured.")
        # initialize analysis parameters
        self.logger.info("Configuring analysis parameters...")
        self.init_analysis(self.cfg.get('analysis', {}))
        self.dump()
        self.logger.info("Analysis parameters configured.")

    def _validate_cfg(self, cfg: dict):
        req_headers = ['selections', 'output']
        if headers_not_in_cfg := [h for h in req_headers if h not in cfg]:
            raise KeyError(f"The following blocks must be present in the YAML configuration: {headers_not_in_cfg}")

    def init_analysis(self, analysis_cfg: dict):
        """Initialize analysis configuration parameters. Override in analyses if necessary."""
        self.logger.info("No analysis configuration to parse.")

    def init_jet_finder(self):
        self.logger.info("Jet finder config found, jets will be loaded into events.")
        self.load_jets = True
        self.logger.info("Configuring jet finder...")
        self.jet_finder = JetFinder.load(self.cfg)
        self.jet_finder.reject_100 = True
        self.jet_finder.dump()
        self.logger.info("Jet finder configured.")

    def init_output(self):
        self.output = Output.load(self.cfg)
        self.hists = self.output.hists
        self.trees = self.output.trees
        self.responses = self.output.responses

    def analyze_events(self):
        self.logger.info("Analyzing events...")
        # slurm_check = is_slurm()
        slurm_check = True
        self.iev = 0
        for e in self.data_source.next_event(disable_bar = slurm_check):
            # build the event only
            self.build_event(e)
            self.iev += 1
            # skip if the event doesn't pass the selection
            if not self.selector.event.selects(self.event):
                continue
            # build the rest of the event and analyze
            self.build_event_objs(e)
            self.analyze_event()
        self.note_time("Events analyzed")

    def analyze_event(self):
        raise NotImplementedError("Must implement an analyze_event method!")

    def build_event(self, event_struct):
        """Build event as Event and store in self.event."""
        self.event = Event(event_struct)
    def build_event_objs(self, event_struct):
        """Build tracks and optionally clusters and jets from associated event."""
        self.tracks = get_selected_tracks(event_struct, self.selector.track)
        if self.load_clusters:
            self.clusters = get_selected_clusters(event_struct, self.selector.cluster)
        if self.load_jets:
            self.jets = self.jet_finder.find_jets(self.tracks)

    def finalize(self):
        """Finalize analysis before saving output. Should be overriden in analyses if necessary."""
        self.logger.info("Nothing to finalize.")
    def save(self):
        self.logger.info(f"Saving output to: {self.output_file}")
        with TFile(self.output_file, "RECREATE") as f:
            self.output.save(f)
        self.logger.info("Output saved.")
    def __getattr__(self, attr):
        if attr == "jets":
            raise AttributeError("Jets not defined; include a `jet_finder` block in your config!")
        if attr == "clusters":
            raise AttributeError("Clusters not defined; include a `clusters` block in the `selections` section of your config! The block can be empty if you want no selections applied.")
        raise AttributeError(f"'{type(self).__name__}' object has no attribute '{attr}'")

    def dump(self):
        cfg = "\n".join([f"\t{param}: {repr(getattr(self, param))}" for param in self._defaults])
        self.logger.info(f"{type(self).__name__} configuration:\n{cfg}", stacklevel = 2)
    def note_time(self, msg):
        self.logger.info(f"{msg}: ------- {self.fmt_time(perf_counter() - self.start_time)} -------", stacklevel = 2)
    def fmt_time(self, seconds):
        m, s = divmod(round(seconds), 60)
        h, m = divmod(m, 60)
        return f"{h:d}h {m:02d}m {s:02d}s"

class AnalysisMCBase(AnalysisBase):
    """A base class for MC analysis tasks.

    Overrides only a few of the AnalysisBase class for MC.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def init_jet_finder(self):
        self.logger.info("Jet finder config found, jets will be loaded into events.")
        self.load_jets = True
        self.logger.info("Configuring jet finders...")
        self.logger.info("Configuring detector-level jet finder...")
        self.jet_finder_det = JetFinder.load(self.cfg)
        self.jet_finder_det.reject_100 = True
        self.jet_finder_det.dump()
        self.logger.info("Configuring generator-level jet finder...")
        self.jet_finder_gen = JetFinder.load(self.cfg)
        self.jet_finder_gen.reject_100 = False
        self.jet_finder_gen.dump()
        self.logger.info("Jet finders configured.")
    def build_event(self, event_struct):
        """Build event as EventMC and store in self.event."""
        self.event = EventMC(event_struct)
        self.weight = self.event.weight
        self.pThat = self.event.pThat
    def build_event_objs(self, event_struct):
        """Build tracks and optionally clusters and jets from associated event."""
        self.tracks = get_selected_tracks_mc(event_struct, self.selector.track)
        self.particles = get_particles(event_struct)
        mapping = self.build_mapping(self.tracks, self.particles)
        self.mapping = self.set_track_part_matching(mapping, self.tracks, self.particles)
        if self.load_jets:
            self.jets_det = self.jet_finder_det.find_jets(self.tracks)
            self.jets_gen = self.jet_finder_gen.find_jets(self.particles)
    def build_mapping(self, dets, gens):
        gen_idx_map = {}
        # build map of gen mcid -> index in gens list
        for gen_idx, gen in enumerate(gens):
            gen_mcid = gen.user_info[alian.TrackInfo]().mcid()
            if gen_mcid in gen_idx_map:
                raise ValueError(f"duplicate MC ID {gen_mcid} in generated particles")
            gen_idx_map[gen_mcid] = gen_idx

        get = gen_idx_map.get # skip the attribute lookup per iteration with local bind
        det_gen_mapping = {}
        for det_idx, det in enumerate(dets):
            det_mcid = det.user_info[alian.TrackInfo]().mcid()
            if det_mcid != -1 and (gen_idx := get(det_mcid)) is not None:
                det_gen_mapping[det_idx] = gen_idx
        return det_gen_mapping
    def set_track_part_matching(self, mapping, tracks, particles):
        mapping_inv = defaultdict(list)
        for det_idx, gen_idx in mapping.items():
            mapping_inv[gen_idx].append(det_idx)

        shared = {gen_idx: det_idxs for gen_idx, det_idxs in mapping_inv.items() if len(det_idxs) > 1}
        if shared:
            self.logger.warning(f"Found {len(shared)} particle matching to multiple tracks, resolving via NN")
            idxs_to_delete = []
            for gen_idx, det_idxs in shared.items():
                self.logger.warning(f"  Particle matches to {len(det_idxs)} tracks")
                min_det_idx = self._get_nearest_track_match(gen_idx, det_idxs)
                idxs_to_delete += [idx for idx in det_idxs if idx != min_det_idx]
            for idx in idxs_to_delete:
                del mapping[idx]

        for match_idx, (det_idx, gen_idx) in enumerate(mapping.items()):
            tracks[det_idx].user_info[alian.TrackInfo]().set_match_idx(match_idx)
            particles[gen_idx].user_info[alian.TrackInfo]().set_match_idx(match_idx)
            tracks[det_idx].user_info[alian.TrackInfo]().set_match(particles[gen_idx])
            particles[gen_idx].user_info[alian.TrackInfo]().set_match(tracks[det_idx])
        return mapping

    def _get_nearest_track_match(self, gen_idx, det_idxs):
        return min(det_idxs, key=lambda det_idx: delta_R(self.particles[gen_idx], self.tracks[det_idx]))

    def _get_jet_matches(self, jets_det, jets_gen, distance):
        pairs = []
        for i, j_det in enumerate(jets_det):
            match = None
            for j, j_gen in enumerate(jets_gen):
                if delta_R(j_det, j_gen) <= distance:
                    if match is not None:
                        match = None
                        break
                    match = j
            if match is None:
                continue
            count = 0
            for j_det_other in jets_det:
                if delta_R(j_det_other, jets_gen[match]) <= distance:
                    count += 1
                    if count > 1:
                        break
            if count == 1:
                pairs.append((i, match))
        return pairs
    def _get_jet_matches_cpp(self, jets_det, jets_gen, distance):
        vec = alian.get_jet_matches(jets_det, jets_gen, distance)
        return [(p.first, p.second) for p in vec]

def add_default_args(parser):
    parser.add_argument('-i', '--input-file',  type = str, required = True,           help = "Input file or file containing a list of input files.")
    parser.add_argument("-o", "--output-file", type = str, default = "analysis.root", help = "File to write the analysis.")
    parser.add_argument('-c', '--config-file', type = str, required = True,           help = "Path to a YAML configuration file describing the cuts to be applied to the data.")
    parser.add_argument("-n", "--nev",         type = int, default = -1,              help = "Number of entries to process, set to -1 for all events.")
    parser.add_argument('-t', '--tree-struct', type = str, default = None,            help = "Path to a YAML file describing the tree structure.")
    parser.add_argument('--lhc-run',           type = int, default = 3,               help = 'LHC Run')

    return parser
