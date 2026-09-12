#!/usr/bin/env python

import argparse
import itertools
import heppyy
from alian.analysis.base import AnalysisBase, add_default_args
from alian.analysis.base.utils import delta_R

fj = heppyy.load_cppyy('fastjet')
alian = heppyy.load_cppyy("alian")

class BasicEEC3D(AnalysisBase):
    _defaults = {
        'pt_min_eec': 1.0,
    }
    def init_analysis(self, analysis_cfg: dict):
        config = self._defaults | analysis_cfg
        for setting, value in config.items():
            setattr(self, setting, value)
        self.eec_trk_selector = fj.SelectorPtMin(self.pt_min_eec)

    def analyze_event(self):
        for jet in self.jets:
            self.do_eec(jet)

    def do_eec(self, jet):
        self.hists["jet_pT"].Fill(jet.pt())
        tracks = self.eec_trk_selector(jet.constituents())
        for p1, p2 in itertools.permutations(tracks, 2):
            ew = p1.pt() * p2.pt() / jet.pt() / jet.pt()
            angle = delta_R(p1, p2)
            q1 = p1.user_info[alian.TrackInfo]().q()
            q2 = p2.user_info[alian.TrackInfo]().q()

            self.hists["eec_T"].Fill(angle, ew, jet.pt())
            if q1 > 0 and q2 > 0:
                self.hists["eec_P"].Fill(angle, ew, jet.pt())
            elif q1 < 0 and q2 < 0:
                self.hists["eec_M"].Fill(angle, ew, jet.pt())
            else:
                self.hists["eec_PM"].Fill(angle, ew, jet.pt())


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Run analysis on ROOT file using YAML configuration.")
    add_default_args(parser)
    args = parser.parse_args()

    ana = BasicEEC3D(args.input_file, args.output_file, args.config_file, args.tree_struct, args.nev, args.lhc_run)
    ana.run()