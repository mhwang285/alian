#!/usr/bin/env python3

import argparse
import itertools

from alian.analysis.base import AnalysisMCBase, add_default_args, delta_R

import heppyy

fj = heppyy.load_cppyy('fastjet')
alian = heppyy.load_cppyy("alian")


class TrackingPerformance(AnalysisMCBase):
    _defaults = {
        'pt_min_eec': 1.0,
    }
    def init_analysis(self, analysis_cfg: dict):
        config = self._defaults | analysis_cfg
        for setting, value in config.items():
            setattr(self, setting, value)
        self.eec_trk_selector = fj.SelectorPtMin(self.pt_min_eec)

    def analyze_event(self):
        for p in self.particles:
            self.hists['track_pT_gen'].Fill(p.pt(), self.weight)
            if p.user_info[alian.TrackInfo]().has_match():
                self.hists['eff_track_pT_gen'].Fill(p.pt(), self.weight)

                match = p.user_info[alian.TrackInfo]().match()
                match_pt = match.pt()
                res = (match_pt - p.pt()) / p.pt() * 100
                self.hists['pT_res'].Fill(p.pt(), res, self.weight)

                self.hists['ch_track_pT_gen'].Fill(p.pt(), self.weight)
                if match.user_info[alian.TrackInfo]().q() == p.user_info[alian.TrackInfo]().q():
                    self.hists['ch_eff_track_pT_gen'].Fill(p.pt(), self.weight)
        for t in self.tracks:
            self.hists['track_pT_det'].Fill(t.pt(), self.weight)
            if t.user_info[alian.TrackInfo]().has_match():
                self.hists['pur_track_pT_det'].Fill(t.pt(), self.weight)


    def finalize(self):
        self.hists['track_pT_det'].Scale(1, "width")
        self.hists['track_pT_gen'].Scale(1, "width")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Run analysis on ROOT file using YAML configuration.")
    parser = add_default_args(parser)

    args = parser.parse_args()

    ana = TrackingPerformance(args.input_file, args.output_file, args.config_file, args.tree_struct, args.nev, args.lhc_run)
    ana.run()
