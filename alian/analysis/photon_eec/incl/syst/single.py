#!/usr/bin/env python3

"""Example script on using the alian framework for Run 3 analysis.

This script is used on the ROOT file (provided as an argument) to
perform the analysis of the events containing isolated photons and
jets. The analysis is configured with a config YAML file, creates
and fills histograms for the cleaned data, and stores them in the
output file.

To be used with alian/config/example/example_analysis.yaml
"""

import argparse
import itertools

import numpy as np
from alian.analysis.base import AnalysisMCBase, JetFinder, add_default_args, delta_R
from cppyy.gbl import std

import heppyy

fj = heppyy.load_cppyy('fastjet')
alian = heppyy.load_cppyy("alian")


class SystSingle(AnalysisMCBase):
    _defaults = {
        'pt_min_eec': 1.0,
        'rng_seed': 5,
        "pThat_scale_det": 1000,
        "pThat_scale_gen": 1000,
    }
    def init_analysis(self, analysis_cfg: dict):
        config = self._defaults | analysis_cfg
        for setting, value in config.items():
            setattr(self, setting, value)
        self.eec_trk_selector = fj.SelectorPtMin(self.pt_min_eec)
        self.rng = np.random.default_rng(seed = self.rng_seed)
        self.jet_finder_rej = JetFinder.load(self.cfg)
        self.jet_finder_rej.reject_100 = True
        self.jet_finder_rej.dump()

    def reject_tracks_pt(self, tracks):
        # 2022o
        edges = [0.15, 0.3, 0.5, 1, 1.3, 2.5, 25, 40,  50, 70, 100, 200]
        unc   = [   2.5, 1.5,  2, 2.5,  3,  2,  2.5, 3.5, 5, 8,   18]
        unc = [1 - (unc_val / 100) for unc_val in unc]
        # 2023 thinned
        # edges = [0.15, 0.3, 1.3,  3,  4,  30,  40,  50,  70, 100, 200]
        # unc =      [1.5, 1.0,  1.5, 1.0, 1.5, 2.0, 2.5, 3.5, 6.0, 19] 
        # 2024 ppref
        # edges = [0.15, 1.6, 2.5, 12,  20,  25,  30,  40,  50,  70, 100, 200]
        # unc =      [1.5, 2.0,  1.5, 2.0, 2.5, 3.0, 4.0, 5.5, 12, 14.5, 24] 

        # ntracks = len(tracks)

        tracks_new = std.vector[fj.PseudoJet]([])
        pts = [part.pt() for part in tracks]
        pt_bin_idxs = np.searchsorted(edges, pts, side='left') - 1
        for i in range(len(pt_bin_idxs)):
            if pt_bin_idxs[i] == len(unc):
                pt_bin_idxs[i] = pt_bin_idxs[i] - 1
        probs = [unc[idx] for idx in pt_bin_idxs]
        results = self.rng.binomial(n = 1, p = probs) == 1
        [tracks_new.push_back(part) for part, res in zip(tracks, results, strict = True) if res]

        # ntracks_new = len(tracks_new)
        # nrej = ntracks - ntracks_new
        # rate = (ntracks - ntracks_new) / ntracks * 100 if ntracks != 0 else 0
        # self.logger.info(f"Out of {ntracks} tracks, {nrej} tracks were rejected, overall rejection rate {rate:.4f}%.")
        # self.note_time("pT-dependent rejection complete")
        return tracks_new


    def analyze_event(self):
        pass_cut = True
        if self.jets_det and (self.jets_det[0].pt() / self.pThat) > self.pThat_scale_det:
            pass_cut = False
        if self.jets_gen and (self.jets_gen[0].pt() / self.pThat) > self.pThat_scale_gen:
            pass_cut = False
        if not pass_cut:
            self.logger.warning(f"Rejecting event: {len(self.jets_det)} det jets, {len(self.jets_gen)} gen jets")
            return

        rej_tracks = self.reject_tracks_pt(self.tracks)
        jets_rej = self.jet_finder_rej.find_jets(rej_tracks)
        [self.hists['track_pT_rej'].Fill(t.pt(), self.weight) for t in rej_tracks]
        [self.hists['track_pT_det'].Fill(t.pt(), self.weight) for t in self.tracks]
        [self.hists['track_pT_gen'].Fill(t.pt(), self.weight) for t in self.particles]

        for j in jets_rej:
            self.do_eec(j, "rej")
        for j in self.jets_det:
            self.do_eec(j, "det")
        for j in self.jets_gen:
            self.do_eec(j, "gen")

    def do_eec(self, jet, suffix):
        self.hists[f'jet_pT_{suffix}'].Fill(jet.pt(), self.weight)
        tracks = self.eec_trk_selector(jet.constituents())

        for p1, p2 in itertools.permutations(tracks, 2):
            ew = p1.pt() * p2.pt() / jet.pt() / jet.pt()
            angle = delta_R(p1, p2)
            q1 = p1.user_info[alian.TrackInfo]().q()
            q2 = p2.user_info[alian.TrackInfo]().q()

            self.hists[f"eec_T_{suffix}"].Fill(jet.pt(), angle, ew * self.weight)
            self.hists[f"eec_Q_{suffix}"].Fill(jet.pt(), angle, ew * self.weight * q1 * q2)
            if q1 > 0 and q2 > 0:
                self.hists[f"eec_P_{suffix}"].Fill(jet.pt(), angle, ew * self.weight)
            elif q1 < 0 and q2 < 0:
                self.hists[f"eec_M_{suffix}"].Fill(jet.pt(), angle, ew * self.weight)
            else:
                self.hists[f"eec_PM_{suffix}"].Fill(jet.pt(), angle, ew * self.weight)

    def finalize(self):
        self.hists['track_pT_rej'].Scale(1, "width")
        self.hists['track_pT_det'].Scale(1, "width")
        self.hists['track_pT_gen'].Scale(1, "width")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Run analysis on ROOT file using YAML configuration.")
    parser = add_default_args(parser)

    args = parser.parse_args()

    ana = SystSingle(args.input_file, args.output_file, args.config_file, args.tree_struct, args.nev, args.lhc_run)
    ana.run()
