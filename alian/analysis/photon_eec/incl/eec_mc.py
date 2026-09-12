#!/usr/bin/env python3

import argparse
import itertools

from alian.analysis.base import AnalysisMCBase, add_default_args, delta_R

import heppyy

fj = heppyy.load_cppyy('fastjet')
alian = heppyy.load_cppyy("alian")


class EECMC(AnalysisMCBase):
    _defaults = {
        'pt_min_eec': 1.0,
        "pThat_scale_det": 1000,
        "pThat_scale_gen": 1000,
    }
    def init_analysis(self, analysis_cfg: dict):
        config = self._defaults | analysis_cfg
        for setting, value in config.items():
            setattr(self, setting, value)
        self.eec_trk_selector = fj.SelectorPtMin(self.pt_min_eec)

    def analyze_event(self):
        pass_cut = True
        if self.jets_det and (self.jets_det[0].pt() / self.pThat) > self.pThat_scale_det:
            pass_cut = False
        if self.jets_gen and (self.jets_gen[0].pt() / self.pThat) > self.pThat_scale_gen:
            pass_cut = False
        if not pass_cut:
            self.logger.warning(f"Rejecting event: {len(self.jets_det)} det jets, {len(self.jets_gen)} gen jets")
            return

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
            if q1 < 0 and q2 < 0:
                self.hists[f"eec_M_{suffix}"].Fill(jet.pt(), angle, ew * self.weight)
            else:
                self.hists[f"eec_PM_{suffix}"].Fill(jet.pt(), angle, ew * self.weight)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser = add_default_args(parser)

    args = parser.parse_args()

    ana = EECMC(args.input_file, args.output_file, args.config_file, args.tree_struct, args.nev, args.lhc_run)
    ana.run()
