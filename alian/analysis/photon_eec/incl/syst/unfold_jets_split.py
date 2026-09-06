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

from alian.analysis.base import AnalysisMCBase, add_default_args

import heppyy

fj = heppyy.load_cppyy('fastjet')
alian = heppyy.load_cppyy("alian")


class UnfoldJetsSplit(AnalysisMCBase):
    _defaults = {
        # 'pt_min_eec': 1.0,
        "jet_matching": 0.6,
        "jet_pT_max_gen": 200,
        "jet_pT_max_det": 200,
        "split": 0.5,
        "rng_seed": 5,
    }
    def init_analysis(self, analysis_cfg: dict):
        config = self._defaults | analysis_cfg
        for setting, value in config.items():
            setattr(self, setting, value)
        self.eec_trk_selector = fj.SelectorPtMin(self.pt_min_eec)
        self.jet_pT_max_det = self.hists["1_jet_pT_det"].GetXaxis().GetXmax()
        self.jet_pT_max_gen = self.hists["1_jet_pT_gen"].GetXaxis().GetXmax()
        self.analyze_event = self.analyze_event_split
        if self.split == 0.5:
            self.is_1 = True
        else:
            self.analyze_event = self.analyze_event_split
            self.rng = np.random.default_rng(seed = self.rng_seed)

    def analyze_event_split(self):
        if self.is_1:
            split = 1
        else:
            split = 2
        self.is_1 = not self.is_1

        for j in self.jets_det:
            self.hists[f'{split}_jet_pT_det'].Fill(j.pt(), self.weight)
            # self.do_eec_det(j)
        for j in self.jets_gen:
            self.hists[f'{split}_jet_pT_gen'].Fill(j.pt(), self.weight)
            # self.do_eec_gen(j)

        # pairs_py = self._get_jet_matches(self.jets_det, self.jets_gen, 0.4*0.6)
        # pairs_cpp = self._get_jet_matches_cpp(self.jets_det, self.jets_gen, 0.4*self.jet_matching)
        # if pairs_py != pairs_cpp:
        #     raise ValueError(f"Pairs don't match:\n\tPython: {pairs_py}\n\tC++: {pairs_cpp}")
        pairs_cpp = alian.get_jet_matches_filtered(self.jets_det, self.jets_gen, 0.4*self.jet_matching,
                                                self.jet_pT_max_det,
                                                self.jet_pT_max_gen)
        pairs_cpp = [(p.first, p.second) for p in pairs_cpp]

        resp = self.responses[f"{split}_jet_pT_unf"]
        for idx_det, idx_gen in pairs_cpp:
            jet_det_matched = self.jets_det[idx_det]
            jet_gen_matched = self.jets_gen[idx_gen]

            self.hists[f'{split}_jet_pT_gen_matched'].Fill(jet_gen_matched.pt(), self.weight)
            self.hists[f'{split}_jet_pT_det_matched'].Fill(jet_det_matched.pt(), self.weight)
            self.hists[f'{split}_jet_pT_resp'].Fill(jet_det_matched.pt(), jet_gen_matched.pt(), self.weight)
            resp.Fill(jet_det_matched.pt(), jet_gen_matched.pt(), self.weight)

        idxs_det_matched = [pair[0] for pair in pairs_cpp]
        idxs_fake = [idx for idx in range(len(self.jets_det)) if idx not in idxs_det_matched]
        idxs_gen_matched = [pair[1] for pair in pairs_cpp]
        idxs_miss = [idx for idx in range(len(self.jets_gen)) if idx not in idxs_gen_matched]

        for idx_fake in idxs_fake:
            jet_fake = self.jets_det[idx_fake]
            resp.Fake(jet_fake.pt(), self.weight)

        for idx_miss in idxs_miss:
            jet_miss = self.jets_gen[idx_miss]
            resp.Miss(jet_miss.pt(), self.weight)

    def finalize(self):
        pass
        # hresponse = self.responses["jet_pT_unf"].Hresponse()
        # hresponse.SetName("hresponse")
        # self.responses["hresponse"] = hresponse

    # def do_eec_det(self, jet):
    #     tracks = self.eec_trk_selector(jet.constituents())

    #     for p1, p2 in itertools.permutations(tracks, 2):
    #         ew = p1.pt() * p2.pt() / jet.pt() / jet.pt()
    #         angle = delta_R(p1, p2)
    #         q1 = p1.user_info[alian.TrackInfo]().q()
    #         q2 = p2.user_info[alian.TrackInfo]().q()
    #         phistar = self.calc_phistar(p1, p2, q1, q2)
    #         delta_eta = p2.eta() - p1.eta()
    #         if np.abs(phistar) < self.phistar_cut and np.abs(delta_eta) < self.eta_cut:
    #             continue

    #         self.hists["eec_T_det"].Fill(jet.pt(), angle, ew * self.weight)
    #         self.hists["eec_Q_det"].Fill(jet.pt(), angle, ew * self.weight * q1 * q2)
    #         if q1 > 0 and q2 > 0:
    #             self.hists["eec_P_det"].Fill(jet.pt(), angle, ew * self.weight)
    #         elif q1 < 0 and q2 < 0:
    #             self.hists["eec_M_det"].Fill(jet.pt(), angle, ew * self.weight)
    #         else:
    #             self.hists["eec_PM_det"].Fill(jet.pt(), angle, ew * self.weight)

    # def do_eec_gen(self, jet):
    #     particles = self.eec_trk_selector(jet.constituents())

    #     for p1, p2 in itertools.permutations(particles, 2):
    #         ew = p1.pt() * p2.pt() / jet.pt() / jet.pt()
    #         angle = delta_R(p1, p2)
    #         q1 = p1.user_info[alian.ParticleInfo]().q()
    #         q2 = p2.user_info[alian.ParticleInfo]().q()

    #         self.hists["eec_T_gen"].Fill(jet.pt(), angle, ew * self.weight)
    #         self.hists["eec_Q_gen"].Fill(jet.pt(), angle, ew * self.weight * q1 * q2)
    #         if q1 > 0 and q2 > 0:
    #             self.hists["eec_P_gen"].Fill(jet.pt(), angle, ew * self.weight)
    #         elif q1 < 0 and q2 < 0:
    #             self.hists["eec_M_gen"].Fill(jet.pt(), angle, ew * self.weight)
    #         else:
    #             self.hists["eec_PM_gen"].Fill(jet.pt(), angle, ew * self.weight)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Run analysis on ROOT file using YAML configuration.")
    parser = add_default_args(parser)

    args = parser.parse_args()

    ana = UnfoldJetsSplit(args.input_file, args.output_file, args.config_file, args.tree_struct, args.nev, args.lhc_run)
    ana.run()
