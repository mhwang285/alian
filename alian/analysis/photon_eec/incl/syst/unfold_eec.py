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
from itertools import combinations
from typing import NamedTuple

from alian.analysis.base import AnalysisMCBase, add_default_args
from alian.analysis.base.utils import delta_R

import heppyy

fj = heppyy.load_cppyy('fastjet')
alian = heppyy.load_cppyy("alian")

class Pair(NamedTuple):
    ctype: str
    RL: float
    ew: float
    jet_pT: float
    hasmatch1: bool
    hasmatch2: bool
    midx1: int
    midx2: int

    @classmethod
    def from_pseudojets(cls, p1, p2, jet_pT):
        q1 = p1.user_info[alian.TrackInfo]().q()
        q2 = p2.user_info[alian.TrackInfo]().q()
        if q1 > 0 and q2 > 0:
            ctype = "P"
        elif q1 < 0 and q2 < 0:
            ctype = "M"
        else:
            ctype = "PM"
        ew = p1.pt() * p2.pt() / jet_pT / jet_pT
        RL = delta_R(p1, p2)
        return cls(ctype, RL, ew, jet_pT,
                   p1.user_info[alian.TrackInfo]().has_match(),
                   p2.user_info[alian.TrackInfo]().has_match(),
                   p1.user_info[alian.TrackInfo]().match_idx(),
                   p2.user_info[alian.TrackInfo]().match_idx(),
                   )
    def is_idx_match_to(self, other) -> bool:
        return self.hasmatch1 and self.hasmatch2 and other.hasmatch1 and other.hasmatch2 \
                and ((self.midx1 == other.midx1 and self.midx2 == other.midx2)  \
                or   (self.midx1 == other.midx2 and self.midx2 == other.midx1))
    def is_charge_match_to(self, other) -> bool:
        return self.ctype == other.ctype
    def is_full_match_to(self, other) -> bool:
        return self.is_charge_match_to(other) and self.is_idx_match_to(other)
    def is_in_range(self, RL_min, RL_max, ew_min, ew_max, jet_pT_min, jet_pT_max):
        return self.RL > RL_min and self.RL < RL_max \
            and self.ew > ew_min and self.ew < ew_max \
            and self.jet_pT > jet_pT_min and self.jet_pT < jet_pT_max

def pair_full_match(pair1, pair2):
    return pair1.is_full_match_to(pair2)
def pair_idx_match(pair1, pair2):
    return pair1.is_idx_match_to(pair2)

class UnfoldEECPairs(AnalysisMCBase):
    _defaults = {
        'pt_min_eec': 1.0,
        'eec_RL_min': 1.0e-2,
        'eec_RL_max': 0.4,
        'eec_ew_min': 1.0e-5,
        'eec_ew_max': 1.0,
        'eec_jet_pT_min': 20,
        'eec_jet_pT_max': 160,
        "pThat_scale_det": 1000,
        "pThat_scale_gen": 1000,
    }
    def init_analysis(self, analysis_cfg: dict):
        config = self._defaults | analysis_cfg
        for setting, value in config.items():
            setattr(self, setting, value)
        self.eec_trk_selector = fj.SelectorPtMin(self.pt_min_eec)
        self.ctypes = ["T", "P", "M", "PM"]
        self.eec_RL_min = self.output._bins["RL"][0]
        self.eec_RL_max = self.output._bins["RL"][-1]
        self.eec_ew_min = self.output._bins["ew"][0]
        self.eec_ew_max = self.output._bins["ew"][-1]
        self.eec_jet_pT_min = self.output._bins["jet_pT"][0]
        self.eec_jet_pT_max = self.output._bins["jet_pT"][-1]
        # self.eec_RL_min = self.hists["eec_T_det"].GetXaxis().GetXmin()
        # self.eec_RL_max = self.hists["eec_T_det"].GetXaxis().GetXmax()
        # self.eec_ew_min = self.hists["eec_T_det"].GetYaxis().GetXmin()
        # self.eec_ew_max = self.hists["eec_T_det"].GetYaxis().GetXmax()
        # self.eec_jet_pT_min = self.hists["eec_T_det"].GetZaxis().GetXmin()
        # self.eec_jet_pT_max = self.hists["eec_T_det"].GetZaxis().GetXmax()

    def is_matched(self, psj1, psj2):
        return psj1.user_info[alian.TrackInfo]().is_matched_to(psj2)

    def analyze_event(self):
        pass_cut = True
        if self.jets_det and (self.jets_det[0].pt() / self.pThat) > self.pThat_scale_det:
            pass_cut = False
        if self.jets_gen and (self.jets_gen[0].pt() / self.pThat) > self.pThat_scale_gen:
            pass_cut = False
        if not pass_cut:
            self.logger.warning(f"Rejecting event: {len(self.jets_det)} det jets, {len(self.jets_gen)} gen jets")
            return

        pairs_det = {ctype: [] for ctype in self.ctypes}
        for jet_det in self.jets_det:
            jet_pT_det = jet_det.pt()
            const_det = self.eec_trk_selector(jet_det.constituents())
            for cdet1, cdet2 in combinations(const_det, 2):
                pair_det = Pair.from_pseudojets(cdet1, cdet2, jet_pT_det)
                ctype_det = pair_det.ctype
                pairs_det["T"].append(pair_det)
                pairs_det[ctype_det].append(pair_det)

        pairs_gen = {ctype: [] for ctype in self.ctypes}
        for jet_gen in self.jets_gen:
            jet_pT_gen = jet_gen.pt()
            const_gen = self.eec_trk_selector(jet_gen.constituents())
            for cgen1, cgen2 in combinations(const_gen, 2):
                pair_gen = Pair.from_pseudojets(cgen1, cgen2, jet_pT_gen)
                ctype_gen = pair_gen.ctype
                pairs_gen["T"].append(pair_gen)
                pairs_gen[ctype_gen].append(pair_gen)

        # for ctype in self.ctypes:
        self.analyze_pairs(pairs_det["T"], pairs_gen["T"], "T", pair_idx_match)
        self.analyze_pairs(pairs_det["P"], pairs_gen["P"], "P", pair_full_match)
        self.analyze_pairs(pairs_det["M"], pairs_gen["M"], "M", pair_full_match)
        self.analyze_pairs(pairs_det["PM"], pairs_gen["PM"], "PM", pair_full_match)

    def analyze_pairs(self, pairs_det, pairs_gen, ctype, match_fcn):
        resp = self.responses[f"eec_{ctype}_unf"]
        matched_idxs_det = []
        matched_idxs_gen = []
        for ipair_det, pair_det in enumerate(pairs_det):
            for ipair_gen, pair_gen in enumerate(pairs_gen):
                if match_fcn(pair_det, pair_gen) \
                    and pair_det.is_in_range(self.eec_RL_min, self.eec_RL_max, self.eec_ew_min, self.eec_ew_max, self.eec_jet_pT_min, self.eec_jet_pT_max) \
                    and pair_gen.is_in_range(self.eec_RL_min, self.eec_RL_max, self.eec_ew_min, self.eec_ew_max, self.eec_jet_pT_min, self.eec_jet_pT_max):
                    matched_idxs_det.append(ipair_det)
                    matched_idxs_gen.append(ipair_gen)

        # self.logger.info("Starting pairs analysis")
        for ipair_det, ipair_gen in zip(matched_idxs_det, matched_idxs_gen):
            # self.logger.info(f"{ipair_det} {ipair_gen}")
            pair_det = pairs_det[ipair_det]
            pair_gen = pairs_gen[ipair_gen]
            resp.Fill(pair_det.RL, pair_det.ew, pair_det.jet_pT, pair_gen.RL, pair_gen.ew, pair_gen.jet_pT, self.weight)
            resp.Fill(pair_det.RL, pair_det.ew, pair_det.jet_pT, pair_gen.RL, pair_gen.ew, pair_gen.jet_pT, self.weight)

        idxs_fake = [idx for idx in range(len(pairs_det)) if idx not in matched_idxs_det]
        idxs_miss = [idx for idx in range(len(pairs_gen)) if idx not in matched_idxs_gen]

        for idx_fake in idxs_fake:
            pair_fake = pairs_det[idx_fake]
            if pair_fake.is_in_range(self.eec_RL_min, self.eec_RL_max, self.eec_ew_min, self.eec_ew_max, self.eec_jet_pT_min, self.eec_jet_pT_max):
                resp.Fake(pair_fake.RL, pair_fake.ew, pair_fake.jet_pT, self.weight)
                resp.Fake(pair_fake.RL, pair_fake.ew, pair_fake.jet_pT, self.weight)

        for idx_miss in idxs_miss:
            pair_miss = pairs_gen[idx_miss]
            if pair_miss.is_in_range(self.eec_RL_min, self.eec_RL_max, self.eec_ew_min, self.eec_ew_max, self.eec_jet_pT_min, self.eec_jet_pT_max):
                resp.Miss(pair_miss.RL, pair_miss.ew, pair_miss.jet_pT, self.weight)
                resp.Miss(pair_miss.RL, pair_miss.ew, pair_miss.jet_pT, self.weight)

        # for pair_det in pairs_det:
        #     if pair_det.is_in_range(self.eec_RL_min, self.eec_RL_max, self.eec_ew_min, self.eec_ew_max, self.eec_jet_pT_min, self.eec_jet_pT_max):
        #         self.hists[f"eec_{ctype}_det"].Fill(pair_det.RL, pair_det.ew, pair_det.jet_pT, self.weight)
        #         self.hists[f"eec_{ctype}_det"].Fill(pair_det.RL, pair_det.ew, pair_det.jet_pT, self.weight)
        # for pair_gen in pairs_gen:
        #     if pair_gen.is_in_range(self.eec_RL_min, self.eec_RL_max, self.eec_ew_min, self.eec_ew_max, self.eec_jet_pT_min, self.eec_jet_pT_max):
        #         self.hists[f"eec_{ctype}_gen"].Fill(pair_gen.RL, pair_gen.ew, pair_gen.jet_pT, self.weight)
        #         self.hists[f"eec_{ctype}_gen"].Fill(pair_gen.RL, pair_gen.ew, pair_gen.jet_pT, self.weight)

    def finalize(self):
        pass
        # hmeasured = self.responses["eec_P_unf"].Hmeasured()
        # hmeasured.SetName("hmeasured_P")
        # self.hists["hmeasured_P"] = hmeasured
        # htruth = self.responses["eec_P_unf"].Htruth()
        # htruth.SetName("htruth_P")
        # self.hists["htruth_P"] = htruth
        # hresponse = self.responses["eec_P_unf"].Hresponse()
        # hresponse.SetName("hresponse_P")
        # self.hists["hresponse_P"] = hresponse

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

    ana = UnfoldEECPairs(args.input_file, args.output_file, args.config_file, args.tree_struct, args.nev, args.lhc_run)
    ana.run()
