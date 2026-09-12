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
from alian.analysis.base import AnalysisMCBase, add_default_args, delta_R

import heppyy

fj = heppyy.load_cppyy('fastjet')
alian = heppyy.load_cppyy("alian")


class SystPairMC(AnalysisMCBase):
    _defaults = {
        'pt_min_eec': 1.0,
        'phistar_cut': 0.01,
        'eta_cut': 0.008,
        "pThat_scale_det": 1000,
        "pThat_scale_gen": 1000,
    }
    def init_analysis(self, analysis_cfg: dict):
        config = self._defaults | analysis_cfg
        for setting, value in config.items():
            setattr(self, setting, value)
        self.eec_trk_selector = fj.SelectorPtMin(self.pt_min_eec)


    def calc_phistar(self, p1, p2, q1, q2):
        R = 1.1 # reference radius for TPC
        Bz = 0.5 # LHC22o is ++ field
        delta_phi = p1.delta_phi_to(p2) # returns phi2 - phi1
        dalpha = - q2 * np.arcsin(0.1499 * Bz * R / p2.pt()) + q1 * np.arcsin(0.1499 * Bz * R / p1.pt())
        delta_phistar = delta_phi + dalpha

        return self.Phi_mpi_pi(delta_phistar)

    def Phi_mpi_pi(self, dphi):
        while dphi >= np.pi:
            dphi -= (np.pi * 2)
        while dphi < -np.pi:
            dphi += (np.pi * 2)
        return dphi

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
            self.do_eec(j, "rej")
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
            if suffix == "rej":
                phistar = self.calc_phistar(p1, p2, q1, q2)
                delta_eta = p2.eta() - p1.eta()
                if np.abs(phistar) < self.phistar_cut and np.abs(delta_eta) < self.eta_cut:
                    continue

            self.hists[f"eec_T_{suffix}"].Fill(jet.pt(), angle, ew * self.weight)
            self.hists[f"eec_Q_{suffix}"].Fill(jet.pt(), angle, ew * self.weight * q1 * q2)
            if q1 > 0 and q2 > 0:
                self.hists[f"eec_P_{suffix}"].Fill(jet.pt(), angle, ew * self.weight)
            elif q1 < 0 and q2 < 0:
                self.hists[f"eec_M_{suffix}"].Fill(jet.pt(), angle, ew * self.weight)
            else:
                self.hists[f"eec_PM_{suffix}"].Fill(jet.pt(), angle, ew * self.weight)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Run analysis on ROOT file using YAML configuration.")
    parser = add_default_args(parser)

    args = parser.parse_args()

    ana = SystPairMC(args.input_file, args.output_file, args.config_file, args.tree_struct, args.nev, args.lhc_run)
    ana.run()
