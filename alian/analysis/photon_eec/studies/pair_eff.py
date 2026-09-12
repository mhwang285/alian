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
from alian.analysis.base import AnalysisMCBase, add_default_args

import heppyy

fj = heppyy.load_cppyy('fastjet')
alian = heppyy.load_cppyy("alian")


class AnalysisExample(AnalysisMCBase):
    _defaults = {
        "pThat_scale_det": 1000,
        "pThat_scale_gen": 1000,
    }
    def init_analysis(self, analysis_cfg: dict):
        config = self._defaults | analysis_cfg
        for setting, value in config.items():
            setattr(self, setting, value)

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

    def calc_kt(self, p1, p2):
        return 0.5*np.sqrt( pow(p1.pt(),2)+pow(p2.pt(),2)+2*p1.pt()*p2.pt()*np.cos(p1.phi()-p2.phi()) )
        # NOTE: this is actually law of cosines but vector causes plus not minus on cos term
        # e.g. consider case when phis are equal, then sum has double magnitude, which is only possible
        # with + cos not - cos (proper angle is 180-theta not theta)

    def analyze_event(self):
        pass_cut = True
        if self.jets_det and (self.jets_det[0].pt() / self.pThat) > self.pThat_scale_det:
            pass_cut = False
        if self.jets_gen and (self.jets_gen[0].pt() / self.pThat) > self.pThat_scale_gen:
            pass_cut = False
        if not pass_cut:
            self.logger.warning(f"Rejecting event: {len(self.jets_det)} det jets, {len(self.jets_gen)} gen jets")
            return

        for p1, p2 in itertools.combinations(self.particles, 2):
            info1 = p1.user_info[alian.TrackInfo]()
            info2 = p2.user_info[alian.TrackInfo]()
            q1 = info1.q()
            q2 = info2.q()
            p1_has_match = info1.has_match()
            p2_has_match = info2.has_match()
            delta_phistar = self.calc_phistar(p1, p2, q1, q2)
            delta_eta = p2.eta() - p1.eta()
            pair_kt = self.calc_kt(p1, p2)
            # pair_kt = (p1.pt() + p2.pt()) / 2

            if q1 > 0 and q2 > 0:
                pair_type = "P"
            elif q1 < 0 and q2 < 0:
                pair_type = "M"
            else:
                pair_type = "PM"

            self.hists["pair_eff_T_gen"].Fill(pair_kt, delta_phistar, delta_eta, self.weight)
            self.hists[f"pair_eff_{pair_type}_gen"].Fill(pair_kt, delta_phistar, delta_eta, self.weight)


            if p1_has_match and p2_has_match:
                self.hists["pair_eff_T_det"].Fill(pair_kt, delta_phistar, delta_eta, self.weight)
                self.hists[f"pair_eff_{pair_type}_det"].Fill(pair_kt, delta_phistar, delta_eta, self.weight)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Run analysis on ROOT file using YAML configuration.")
    parser = add_default_args(parser)

    args = parser.parse_args()

    ana = AnalysisExample(args.input_file, args.output_file, args.config_file, args.tree_struct, args.nev, args.lhc_run)
    ana.run()
