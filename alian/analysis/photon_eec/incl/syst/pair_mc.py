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
        for j in self.jets_det:
            self.hists['jet_pT_det'].Fill(j.pt(), self.weight)
            self.do_eec_det(j)
        for j in self.jets_gen:
            self.hists['jet_pT_gen'].Fill(j.pt(), self.weight)
            self.do_eec_gen(j)

    def do_eec_det(self, jet):
        tracks = self.eec_trk_selector(jet.constituents())

        for p1, p2 in itertools.permutations(tracks, 2):
            ew = p1.pt() * p2.pt() / jet.pt() / jet.pt()
            angle = delta_R(p1, p2)
            q1 = p1.user_info[alian.TrackInfo]().q()
            q2 = p2.user_info[alian.TrackInfo]().q()
            phistar = self.calc_phistar(p1, p2, q1, q2)
            delta_eta = p2.eta() - p1.eta()
            if np.abs(phistar) < self.phistar_cut and np.abs(delta_eta) < self.eta_cut:
                continue

            self.hists["eec_T_det"].Fill(jet.pt(), angle, ew * self.weight)
            self.hists["eec_Q_det"].Fill(jet.pt(), angle, ew * self.weight * q1 * q2)
            if q1 > 0 and q2 > 0:
                self.hists["eec_P_det"].Fill(jet.pt(), angle, ew * self.weight)
            elif q1 < 0 and q2 < 0:
                self.hists["eec_M_det"].Fill(jet.pt(), angle, ew * self.weight)
            else:
                self.hists["eec_PM_det"].Fill(jet.pt(), angle, ew * self.weight)

    def do_eec_gen(self, jet):
        particles = self.eec_trk_selector(jet.constituents())

        for p1, p2 in itertools.permutations(particles, 2):
            ew = p1.pt() * p2.pt() / jet.pt() / jet.pt()
            angle = delta_R(p1, p2)
            q1 = p1.user_info[alian.TrackInfo]().q()
            q2 = p2.user_info[alian.TrackInfo]().q()

            self.hists["eec_T_gen"].Fill(jet.pt(), angle, ew * self.weight)
            self.hists["eec_Q_gen"].Fill(jet.pt(), angle, ew * self.weight * q1 * q2)
            if q1 > 0 and q2 > 0:
                self.hists["eec_P_gen"].Fill(jet.pt(), angle, ew * self.weight)
            elif q1 < 0 and q2 < 0:
                self.hists["eec_M_gen"].Fill(jet.pt(), angle, ew * self.weight)
            else:
                self.hists["eec_PM_gen"].Fill(jet.pt(), angle, ew * self.weight)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Run analysis on ROOT file using YAML configuration.")
    parser = add_default_args(parser)

    args = parser.parse_args()

    ana = SystPairMC(args.input_file, args.output_file, args.config_file, args.tree_struct, args.nev, args.lhc_run)
    ana.run()
