import unittest
import numpy as np
from nipype.interfaces.base import CommandLine
from QUIT.interfaces.core import NewImage, Diff
from QUIT.interfaces.mt import Lorentzian, LorentzianSim, Lineshape, qMT, qMTSim, ZSpec

vb = True
CommandLine.terminal_output = 'allatonce'


class MT(unittest.TestCase):
    def test_lorentzian(self):
        sequence = {'MTSat': {'pulse': {'p1': 0.4,
                                        'p2': 0.3,
                                        'bandwidth': 0.39},
                              'TR': 4,
                              'Trf': 0.02,
                              'FA': 5,
                              'sat_f0': np.linspace(-5, 5, 21).squeeze().tolist(),
                              'sat_angle': np.repeat(180.0, 21).squeeze().tolist()}}
        pools = [{'name': 'DS',
                  'df0': [0, -2.5, 2.5],
                  'fwhm': [1.0, 1.e-6, 3.0],
                  'A': [0.2, 1.e-3, 1.0],
                  'use_bandwidth': True}]
        lorentz_file = 'lorentz_sim.nii.gz'
        img_sz = [32, 32, 32]
        noise = 0.001

        NewImage(out_file='f0.nii.gz', verbose=vb, img_size=img_sz,
                 grad_dim=0, grad_vals=(-0.5, 0.5)).run()
        NewImage(out_file='fwhm.nii.gz', verbose=vb, img_size=img_sz,
                 grad_dim=1, grad_vals=(0.5, 2.5)).run()
        NewImage(out_file='A.nii.gz', verbose=vb, img_size=img_sz,
                 grad_dim=2, grad_vals=(0.5, 1)).run()

        LorentzianSim(sequence=sequence, pools=pools, in_file=lorentz_file,
                      noise=noise, verbose=vb,
                      DS_f0='f0.nii.gz',
                      DS_fwhm='fwhm.nii.gz',
                      DS_A='A.nii.gz').run()
        Lorentzian(sequence=sequence, pools=pools,
                   in_file=lorentz_file, verbose=vb).run()

        diff_f0 = Diff(in_file='LTZ_DS_f0.nii.gz', baseline='f0.nii.gz',
                       noise=noise, verbose=vb).run()
        diff_fwhm = Diff(in_file='LTZ_DS_fwhm.nii.gz', baseline='fwhm.nii.gz',
                         noise=noise, verbose=vb).run()
        diff_A = Diff(in_file='LTZ_DS_A.nii.gz', baseline='A.nii.gz',
                      noise=noise, verbose=vb).run()
        self.assertLessEqual(diff_f0.outputs.out_diff, 60)
        self.assertLessEqual(diff_fwhm.outputs.out_diff, 20)
        self.assertLessEqual(diff_A.outputs.out_diff, 25)

    def test_lorentzian2(self):
        sat_f0 = [*np.linspace(-40, 40,
                               12).squeeze().tolist(), -0.5, -0.25, 0, 0.25, 0.5]
        sequence = {'MTSat': {'pulse': {'p1': 0.4,
                                        'p2': 0.3,
                                        'bandwidth': 0.39},
                              'TR': 4,
                              'Trf': 0.02,
                              'FA': 5,
                              'sat_f0': sat_f0,
                              'sat_angle': np.repeat(180.0, 17).squeeze().tolist()}}
        pools = [{'name': 'DS',
                  'df0': [0, -2.5, 2.5],
                  'fwhm': [1.0, 1.e-6, 3.0],
                  'A': [0.2, 1.e-3, 1.0],
                  'use_bandwidth': True},
                 {'name': 'MT',
                  'df0': [-2.5, -5.0, -0.5],
                  'fwhm': [50.0, 35.0, 200.0],
                  'A': [0.3, 1.e-3, 1.0]}]

        lorentz_file = 'lorentz_sim.nii.gz'
        img_sz = [32, 32, 32]
        noise = 0.001

        NewImage(out_file='PD.nii.gz', verbose=vb, img_size=img_sz,
                 fill=1.0).run()
        NewImage(out_file='f0.nii.gz', verbose=vb, img_size=img_sz,
                 grad_dim=0, grad_vals=(-0.25, 0.25)).run()

        NewImage(out_file='fwhm.nii.gz', verbose=vb, img_size=img_sz,
                 fill=1.8).run()
        NewImage(out_file='A.nii.gz', verbose=vb, img_size=img_sz,
                 grad_dim=2, grad_vals=(0.3, 0.4)).run()

        NewImage(out_file='MTf0.nii.gz', verbose=vb, img_size=img_sz,
                 fill=-2.2).run()
        NewImage(out_file='MTfwhm.nii.gz', verbose=vb, img_size=img_sz,
                 fill=90).run()
        NewImage(out_file='MTA.nii.gz', verbose=vb, img_size=img_sz,
                 fill=0.4).run()

        LorentzianSim(sequence=sequence, pools=pools, in_file=lorentz_file,
                      noise=noise, verbose=vb,
                      DS_f0='f0.nii.gz',
                      DS_fwhm='fwhm.nii.gz',
                      DS_A='A.nii.gz',
                      MT_fwhm='MTfwhm.nii.gz',
                      MT_f0='f0.nii.gz',
                      MT_A='MTA.nii.gz').run()
        Lorentzian(sequence=sequence, pools=pools,
                   in_file=lorentz_file, verbose=vb).run()

        diff_fwhm = Diff(in_file='LTZ_DS_fwhm.nii.gz', baseline='fwhm.nii.gz',
                         noise=noise, verbose=vb).run()
        diff_A = Diff(in_file='LTZ_DS_A.nii.gz', baseline='A.nii.gz',
                      noise=noise, verbose=vb).run()

        diff_MTfwhm = Diff(in_file='LTZ_MT_fwhm.nii.gz', baseline='MTfwhm.nii.gz',
                           noise=noise, verbose=vb).run()
        diff_MTA = Diff(in_file='LTZ_MT_A.nii.gz', baseline='MTA.nii.gz',
                        noise=noise, verbose=vb).run()

        self.assertLessEqual(diff_fwhm.outputs.out_diff, 15)
        self.assertLessEqual(diff_A.outputs.out_diff, 300)
        self.assertLessEqual(diff_MTfwhm.outputs.out_diff, 15)
        self.assertLessEqual(diff_MTA.outputs.out_diff, 300)

    def test_qMT(self):
        qmt = {'MTSat': {'sat_f0': [1000, 3000, 5000, 7000, 9000, 1000, 3000, 5000, 7000, 9000],
                         'sat_angle': [360, 360, 360, 360, 360, 720, 720, 720, 720, 720],
                         'bw': 100,
                         'FA': 5,
                         'TR': 0.03,
                         'Trf': 0.015,
                         'pulse': {'p1': 0.416, 'p2': 0.295, 'bandwidth': 100}},
               }
        qmt_file = 'qmt_sim.nii.gz'
        t1app = 'qmt_t1app.nii.gz'
        img_sz = [32, 32, 32]
        noise = 0.001

        lineshape_file = '_qmt_lineshape.json'

        Lineshape(out_file=lineshape_file, lineshape='SuperLorentzian',
                  frq_start=500, frq_space=500, frq_count=150).run()

        NewImage(out_file='gM0_a.nii.gz', verbose=vb, img_size=img_sz,
                 fill=1.0).run()
        NewImage(out_file='T1_a_over_T2_a.nii.gz', verbose=vb, img_size=img_sz,
                 grad_dim=0, grad_vals=(5, 15)).run()
        NewImage(out_file='f_over_R_af.nii.gz', verbose=vb, img_size=img_sz,
                 fill=0.1).run()
        NewImage(out_file='T2_b.nii.gz', verbose=vb, img_size=img_sz,
                 fill=12e-6).run()
        NewImage(out_file='RM0a.nii.gz', verbose=vb, img_size=img_sz,
                 grad_dim=1, grad_vals=(1.0, 5.0)).run()
        NewImage(out_file=t1app, verbose=vb, img_size=img_sz,
                 grad_dim=2, grad_vals=(0.5, 1.5)).run()

        qMTSim(sequence=qmt, in_file=qmt_file, t1_map=t1app, lineshape=lineshape_file,
               noise=noise, verbose=vb,
               RM0a='RM0a.nii.gz', f_over_R_af='f_over_R_af.nii.gz', T2_b='T2_b.nii.gz', T1_a_over_T2_a='T1_a_over_T2_a.nii.gz', gM0_a='gM0_a.nii.gz').run()
        qMT(sequence=qmt, in_file=qmt_file, t1_map=t1app,
            lineshape=lineshape_file, verbose=vb).run()

        diff_PD = Diff(in_file='QMT_PD.nii.gz', baseline='gM0_a.nii.gz',
                       noise=noise, verbose=vb).run()
        self.assertLessEqual(diff_PD.outputs.out_diff, 32)

    def test_ZSpec(self):
        NewImage(out_file='zspec_linear.nii.gz', verbose=vb, img_size=[8, 8, 8, 4],
                 grad_dim=3, grad_vals=(-3, 3)).run()
        NewImage(out_file='zspec_zero.nii.gz', verbose=vb,
                 img_size=[8, 8, 8], fill=0).run()
        ZSpec(in_file='zspec_linear.nii.gz',
              in_freqs=[-3, -1, 1, 3],
              out_freqs=[0],
              verbose=vb).run()
        diff_zero = Diff(in_file='zspec_linear_interp.nii.gz', abs_diff=True,
                         baseline='zspec_zero.nii.gz', noise=1, verbose=vb).run()
        self.assertLessEqual(diff_zero.outputs.out_diff, 0.01)


if __name__ == '__main__':
    unittest.main()
