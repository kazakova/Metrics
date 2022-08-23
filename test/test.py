import unittest
from tempfile import TemporaryDirectory
import pandas as pd
import numpy as np
import argparse
import os
import sys
sys.path.insert(0, os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'QRePS'))
from stattest_metrics import QRePS

class QRePSResultTest(unittest.TestCase):
    def setUp(self):
        res_dir = os.path.dirname(os.path.abspath(__file__))
        self.DBTRG_q = pd.read_csv(os.path.join(res_dir, 'Quant_res_DBTRG_I,DBTRG_K.tsv'), sep = '\t', index_col = 0)
        self.DBTRG_m = pd.read_csv(os.path.join(res_dir, 'metrics_DBTRG_I,DBTRG_K.tsv'), sep = '\t')
        self.A172_q = pd.read_csv(os.path.join(res_dir, 'Quant_res_A172_I,A172_K.tsv'), sep = '\t', index_col = 0)
        self.A172_m = pd.read_csv(os.path.join(res_dir, 'metrics_A172_I,A172_K.tsv'), sep = '\t')

    def test_qreps(self):
        true_q = [self.DBTRG_q, self.A172_q]
        true_m = [self.DBTRG_m, self.A172_m]
        datadir =  os.path.join(
                os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'example')

        with TemporaryDirectory() as tmpdir:
            args = argparse.Namespace(sample_file=os.path.join(datadir, 'a172_dbtrg_sample.csv'), labels=['DBTRG_I,DBTRG_K','A172_I,A172_K'], input_dir=datadir, output_dir=tmpdir, imputation='kNN', thresholds='dynamic', regulation='UP', pattern = '_protein_groups.tsv', species = '9606', fold_change=2, alpha=0.01)
            self.quants, self.gos, self.metrics = QRePS(args)
        for i in range(2):
            print(i)
            q = self.quants[i]
            m = self.metrics[i]
            self.assertEqual(true_q[i].index.tolist(), q.index.tolist())
            n_proteins = q.shape[0]
            n_err = ((true_q[i]['log2(fold_change)'] - q['log2(fold_change)']).abs() > 1e-5).sum()
            self.assertLess(n_err, n_proteins//100)
#             for j in q.index:
#                 if abs(true_q[i].loc[j, 'log2(fold_change)'] - q.loc[j, 'log2(fold_change)']) > 1e-5:
#                     print(j)
#                 self.assertAlmostEqual()
#             np.testing.assert_almost_equal(q['log2(fold_change)'].values, true_q[i]['log2(fold_change)'].values, 5)
#             np.testing.assert_almost_equal(q['-log10(fdr_BH)'].values, true_q[i]['-log10(fdr_BH)'].values, 5)
            np.testing.assert_almost_equal(m.iloc[0].values, true_m[i].iloc[0].values, 5)