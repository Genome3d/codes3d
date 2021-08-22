import unittest
import pandas as pd
from sqlalchemy import create_engine
import os
import codes3d
import snps
import create_test_data as test


class VariantTest(unittest.TestCase):
    def setUp(self):
        self.output_dir = test.output_dir
        self.snp_list = test.snp_list
        self.snp_rsid_df = test.snp_rsid_df
        self.snp_fragment_df = test.snp_fragment_df
        self.c = codes3d.CODES3D(test.config_fp)
        self.db = create_engine(self.c.postgres_commons, echo=False)
        self.enzymes = test.enzymes
        self.logger = codes3d.Logger(
            #logfile=os.path.join(self.output_dir, 'codes3d.log'),
            verbose=False)

    # Hi-C libraries
    # eQTL tissues
    def test_process_snp_rsid(self):
        df = snps.process_rs_df(
            self.snp_list, self.db
        ).sort_values(by=['id'], ascending=True).reset_index(
            drop=True)
        df[['locus', 'id']] = df[['locus', 'id']].astype('int64')
        self.assertEqual(
            df
            .equals(test.snp_rsid_df),
            True)

    def test_get_snp_fragments(self):
        df = snps.get_snp_fragments(
            self.snp_rsid_df,
            self.enzymes,
            self.db).sort_values(
            by=['id', 'chr', 'frag_id', 'enzyme']
        ).reset_index(drop=True)
        self.assertEqual(
            df.equals(test.snp_fragment_df),
            True
        )

    def test_merged_snps(self):
        df, _, _, _ = snps.process_snp_input(
            [test.merged_snps_fp], self.output_dir, self.c.postgres_commons,
            self.c.rs_merge_arch_fp, self.logger)
        self.assertFalse(df.empty)

    def tearDown(self):
        self.snp_df = None
        if os.path.exists(self.output_dir):
            os.removedirs(self.output_dir)


if __name__ == '__main__':
    unittest.main()
