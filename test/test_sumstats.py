"""
Unit tests for the sumstats module.

This module contains unit tests for the functions in the sumstats.py module,
ensuring correctness and robustness.

Usage:
    Run this script with a test runner like pytest or nose.

Author: [Your Name]
"""

import os
from unittest import TestCase

import numpy as np
import pandas as pd
from nose.plugins.attrib import attr
from numpy.testing import (
    assert_allclose,
    assert_almost_equal,
    assert_array_almost_equal,
    assert_array_equal,
)
from pandas.testing import assert_series_equal

import ldscore.parse as ps
import ldscore.sumstats as s
from ldsc import parser

# Constants
TEST_DIR = os.path.dirname(__file__)
NUM_REPETITIONS = 500
s.NUM_CHROMOSOMES = 2  # Mocking chromosomes for testing purposes


class MockLogger:
    """
    Mock logger class for capturing log outputs during testing.
    """

    def log(self, message: str) -> None:
        # For debugging purposes, you can print the message
        # print(message)
        pass


logger = MockLogger()
args = parser.parse_args("")


def get_attr(attr: str):
    """
    Helper function to get an attribute from an object.

    Args:
        attr (str): Attribute name.

    Returns:
        Callable: Function that retrieves the attribute from an object.
    """
    return lambda obj: getattr(obj, attr, float("nan"))


class TestSumstatsFunctions(TestCase):
    """
    Unit tests for individual functions in sumstats.py.
    """

    def test_check_ld_condition_number(self):
        """
        Test the check_ld_condition_number function.
        """
        ld_matrix = np.ones((2, 2))
        ld_matrix[1, 1] += 1e-5
        args.invert_anyway = False
        with self.assertRaises(ValueError):
            s.check_ld_condition_number(args, logger, ld_matrix)
        args.invert_anyway = True
        # Should not raise an error
        s.check_ld_condition_number(args, logger, ld_matrix)

    def test_check_variance(self):
        """
        Test the check_variance function for removing zero-variance LD Scores.
        """
        ld_scores = pd.DataFrame({"SNP": ["a", "b", "c"], "LD1": np.ones(3), "LD2": np.arange(3)})
        m_annot = np.array([[1, 2]])
        m_annot_updated, ld_scores_updated, novar_cols = s.check_variance(logger, m_annot, ld_scores)
        self.assertEqual(m_annot_updated.shape, (1, 1))
        assert_array_equal(m_annot_updated, [[2]])
        assert_allclose(ld_scores_updated.iloc[:, 1].values, [0, 1, 2])
        assert_array_equal(novar_cols.values, [True, False])

    def test_align_alleles(self):
        """
        Test the align_alleles function for aligning Z-scores based on allele orientation.
        """
        z_scores = pd.Series(np.ones(6))
        alleles = pd.Series(["ACAC", "TGTG", "GTGT", "AGCT", "AGTC", "TCTC"])
        aligned_z_scores = s.align_alleles(z_scores, alleles)
        expected_z_scores = pd.Series([1.0, 1, 1, -1, 1, 1])
        assert_series_equal(aligned_z_scores, expected_z_scores)

    def test_filter_alleles(self):
        """
        Test the filter_alleles function for identifying valid SNPs.
        """
        alleles = pd.Series(["ATAT", "ATAG", "DIID", "ACAC"])
        valid_indices = s.filter_alleles(alleles)
        expected_indices = pd.Series([False, False, False, True])
        assert_series_equal(valid_indices, expected_indices)

    def test_read_annotation_matrix(self):
        """
        Test reading the annotation matrix from files.
        """
        ref_ld_chr = None
        ref_ld = os.path.join(TEST_DIR, "annot_test/test")
        overlap_matrix, m_tot = s.read_chr_split_files(
            ref_ld_chr, ref_ld, logger, "annot matrix", ps.annot, frqfile=None
        )
        assert_array_equal(overlap_matrix, np.array([[1, 0, 0], [0, 2, 2], [0, 2, 2]]))
        assert_array_equal(m_tot, np.array(3))

        frqfile = os.path.join(TEST_DIR, "annot_test/test1")
        overlap_matrix, m_tot = s.read_chr_split_files(
            ref_ld_chr, ref_ld, logger, "annot matrix", ps.annot, frqfile=frqfile
        )
        assert_array_equal(overlap_matrix, np.array([[1, 0, 0], [0, 1, 1], [0, 1, 1]]))
        assert_array_equal(m_tot, np.array(2))

    def test_valid_snps(self):
        """
        Test the VALID_SNPS set for correctness.
        """
        expected_valid_snps = {"AC", "AG", "CA", "CT", "GA", "GT", "TC", "TG"}
        self.assertEqual(expected_valid_snps, s.VALID_SNPS)

    def test_bases(self):
        """
        Test the BASES set for correctness.
        """
        expected_bases = {"A", "T", "G", "C"}
        self.assertEqual(expected_bases, set(s.BASES))

    def test_complement(self):
        """
        Test the COMPLEMENT dictionary for correctness.
        """
        expected_complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
        self.assertEqual(expected_complement, s.COMPLEMENT)

    def test_warn_if_few_snps(self):
        """
        Test that the warn_if_few_snps function executes without error.
        """
        s.warn_if_few_snps(logger, pd.DataFrame({"SNP": [1]}))

    def test_match_alleles(self):
        """
        Test the MATCH_ALLELES set for correctness.
        """
        expected_match_alleles = {
            "ACAC",
            "ACCA",
            "ACGT",
            "ACTG",
            "AGAG",
            "AGCT",
            "AGGA",
            "AGTC",
            "CAAC",
            "CACA",
            "CAGT",
            "CATG",
            "CTAG",
            "CTCT",
            "CTGA",
            "CTTC",
            "GAAG",
            "GACT",
            "GAGA",
            "GATC",
            "GTAC",
            "GTCA",
            "GTGT",
            "GTTG",
            "TCAG",
            "TCCT",
            "TCGA",
            "TCTC",
            "TGAC",
            "TGCA",
            "TGGT",
            "TGTG",
        }
        self.assertEqual(expected_match_alleles, s.MATCH_ALLELES)

    def test_flip_alleles(self):
        """
        Test the FLIP_ALLELES dictionary for correctness.
        """
        expected_flip_alleles = {
            "ACAC": False,
            "ACCA": True,
            "ACGT": True,
            "ACTG": False,
            "AGAG": False,
            "AGCT": True,
            "AGGA": True,
            "AGTC": False,
            "CAAC": True,
            "CACA": False,
            "CAGT": False,
            "CATG": True,
            "CTAG": True,
            "CTCT": False,
            "CTGA": False,
            "CTTC": True,
            "GAAG": True,
            "GACT": False,
            "GAGA": False,
            "GATC": True,
            "GTAC": True,
            "GTCA": False,
            "GTGT": False,
            "GTTG": True,
            "TCAG": False,
            "TCCT": True,
            "TCGA": True,
            "TCTC": False,
            "TGAC": False,
            "TGCA": True,
            "TGGT": True,
            "TGTG": False,
        }
        self.assertEqual(expected_flip_alleles, s.FLIP_ALLELES)

    def test_strand_ambiguous(self):
        """
        Test the STRAND_AMBIGUOUS dictionary for correctness.
        """
        expected_strand_ambiguous = {
            "AC": False,
            "AG": False,
            "AT": True,
            "CA": False,
            "CG": True,
            "CT": False,
            "GA": False,
            "GC": True,
            "GT": False,
            "TA": True,
            "TC": False,
            "TG": False,
        }
        self.assertEqual(expected_strand_ambiguous, s.STRAND_AMBIGUOUS)


@attr("rg")
@attr("slow")
class TestGeneticCorrelationStatistical(TestCase):
    """
    Statistical tests for genetic correlation estimation.
    """

    @classmethod
    def setUpClass(cls):
        """
        Set up test cases by running genetic correlation estimation.
        """
        args = parser.parse_args("")
        args.ref_ld = os.path.join(TEST_DIR, "simulate_test/ldscore/twold_onefile")
        args.w_ld = os.path.join(TEST_DIR, "simulate_test/ldscore/w")
        args.rg = ",".join([os.path.join(TEST_DIR, f"simulate_test/sumstats/{i}") for i in range(NUM_REPETITIONS)])
        args.out = os.path.join(TEST_DIR, "simulate_test/1")
        cls.rg_results = s.estimate_genetic_correlation(args, logger)
        args.intercept_gencov = ",".join(["0"] * NUM_REPETITIONS)
        args.intercept_h2 = ",".join(["1"] * NUM_REPETITIONS)
        cls.rg_results_no_intercept = s.estimate_genetic_correlation(args, logger)

    def test_rg_ratio(self):
        """
        Test that the mean rg_ratio is close to 0.
        """
        rg_ratios = [get_attr("rg_ratio")(rg) for rg in self.rg_results]
        mean_rg_ratio = np.nanmean(rg_ratios)
        self.assertAlmostEqual(mean_rg_ratio, 0, delta=0.02)

    def test_rg_ratio_no_intercept(self):
        """
        Test that the mean rg_ratio without intercept is close to 0.
        """
        rg_ratios = [get_attr("rg_ratio")(rg) for rg in self.rg_results_no_intercept]
        mean_rg_ratio = np.nanmean(rg_ratios)
        self.assertAlmostEqual(mean_rg_ratio, 0, delta=0.02)

    def test_rg_se(self):
        """
        Test that the standard error of rg matches the standard deviation of rg_ratio.
        """
        rg_ratios = [get_attr("rg_ratio")(rg) for rg in self.rg_results]
        rg_ses = [get_attr("rg_se")(rg) for rg in self.rg_results]
        self.assertAlmostEqual(np.nanmean(rg_ses), np.nanstd(rg_ratios), delta=0.02)

    def test_rg_se_no_intercept(self):
        """
        Test that the standard error of rg without intercept matches the standard deviation of rg_ratio.
        """
        rg_ratios = [get_attr("rg_ratio")(rg) for rg in self.rg_results_no_intercept]
        rg_ses = [get_attr("rg_se")(rg) for rg in self.rg_results_no_intercept]
        self.assertAlmostEqual(np.nanmean(rg_ses), np.nanstd(rg_ratios), delta=0.02)

    # Additional tests for genetic covariance and other statistics can be added here.


@attr("h2")
@attr("slow")
class TestHeritabilityStatistical(TestCase):
    """
    Statistical tests for heritability estimation.
    """

    @classmethod
    def setUpClass(cls):
        """
        Set up test cases by running heritability estimation.
        """
        args = parser.parse_args("")
        args.ref_ld = os.path.join(TEST_DIR, "simulate_test/ldscore/twold_onefile")
        args.w_ld = os.path.join(TEST_DIR, "simulate_test/ldscore/w")
        args.chisq_max = 99999
        cls.h2_results = []
        cls.h2_results_no_intercept = []
        for i in range(NUM_REPETITIONS):
            args.intercept_h2 = None
            args.h2 = os.path.join(TEST_DIR, f"simulate_test/sumstats/{i}")
            args.out = os.path.join(TEST_DIR, "simulate_test/1")
            h2 = s.estimate_heritability(args, logger)
            cls.h2_results.append(h2)
            args.intercept_h2 = 1
            h2_no_intercept = s.estimate_heritability(args, logger)
            cls.h2_results_no_intercept.append(h2_no_intercept)

    def test_total_heritability(self):
        """
        Test that the mean total heritability estimate is close to 0.9.
        """
        total_h2 = [get_attr("tot")(h2) for h2 in self.h2_results]
        mean_total_h2 = np.nanmean(total_h2)
        self.assertAlmostEqual(mean_total_h2, 0.9, delta=0.05)

    def test_total_heritability_no_intercept(self):
        """
        Test that the mean total heritability estimate without intercept is close to 0.9.
        """
        total_h2 = [get_attr("tot")(h2) for h2 in self.h2_results_no_intercept]
        mean_total_h2 = np.nanmean(total_h2)
        self.assertAlmostEqual(mean_total_h2, 0.9, delta=0.05)

    def test_total_heritability_se(self):
        """
        Test that the standard error of total heritability matches the standard deviation of estimates.
        """
        total_h2 = [get_attr("tot")(h2) for h2 in self.h2_results]
        total_h2_se = [get_attr("tot_se")(h2) for h2 in self.h2_results]
        self.assertAlmostEqual(np.nanmean(total_h2_se), np.nanstd(total_h2), delta=0.05)

    def test_total_heritability_se_no_intercept(self):
        """
        Test that the standard error of total heritability without intercept matches the standard deviation of estimates.
        """
        total_h2 = [get_attr("tot")(h2) for h2 in self.h2_results_no_intercept]
        total_h2_se = [get_attr("tot_se")(h2) for h2 in self.h2_results_no_intercept]
        self.assertAlmostEqual(np.nanmean(total_h2_se), np.nanstd(total_h2), delta=0.05)

    # Additional tests for category-specific heritability and other statistics can be added here.


class TestEstimateFunctions(TestCase):
    """
    Tests for the estimate_h2 and estimate_rg functions.
    """

    def test_estimate_h2_with_M(self):
        """
        Test estimate_h2 function with provided M values.
        """
        args = parser.parse_args("")
        args.ref_ld = os.path.join(TEST_DIR, "simulate_test/ldscore/oneld_onefile")
        args.w_ld = os.path.join(TEST_DIR, "simulate_test/ldscore/w")
        args.h2 = os.path.join(TEST_DIR, "simulate_test/sumstats/1")
        args.out = os.path.join(TEST_DIR, "simulate_test/1")
        args.print_cov = True
        args.print_delete_vals = True
        h2_result = s.estimate_heritability(args, logger)
        with open(os.path.join(TEST_DIR, "simulate_test/ldscore/oneld_onefile.l2.M_5_50"), "r") as f:
            m_value = f.read().strip()
        args.M = m_value
        h2_result_with_M = s.estimate_heritability(args, logger)
        assert_array_almost_equal(h2_result.tot, h2_result_with_M.tot)
        assert_array_almost_equal(h2_result.tot_se, h2_result_with_M.tot_se)

    def test_estimate_rg_with_M(self):
        """
        Test estimate_rg function with provided M values.
        """
        args = parser.parse_args("")
        args.ref_ld = os.path.join(TEST_DIR, "simulate_test/ldscore/oneld_onefile")
        args.w_ld = os.path.join(TEST_DIR, "simulate_test/ldscore/w")
        args.rg = ",".join([os.path.join(TEST_DIR, "simulate_test/sumstats/1") for _ in range(2)])
        args.out = os.path.join(TEST_DIR, "simulate_test/1")
        rg_result = s.estimate_genetic_correlation(args, logger)[0]
        with open(os.path.join(TEST_DIR, "simulate_test/ldscore/oneld_onefile.l2.M_5_50"), "r") as f:
            m_value = f.read().strip()
        args.M = m_value
        rg_result_with_M = s.estimate_genetic_correlation(args, logger)[0]
        assert_array_almost_equal(rg_result.rg_ratio, rg_result_with_M.rg_ratio)
        assert_array_almost_equal(rg_result.rg_se, rg_result_with_M.rg_se)

    def test_no_check_alleles(self):
        """
        Test estimate_rg function with the no_check_alleles option.
        """
        args = parser.parse_args("")
        args.ref_ld = os.path.join(TEST_DIR, "simulate_test/ldscore/oneld_onefile")
        args.w_ld = os.path.join(TEST_DIR, "simulate_test/ldscore/w")
        args.rg = ",".join([os.path.join(TEST_DIR, "simulate_test/sumstats/1") for _ in range(2)])
        args.out = os.path.join(TEST_DIR, "simulate_test/1")
        rg_result = s.estimate_genetic_correlation(args, logger)[0]
        args.no_check_alleles = True
        rg_result_no_check = s.estimate_genetic_correlation(args, logger)[0]
        self.assertEqual(rg_result.rg_ratio, rg_result_no_check.rg_ratio)
        assert_almost_equal(rg_result.gencov.tot, rg_result_no_check.gencov.tot)

    def test_two_step_h2(self):
        """
        Test estimate_heritability with different two-step estimator cutoffs.
        """
        args = parser.parse_args("")
        args.ref_ld = os.path.join(TEST_DIR, "simulate_test/ldscore/oneld_onefile")
        args.w_ld = os.path.join(TEST_DIR, "simulate_test/ldscore/w")
        args.h2 = os.path.join(TEST_DIR, "simulate_test/sumstats/1")
        args.out = os.path.join(TEST_DIR, "simulate_test/1")
        args.chisq_max = 9999999
        args.two_step = 999
        h2_result = s.estimate_heritability(args, logger)
        args.two_step = 99999
        h2_result_large_cutoff = s.estimate_heritability(args, logger)
        assert_allclose(h2_result.tot, h2_result_large_cutoff.tot, atol=1e-5)

    def test_two_step_rg(self):
        """
        Test estimate_genetic_correlation with different two-step estimator cutoffs.
        """
        args = parser.parse_args("")
        args.ref_ld_chr = os.path.join(TEST_DIR, "simulate_test/ldscore/oneld_onefile")
        args.w_ld = os.path.join(TEST_DIR, "simulate_test/ldscore/w")
        args.rg = ",".join([os.path.join(TEST_DIR, "simulate_test/sumstats/1") for _ in range(2)])
        args.out = os.path.join(TEST_DIR, "simulate_test/rg")
        args.two_step = 999
        rg_result = s.estimate_genetic_correlation(args, logger)[0]
        args.two_step = 99999
        rg_result_large_cutoff = s.estimate_genetic_correlation(args, logger)[0]
        assert_allclose(rg_result.rg_ratio, rg_result_large_cutoff.rg_ratio, atol=1e-5)
        assert_allclose(rg_result.gencov.tot, rg_result_large_cutoff.gencov.tot, atol=1e-5)
