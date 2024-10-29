"""
Unit Tests for LD Score Calculation Module.

This module contains unit tests for the ldscore.py module, specifically testing
the functions and classes related to LD score calculation using PLINK .bed files.

Tests:
    - test_get_block_lefts: Tests the get_block_lefts function.
    - test_block_left_to_right: Tests the block_left_to_right function.
    - TestPlinkBEDFile: Unit tests for the PlinkBEDFile class.

Note:
    Ensure that the test data files are located in the 'test/plink_test' directory.

"""

from unittest import TestCase

import bitarray as ba
import numpy as np

import ldscore.ldscore as ld
import ldscore.parse as ps


def test_get_block_lefts():
    """
    Test the get_block_lefts function with various inputs.
    """
    test_cases = [
        (np.arange(1, 6), 5, np.zeros(5, dtype=int)),
        (np.arange(1, 6), 0, np.arange(5)),
        (np.array([1, 4, 6, 7, 7, 8]), 2, np.array([0, 1, 1, 2, 2, 2])),
    ]
    for coords, max_dist, expected in test_cases:
        result = ld.get_block_lefts(coords, max_dist)
        assert np.array_equal(result, expected), f"Failed for coords={coords}, max_dist={max_dist}"


def test_block_left_to_right():
    """
    Test the block_left_to_right function with various inputs.
    """
    test_cases = [
        (np.array([0, 0, 0, 0, 0]), np.array([5, 5, 5, 5, 5])),
        (np.array([0, 1, 2, 3, 4, 5]), np.array([1, 2, 3, 4, 5, 6])),
        (np.array([0, 0, 2, 2]), np.array([2, 2, 4, 4])),
    ]
    for block_left, expected in test_cases:
        result = ld.block_left_to_right(block_left)
        assert np.array_equal(result, expected), f"Failed for block_left={block_left}"


class TestPlinkBEDFile(TestCase):
    """
    Unit tests for the PlinkBEDFile class in ldscore.py.
    """

    def setUp(self):
        """
        Set up the test environment by initializing necessary variables.
        """
        self.m = 8  # Total number of SNPs in test data
        self.n = 5  # Total number of individuals in test data
        self.bim = ps.PlinkBIMFile("test/plink_test/plink.bim")

    def test_bed_initialization(self):
        """
        Test the initialization of the PlinkBEDFile class.
        """
        bed = ld.PlinkBEDFile("test/plink_test/plink.bed", self.n, self.bim)
        # After filtering monomorphic SNPs, m should be 4
        self.assertEqual(bed.m, 4, "Number of SNPs after filtering should be 4.")
        # No individuals should be removed
        self.assertEqual(bed.n, self.n, "Number of individuals should remain unchanged.")
        # Check the length of the genotype bitarray
        expected_length = 2 * bed.m * bed.nru
        self.assertEqual(len(bed.geno), expected_length, "Genotype bitarray length mismatch.")
        # Check allele frequencies
        expected_freq = np.array([0.6, 0.6, 0.625, 0.625])
        np.testing.assert_array_almost_equal(
            bed.freq, expected_freq, err_msg="Allele frequencies do not match expected values."
        )

    def test_filter_snps(self):
        """
        Test SNP filtering in PlinkBEDFile.
        """
        keep_snps = [1, 4]
        bed = ld.PlinkBEDFile("test/plink_test/plink.bed", self.n, self.bim, keep_snps=np.array(keep_snps))
        # Only SNP index 1 should remain after filtering (since SNP at index 4 is monomorphic)
        self.assertEqual(bed.m, 1, "Number of SNPs after filtering should be 1.")
        self.assertEqual(bed.n, self.n, "Number of individuals should remain unchanged.")
        # Test the genotype bitarray (cannot test pad bits)
        expected_bits = ba.bitarray("0001011111")
        self.assertEqual(
            bed.geno[0:10],
            expected_bits,
            "Genotype bitarray does not match expected values after SNP filtering.",
        )

    def test_filter_individuals(self):
        """
        Test individual filtering in PlinkBEDFile.
        """
        keep_indivs = [0, 1]
        bed = ld.PlinkBEDFile("test/plink_test/plink.bed", self.n, self.bim, keep_indivs=np.array(keep_indivs))
        self.assertEqual(bed.m, 2, "Number of SNPs should be 2 after filtering monomorphic SNPs.")
        self.assertEqual(bed.n, 2, "Number of individuals after filtering should be 2.")
        # Test the genotype bitarray (cannot test pad bits)
        expected_bits_snp1 = ba.bitarray("0001")
        expected_bits_snp2 = ba.bitarray("0001")
        self.assertEqual(
            bed.geno[0:4],
            expected_bits_snp1,
            "Genotype bitarray for SNP 1 does not match expected values after individual filtering.",
        )
        self.assertEqual(
            bed.geno[8:12],
            expected_bits_snp2,
            "Genotype bitarray for SNP 2 does not match expected values after individual filtering.",
        )

    def test_filter_individuals_and_snps(self):
        """
        Test simultaneous SNP and individual filtering in PlinkBEDFile.
        """
        keep_indivs = [0, 1]
        keep_snps = [1, 5]
        bed = ld.PlinkBEDFile(
            "test/plink_test/plink.bed",
            self.n,
            self.bim,
            keep_snps=np.array(keep_snps),
            keep_indivs=np.array(keep_indivs),
        )
        # Only SNP at index 1 should remain after filtering
        self.assertEqual(bed.m, 1, "Number of SNPs after filtering should be 1.")
        self.assertEqual(bed.n, 2, "Number of individuals after filtering should be 2.")
        expected_bits = ba.bitarray("0001")
        self.assertEqual(
            bed.geno[0:4],
            expected_bits,
            "Genotype bitarray does not match expected values after filtering.",
        )

    def test_bad_filename(self):
        """
        Test error handling when an incorrect filename is provided.
        """
        with self.assertRaises(ValueError):
            ld.PlinkBEDFile("test/plink_test/plink.bim", self.n, self.bim)

    def test_next_snps_errors(self):
        """
        Test error handling in the next_snps method.
        """
        bed = ld.PlinkBEDFile("test/plink_test/plink.bed", self.n, self.bim)
        with self.assertRaises(ValueError):
            bed.next_snps(0)
        with self.assertRaises(ValueError):
            bed.next_snps(5)  # Requesting more SNPs than available

    def test_next_snps(self):
        """
        Test the next_snps method for retrieving SNPs.
        """
        for b in [1, 2, 3]:
            bed = ld.PlinkBEDFile("test/plink_test/plink.bed", self.n, self.bim)
            x = bed.next_snps(b)
            self.assertEqual(x.shape, (self.n, b), f"Shape of SNP matrix should be ({self.n}, {b}).")
            np.testing.assert_array_almost_equal(
                np.mean(x, axis=0),
                np.zeros(b),
                decimal=2,
                err_msg="Mean of SNP matrix columns should be approximately zero.",
            )
            np.testing.assert_array_almost_equal(
                np.std(x, axis=0),
                np.ones(b),
                decimal=2,
                err_msg="Standard deviation of SNP matrix columns should be approximately one.",
            )

    def test_next_snps_minor_ref(self):
        """
        Test the next_snps method with minor allele as the reference.
        """
        b = 4
        bed = ld.PlinkBEDFile("test/plink_test/plink.bed", self.n, self.bim)
        x = bed.next_snps(b)
        bed._current_snp -= b  # Reset the current SNP index
        y = bed.next_snps(b, minor_ref=True)
        np.testing.assert_array_almost_equal(
            x, -y, decimal=5, err_msg="SNP matrices should be negatives of each other."
        )
