"""
LD Score Calculation Module.

This module provides classes and functions for calculating linkage disequilibrium (LD) scores,
which are useful in genetic studies for understanding the correlation structure of genetic variants.

Classes:
    GenotypeArrayInMemory: Base class for genotype data handling in memory.
    PlinkBEDFile: Class for handling PLINK .bed genotype files.

Functions:
    get_block_lefts(coords, max_dist): Compute indices of leftmost SNPs within a specified distance.
    block_left_to_right(block_left): Convert block left indices to block right indices.

(c) 2015 Brendan Bulik-Sullivan and Hilary Finucane
(c) 2024 Thomas Reimonn
"""

from typing import Callable, Optional, Tuple

import bitarray as ba
import numpy as np


def get_block_lefts(coords: np.ndarray, max_dist: float) -> np.ndarray:
    """
    Compute indices of the leftmost SNPs within a specified maximum distance.

    Args:
        coords (np.ndarray): Array of genomic coordinates (must be sorted).
        max_dist (float): Maximum distance between SNPs to be included in the same window.

    Returns:
        np.ndarray: Array where each element is the index of the leftmost SNP included
                    in the LD score calculation for the corresponding SNP.

    Raises:
        ValueError: If coords is not a one-dimensional array.
    """
    if coords.ndim != 1:
        raise ValueError("coords must be a one-dimensional array.")
    M = len(coords)
    block_left = np.zeros(M, dtype=int)
    j = 0
    for i in range(M):
        while j < M and abs(coords[j] - coords[i]) > max_dist:
            j += 1
        block_left[i] = j
    return block_left


def block_left_to_right(block_left: np.ndarray) -> np.ndarray:
    """
    Convert block left indices to block right indices.

    Args:
        block_left (np.ndarray): Array of block left indices.

    Returns:
        np.ndarray: Array where each element is the index of the rightmost SNP included
                    in the LD score calculation for the corresponding SNP.
    """
    M = len(block_left)
    block_right = np.zeros(M, dtype=int)
    j = 0
    for i in range(M):
        while j < M and block_left[j] <= i:
            j += 1
        block_right[i] = j
    return block_right


class GenotypeArrayInMemory:
    """
    Base class for genotype data handling in memory.

    This class provides methods to read genotype data, filter SNPs and individuals,
    and compute LD scores.

    Attributes:
        m (int): Number of SNPs.
        n (int): Number of individuals.
        df (np.ndarray): SNP metadata array (e.g., chromosome, SNP ID, base pair position).
        colnames (list): Column names for the SNP metadata.
        maf_min (float): Minimum minor allele frequency for filtering.
        geno (bitarray.bitarray): Bitarray representing genotype data.
        kept_snps (list): Indices of SNPs kept after filtering.
        freq (np.ndarray): Allele frequencies of the kept SNPs.
        maf (np.ndarray): Minor allele frequencies of the kept SNPs.
        sqrtpq (np.ndarray): Square root of p * q for each SNP.
    """

    def __init__(
        self,
        fname: str,
        n: int,
        snp_list,
        keep_snps: Optional[np.ndarray] = None,
        keep_indivs: Optional[np.ndarray] = None,
        maf_min: Optional[float] = None,
    ) -> None:
        """
        Initialize the GenotypeArrayInMemory object.

        Args:
            fname (str): Filename of the genotype data.
            n (int): Number of individuals.
            snp_list: SNP list object containing SNP metadata.
            keep_snps (Optional[np.ndarray]): Indices of SNPs to keep.
            keep_indivs (Optional[np.ndarray]): Indices of individuals to keep.
            maf_min (Optional[float]): Minimum minor allele frequency for filtering.

        Raises:
            ValueError: If filtering results in zero individuals or SNPs remaining.
        """
        self.m = len(snp_list.IDList)
        self.n = n
        self.keep_snps = keep_snps
        self.keep_indivs = keep_indivs
        self.df = np.array(snp_list.df[["CHR", "SNP", "BP", "CM"]])
        self.colnames = ["CHR", "SNP", "BP", "CM"]
        self.maf_min = maf_min if maf_min is not None else 0.0
        self._current_snp = 0

        self.nru, self.geno = self._read(fname, self.m, n)

        # Filter individuals
        if self.keep_indivs is not None:
            self.geno, self.m, self.n = self._filter_indivs(self.geno, self.keep_indivs, self.m, self.n)
            if self.n == 0:
                raise ValueError("After filtering, no individuals remain.")
            else:
                print(f"After filtering, {self.n} individuals remain.")

        # Filter SNPs
        self.geno, self.m, self.n, self.kept_snps, self.freq = self._filter_snps_maf(
            self.geno, self.m, self.n, self.maf_min, self.keep_snps
        )
        if self.m == 0:
            raise ValueError("After filtering, no SNPs remain.")
        else:
            print(f"After filtering, {self.m} SNPs remain.")

        self.df = self.df[self.kept_snps, :]
        self.maf = np.minimum(self.freq, 1.0 - self.freq)
        self.sqrtpq = np.sqrt(self.freq * (1.0 - self.freq))
        self.df = np.column_stack((self.df, self.maf))
        self.colnames.append("MAF")

    def _read(self, fname: str, m: int, n: int) -> Tuple[int, ba.bitarray]:
        """
        Read genotype data from a file.

        Args:
            fname (str): Filename of the genotype data.
            m (int): Number of SNPs.
            n (int): Number of individuals.

        Returns:
            Tuple[int, ba.bitarray]: Tuple containing the number of units (nru) and
                                     the genotype bitarray.

        Raises:
            NotImplementedError: Must be implemented in subclasses.
        """
        raise NotImplementedError("Subclasses must implement the _read method.")

    def _filter_indivs(
        self, geno: ba.bitarray, keep_indivs: np.ndarray, m: int, n: int
    ) -> Tuple[ba.bitarray, int, int]:
        """
        Filter individuals from the genotype data.

        Args:
            geno (ba.bitarray): Genotype bitarray.
            keep_indivs (np.ndarray): Indices of individuals to keep.
            m (int): Number of SNPs.
            n (int): Number of individuals.

        Returns:
            Tuple[ba.bitarray, int, int]: Tuple containing the filtered genotype bitarray,
                                          number of SNPs, and new number of individuals.

        Raises:
            NotImplementedError: Must be implemented in subclasses.
        """
        raise NotImplementedError("Subclasses must implement the _filter_indivs method.")

    def _filter_snps_maf(
        self,
        geno: ba.bitarray,
        m: int,
        n: int,
        maf_min: float,
        keep_snps: Optional[np.ndarray],
    ) -> Tuple[ba.bitarray, int, int, list, np.ndarray]:
        """
        Filter SNPs based on minor allele frequency (MAF) and SNP indices.

        Args:
            geno (ba.bitarray): Genotype bitarray.
            m (int): Number of SNPs.
            n (int): Number of individuals.
            maf_min (float): Minimum minor allele frequency.
            keep_snps (Optional[np.ndarray]): Indices of SNPs to keep.

        Returns:
            Tuple containing:
                - ba.bitarray: Filtered genotype bitarray.
                - int: Number of polymorphic SNPs.
                - int: Number of individuals.
                - list: Indices of kept SNPs.
                - np.ndarray: Allele frequencies of kept SNPs.
        """
        raise NotImplementedError("Subclasses must implement the _filter_snps_maf method.")

    def ld_score_var_blocks(self, block_left: np.ndarray, c: int, annot: Optional[np.ndarray] = None) -> np.ndarray:
        """
        Compute an unbiased estimate of LD scores using variable block sizes.

        Args:
            block_left (np.ndarray): Array of block left indices.
            c (int): Chunk size.
            annot (Optional[np.ndarray]): SNP annotations (shape: (m, n_a)).

        Returns:
            np.ndarray: LD scores (shape: (m, n_a)).
        """
        func = lambda x: self._l2_unbiased(x, self.n)
        snp_getter = self.next_snps
        return self._cor_sum_var_blocks(block_left, c, func, snp_getter, annot)

    def ld_score_block_jackknife(
        self, block_left: np.ndarray, c: int, annot: Optional[np.ndarray] = None, jn: int = 10
    ) -> np.ndarray:
        """
        Compute LD scores using block jackknife.

        Args:
            block_left (np.ndarray): Array of block left indices.
            c (int): Chunk size.
            annot (Optional[np.ndarray]): SNP annotations.
            jn (int): Number of jackknife blocks.

        Returns:
            np.ndarray: LD scores with jackknife variance estimates.
        """
        func = lambda x: np.square(x)
        snp_getter = self.next_snps
        return self._cor_sum_block_jackknife(block_left, c, func, snp_getter, annot, jn)

    @staticmethod
    def _l2_unbiased(x: np.ndarray, n: int) -> np.ndarray:
        """
        Compute an unbiased estimate of squared correlation coefficients.

        Args:
            x (np.ndarray): Correlation coefficients.
            n (int): Number of individuals.

        Returns:
            np.ndarray: Unbiased estimate of squared correlation coefficients.

        Notes:
            The unbiased estimator is calculated as:
                l2_unbiased = x^2 - (1 - x^2) / (n - 2)
        """
        denom = n - 2 if n > 2 else n  # Allow n < 2 for testing purposes
        sq = np.square(x)
        return sq - (1 - sq) / denom

    def _cor_sum_var_blocks(
        self,
        block_left: np.ndarray,
        c: int,
        func: Callable[[np.ndarray], np.ndarray],
        snp_getter: Callable[[int], np.ndarray],
        annot: Optional[np.ndarray] = None,
    ) -> np.ndarray:
        """
        General method for calculating sums of transformed Pearson correlation coefficients.

        Args:
            block_left (np.ndarray): Array of block left indices.
            c (int): Chunk size.
            func (Callable[[np.ndarray], np.ndarray]): Function to apply to the correlation matrix.
            snp_getter (Callable[[int], np.ndarray]): Function to retrieve SNPs.
            annot (Optional[np.ndarray]): SNP annotations (shape: (m, n_a)).

        Returns:
            np.ndarray: Summed values after applying the function and weighting by annotations.
        """
        m, n = self.m, self.n
        if annot is None:
            annot = np.ones((m, 1))
        else:
            if annot.shape[0] != m:
                raise ValueError("Incorrect number of SNPs in annotations.")

        n_a = annot.shape[1]  # Number of annotations
        cor_sum = np.zeros((m, n_a))
        block_sizes = np.array(np.arange(m) - block_left)
        block_sizes = np.ceil(block_sizes / c) * c

        b = np.nonzero(block_left > 0)[0]
        b = b[0] if b.size > 0 else m
        b = int(np.ceil(b / c) * c)
        if b > m:
            c = 1
            b = m

        l_a = 0  # Index of leftmost SNP in matrix A
        A = snp_getter(b)
        rfunc_ab = np.zeros((b, c))
        rfunc_bb = np.zeros((c, c))

        # Process chunks inside the block
        for l_b in range(0, b, c):
            B = A[:, l_b : l_b + c]
            np.dot(A.T, B / n, out=rfunc_ab)
            rfunc_ab = func(rfunc_ab)
            cor_sum[l_a : l_a + b, :] += rfunc_ab @ annot[l_b : l_b + c, :]

        # Process chunks to the right of the block
        b0 = b
        md = int(c * np.floor(m / c))
        end = md + 1 if md != m else md

        for l_b in range(b0, end, c):
            old_b = b
            b = int(block_sizes[l_b])
            if l_b > b0 and b > 0:
                A = np.hstack((A[:, old_b - b + c : old_b], B))
                l_a += old_b - b + c
            elif l_b == b0 and b > 0:
                A = A[:, b0 - b : b0]
                l_a = b0 - b
            elif b == 0:
                A = np.empty((n, 0))
                l_a = l_b
            if l_b == md:
                c = m - md
                rfunc_ab = np.zeros((b, c))
                rfunc_bb = np.zeros((c, c))
            if b != old_b:
                rfunc_ab = np.zeros((b, c))

            B = snp_getter(c)
            if np.all(annot[l_a : l_a + b, :] == 0) and np.all(annot[l_b : l_b + c, :] == 0):
                continue

            np.dot(A.T, B / n, out=rfunc_ab)
            rfunc_ab = func(rfunc_ab)
            cor_sum[l_a : l_a + b, :] += rfunc_ab @ annot[l_b : l_b + c, :]
            cor_sum[l_b : l_b + c, :] += (annot[l_a : l_a + b, :].T @ rfunc_ab).T
            np.dot(B.T, B / n, out=rfunc_bb)
            rfunc_bb = func(rfunc_bb)
            cor_sum[l_b : l_b + c, :] += rfunc_bb @ annot[l_b : l_b + c, :]

        return cor_sum

    def next_snps(self, b: int, minor_ref: Optional[bool] = None) -> np.ndarray:
        """
        Retrieve the next b SNPs from the genotype data.

        Args:
            b (int): Number of SNPs to retrieve.
            minor_ref (Optional[bool]): Whether to flip reference alleles to the minor allele.

        Returns:
            np.ndarray: Matrix of normalized genotypes (shape: (n, b)).

        Raises:
            ValueError: If b is not a positive integer or if insufficient SNPs remain.
        """
        raise NotImplementedError("Subclasses must implement the next_snps method.")


class PlinkBEDFile(GenotypeArrayInMemory):
    """
    Class for handling PLINK .bed genotype files.

    This class provides methods to read PLINK .bed files, filter data, and compute LD scores.
    """

    def __init__(
        self,
        fname: str,
        n: int,
        snp_list,
        keep_snps: Optional[np.ndarray] = None,
        keep_indivs: Optional[np.ndarray] = None,
        maf_min: Optional[float] = None,
    ) -> None:
        """
        Initialize the PlinkBEDFile object.

        Args:
            fname (str): Filename of the .bed file.
            n (int): Number of individuals.
            snp_list: SNP list object containing SNP metadata.
            keep_snps (Optional[np.ndarray]): Indices of SNPs to keep.
            keep_indivs (Optional[np.ndarray]): Indices of individuals to keep.
            maf_min (Optional[float]): Minimum minor allele frequency for filtering.
        """
        self._bedcode = {
            2: ba.bitarray("11"),
            9: ba.bitarray("10"),
            1: ba.bitarray("01"),
            0: ba.bitarray("00"),
        }
        super().__init__(
            fname,
            n,
            snp_list,
            keep_snps=keep_snps,
            keep_indivs=keep_indivs,
            maf_min=maf_min,
        )

    def _read(self, fname: str, m: int, n: int) -> Tuple[int, ba.bitarray]:
        """
        Read genotype data from a PLINK .bed file.

        Args:
            fname (str): Filename of the .bed file.
            m (int): Number of SNPs.
            n (int): Number of individuals.

        Returns:
            Tuple[int, ba.bitarray]: Number of units (nru) and genotype bitarray.

        Raises:
            ValueError: If the file format is incorrect or the magic number is unrecognized.
            IOError: If the .bed file is not in SNP-major mode.
        """
        if not fname.endswith(".bed"):
            raise ValueError("Filename must end with '.bed'.")

        with open(fname, "rb") as fh:
            magic_number = ba.bitarray(endian="little")
            magic_number.fromfile(fh, 2)
            bed_mode = ba.bitarray(endian="little")
            bed_mode.fromfile(fh, 1)
            e = (4 - n % 4) if n % 4 != 0 else 0
            nru = n + e

            # Check magic number
            if magic_number != ba.bitarray("0011011011011000"):
                raise IOError("Unrecognized magic number in PLINK .bed file.")

            if bed_mode != ba.bitarray("10000000"):
                raise IOError("PLINK .bed file must be in default SNP-major mode.")

            # Read genotype data
            geno = ba.bitarray(endian="little")
            geno.fromfile(fh)
            self._test_length(geno, m, nru)
            return nru, geno

    @staticmethod
    def _test_length(geno: ba.bitarray, m: int, nru: int) -> None:
        """
        Verify the length of the genotype bitarray.

        Args:
            geno (ba.bitarray): Genotype bitarray.
            m (int): Number of SNPs.
            nru (int): Number of units (number of individuals plus padding).

        Raises:
            IOError: If the actual length does not match the expected length.
        """
        expected_len = 2 * m * nru
        actual_len = len(geno)
        if actual_len != expected_len:
            raise IOError(f"PLINK .bed file has {actual_len} bits; expected {expected_len} bits.")

    def _filter_indivs(
        self, geno: ba.bitarray, keep_indivs: np.ndarray, m: int, n: int
    ) -> Tuple[ba.bitarray, int, int]:
        """
        Filter individuals from the genotype data.

        Args:
            geno (ba.bitarray): Genotype bitarray.
            keep_indivs (np.ndarray): Indices of individuals to keep.
            m (int): Number of SNPs.
            n (int): Number of individuals.

        Returns:
            Tuple[ba.bitarray, int, int]: Filtered genotype bitarray, number of SNPs, and new n.

        Raises:
            ValueError: If keep_indivs indices are out of bounds.
        """
        if np.any(keep_indivs >= n):
            raise ValueError("keep_indivs indices out of bounds.")

        n_new = len(keep_indivs)
        e = (4 - n_new % 4) if n_new % 4 != 0 else 0
        nru_new = n_new + e
        nru = self.nru
        z = ba.bitarray(m * 2 * nru_new, endian="little")
        z.setall(0)
        for idx, i in enumerate(keep_indivs):
            z[2 * idx :: 2 * nru_new] = geno[2 * i :: 2 * nru]
            z[2 * idx + 1 :: 2 * nru_new] = geno[2 * i + 1 :: 2 * nru]
        self.nru = nru_new
        return z, m, n_new

    def _filter_snps_maf(
        self,
        geno: ba.bitarray,
        m: int,
        n: int,
        maf_min: float,
        keep_snps: Optional[np.ndarray],
    ) -> Tuple[ba.bitarray, int, int, list, np.ndarray]:
        """
        Filter SNPs based on MAF and specified SNP indices.

        Args:
            geno (ba.bitarray): Genotype bitarray.
            m (int): Number of SNPs.
            n (int): Number of individuals.
            maf_min (float): Minimum minor allele frequency.
            keep_snps (Optional[np.ndarray]): Indices of SNPs to keep.

        Returns:
            Tuple containing:
                - ba.bitarray: Filtered genotype bitarray.
                - int: Number of polymorphic SNPs.
                - int: Number of individuals.
                - list: Indices of kept SNPs.
                - np.ndarray: Allele frequencies of kept SNPs.
        """
        nru = self.nru
        m_poly = 0
        filtered_geno = ba.bitarray(endian="little")
        if keep_snps is None:
            keep_snps = range(m)
        kept_snps = []
        freq = []
        for idx, j in enumerate(keep_snps):
            z = geno[2 * nru * j : 2 * nru * (j + 1)]
            A = z[0::2]
            B = z[1::2]
            a = A.count()
            b = B.count()
            c = (A & B).count()
            major_ct = b + c
            n_nomiss = n - a + c
            f = major_ct / (2 * n_nomiss) if n_nomiss > 0 else 0
            het_miss_ct = a + b - 2 * c
            if min(f, 1 - f) > maf_min and het_miss_ct < n:
                freq.append(f)
                filtered_geno += z
                m_poly += 1
                kept_snps.append(j)
        return filtered_geno, m_poly, n, kept_snps, np.array(freq)

    def next_snps(self, b: int, minor_ref: Optional[bool] = None) -> np.ndarray:
        """
        Retrieve the next b SNPs from the genotype data.

        Args:
            b (int): Number of SNPs to retrieve.
            minor_ref (Optional[bool]): Whether to flip reference alleles to the minor allele.

        Returns:
            np.ndarray: Matrix of normalized genotypes (shape: (n, b)).

        Raises:
            ValueError: If b is not a positive integer or if insufficient SNPs remain.
        """
        if not isinstance(b, int) or b <= 0:
            raise ValueError("b must be a positive integer.")
        if self._current_snp + b > self.m:
            remaining = self.m - self._current_snp
            raise ValueError(f"{b} SNPs requested; only {remaining} SNPs remain.")

        c = self._current_snp
        n = self.n
        nru = self.nru
        slice_start = 2 * c * nru
        slice_end = 2 * (c + b) * nru
        geno_slice = self.geno[slice_start:slice_end]
        X = np.array(list(geno_slice.decode(self._bedcode)), dtype="float64").reshape((b, nru)).T
        X = X[:n, :]
        Y = np.zeros_like(X)
        for j in range(b):
            snp = X[:, j]
            valid_idx = snp != 9
            avg = np.mean(snp[valid_idx])
            snp[~valid_idx] = avg
            denom = np.std(snp)
            if denom == 0:
                denom = 1.0
            if minor_ref is not None and self.freq[self._current_snp + j] > 0.5:
                denom *= -1
            Y[:, j] = (snp - avg) / denom
        self._current_snp += b
        return Y
