"""
Module for handling summary statistics and performing LD Score regression.

This module deals with loading the necessary data for LD Score regression from files,
checking the validity of the inputs, and performing the regression analysis.

(c) 2014 Brendan Bulik-Sullivan and Hilary Finucane
(c) 2024 Thomas Reimonn
"""

import copy
import itertools
import os
import traceback
from typing import Any, Callable, Dict, List, Optional, Set, Tuple, Union

import numpy as np
import pandas as pd
from scipy import stats

from . import parse as ps
from . import regressions as reg

# Constants
NUM_CHROMOSOMES = 22

# Complementary DNA bases
COMPLEMENT: Dict[str, str] = {"A": "T", "T": "A", "C": "G", "G": "C"}

# DNA bases
BASES: List[str] = list(COMPLEMENT.keys())

# Dictionary indicating whether a SNP is strand ambiguous
STRAND_AMBIGUOUS: Dict[str, bool] = {
    "".join(alleles): alleles[0] == COMPLEMENT[alleles[1]]
    for alleles in itertools.product(BASES, BASES)
    if alleles[0] != alleles[1]
}

# Set of valid SNPs (pairs of alleles that are not strand ambiguous)
VALID_SNPS: Set[str] = {
    "".join(alleles)
    for alleles in itertools.product(BASES, BASES)
    if alleles[0] != alleles[1] and not STRAND_AMBIGUOUS["".join(alleles)]
}

# Set of allele combinations indicating matching alleles (allowing for strand or reference allele flip)
MATCH_ALLELES: Set[str] = {
    "".join(alleles1 + alleles2)
    for alleles1, alleles2 in itertools.product(VALID_SNPS, repeat=2)
    if (
        # Strand and reference match
        (alleles1[0] == alleles2[0] and alleles1[1] == alleles2[1])
        or
        # Reference match, strand flip
        (alleles1[0] == COMPLEMENT[alleles2[0]] and alleles1[1] == COMPLEMENT[alleles2[1]])
        or
        # Reference flip, strand match
        (alleles1[0] == alleles2[1] and alleles1[1] == alleles2[0])
        or
        # Reference flip, strand flip
        (alleles1[0] == COMPLEMENT[alleles2[1]] and alleles1[1] == COMPLEMENT[alleles2[0]])
    )
}

# Dictionary indicating whether SNP1 has the same alleles as SNP2 with reference allele flip
FLIP_ALLELES: Dict[str, bool] = {
    "".join(alleles1 + alleles2): (
        # Strand match with reference flip
        (alleles1[0] == alleles2[1] and alleles1[1] == alleles2[0])
        or
        # Strand flip with reference flip
        (alleles1[0] == COMPLEMENT[alleles2[1]] and alleles1[1] == COMPLEMENT[alleles2[0]])
    )
    for alleles1, alleles2 in itertools.product(VALID_SNPS, repeat=2)
    if "".join(alleles1 + alleles2) in MATCH_ALLELES
}


def split_paths(path_string: str) -> List[str]:
    """
    Split a comma-separated string of file paths into a list,
    expanding user (~) and environment variables.

    Args:
        path_string (str): Comma-separated file paths.

    Returns:
        List[str]: List of expanded file paths.
    """
    path_list = path_string.split(",")
    path_list = [os.path.expanduser(os.path.expandvars(path)) for path in path_list]
    return path_list


def select_and_log(
    data: pd.DataFrame, indices: Union[pd.Series, np.ndarray], logger: Any, message: str
) -> pd.DataFrame:
    """
    Filter the DataFrame to rows where indices are True, and log the number of SNPs remaining.

    Args:
        data (pd.DataFrame): The DataFrame to filter.
        indices (Union[pd.Series, np.ndarray]): Boolean array indicating rows to keep.
        logger (Any): Logger object for logging messages.
        message (str): Message template with a placeholder {N} for the number of SNPs.

    Returns:
        pd.DataFrame: Filtered DataFrame.

    Raises:
        ValueError: If no SNPs remain after filtering.
    """
    num_remaining = indices.sum()
    if num_remaining == 0:
        raise ValueError(message.format(N=0))
    else:
        filtered_data = data[indices]
        logger.log(message.format(N=num_remaining))
    return filtered_data


def smart_merge(data1: pd.DataFrame, data2: pd.DataFrame) -> pd.DataFrame:
    """
    Merge two DataFrames on 'SNP' column. If the 'SNP' columns are equal,
    use concat for efficiency.

    Args:
        data1 (pd.DataFrame): First DataFrame.
        data2 (pd.DataFrame): Second DataFrame.

    Returns:
        pd.DataFrame: Merged DataFrame.
    """
    if len(data1) == len(data2) and (data1.index == data2.index).all() and (data1["SNP"] == data2["SNP"]).all():
        data1 = data1.reset_index(drop=True)
        data2 = data2.reset_index(drop=True).drop(["SNP"], axis=1)
        merged_data = pd.concat([data1, data2], axis=1)
    else:
        merged_data = pd.merge(data1, data2, how="inner", on="SNP")
    return merged_data


def read_reference_ld_scores(args: Any, logger: Any) -> pd.DataFrame:
    """
    Read reference LD Scores from files specified in args.

    Args:
        args (argparse.Namespace): Command-line arguments.
        logger (Any): Logger object.

    Returns:
        pd.DataFrame: Reference LD Scores DataFrame.

    Raises:
        ValueError: If there is an error reading the LD Scores.
    """
    try:
        if args.ref_ld:
            logger.log(f"Reading reference panel LD Scores from {args.ref_ld} ...")
            ref_ld = ps.ldscore_fromlist(split_paths(args.ref_ld))
        elif args.ref_ld_chr:
            pattern = ps.sub_chr(args.ref_ld_chr, "[1-22]")
            logger.log(f"Reading reference panel LD Scores from {pattern} ...")
            ref_ld = ps.ldscore_fromlist(split_paths(args.ref_ld_chr))
        else:
            raise ValueError("No reference LD Scores provided.")
    except Exception as e:
        logger.log("Error parsing reference LD Scores.")
        raise e

    logger.log(f"Read reference panel LD Scores for {len(ref_ld)} SNPs.")
    return ref_ld


def read_annotation_matrix(args: Any, logger: Any) -> Tuple[pd.DataFrame, np.ndarray]:
    """
    Read annotation matrix from files specified in args.

    Args:
        args (argparse.Namespace): Command-line arguments.
        logger (Any): Logger object.

    Returns:
        Tuple[pd.DataFrame, np.ndarray]: Annotation matrix and M_tot array.

    Raises:
        ValueError: If there is an error reading the annotation matrix.
    """
    try:
        if args.ref_ld:
            annot_matrix, m_tot = ps.annot(split_paths(args.ref_ld), frqfile=args.frqfile)
        elif args.ref_ld_chr:
            annot_matrix, m_tot = ps.annot(split_paths(args.ref_ld_chr), frqfile=args.frqfile_chr)
        else:
            raise ValueError("No reference LD Scores provided for annotation.")
    except Exception as e:
        logger.log("Error parsing .annot file.")
        raise e

    return annot_matrix, m_tot


def read_m(args: Any, logger: Any, num_annotations: int) -> np.ndarray:
    """
    Read M (--M, --M-file, etc.) values.

    Args:
        args (argparse.Namespace): Command-line arguments.
        logger (Any): Logger object.
        num_annotations (int): Number of annotations.

    Returns:
        np.ndarray: M_annot array.

    Raises:
        ValueError: If M cannot be parsed or dimensions mismatch.
    """
    if args.M:
        try:
            m_annot_list = [float(x) for x in split_paths(args.M)]
        except ValueError as e:
            raise ValueError(f"Could not cast --M to float: {str(e)}")
    else:
        if args.ref_ld:
            m_annot_list = ps.M_fromlist(split_paths(args.ref_ld), common=(not args.not_M_5_50))
        elif args.ref_ld_chr:
            m_annot_list = ps.M_fromlist(split_paths(args.ref_ld_chr), common=(not args.not_M_5_50))
        else:
            raise ValueError("No reference LD Scores provided for M.")

    try:
        m_annot_array = np.array(m_annot_list).reshape((1, num_annotations))
    except ValueError as e:
        raise ValueError(f"# terms in --M must match # of LD Scores in --ref-ld.\n{str(e)}")

    return m_annot_array


def read_regression_weight_ld_scores(args: Any, logger: Any) -> pd.DataFrame:
    """
    Read regression weight LD Scores from files specified in args.

    Args:
        args (argparse.Namespace): Command-line arguments.
        logger (Any): Logger object.

    Returns:
        pd.DataFrame: Regression weight LD Scores DataFrame.

    Raises:
        ValueError: If there is an error reading the LD Scores.
    """
    if (args.w_ld and "," in args.w_ld) or (args.w_ld_chr and "," in args.w_ld_chr):
        raise ValueError("--w-ld must point to a single fileset (no commas allowed).")
    try:
        if args.w_ld:
            logger.log(f"Reading regression weight LD Scores from {args.w_ld} ...")
            w_ld = ps.ldscore_fromlist(split_paths(args.w_ld))
        elif args.w_ld_chr:
            pattern = ps.sub_chr(args.w_ld_chr, "[1-22]")
            logger.log(f"Reading regression weight LD Scores from {pattern} ...")
            w_ld = ps.ldscore_fromlist(split_paths(args.w_ld_chr))
        else:
            raise ValueError("No regression weight LD Scores provided.")
    except Exception as e:
        logger.log("Error parsing regression weight LD Scores.")
        raise e

    if len(w_ld.columns) != 2:
        raise ValueError("--w-ld may only have one LD Score column.")
    w_ld.columns = ["SNP", "LD_weights"]  # Prevent column name conflicts with ref_ld
    logger.log(f"Read regression weight LD Scores for {len(w_ld)} SNPs.")
    return w_ld


def read_summary_statistics(
    args: Any, logger: Any, filepath: str, alleles: bool = False, dropna: bool = False
) -> pd.DataFrame:
    """
    Parse summary statistics from the specified file.

    Args:
        args (argparse.Namespace): Command-line arguments.
        logger (Any): Logger object.
        filepath (str): Path to the summary statistics file.
        alleles (bool, optional): Whether to include allele columns. Defaults to False.
        dropna (bool, optional): Whether to drop rows with NA values. Defaults to False.

    Returns:
        pd.DataFrame: Summary statistics DataFrame.

    Raises:
        ValueError: If there is an error reading the summary statistics.
    """
    logger.log(f"Reading summary statistics from {filepath} ...")
    try:
        sumstats = ps.sumstats(filepath, alleles=alleles, dropna=dropna)
    except Exception as e:
        logger.log("Error parsing summary statistics.")
        raise e

    logger.log(f"Read summary statistics for {len(sumstats)} SNPs.")
    initial_len = len(sumstats)
    sumstats = sumstats.drop_duplicates(subset="SNP")
    if initial_len > len(sumstats):
        logger.log(f"Dropped {initial_len - len(sumstats)} SNPs with duplicated rs numbers.")

    return sumstats


def check_ld_condition_number(args: Any, logger: Any, ref_ld: pd.DataFrame) -> None:
    """
    Check the condition number of the LD Score matrix to ensure it is well-conditioned.

    Args:
        args (argparse.Namespace): Command-line arguments.
        logger (Any): Logger object.
        ref_ld (pd.DataFrame): Reference LD Scores DataFrame.

    Raises:
        ValueError: If the condition number is too high and inversion is not forced.
    """
    if ref_ld.shape[1] >= 2:
        condition_number = int(np.linalg.cond(ref_ld))
        if condition_number > 100000:
            if args.invert_anyway:
                warning_msg = (
                    f"WARNING: LD Score matrix condition number is {condition_number}. "
                    "Inverting anyway because the --invert-anyway flag is set."
                )
                logger.log(warning_msg)
            else:
                error_msg = (
                    f"ERROR: LD Score matrix condition number is {condition_number}. "
                    "Remove collinear LD Scores or use the --invert-anyway flag."
                )
                raise ValueError(error_msg)


def check_variance(
    logger: Any, m_annot: np.ndarray, ref_ld: pd.DataFrame
) -> Tuple[np.ndarray, pd.DataFrame, np.ndarray]:
    """
    Remove zero-variance LD Scores from the data.

    Args:
        logger (Any): Logger object.
        m_annot (np.ndarray): M_annot array.
        ref_ld (pd.DataFrame): Reference LD Scores DataFrame.

    Returns:
        Tuple[np.ndarray, pd.DataFrame, np.ndarray]: Updated M_annot, ref_ld, and boolean array of columns removed.

    Raises:
        ValueError: If all LD Scores have zero variance.
    """
    variance_zero = ref_ld.iloc[:, 1:].var() == 0  # Exclude 'SNP' column
    if variance_zero.all():
        raise ValueError("All LD Scores have zero variance.")
    else:
        logger.log("Removing partitioned LD Scores with zero variance.")
        columns_to_keep = np.array([True] + list(~variance_zero))  # Include 'SNP' column
        m_annot_columns = np.array(~variance_zero)
        ref_ld = ref_ld.iloc[:, columns_to_keep]
        m_annot = m_annot[:, m_annot_columns]

    return m_annot, ref_ld, variance_zero


def warn_if_few_snps(logger: Any, sumstats: pd.DataFrame) -> None:
    """
    Log a warning if the number of SNPs is less than 200,000.

    Args:
        logger (Any): Logger object.
        sumstats (pd.DataFrame): Summary statistics DataFrame.
    """
    if len(sumstats) < 200000:
        logger.log("WARNING: number of SNPs less than 200k; this is almost always bad.")


def print_covariance_matrix(ldscore_reg: Any, filepath: str, logger: Any) -> None:
    """
    Print the covariance matrix of the regression coefficients to a file.

    Args:
        ldscore_reg (Any): LD Score regression result object.
        filepath (str): Output file path.
        logger (Any): Logger object.
    """
    logger.log(f"Printing covariance matrix of the estimates to {filepath}.")
    np.savetxt(filepath, ldscore_reg.coef_cov)


def print_delete_values(ldscore_reg: Any, filepath: str, logger: Any) -> None:
    """
    Print block jackknife delete values to a file.

    Args:
        ldscore_reg (Any): LD Score regression result object.
        filepath (str): Output file path.
        logger (Any): Logger object.
    """
    logger.log(f"Printing block jackknife delete values to {filepath}.")
    np.savetxt(filepath, ldscore_reg.tot_delete_values)


def print_partitioned_delete_values(ldscore_reg: Any, filepath: str, logger: Any) -> None:
    """
    Print partitioned block jackknife delete values to a file.

    Args:
        ldscore_reg (Any): LD Score regression result object.
        filepath (str): Output file path.
        logger (Any): Logger object.
    """
    logger.log(f"Printing partitioned block jackknife delete values to {filepath}.")
    np.savetxt(filepath, ldscore_reg.part_delete_values)


def merge_and_log(ld_scores: pd.DataFrame, sumstats: pd.DataFrame, description: str, logger: Any) -> pd.DataFrame:
    """
    Merge LD Scores and summary statistics, and log the number of SNPs remaining.

    Args:
        ld_scores (pd.DataFrame): LD Scores DataFrame.
        sumstats (pd.DataFrame): Summary statistics DataFrame.
        description (str): Description of the LD Scores (e.g., "reference panel LD").
        logger (Any): Logger object.

    Returns:
        pd.DataFrame: Merged DataFrame.

    Raises:
        ValueError: If no SNPs remain after merging.
    """
    merged_data = smart_merge(ld_scores, sumstats)
    num_snps = len(merged_data)
    if num_snps == 0:
        raise ValueError(f"After merging with {description}, {num_snps} SNPs remain.")
    else:
        logger.log(f"After merging with {description}, {num_snps} SNPs remain.")

    return merged_data


def read_chr_split_files(
    chr_arg: Optional[str], not_chr_arg: Optional[str], logger: Any, noun: str, parsefunc: Callable, **kwargs
) -> Any:
    """
    Read files split across chromosomes (e.g., annot, ref_ld, w_ld).

    Args:
        chr_arg (Optional[str]): Comma-separated file paths with chromosome placeholders.
        not_chr_arg (Optional[str]): Comma-separated file paths without chromosome placeholders.
        logger (Any): Logger object.
        noun (str): Description of the data being read (e.g., "annot matrix").
        parsefunc (Callable): Function to parse the files.
        **kwargs: Additional keyword arguments to pass to parsefunc.

    Returns:
        Any: Parsed data from the files.

    Raises:
        ValueError: If there is an error parsing the files.
    """
    try:
        if not_chr_arg:
            logger.log(f"Reading {noun} from {not_chr_arg} ... ({parsefunc.__name__})")
            out = parsefunc(split_paths(not_chr_arg), **kwargs)
        elif chr_arg:
            pattern = ps.sub_chr(chr_arg, "[1-22]")
            logger.log(f"Reading {noun} from {pattern} ... ({parsefunc.__name__})")
            out = parsefunc(split_paths(chr_arg), **kwargs)
        else:
            raise ValueError(f"No files specified for {noun}.")
    except ValueError as e:
        logger.log(f"Error parsing {noun}.")
        raise e

    return out


def read_ld_and_sumstats(
    args: Any, logger: Any, filepath: str, alleles: bool = False, dropna: bool = True
) -> Tuple[np.ndarray, str, List[str], pd.DataFrame, np.ndarray]:
    """
    Read LD Scores and summary statistics, and prepare for regression.

    Args:
        args (argparse.Namespace): Command-line arguments.
        logger (Any): Logger object.
        filepath (str): Path to the summary statistics file.
        alleles (bool, optional): Whether to include allele columns. Defaults to False.
        dropna (bool, optional): Whether to drop rows with NA values. Defaults to True.

    Returns:
        Tuple[np.ndarray, str, List[str], pd.DataFrame, np.ndarray]: M_annot, w_ld_cname,
        ref_ld_cnames, sumstats, and novar_cols.

    Raises:
        ValueError: If there is an error in data preparation.
    """
    sumstats = read_summary_statistics(args, logger, filepath, alleles=alleles, dropna=dropna)
    ref_ld = read_reference_ld_scores(args, logger)
    num_annotations = len(ref_ld.columns) - 1  # Exclude 'SNP' column
    m_annot = read_m(args, logger, num_annotations)
    m_annot, ref_ld, novar_cols = check_variance(logger, m_annot, ref_ld)
    w_ld = read_regression_weight_ld_scores(args, logger)
    sumstats = merge_and_log(ref_ld, sumstats, "reference panel LD", logger)
    sumstats = merge_and_log(sumstats, w_ld, "regression SNP LD", logger)
    w_ld_cname = sumstats.columns[-1]
    ref_ld_cnames = ref_ld.columns[1:].tolist()
    return m_annot, w_ld_cname, ref_ld_cnames, sumstats, novar_cols


# [The rest of the module remains the same as previously provided]


def estimate_cell_type_specific_heritability(args: Any, logger: Any) -> None:
    """
    Perform cell-type-specific heritability analysis.

    Args:
        args (argparse.Namespace): Command-line arguments.
        logger (Any): Logger object.
    """
    args = copy.deepcopy(args)
    if args.intercept_h2 is not None:
        args.intercept_h2 = float(args.intercept_h2)
    if args.no_intercept:
        args.intercept_h2 = 1

    m_annot_all_regr, w_ld_cname, ref_ld_cnames_all_regr, sumstats, novar_cols = read_ld_and_sumstats(
        args, logger, args.h2_cts
    )
    m_tot = np.sum(m_annot_all_regr)
    check_ld_condition_number(args, logger, sumstats[ref_ld_cnames_all_regr])
    warn_if_few_snps(logger, sumstats)
    num_snps = len(sumstats)
    num_blocks = min(num_snps, args.n_blocks)
    if args.chisq_max is None:
        chisq_max = max(0.001 * sumstats.N.max(), 80)
    else:
        chisq_max = args.chisq_max

    valid_indices = (sumstats.Z**2 < chisq_max).values
    sumstats = sumstats.loc[valid_indices, :]
    logger.log(
        f"Removed {num_snps - valid_indices.sum()} SNPs with chi^2 > {chisq_max} "
        f"({valid_indices.sum()} SNPs remain)"
    )
    num_snps = valid_indices.sum()
    ref_ld_all_regr = sumstats[ref_ld_cnames_all_regr].values.reshape((num_snps, -1))
    chisq = sumstats.Z**2
    keep_snps = sumstats[["SNP"]]

    def reshape_array(x: pd.Series) -> np.ndarray:
        return x.values.reshape((num_snps, 1))

    results_columns = [
        "Name",
        "Coefficient",
        "Coefficient_std_error",
        "Coefficient_P_value",
    ]
    results_data = []
    with open(args.ref_ld_chr_cts) as f:
        for line in f:
            name, ct_ld_chr = line.strip().split()
            ref_ld_cts_allsnps = ps.ldscore_fromlist(split_paths(ct_ld_chr), n_chr=NUM_CHROMOSOMES)
            logger.log("Performing regression.")
            ref_ld_cts = pd.merge(keep_snps, ref_ld_cts_allsnps, on="SNP", how="left").iloc[:, 1:].values
            if np.any(np.isnan(ref_ld_cts)):
                raise ValueError(
                    "Missing some LD scores from cts files. "
                    "Are you sure all SNPs in ref-ld-chr are also in ref-ld-chr-cts?"
                )

            ref_ld = np.hstack([ref_ld_cts, ref_ld_all_regr])
            m_cts = ps.M_fromlist(split_paths(ct_ld_chr), n_chr=NUM_CHROMOSOMES, common=(not args.not_M_5_50))
            m_annot = np.hstack([m_cts, m_annot_all_regr])
            hsqhat = reg.Hsq(
                chisq.values.reshape((num_snps, 1)),
                ref_ld,
                reshape_array(sumstats[w_ld_cname]),
                reshape_array(sumstats.N),
                m_annot,
                n_blocks=num_blocks,
                intercept=args.intercept_h2,
                twostep=None,
                old_weights=True,
            )
            coef, coef_se = hsqhat.coef[0], hsqhat.coef_se[0]
            p_value = stats.norm.sf(abs(coef / coef_se)) * 2  # Two-tailed p-value
            results_data.append((name, coef, coef_se, p_value))
            if args.print_all_cts:
                for idx in range(1, len(ct_ld_chr.split(","))):
                    coef, coef_se = hsqhat.coef[idx], hsqhat.coef_se[idx]
                    p_value = stats.norm.sf(abs(coef / coef_se)) * 2
                    results_data.append((f"{name}_{idx}", coef, coef_se, p_value))

    df_results = pd.DataFrame(data=results_data, columns=results_columns)
    df_results.sort_values(by="Coefficient_P_value", inplace=True)
    output_filepath = f"{args.out}.cell_type_results.txt"
    df_results.to_csv(output_filepath, sep="\t", index=False)
    logger.log(f"Results printed to {output_filepath}")


def estimate_heritability(args: Any, logger: Any) -> Any:
    """
    Estimate heritability and partitioned heritability.

    Args:
        args (argparse.Namespace): Command-line arguments.
        logger (Any): Logger object.

    Returns:
        Any: Heritability estimation result object.
    """
    args = copy.deepcopy(args)
    if args.samp_prev is not None and args.pop_prev is not None:
        args.samp_prev, args.pop_prev = list(map(float, [args.samp_prev, args.pop_prev]))
    if args.intercept_h2 is not None:
        args.intercept_h2 = float(args.intercept_h2)
    if args.no_intercept:
        args.intercept_h2 = 1
    m_annot, w_ld_cname, ref_ld_cnames, sumstats, novar_cols = read_ld_and_sumstats(args, logger, args.h2)
    ref_ld = sumstats[ref_ld_cnames].values
    check_ld_condition_number(args, logger, sumstats[ref_ld_cnames])
    warn_if_few_snps(logger, sumstats)
    num_snps = len(sumstats)
    num_blocks = min(num_snps, args.n_blocks)
    num_annotations = len(ref_ld_cnames)
    chisq_max = args.chisq_max
    old_weights = False
    if num_annotations == 1:
        if args.two_step is None and args.intercept_h2 is None:
            args.two_step = 30
    else:
        old_weights = True
        if args.chisq_max is None:
            chisq_max = max(0.001 * sumstats.N.max(), 80)

    def reshape_array(x: pd.Series) -> np.ndarray:
        return x.values.reshape((num_snps, 1))

    chisq = sumstats.Z**2
    if chisq_max is not None:
        valid_indices = (chisq < chisq_max).values
        sumstats = sumstats.iloc[valid_indices, :]
        logger.log(
            f"Removed {num_snps - valid_indices.sum()} SNPs with chi^2 > {chisq_max} "
            f"({valid_indices.sum()} SNPs remain)"
        )
        num_snps = valid_indices.sum()
        ref_ld = sumstats[ref_ld_cnames].values
        chisq = chisq[valid_indices]

    if args.two_step is not None:
        logger.log(f"Using two-step estimator with cutoff at {args.two_step}.")

    hsqhat = reg.Hsq(
        chisq.values.reshape((num_snps, 1)),
        ref_ld,
        reshape_array(sumstats[w_ld_cname]),
        reshape_array(sumstats.N),
        m_annot,
        n_blocks=num_blocks,
        intercept=args.intercept_h2,
        twostep=args.two_step,
        old_weights=old_weights,
    )

    if args.print_cov:
        print_covariance_matrix(hsqhat, args.out + ".cov", logger)
    if args.print_delete_vals:
        print_delete_values(hsqhat, args.out + ".delete", logger)
        print_partitioned_delete_values(hsqhat, args.out + ".part_delete", logger)

    logger.log(hsqhat.summary(ref_ld_cnames, P=args.samp_prev, K=args.pop_prev, overlap=args.overlap_annot))
    if args.overlap_annot:
        annot_matrix, m_tot = read_annotation_matrix(args, logger)
        df_results = hsqhat.overlap_output(ref_ld_cnames, annot_matrix, m_annot, m_tot, args.print_coefficients)
        df_results.to_csv(args.out + ".results", sep="\t", index=False)
        logger.log(f"Results printed to {args.out}.results")

    return hsqhat


def estimate_genetic_correlation(args: Any, logger: Any) -> List[Any]:
    """
    Estimate genetic correlation (rg) between trait 1 and a list of other traits.

    Args:
        args (argparse.Namespace): Command-line arguments.
        logger (Any): Logger object.

    Returns:
        List[Any]: List of genetic correlation estimation results.
    """
    args = copy.deepcopy(args)
    rg_paths, rg_files = parse_rg(args.rg)
    num_phenotypes = len(rg_paths)
    args.intercept_h2, args.intercept_gencov, args.samp_prev, args.pop_prev = map(
        lambda x: split_or_none(x, num_phenotypes),
        (args.intercept_h2, args.intercept_gencov, args.samp_prev, args.pop_prev),
    )
    for arg_values, arg_name in [
        (args.intercept_h2, "--intercept-h2"),
        (args.intercept_gencov, "--intercept-gencov"),
        (args.samp_prev, "--samp-prev"),
        (args.pop_prev, "--pop-prev"),
    ]:
        check_arg_length(arg_values, num_phenotypes, arg_name)

    if args.no_intercept:
        args.intercept_h2 = [1] * num_phenotypes
        args.intercept_gencov = [0] * num_phenotypes
    p1 = rg_paths[0]
    out_prefix = args.out + rg_files[0]
    m_annot, w_ld_cname, ref_ld_cnames, sumstats, _ = read_ld_and_sumstats(args, logger, p1, alleles=True, dropna=True)
    rg_results = []
    num_annotations = m_annot.shape[1]
    if num_annotations == 1 and args.two_step is None and args.intercept_h2[0] is None:
        args.two_step = 30
    if args.two_step is not None:
        logger.log(f"Using two-step estimator with cutoff at {args.two_step}.")

    for i, p2 in enumerate(rg_paths[1:]):
        logger.log(f"Computing rg for phenotype {i + 2}/{num_phenotypes}")
        try:
            loop_data = read_other_sumstats(args, logger, p2, sumstats, ref_ld_cnames)
            rghat = compute_rg(loop_data, args, logger, m_annot, ref_ld_cnames, w_ld_cname, i)
            rg_results.append(rghat)
            print_genetic_correlation(args, logger, rghat, ref_ld_cnames, i, rg_paths, i == 0)
            out_prefix_loop = out_prefix + "_" + rg_files[i + 1]
            if args.print_cov:
                print_rg_covariance(rghat, out_prefix_loop, logger)
            if args.print_delete_vals:
                print_rg_delete_values(rghat, out_prefix_loop, logger)
        except Exception as e:
            error_msg = f"ERROR computing rg for phenotype {i + 2}/{num_phenotypes}, from file {rg_paths[i + 1]}."
            logger.log(error_msg)
            logger.log(f"Exception: {e}")
            logger.log(f"Traceback: {traceback.format_exc()}")
            if len(rg_results) <= i:
                rg_results.append(None)

    logger.log("\nSummary of Genetic Correlation Results\n" + get_rg_table(rg_paths, rg_results, args))
    return rg_results


def read_other_sumstats(
    args: Any, logger: Any, filepath: str, sumstats: pd.DataFrame, ref_ld_cnames: List[str]
) -> pd.DataFrame:
    """
    Read and merge summary statistics for another phenotype.

    Args:
        args (argparse.Namespace): Command-line arguments.
        logger (Any): Logger object.
        filepath (str): Path to the summary statistics file for the other phenotype.
        sumstats (pd.DataFrame): Summary statistics DataFrame for the first phenotype.
        ref_ld_cnames (List[str]): List of reference LD Score column names.

    Returns:
        pd.DataFrame: Merged DataFrame with summary statistics for both phenotypes.
    """
    other_sumstats = read_summary_statistics(args, logger, filepath, alleles=True, dropna=False)
    merged_sumstats = merge_sumstats(sumstats, other_sumstats, logger)
    merged_sumstats = merged_sumstats.dropna(how="any")
    alleles = merged_sumstats.A1 + merged_sumstats.A2 + merged_sumstats.A1x + merged_sumstats.A2x
    if not args.no_check_alleles:
        valid_indices = filter_alleles(alleles)
        merged_sumstats = select_and_log(merged_sumstats, valid_indices, logger, "{N} SNPs with valid alleles.")
        merged_sumstats["Z2"] = align_alleles(merged_sumstats.Z2, alleles)

    merged_sumstats = merged_sumstats.drop(["A1", "A1x", "A2", "A2x"], axis=1)
    check_ld_condition_number(args, logger, merged_sumstats[ref_ld_cnames])
    warn_if_few_snps(logger, merged_sumstats)
    return merged_sumstats


def get_rg_table(rg_paths: List[str], rg_results: List[Any], args: Any) -> str:
    """
    Generate a table of genetic correlations.

    Args:
        rg_paths (List[str]): List of summary statistics file paths.
        rg_results (List[Any]): List of genetic correlation results.
        args (argparse.Namespace): Command-line arguments.

    Returns:
        str: Formatted table as a string.
    """

    def get_attribute(obj: Any, attr: str) -> Any:
        return getattr(obj, attr, "NA")

    data = {
        "p1": [rg_paths[0]] * (len(rg_paths) - 1),
        "p2": rg_paths[1:],
        "rg": [get_attribute(rg, "rg_ratio") for rg in rg_results],
        "se": [get_attribute(rg, "rg_se") for rg in rg_results],
        "z": [get_attribute(rg, "z") for rg in rg_results],
        "p": [get_attribute(rg, "p") for rg in rg_results],
    }

    if (
        args.samp_prev is not None
        and args.pop_prev is not None
        and all(i is not None for i in args.samp_prev)
        and all(i is not None for i in args.pop_prev)
    ):
        c = [reg.h2_obs_to_liab(1, sp, pp) for sp, pp in zip(args.samp_prev[1:], args.pop_prev[1:])]
        h2_liab = [c_i * get_attribute(get_attribute(rg, "hsq2"), "tot") for c_i, rg in zip(c, rg_results)]
        h2_liab_se = [c_i * get_attribute(get_attribute(rg, "hsq2"), "tot_se") for c_i, rg in zip(c, rg_results)]
        data["h2_liab"] = h2_liab
        data["h2_liab_se"] = h2_liab_se
    else:
        data["h2_obs"] = [get_attribute(get_attribute(rg, "hsq2"), "tot") for rg in rg_results]
        data["h2_obs_se"] = [get_attribute(get_attribute(rg, "hsq2"), "tot_se") for rg in rg_results]

    data["h2_int"] = [get_attribute(get_attribute(rg, "hsq2"), "intercept") for rg in rg_results]
    data["h2_int_se"] = [get_attribute(get_attribute(rg, "hsq2"), "intercept_se") for rg in rg_results]
    data["gcov_int"] = [get_attribute(get_attribute(rg, "gencov"), "intercept") for rg in rg_results]
    data["gcov_int_se"] = [get_attribute(get_attribute(rg, "gencov"), "intercept_se") for rg in rg_results]

    df = pd.DataFrame(data)
    return df.to_string(header=True, index=False) + "\n"


def print_genetic_correlation(
    args: Any, logger: Any, rghat: Any, ref_ld_cnames: List[str], index: int, rg_paths: List[str], print_hsq1: bool
) -> None:
    """
    Print genetic correlation results.

    Args:
        args (argparse.Namespace): Command-line arguments.
        logger (Any): Logger object.
        rghat (Any): Genetic correlation estimation result.
        ref_ld_cnames (List[str]): List of reference LD Score column names.
        index (int): Index of the phenotype.
        rg_paths (List[str]): List of summary statistics file paths.
        print_hsq1 (bool): Whether to print heritability of the first phenotype.
    """

    def header_line(title: str) -> str:
        return title + "\n" + "-" * len(title)

    P = [args.samp_prev[0], args.samp_prev[index + 1]]
    K = [args.pop_prev[0], args.pop_prev[index + 1]]
    if args.samp_prev is None and args.pop_prev is None:
        args.samp_prev = [None, None]
        args.pop_prev = [None, None]
    if print_hsq1:
        logger.log(header_line("\nHeritability of phenotype 1"))
        logger.log(rghat.hsq1.summary(ref_ld_cnames, P=P[0], K=K[0]))

    logger.log(header_line(f"\nHeritability of phenotype {index + 2}/{len(rg_paths)}"))
    logger.log(rghat.hsq2.summary(ref_ld_cnames, P=P[1], K=K[1]))
    logger.log(header_line("\nGenetic Covariance"))
    logger.log(rghat.gencov.summary(ref_ld_cnames, P=P, K=K))
    logger.log(header_line("\nGenetic Correlation"))
    logger.log(rghat.summary() + "\n")


def merge_sumstats(sumstats1: pd.DataFrame, sumstats2: pd.DataFrame, logger: Any) -> pd.DataFrame:
    """
    Merge two sets of summary statistics.

    Args:
        sumstats1 (pd.DataFrame): Summary statistics DataFrame for the first phenotype.
        sumstats2 (pd.DataFrame): Summary statistics DataFrame for the second phenotype.
        logger (Any): Logger object.

    Returns:
        pd.DataFrame: Merged DataFrame.
    """
    sumstats1 = sumstats1.rename(columns={"N": "N1", "Z": "Z1"})
    sumstats2 = sumstats2.rename(columns={"A1": "A1x", "A2": "A2x", "N": "N2", "Z": "Z2"})
    merged_sumstats = merge_and_log(sumstats1, sumstats2, "summary statistics", logger)
    return merged_sumstats


def filter_alleles(alleles: pd.Series) -> pd.Series:
    """
    Filter out SNPs with invalid alleles (mismatched alleles, non-SNPs, strand ambiguous).

    Args:
        alleles (pd.Series): Series of concatenated allele strings.

    Returns:
        pd.Series: Boolean Series indicating valid SNPs.
    """
    return alleles.apply(lambda x: x in MATCH_ALLELES)


def align_alleles(z_scores: pd.Series, alleles: pd.Series) -> pd.Series:
    """
    Align Z-scores to the same choice of reference allele (allowing for strand flip).

    Args:
        z_scores (pd.Series): Series of Z-scores to align.
        alleles (pd.Series): Series of concatenated allele strings.

    Returns:
        pd.Series: Aligned Z-scores.

    Raises:
        KeyError: If alleles are incompatible between summary statistics files.
    """
    try:
        alignment_factors = alleles.apply(lambda x: (-1) ** FLIP_ALLELES[x])
        aligned_z_scores = z_scores * alignment_factors
    except KeyError as e:
        msg = f"Incompatible alleles in .sumstats files: {e.args}. "
        msg += "Did you forget to use --merge-alleles with munge_sumstats.py?"
        raise KeyError(msg)
    return aligned_z_scores


def compute_rg(
    sumstats: pd.DataFrame,
    args: Any,
    logger: Any,
    m_annot: np.ndarray,
    ref_ld_cnames: List[str],
    w_ld_cname: str,
    index: int,
) -> Any:
    """
    Compute genetic correlation.

    Args:
        sumstats (pd.DataFrame): Summary statistics DataFrame.
        args (argparse.Namespace): Command-line arguments.
        logger (Any): Logger object.
        m_annot (np.ndarray): M_annot array.
        ref_ld_cnames (List[str]): List of reference LD Score column names.
        w_ld_cname (str): Name of the weight LD Score column.
        index (int): Index of the phenotype.

    Returns:
        Any: Genetic correlation estimation result.
    """
    num_snps = len(sumstats)
    chisq_max = args.chisq_max
    if chisq_max is not None:
        valid_indices = (sumstats.Z1**2 * sumstats.Z2**2 < chisq_max**2).values
        num_snps = valid_indices.sum()
        sumstats = sumstats[valid_indices]

    num_blocks = min(args.n_blocks, num_snps)
    ref_ld = sumstats[ref_ld_cnames].values
    intercepts = [
        args.intercept_h2[0],
        args.intercept_h2[index + 1],
        args.intercept_gencov[index + 1],
    ]

    def reshape_array(x: pd.Series) -> np.ndarray:
        return x.values.reshape((num_snps, 1))

    rghat = reg.RG(
        reshape_array(sumstats.Z1),
        reshape_array(sumstats.Z2),
        ref_ld,
        reshape_array(sumstats[w_ld_cname]),
        reshape_array(sumstats.N1),
        reshape_array(sumstats.N2),
        m_annot,
        intercept_hsq1=intercepts[0],
        intercept_hsq2=intercepts[1],
        intercept_gencov=intercepts[2],
        n_blocks=num_blocks,
        twostep=args.two_step,
    )

    return rghat


def parse_rg(rg: str) -> Tuple[List[str], List[str]]:
    """
    Parse the --rg argument into file paths and file names.

    Args:
        rg (str): Comma-separated string of summary statistics file paths.

    Returns:
        Tuple[List[str], List[str]]: List of file paths and list of file names.

    Raises:
        ValueError: If fewer than two phenotypes are provided.
    """
    rg_paths = split_paths(rg)
    rg_files = [os.path.basename(path) for path in rg_paths]
    if len(rg_paths) < 2:
        raise ValueError("Must specify at least two phenotypes for rg estimation.")

    return rg_paths, rg_files


def print_rg_delete_values(rg_result: Any, filepath: str, logger: Any) -> None:
    """
    Print block jackknife delete values for genetic correlation estimation.

    Args:
        rg_result (Any): Genetic correlation estimation result.
        filepath (str): Output file path prefix.
        logger (Any): Logger object.
    """
    print_delete_values(rg_result.hsq1, filepath + ".hsq1.delete", logger)
    print_delete_values(rg_result.hsq2, filepath + ".hsq2.delete", logger)
    print_delete_values(rg_result.gencov, filepath + ".gencov.delete", logger)


def print_rg_covariance(rg_result: Any, filepath: str, logger: Any) -> None:
    """
    Print covariance matrices for genetic correlation estimation.

    Args:
        rg_result (Any): Genetic correlation estimation result.
        filepath (str): Output file path prefix.
        logger (Any): Logger object.
    """
    print_covariance_matrix(rg_result.hsq1, filepath + ".hsq1.cov", logger)
    print_covariance_matrix(rg_result.hsq2, filepath + ".hsq2.cov", logger)
    print_covariance_matrix(rg_result.gencov, filepath + ".gencov.cov", logger)


def split_or_none(value: Optional[str], n: int) -> List[Optional[float]]:
    """
    Split a comma-separated string into a list, or return a list of None.

    Args:
        value (Optional[str]): Comma-separated string or None.
        n (int): Length of the list to return.

    Returns:
        List[Optional[float]]: List of values or None.
    """
    if value is not None:
        value_list = [float(x) if x != "N" else None for x in value.split(",")]
    else:
        value_list = [None] * n
    return value_list


def check_arg_length(arg_list: List[Any], expected_length: int, arg_name: str) -> None:
    """
    Check that a list has the expected length.

    Args:
        arg_list (List[Any]): List of arguments.
        expected_length (int): Expected length of the list.
        arg_name (str): Name of the argument for error messages.

    Raises:
        ValueError: If the length of the list does not match the expected length.
    """
    if len(arg_list) != expected_length:
        raise ValueError(f"{arg_name} must have the same number of arguments as --rg/--h2.")
