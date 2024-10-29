"""
Iteratively Re-weighted Least Squares (IRWLS) module.

This module provides the IRWLS class, which implements iteratively re-weighted
least squares for regression analysis, including methods for jackknife variance
estimation.

(c) 2015 Brendan Bulik-Sullivan and Hilary Finucane
(c) 2024 Thomas Reimonn
"""

from typing import Callable, Optional, Union

import numpy as np

from . import jackknife as jk


class IRWLS:
    """
    Iteratively Re-weighted Least Squares (IRWLS) estimator.

    This class implements the IRWLS algorithm for estimating regression coefficients,
    allowing for heteroscedasticity or other forms of non-constant variance in the
    residuals. It also provides jackknife variance estimation using block jackknife.

    Attributes:
        est (np.ndarray): Estimated regression coefficients (shape: (n_features, 1)).
        jknife_est (np.ndarray): Jackknife estimates of the regression coefficients
            (shape: (n_features, 1)).
        jknife_var (np.ndarray): Variance of the jackknife estimates (shape: (n_features, 1)).
        jknife_se (np.ndarray): Standard errors of the jackknife estimates
            (shape: (n_features, 1)).
        jknife_cov (np.ndarray): Covariance matrix of the jackknife estimates
            (shape: (n_features, n_features)).
        delete_values (np.ndarray): Jackknife delete values (shape: (n_blocks, n_features)).
        separators (Optional[np.ndarray]): Block boundaries for jackknife
            (shape: (n_blocks + 1,)).

    """

    def __init__(
        self,
        X: np.ndarray,
        y: np.ndarray,
        update_func: Callable[[tuple], np.ndarray],
        n_blocks: int,
        w: Optional[np.ndarray] = None,
        slow: bool = False,
        separators: Optional[np.ndarray] = None,
        max_iter: int = 2,
    ) -> None:
        """
        Initialize the IRWLS estimator.

        Args:
            X (np.ndarray): Independent variables (shape: (n_samples, n_features)).
            y (np.ndarray): Dependent variable (shape: (n_samples,) or (n_samples, 1)).
            update_func (Callable[[tuple], np.ndarray]): Function to update weights.
                Should take the output of np.linalg.lstsq and return new weights
                (shape: (n_samples, 1)).
            n_blocks (int): Number of jackknife blocks for variance estimation.
            w (Optional[np.ndarray]): Initial regression weights (shape: (n_samples,) or
                (n_samples, 1)). Defaults to ones if None.
            slow (bool): Whether to use the slow block jackknife method (for testing).
            separators (Optional[np.ndarray]): Optional block boundaries for jackknife.
            max_iter (int): Maximum number of iterations for the IRWLS algorithm.

        Raises:
            ValueError: If input arrays have incompatible shapes.
        """
        n_samples, _ = X.shape
        y = y.reshape(-1, 1)
        if w is None:
            w = np.ones((n_samples, 1))
        else:
            w = w.reshape(-1, 1)

        if w.shape != (n_samples, 1):
            raise ValueError(f"w has shape {w.shape}. Expected shape: ({n_samples}, 1).")

        jknife = self.irwls(
            X,
            y,
            update_func,
            n_blocks,
            w,
            slow=slow,
            separators=separators,
            max_iter=max_iter,
        )
        self.est = jknife.est
        self.jknife_est = jknife.jknife_est
        self.jknife_var = jknife.jknife_var
        self.jknife_se = jknife.jknife_se
        self.jknife_cov = jknife.jknife_cov
        self.delete_values = jknife.delete_values
        self.separators = jknife.separators

    @classmethod
    def irwls(
        cls,
        X: np.ndarray,
        y: np.ndarray,
        update_func: Callable[[tuple], np.ndarray],
        n_blocks: int,
        w: np.ndarray,
        slow: bool = False,
        separators: Optional[np.ndarray] = None,
        max_iter: int = 2,
    ) -> Union[jk.LstsqJackknifeFast, jk.LstsqJackknifeSlow]:
        """
        Perform Iteratively Re-weighted Least Squares (IRWLS).

        Args:
            X (np.ndarray): Independent variables (shape: (n_samples, n_features)).
            y (np.ndarray): Dependent variable (shape: (n_samples, 1)).
            update_func (Callable[[tuple], np.ndarray]): Function to update weights.
            n_blocks (int): Number of jackknife blocks.
            w (np.ndarray): Initial regression weights (shape: (n_samples, 1)).
            slow (bool): Whether to use the slow block jackknife method.
            separators (Optional[np.ndarray]): Optional block boundaries.
            max_iter (int): Maximum number of iterations for the IRWLS algorithm.

        Returns:
            Union[jk.LstsqJackknifeFast, jk.LstsqJackknifeSlow]: Jackknife regression object
            with final IRWLS weights.

        Raises:
            ValueError: If input arrays have incompatible shapes or weights are invalid.
        """
        n_samples, _ = X.shape
        y = y.reshape(-1, 1)
        w = w.reshape(-1, 1)

        if y.shape != (n_samples, 1):
            raise ValueError(f"y has shape {y.shape}. Expected shape: ({n_samples}, 1).")
        if w.shape != (n_samples, 1):
            raise ValueError(f"w has shape {w.shape}. Expected shape: ({n_samples}, 1).")

        # Initialize weights
        w_sqrt = np.sqrt(w)

        # Iteratively update weights
        for iteration in range(max_iter):
            coef = cls.wls(X, y, w_sqrt)
            new_w = np.sqrt(update_func(coef))
            if new_w.shape != w_sqrt.shape:
                raise ValueError(f"New weights have shape {new_w.shape}, expected {w_sqrt.shape}.")
            w_sqrt = new_w

        # Weight the data
        X_weighted = cls._weight(X, w_sqrt)
        y_weighted = cls._weight(y, w_sqrt)

        # Perform jackknife estimation
        if slow:
            jknife = jk.LstsqJackknifeSlow(X_weighted, y_weighted, n_blocks, separators=separators)
        else:
            jknife = jk.LstsqJackknifeFast(X_weighted, y_weighted, n_blocks, separators=separators)

        return jknife

    @classmethod
    def wls(
        cls,
        X: np.ndarray,
        y: np.ndarray,
        w_sqrt: np.ndarray,
    ) -> tuple:
        """
        Perform Weighted Least Squares regression.

        Args:
            X (np.ndarray): Independent variables (shape: (n_samples, n_features)).
            y (np.ndarray): Dependent variable (shape: (n_samples, 1)).
            w_sqrt (np.ndarray): Square root of weights (shape: (n_samples, 1)).

        Returns:
            tuple: Output of np.linalg.lstsq (coefficients, residuals, rank, singular values).

        Raises:
            ValueError: If input arrays have incompatible shapes.
        """
        n_samples, _ = X.shape
        y = y.reshape(-1, 1)
        w_sqrt = w_sqrt.reshape(-1, 1)

        if y.shape != (n_samples, 1):
            raise ValueError(f"y has shape {y.shape}. Expected shape: ({n_samples}, 1).")
        if w_sqrt.shape != (n_samples, 1):
            raise ValueError(f"w_sqrt has shape {w_sqrt.shape}. Expected shape: ({n_samples}, 1).")

        # Weight the data
        X_weighted = cls._weight(X, w_sqrt)
        y_weighted = cls._weight(y, w_sqrt)

        # Perform least squares regression
        coef = np.linalg.lstsq(X_weighted, y_weighted, rcond=None)

        return coef

    @staticmethod
    def _weight(
        X: np.ndarray,
        w_sqrt: np.ndarray,
    ) -> np.ndarray:
        """
        Weight the data matrix X by w_sqrt.

        Args:
            X (np.ndarray): Data matrix (shape: (n_samples, n_features) or (n_samples, 1)).
            w_sqrt (np.ndarray): Square root of weights (shape: (n_samples, 1)).

        Returns:
            np.ndarray: Weighted data matrix (shape: (n_samples, n_features) or (n_samples, 1)).

        Raises:
            ValueError: If weights contain non-positive values or shapes are incompatible.
        """
        if np.any(w_sqrt <= 0):
            raise ValueError("Weights must be positive.")

        n_samples = X.shape[0]
        if w_sqrt.shape != (n_samples, 1):
            raise ValueError(f"w_sqrt has shape {w_sqrt.shape}. Expected shape: ({n_samples}, 1).")

        # Normalize weights to have sum 1
        w_normalized = w_sqrt / np.sum(w_sqrt)

        # Multiply each row of X by the corresponding weight
        X_weighted = X * w_normalized

        return X_weighted
