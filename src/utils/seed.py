"""
Seed Management for Reproducibility

This module provides unified seed setting across all random number generators
used in the project: Python's random, NumPy, and PyTorch.
"""

import random
import numpy as np
import torch
from typing import Optional


def set_seed(seed: int = 42, deterministic: bool = True) -> None:
    """
    Set random seed for reproducibility across all libraries.

    Args:
        seed: Random seed value
        deterministic: If True, use deterministic algorithms in PyTorch
                      (may impact performance)
    """
    # Python's built-in random
    random.seed(seed)

    # NumPy
    np.random.seed(seed)

    # PyTorch CPU
    torch.manual_seed(seed)

    # PyTorch GPU (all GPUs)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)

        # Deterministic algorithms
        if deterministic:
            torch.backends.cudnn.deterministic = True
            torch.backends.cudnn.benchmark = False


def get_device(prefer_cuda: bool = True) -> torch.device:
    """
    Get the appropriate device for computation.

    Args:
        prefer_cuda: If True, use CUDA when available

    Returns:
        torch.device: The device to use for computation
    """
    if prefer_cuda and torch.cuda.is_available():
        return torch.device('cuda')
    return torch.device('cpu')


def seed_worker(worker_id: int) -> None:
    """
    Worker initialization function for DataLoader.

    Use this with DataLoader's worker_init_fn to ensure reproducibility
    when using multiple workers.

    Args:
        worker_id: Worker ID assigned by DataLoader

    Example:
        >>> from torch.utils.data import DataLoader
        >>> loader = DataLoader(dataset, num_workers=4, worker_init_fn=seed_worker)
    """
    worker_seed = torch.initial_seed() % 2**32
    np.random.seed(worker_seed)
    random.seed(worker_seed)


class SeedContext:
    """
    Context manager for temporary seed setting.

    Useful for ensuring reproducibility in specific code blocks
    while not affecting the global state.

    Example:
        >>> with SeedContext(42):
        ...     random_value = torch.randn(3)
    """

    def __init__(self, seed: int):
        self.seed = seed
        self._rng_state = None
        self._np_state = None
        self._torch_state = None
        self._cuda_state = None

    def __enter__(self):
        # Save current states
        self._rng_state = random.getstate()
        self._np_state = np.random.get_state()
        self._torch_state = torch.get_rng_state()
        if torch.cuda.is_available():
            self._cuda_state = torch.cuda.get_rng_state_all()

        # Set new seed
        set_seed(self.seed, deterministic=False)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        # Restore states
        random.setstate(self._rng_state)
        np.random.set_state(self._np_state)
        torch.set_rng_state(self._torch_state)
        if torch.cuda.is_available() and self._cuda_state is not None:
            torch.cuda.set_rng_state_all(self._cuda_state)
        return False
