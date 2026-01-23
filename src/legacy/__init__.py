"""
Legacy Code for Backward Compatibility

This package contains original implementations preserved for:
- Reproducing published results
- Comparing with new standardized implementations
- Historical reference

WARNING: These modules are not maintained and may have inconsistencies
with the current standardized pipeline. Use src/train.py and
src/data_pipeline.py for new experiments.
"""

from .legacy_bipartite_builder import build_bipartite_graph as legacy_build_bipartite
from .legacy_trainer import train as legacy_train

__all__ = ['legacy_build_bipartite', 'legacy_train']
