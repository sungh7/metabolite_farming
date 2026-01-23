"""
Centralized Configuration for Ethylene-Isoflavonoid GNN Project

This module provides a single source of truth for all configuration parameters
used across the data pipeline, training, and evaluation scripts.
"""

from pathlib import Path

# =============================================================================
# DATA PATHS
# =============================================================================
DATA_PATHS = {
    # Raw data
    'raw_dir': Path('data/raw'),
    'kegg_dir': Path('data/kegg'),

    # Processed data
    'processed_dir': Path('data/processed'),
    'base_graph': Path('data/processed/strict_graph.pt'),
    'enhanced_graph': Path('data/processed/enhanced_bipartite_graph.pt'),
    'bipartite_graph_v2': Path('data/processed/strict_bipartite_v2.pt'),

    # Mappings
    'enzyme_mapping': Path('data/processed/enzyme_string_mapping.csv'),
    'tf_mapping': Path('data/processed/tf_string_mapping.csv'),

    # KEGG specific
    'kegg_gene_ec': Path('data/kegg/gene_ec_mapping.tsv'),
    'kegg_uniprot': Path('data/kegg/kegg_uniprot_mapping.csv'),
    'kegg_metabolites': Path('data/kegg/metabolites.csv'),
    'kegg_edges': Path('data/kegg/full_enzyme_metabolite_edges.tsv'),

    # Models
    'model_dir': Path('data/models'),

    # Results
    'results_dir': Path('results'),
    'gnn_results': Path('results/gnn'),

    # V3 specific
    'graph_v3': Path('data/processed/graph_v3.pt'),
}

# =============================================================================
# GRAPH BUILDING CONFIGURATION
# =============================================================================
GRAPH_CONFIG = {
    # STRING-DB threshold for PPI edges
    'ppi_threshold': 700,

    # Use strict mode (exclude text-mining evidence)
    'strict_mode': True,

    # Node feature dimensionality
    'feature_dim': 64,

    # Currency metabolites to exclude from graph
    'currency_metabolites': {
        'C00001', 'C00002', 'C00003', 'C00004', 'C00005', 'C00006',
        'C00008', 'C00009', 'C00010', 'C00011', 'C00013', 'C00014',
        'C00020', 'C00027', 'C00044', 'C00080'
    },

    # Edge tier weights (for enhanced graph)
    'tier_r_weight': 1.0,  # Reaction-grounded edges
    'tier_p_weight': 0.5,  # Pathway-supported edges

    # Maximum context metabolites from KEGG
    'max_context_metabolites': 300,
}

# =============================================================================
# TRAINING CONFIGURATION
# =============================================================================
TRAINING_CONFIG = {
    # Model architecture
    'hidden_channels': 64,
    'out_channels': 64,
    'num_heads': 4,
    'num_layers': 2,

    # Optimization
    'learning_rate': 0.01,
    'weight_decay': 1e-5,
    'epochs': 50,
    'patience': 10,  # Early stopping patience

    # Gradient clipping
    'max_grad_norm': 1.0,

    # Batch size (if using mini-batch training)
    'batch_size': 1024,
}

# =============================================================================
# DATA SPLIT CONFIGURATION
# =============================================================================
SPLIT_CONFIG = {
    # Default split strategy
    # Options: 'node_split', 'edge_split', 'random_link_split'
    'strategy': 'node_split',

    # Split ratios (must sum to 1.0)
    'train_ratio': 0.70,
    'val_ratio': 0.15,
    'test_ratio': 0.15,

    # Whether to use disjoint node splits (for inductive evaluation)
    'disjoint': True,
}

# =============================================================================
# NEGATIVE SAMPLING CONFIGURATION
# =============================================================================
NEGATIVE_SAMPLING_CONFIG = {
    # Sampling strategy
    # Options: 'random', 'hard', 'mixed', 'ec_class'
    # 'ec_class' is recommended: avoids false negatives from same EC class
    'strategy': 'ec_class',

    # For EC-class sampling: fallback when EC info is unavailable
    # Options: 'random', 'hard'
    'ec_class_fallback': 'random',

    # For hard negative sampling: offset range for nearby metabolites
    # (legacy, kept for backward compatibility)
    'hard_offset_min': 1,
    'hard_offset_max': 5,

    # For mixed sampling: ratio of hard negatives
    'hard_ratio': 0.5,

    # Negative to positive ratio
    'neg_ratio': 1.0,
}

# =============================================================================
# EVALUATION CONFIGURATION
# =============================================================================
EVALUATION_CONFIG = {
    # Hits@K values to compute
    'hits_at_k': [1, 3, 10, 20, 50],

    # Number of random seeds for multi-seed evaluation
    'seeds': [42, 123, 456, 789, 1024],

    # Layer configurations to test
    'layer_configs': [2, 3],
}

# =============================================================================
# REPRODUCIBILITY
# =============================================================================
SEED_CONFIG = {
    # Default random seed
    'default_seed': 42,

    # Whether to use deterministic algorithms
    'deterministic': True,
}


# =============================================================================
# V3 CONFIGURATION (학술적으로 방어 가능한 설계)
# =============================================================================
V3_CONFIG = {
    # --- Edge Types ---
    # Layer 1: All edge types for local context
    'layer1_edge_types': [
        ('Enzyme', 'interacts', 'Enzyme'),
        ('Enzyme', 'interacts', 'Protein'),
        ('Protein', 'interacts', 'Enzyme'),
        ('Enzyme', 'catalyzes_R', 'Metabolite'),
        ('Metabolite', 'rev_catalyzes_R', 'Enzyme'),
        ('Metabolite', 'rxn_neighbor', 'Metabolite'),
        ('TF', 'associates', 'Enzyme'),
        ('Enzyme', 'rev_associates', 'TF'),
    ],
    # Layer 2+: Only strong relationships (exclude weak associations)
    # Rationale: rxn_neighbor/TF are weak associations; multi-hop propagation causes noise
    'layer2_edge_types': [
        ('Enzyme', 'interacts', 'Enzyme'),
        ('Enzyme', 'interacts', 'Protein'),
        ('Protein', 'interacts', 'Enzyme'),
        ('Enzyme', 'catalyzes_R', 'Metabolite'),
        ('Metabolite', 'rev_catalyzes_R', 'Enzyme'),
    ],

    # --- rxn_neighbor Configuration ---
    # Maximum reaction distance for Metabolite↔Metabolite edges
    'rxn_neighbor_max_dist': 2,
    # Dropout rate for rxn_neighbor edges (information leakage control)
    'rxn_neighbor_dropout': 0.3,
    # Dropout mode: 'per_run_fixed' (reproducibility) or 'per_epoch' (regularization)
    'rxn_neighbor_dropout_mode': 'per_run_fixed',

    # --- Node Features ---
    # Feature dimension breakdown
    'feature_dim': 64,
    # Metabolite: learnable (61) + omics (3: log2fc, pvalue, n_pathways)
    'met_learnable_dim': 61,
    'met_omics_dim': 3,
    # Enzyme: learnable (61) + proteomics (3: log2fc, pvalue, abundance)
    'enz_learnable_dim': 61,
    'enz_omics_dim': 3,
    # TF: learnable (58) + domain one-hot (6: ERF, WRKY, MYB, bHLH, zinc_finger, other)
    'tf_learnable_dim': 58,
    'tf_domain_dim': 6,

    # TF Domain categories for one-hot encoding
    'tf_domains': ['ERF', 'WRKY', 'MYB', 'bHLH', 'zinc_finger', 'other'],

    # --- Evaluation ---
    # Candidate set sizes
    'full_candidate_size': 306,  # All metabolites
    'experimental_candidate_size': 10,  # MTBLS531 experimental metabolites

    # --- Data Paths ---
    'graph_v3': 'data/processed/graph_v3.pt',

    # --- Training ---
    'learning_rate': 0.01,
    'weight_decay': 1e-5,
    'epochs': 50,
    'patience': 10,
}


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================
def get_config(name: str) -> dict:
    """Get a specific configuration dictionary by name."""
    configs = {
        'data_paths': DATA_PATHS,
        'graph': GRAPH_CONFIG,
        'training': TRAINING_CONFIG,
        'split': SPLIT_CONFIG,
        'negative_sampling': NEGATIVE_SAMPLING_CONFIG,
        'evaluation': EVALUATION_CONFIG,
        'seed': SEED_CONFIG,
    }
    if name not in configs:
        raise ValueError(f"Unknown config: {name}. Available: {list(configs.keys())}")
    return configs[name]


def ensure_directories():
    """Create all necessary directories."""
    for key, path in DATA_PATHS.items():
        if key.endswith('_dir'):
            path.mkdir(parents=True, exist_ok=True)


# Create directories on import
ensure_directories()
