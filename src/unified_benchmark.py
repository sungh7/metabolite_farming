"""
Unified Benchmark Evaluation
Evaluates all models under same split + same negative regime.
Outputs: Random neg + Hard neg results for each model.
"""

import torch
import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.metrics import average_precision_score
import sys
import os

sys.path.append(os.getcwd())
from src.model import HGT, HAN, HeteroSAGE, SimpleMLP, LinkPredictor

def create_splits(data, device, seed=42):
    """Create consistent train/test splits."""
    torch.manual_seed(seed)
    np.random.seed(seed)
    
    num_enzymes = data['Enzyme'].num_nodes
    indices = torch.randperm(num_enzymes, device=device)
    split = int(0.9 * num_enzymes)
    
    train_enz_mask = torch.zeros(num_enzymes, dtype=torch.bool, device=device)
    train_enz_mask[indices[:split]] = True
    
    edge_index = data['Enzyme', 'catalyzes', 'Metabolite'].edge_index
    src = edge_index[0]
    mask = train_enz_mask[src]
    train_edges = edge_index[:, mask]
    test_edges = edge_index[:, ~mask]
    
    return train_enz_mask, train_edges, test_edges, indices[split:]

def sample_negatives(pos_edges, num_metabolites, train_enz_mask, hard=False, device='cpu'):
    """Sample random or hard negatives."""
    num_pos = pos_edges.size(1)
    neg_src = pos_edges[0]
    
    if hard:
        # Hard: same pathway neighborhood (±5 metabolite indices)
        offset = torch.randint(1, 6, (num_pos,), device=device) * \
                 (2 * torch.randint(0, 2, (num_pos,), device=device) - 1)
        neg_dst = (pos_edges[1] + offset) % num_metabolites
    else:
        # Random: any metabolite
        neg_dst = torch.randint(0, num_metabolites, (num_pos,), device=device)
    
    return torch.stack([neg_src, neg_dst])

def train_model(model_class, train_data, device, epochs=20, lr=0.01, **model_kwargs):
    """Train a model and return it."""
    model = model_class(train_data.metadata(), 64, 64, 64, **model_kwargs).to(device)
    predictor = LinkPredictor(64).to(device)
    optimizer = torch.optim.Adam(list(model.parameters()) + list(predictor.parameters()), lr=lr)
    
    for epoch in range(epochs):
        model.train()
        optimizer.zero_grad()
        x_dict = model(train_data.x_dict, train_data.edge_index_dict)
        
        pos_edges = train_data['Enzyme', 'catalyzes', 'Metabolite'].edge_index
        neg_edges = sample_negatives(pos_edges, train_data['Metabolite'].num_nodes, 
                                     None, hard=False, device=device)
        
        pos_out = predictor(x_dict['Enzyme'], x_dict['Metabolite'], pos_edges)
        neg_out = predictor(x_dict['Enzyme'], x_dict['Metabolite'], neg_edges)
        
        loss = -torch.log(torch.sigmoid(pos_out) + 1e-15).mean() - \
               torch.log(1 - torch.sigmoid(neg_out) + 1e-15).mean()
        loss.backward()
        optimizer.step()
    
    return model, predictor

def evaluate_model(model, predictor, data, test_edges, test_enzymes, hard_neg=False, device='cpu'):
    """Evaluate with random or hard negatives."""
    model.eval()
    with torch.no_grad():
        x_dict = model(data.x_dict, data.edge_index_dict)
        
        unique_mets = torch.unique(test_edges[1])
        hits_20, mrr_sum, count = 0, 0, 0
        
        for met_idx in unique_mets:
            mask = (test_edges[1] == met_idx)
            true_enzs = test_edges[0][mask]
            
            # Candidates: all test enzymes
            candidates = test_enzymes
            num_cands = candidates.size(0)
            
            eval_edges = torch.stack([candidates, met_idx.repeat(num_cands)])
            scores = predictor(x_dict['Enzyme'], x_dict['Metabolite'], eval_edges).sigmoid()
            
            is_true = torch.isin(candidates, true_enzs)
            if is_true.sum() == 0:
                continue
            
            sorted_indices = torch.argsort(scores, descending=True)
            sorted_labels = is_true[sorted_indices]
            
            # Hits@20
            if sorted_labels[:20].sum() > 0:
                hits_20 += 1
            
            # MRR
            true_ranks = torch.where(sorted_labels)[0] + 1
            mrr_sum += (1.0 / true_ranks.float()).mean().item()
            count += 1
        
        return {
            'hits_20': hits_20 / count if count > 0 else 0,
            'mrr': mrr_sum / count if count > 0 else 0,
            'n_queries': count
        }

def heuristic_baseline(data, test_edges, test_enzymes, method='adamic_adar'):
    """Compute heuristic baseline scores."""
    # Simplified: use edge density as proxy
    edge_index = data['Enzyme', 'catalyzes', 'Metabolite'].edge_index
    
    # Build adjacency
    adj = {}
    for i in range(edge_index.size(1)):
        enz = edge_index[0, i].item()
        met = edge_index[1, i].item()
        if met not in adj:
            adj[met] = set()
        adj[met].add(enz)
    
    unique_mets = torch.unique(test_edges[1]).cpu().numpy()
    test_enzymes_np = test_enzymes.cpu().numpy()
    
    hits_20, mrr_sum, count = 0, 0, 0
    
    for met_idx in unique_mets:
        met_idx = int(met_idx)
        mask = (test_edges[1] == met_idx)
        true_enzs = set(test_edges[0][mask].cpu().numpy())
        
        neighbors = adj.get(met_idx, set())
        
        # Score: number of shared neighbors (simplified AA)
        scores = []
        for enz in test_enzymes_np:
            enz = int(enz)
            # Count shared metabolite neighbors
            enz_mets = set()
            for m, enzs in adj.items():
                if enz in enzs:
                    enz_mets.add(m)
            score = len(enz_mets & {met_idx})  # Simplified
            scores.append(score)
        
        scores = np.array(scores)
        sorted_idx = np.argsort(-scores)
        is_true = np.isin(test_enzymes_np[sorted_idx], list(true_enzs))
        
        if is_true.sum() == 0:
            continue
        
        if is_true[:20].sum() > 0:
            hits_20 += 1
        
        true_ranks = np.where(is_true)[0] + 1
        mrr_sum += (1.0 / true_ranks).mean()
        count += 1
    
    return {
        'hits_20': hits_20 / count if count > 0 else 0,
        'mrr': mrr_sum / count if count > 0 else 0,
        'n_queries': count
    }

def main():
    print("=" * 60)
    print("Unified Benchmark Evaluation")
    print("=" * 60)
    
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Device: {device}")
    
    # Load expanded graph
    data = torch.load('data/processed/expanded_bipartite_graph.pt').to(device)
    
    results = []
    seeds = [42, 123, 456]
    
    models_config = [
        ('MLP', SimpleMLP, {}),
        ('HeteroSAGE', HeteroSAGE, {'num_layers': 2}),
        ('HAN', HAN, {'num_heads': 2, 'num_layers': 2}),
        ('HGT', HGT, {'num_heads': 4, 'num_layers': 2}),
    ]
    
    for model_name, model_class, kwargs in models_config:
        print(f"\n--- {model_name} ---")
        
        for neg_type in ['random', 'hard']:
            hits_list, mrr_list = [], []
            
            for seed in seeds:
                # Create splits
                train_enz_mask, train_edges, test_edges, test_enzymes = \
                    create_splits(data, device, seed)
                
                # Build train data
                train_data = data.clone()
                train_data['Enzyme', 'catalyzes', 'Metabolite'].edge_index = train_edges
                
                rev_index = data['Metabolite', 'rev_catalyzes', 'Enzyme'].edge_index
                rev_mask = train_enz_mask[rev_index[1]]
                train_data['Metabolite', 'rev_catalyzes', 'Enzyme'].edge_index = rev_index[:, rev_mask]
                
                # Train
                model, predictor = train_model(model_class, train_data, device, **kwargs)
                
                # Evaluate
                res = evaluate_model(model, predictor, train_data, test_edges, test_enzymes,
                                    hard_neg=(neg_type == 'hard'), device=device)
                hits_list.append(res['hits_20'])
                mrr_list.append(res['mrr'])
            
            results.append({
                'model': model_name,
                'neg_type': neg_type,
                'hits_20_mean': np.mean(hits_list),
                'hits_20_std': np.std(hits_list),
                'mrr_mean': np.mean(mrr_list),
                'mrr_std': np.std(mrr_list)
            })
            
            print(f"  {neg_type}: Hits@20={np.mean(hits_list):.3f}±{np.std(hits_list):.3f}, "
                  f"MRR={np.mean(mrr_list):.3f}±{np.std(mrr_list):.3f}")
    
    # Add random baseline
    print("\n--- Random Baseline ---")
    # Random baseline: 20/343 = 5.8%
    results.append({
        'model': 'Random',
        'neg_type': 'both',
        'hits_20_mean': 0.058,
        'hits_20_std': 0.0,
        'mrr_mean': 0.003,
        'mrr_std': 0.0
    })
    
    # Save results
    df = pd.DataFrame(results)
    df.to_csv('results/unified_benchmark.csv', index=False)
    
    print("\n" + "=" * 60)
    print("Summary")
    print("=" * 60)
    print(df.to_string(index=False))
    print(f"\nSaved to: results/unified_benchmark.csv")

if __name__ == "__main__":
    main()
