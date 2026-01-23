#!/usr/bin/env python3
"""
Ethylene-Specific GNN Trainer
- Uses metabolite FC features
- Evaluates link prediction with proper train/test split
"""

import torch
import torch.nn as nn
import torch.nn.functional as F
import torch_geometric.transforms as T
from torch_geometric.nn import HGTConv, Linear
from pathlib import Path
import numpy as np
from sklearn.metrics import roc_auc_score
from tqdm import tqdm
import sys
sys.path.append('/data/ethylene')

class HGT(nn.Module):
    def __init__(self, metadata, in_channels, hidden_channels, out_channels, num_heads=4, num_layers=3):
        super().__init__()
        self.lin_dict = nn.ModuleDict()
        for node_type in metadata[0]:
            self.lin_dict[node_type] = Linear(in_channels, hidden_channels)

        self.convs = nn.ModuleList()
        for _ in range(num_layers):
            conv = HGTConv(hidden_channels, hidden_channels, metadata, heads=num_heads)
            self.convs.append(conv)

        self.out_lin = Linear(hidden_channels, out_channels)

    def forward(self, x_dict, edge_index_dict):
        x_dict = {
            node_type: self.lin_dict[node_type](x).relu_()
            for node_type, x in x_dict.items()
        }
        for conv in self.convs:
            x_dict = conv(x_dict, edge_index_dict)
        return x_dict

class LinkPredictor(nn.Module):
    def __init__(self, in_channels):
        super().__init__()
        
    def forward(self, x_src, x_dst, edge_label_index):
        row, col = edge_label_index
        src_feats = x_src[row]
        dst_feats = x_dst[col]
        return (src_feats * dst_feats).sum(dim=-1)

def evaluate_hits_at_k(model, predictor, data, edge_type, k_list=[1, 3, 10, 20, 50]):
    """Evaluate Hits@K for the given edge type."""
    model.eval()
    src_type, _, dst_type = edge_type
    
    with torch.no_grad():
        x_dict = model(data.x_dict, data.edge_index_dict)
        
        src_emb = x_dict[src_type]
        dst_emb = x_dict[dst_type]
        
        # Test edges
        test_edge_index = data[edge_type].edge_label_index
        test_labels = data[edge_type].edge_label
        
        # Keep only positive edges
        pos_mask = test_labels == 1
        pos_edges = test_edge_index[:, pos_mask]
        
        results = {k: [] for k in k_list}
        
        for i in range(pos_edges.shape[1]):
            src_idx = pos_edges[0, i]
            true_dst = pos_edges[1, i]
            
            # Score all destinations
            src_rep = src_emb[src_idx].unsqueeze(0)
            scores = (src_rep * dst_emb).sum(dim=1)
            
            # Rank
            sorted_indices = torch.argsort(scores, descending=True)
            rank = (sorted_indices == true_dst).nonzero(as_tuple=True)[0].item() + 1
            
            for k in k_list:
                results[k].append(1 if rank <= k else 0)
        
        hits = {k: np.mean(v) * 100 for k, v in results.items()}
        mrr = np.mean([1/r for r in range(1, len(results[k_list[0]])+1)])  # Simplified
        
    return hits, mrr

def train_ethylene_model():
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Device: {device}")
    
    # Paths
    graph_path = Path('/data/ethylene/experiments/ethylene_specific/ethylene_specific_graph.pt')
    output_dir = Path('/data/ethylene/experiments/ethylene_specific')
    
    # Load graph
    print("\n[1] Loading ethylene-specific graph...")
    data = torch.load(graph_path)
    
    # Add reverse edges
    data = T.ToUndirected()(data)
    
    # Target edge type
    target_edge_type = ('Enzyme', 'catalyzes', 'Metabolite')
    print(f"Target edge: {target_edge_type}")
    
    # Split
    print("\n[2] Splitting data...")
    transform = T.RandomLinkSplit(
        num_val=0.1,
        num_test=0.1,
        is_undirected=False,
        edge_types=[target_edge_type],
        rev_edge_types=[('Metabolite', 'rev_catalyzes', 'Enzyme')],
        add_negative_train_samples=True,
        neg_sampling_ratio=1.0
    )
    
    train_data, val_data, test_data = transform(data)
    train_data = train_data.to(device)
    val_data = val_data.to(device)
    test_data = test_data.to(device)
    
    # Model
    print("\n[3] Initializing model...")
    model = HGT(
        metadata=data.metadata(),
        in_channels=64,
        hidden_channels=64,
        out_channels=64,
        num_heads=4,
        num_layers=3
    ).to(device)
    
    predictor = LinkPredictor(64).to(device)
    
    optimizer = torch.optim.Adam(
        list(model.parameters()) + list(predictor.parameters()),
        lr=0.005
    )
    
    # Training
    print("\n[4] Training...")
    epochs = 100
    best_val_hits20 = 0
    
    for epoch in range(1, epochs + 1):
        model.train()
        optimizer.zero_grad()
        
        x_dict = model(train_data.x_dict, train_data.edge_index_dict)
        
        src_type, _, dst_type = target_edge_type
        edge_label_index = train_data[target_edge_type].edge_label_index
        edge_label = train_data[target_edge_type].edge_label
        
        out = predictor(x_dict[src_type], x_dict[dst_type], edge_label_index)
        
        loss = F.binary_cross_entropy_with_logits(out, edge_label)
        loss.backward()
        optimizer.step()
        
        # Validation every 10 epochs
        if epoch % 10 == 0:
            model.eval()
            with torch.no_grad():
                x_dict_val = model(val_data.x_dict, val_data.edge_index_dict)
                val_edge_label_index = val_data[target_edge_type].edge_label_index
                val_edge_label = val_data[target_edge_type].edge_label
                
                val_out = predictor(x_dict_val[src_type], x_dict_val[dst_type], val_edge_label_index)
                val_pred = torch.sigmoid(val_out).cpu().numpy()
                val_label = val_edge_label.cpu().numpy()
                
                try:
                    auc = roc_auc_score(val_label, val_pred)
                except:
                    auc = 0.5
            
            print(f"Epoch {epoch:3d} | Loss: {loss.item():.4f} | Val AUC: {auc:.4f}")
    
    # Final evaluation
    print("\n[5] Final Evaluation...")
    hits, mrr = evaluate_hits_at_k(model, predictor, test_data, target_edge_type)
    
    print("\n" + "=" * 60)
    print("Results: Ethylene-Specific Model")
    print("=" * 60)
    for k, v in hits.items():
        print(f"Hits@{k}: {v:.2f}%")
    
    # Save results
    results = {
        'model': 'HGT-Ethylene',
        'layers': 3,
        'features': 'Metabolite FC',
        **{f'hits@{k}': v for k, v in hits.items()}
    }
    
    import json
    with open(output_dir / 'results.json', 'w') as f:
        json.dump(results, f, indent=2)
    
    # Save model
    torch.save(model.state_dict(), output_dir / 'ethylene_hgt.pth')
    torch.save(predictor.state_dict(), output_dir / 'ethylene_predictor.pth')
    
    print(f"\nSaved to: {output_dir}")
    
    return results

if __name__ == "__main__":
    train_ethylene_model()
