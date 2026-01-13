import torch
import torch_geometric.transforms as T
from torch_geometric.data import HeteroData
from src.model import SimpleMLP, HeteroSAGE, LinkPredictor, HGT, HAN
import torch.nn.functional as F
from sklearn.metrics import average_precision_score, roc_auc_score
import numpy as np
import os
import argparse

def train_baseline(model_type, graph_path):
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Training {model_type} on {graph_path} using {device}")
    
    # 1. Load Data
    data = torch.load(graph_path)
    data = data.to(device)
    
    # 2. Split (Same Seed as Refined Trainer)
    torch.manual_seed(42)
    np.random.seed(42)
    
    num_enzymes = data['Enzyme'].num_nodes
    indices = torch.randperm(num_enzymes, device=device)
    split = int(0.9 * num_enzymes)
    train_enz_mask = torch.zeros(num_enzymes, dtype=torch.bool, device=device)
    train_enz_mask[indices[:split]] = True
    test_enz_mask = ~train_enz_mask
    
    # Create Train Graph (Remove edges connected to Test Enzymes)
    edge_index = data['Enzyme', 'catalyzes', 'Metabolite'].edge_index
    src = edge_index[0]
    mask = train_enz_mask[src]
    train_edges = edge_index[:, mask]
    test_edges = edge_index[:, ~mask]
    
    train_data = data.clone()
    train_data['Enzyme', 'catalyzes', 'Metabolite'].edge_index = train_edges
    
    rev_index = data['Metabolite', 'rev_catalyzes', 'Enzyme'].edge_index
    rev_dst = rev_index[1] 
    rev_mask = train_enz_mask[rev_dst]
    train_data['Metabolite', 'rev_catalyzes', 'Enzyme'].edge_index = rev_index[:, rev_mask]
    
    # 3. Model Logic
    hidden_dim = 64
    out_dim = 64
    
    if model_type == 'MLP':
        model = SimpleMLP(train_data.metadata(), 64, hidden_dim, out_dim).to(device)
    elif model_type == 'SAGE':
        model = HeteroSAGE(train_data.metadata(), 64, hidden_dim, out_dim).to(device)
    elif model_type == 'HAN':
        model = HAN(train_data.metadata(), 64, hidden_dim, out_dim, num_heads=2).to(device)
    else:
        raise ValueError("Unknown model type")
        
    predictor = LinkPredictor(64).to(device)
    optimizer = torch.optim.Adam(list(model.parameters()) + list(predictor.parameters()), lr=0.01)
    
    # 4. Training Loop
    epochs = 20
    for epoch in range(1, epochs + 1):
        model.train()
        optimizer.zero_grad()
        
        # Forward
        x_dict = model(train_data.x_dict, train_data.edge_index_dict)
        
        # Pos
        pos_edge_index = train_data['Enzyme', 'catalyzes', 'Metabolite'].edge_index
        
        # Neg
        num_pos = pos_edge_index.size(1)
        neg_src = indices[:split][torch.randint(0, split, (num_pos,), device=device)]
        neg_dst = torch.randint(0, train_data['Metabolite'].num_nodes, (num_pos,), device=device)
        
        pos_out = predictor(x_dict['Enzyme'], x_dict['Metabolite'], pos_edge_index)
        neg_out = predictor(x_dict['Enzyme'], x_dict['Metabolite'], torch.stack([neg_src, neg_dst]))
        
        loss = -torch.log(torch.sigmoid(pos_out) + 1e-15).mean() - torch.log(1 - torch.sigmoid(neg_out) + 1e-15).mean()
        loss.backward()
        optimizer.step()
        
        if epoch % 5 == 0:
            stats = evaluate(model, predictor, train_data, test_edges, indices[split:], device)
            print(f"Epoch {epoch:02d}, Loss: {loss.item():.4f} | Test Hits@20: {stats['Hits@20']:.4f}, MAP: {stats['MAP']:.4f}")

def evaluate(model, predictor, data, test_edges, test_enzymes, device):
    model.eval()
    with torch.no_grad():
        x_dict = model(data.x_dict, data.edge_index_dict)
        enz_emb = x_dict['Enzyme']
        met_emb = x_dict['Metabolite']
        
        test_src = test_edges[0] 
        test_dst = test_edges[1]
        unique_mets = torch.unique(test_dst)
        
        hits_20 = 0
        ap_sum = 0
        count = 0
        
        for met_idx in unique_mets:
            mask = (test_dst == met_idx)
            true_enzs = test_src[mask]
            
            candidate_enzs = test_enzymes
            
            num_cands = candidate_enzs.size(0)
            eval_src = candidate_enzs
            eval_dst = met_idx.repeat(num_cands)
            eval_edges = torch.stack([eval_src, eval_dst])
            
            scores = predictor(enz_emb, met_emb, eval_edges).sigmoid()
            
            is_true = torch.isin(candidate_enzs, true_enzs)
            
            if is_true.sum() == 0: continue
            
            sorted_indices = torch.argsort(scores, descending=True)
            sorted_labels = is_true[sorted_indices]
            
            if sorted_labels[:20].sum() > 0:
                hits_20 += 1
                
            ap = average_precision_score(is_true.cpu().numpy(), scores.cpu().numpy())
            ap_sum += ap
            count += 1
            
        if count == 0: return {'Hits@20': 0, 'MAP': 0}
        return {'Hits@20': hits_20/count, 'MAP': ap_sum/count}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--model', type=str, required=True, choices=['MLP', 'SAGE', 'HAN'])
    parser.add_argument('--graph', type=str, default='data/processed/strict_bipartite.pt')
    args = parser.parse_args()
    
    train_baseline(args.model, args.graph)
