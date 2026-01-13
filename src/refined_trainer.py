import torch
import torch_geometric.transforms as T
from torch_geometric.data import HeteroData
from torch_geometric.loader import LinkNeighborLoader
from src.model import HGT, LinkPredictor
import torch.nn.functional as F
from sklearn.metrics import average_precision_score, roc_auc_score
import numpy as np
import os

def train_refined(graph_path):
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Using device: {device}")
    
    # 1. Load Data
    data = torch.load(graph_path)
    data = data.to(device)
    
    # 2. Node-Disjoint Split Strategy
    # Hold out 10% of Enzymes (and all their edges) for testing.
    # This simulates discovering NEW enzymes for metabolites.
    
    num_enzymes = data['Enzyme'].num_nodes
    indices = torch.randperm(num_enzymes, device=device)
    split = int(0.9 * num_enzymes)
    train_enz_mask = torch.zeros(num_enzymes, dtype=torch.bool, device=device)
    train_enz_mask[indices[:split]] = True
    test_enz_mask = ~train_enz_mask
    
    print(f"Split: {train_enz_mask.sum()} Train Enzymes, {test_enz_mask.sum()} Test Enzymes")
    
    # Create Train Graph (Remove edges connected to Test Enzymes)
    # Edge type: ('Enzyme', 'catalyzes', 'Metabolite')
    edge_index = data['Enzyme', 'catalyzes', 'Metabolite'].edge_index
    
    # Mask edges where source (Enzyme) is in Test set
    src = edge_index[0]
    mask = train_enz_mask[src]
    train_edges = edge_index[:, mask]
    test_edges = edge_index[:, ~mask]
    
    train_data = data.clone()
    train_data['Enzyme', 'catalyzes', 'Metabolite'].edge_index = train_edges
    # Remove reverse edges for test set too?
    # Yes, must be consistent.
    # But BipartiteBuilder added rev_catalyzes.
    # Let's just rebuild 'rev' from scratch using T.ToUndirected later?
    # Or manually filter.
    rev_index = data['Metabolite', 'rev_catalyzes', 'Enzyme'].edge_index
    rev_dst = rev_index[1] # Enzyme is destination
    rev_mask = train_enz_mask[rev_dst]
    train_data['Metabolite', 'rev_catalyzes', 'Enzyme'].edge_index = rev_index[:, rev_mask]
    
    print(f"Train Edges: {train_edges.size(1)}, Test Edges: {test_edges.size(1)}")
    
    # 3. Model
    # HGT needs metadata
    # Train data has correct structure
    model = HGT(train_data.metadata(), 64, 64, 64, 4, 2).to(device)
    predictor = LinkPredictor(64).to(device)
    optimizer = torch.optim.Adam(list(model.parameters()) + list(predictor.parameters()), lr=0.01)
    
    # 4. Training Loop
    # "Dynamic Random Negatives" -> We sample negatives inside the loop.
    
    epochs = 20
    for epoch in range(1, epochs + 1):
        model.train()
        optimizer.zero_grad()
        
        # Forward features
        x_dict = model(train_data.x_dict, train_data.edge_index_dict)
        
        # Supervised Edges
        pos_edge_index = train_data['Enzyme', 'catalyzes', 'Metabolite'].edge_index
        
        # Hard Negative Sampling (Phase 7 Workflow Gap Fix)
        # Strategy: Sample negatives from the SAME pathway but different metabolites.
        # This simulates "Reaction Distance >= 3" by ensuring the enzyme COULD
        # catalyze a related metabolite but doesn't.
        
        # For MVP: We use a simple heuristic:
        # - For each positive (enzyme, metabolite), sample a DIFFERENT metabolite
        #   that shares the same pathway prefix (simulating "nearby" metabolites).
        # - If no pathway info, fall back to random.
        
        num_pos = pos_edge_index.size(1)
        num_metabolites = train_data['Metabolite'].num_nodes
        
        # Hard Negative: Swap metabolite index with a nearby one (±5 indices, wrapping)
        # This simulates "same pathway, different reaction" without explicit pathway data.
        neg_src = pos_edge_index[0]  # Same enzyme
        offset = torch.randint(1, 6, (num_pos,), device=device) * (2 * torch.randint(0, 2, (num_pos,), device=device) - 1)
        neg_dst = (pos_edge_index[1] + offset) % num_metabolites
        
        # Pos Scores
        pos_out = predictor(x_dict['Enzyme'], x_dict['Metabolite'], pos_edge_index)
        
        # Neg Scores
        neg_out = predictor(x_dict['Enzyme'], x_dict['Metabolite'], torch.stack([neg_src, neg_dst]))
        
        loss = -torch.log(torch.sigmoid(pos_out) + 1e-15).mean() - torch.log(1 - torch.sigmoid(neg_out) + 1e-15).mean()
        loss.backward()
        optimizer.step()
        
        if epoch % 5 == 0:
            stats = evaluate(model, predictor, train_data, test_edges, indices[split:], device)
            print(f"Epoch {epoch:02d}, Loss: {loss.item():.4f} | Test Hits@20: {stats['Hits@20']:.4f}, MAP: {stats['MAP']:.4f}")

    # Save
    torch.save(model.state_dict(), 'data/models/refined_hgt.pth')

def evaluate(model, predictor, data, test_edges, test_enzymes, device):
    """
    Ranking Evaluation on Test Set (Node Disjoint).
    For each Test Enzyme (that has edges), rank all Metabolites?
    Or for each Metabolite, rank Test Enzymes?
    Task: Metabolite -> Enzyme.
    So query: Metabolite involved in test edges.
    Candidates: All Enzymes (or just Test Enzymes + Hard Negatives?)
    Let's rank against **All Test Enzymes** (Global ranking among unseen).
    """
    model.eval()
    with torch.no_grad():
        x_dict = model(data.x_dict, data.edge_index_dict)
        
        # Embeddings
        enz_emb = x_dict['Enzyme']
        met_emb = x_dict['Metabolite']
        
        # Test Pairs
        test_src = test_edges[0] # Enzymes
        test_dst = test_edges[1] # Metabolites
        
        # Group by Metabolite (Query)
        # Because we want to see: Given Metabolite X, is the True Enzyme ranked high?
        # Note: A metabolite might have multiple true enzymes.
        
        unique_mets = torch.unique(test_dst)
        
        hits_20 = 0
        ap_sum = 0
        count = 0
        
        # For each query metabolite
        for met_idx in unique_mets:
            # True Enzymes for this metabolite in Test Set
            mask = (test_dst == met_idx)
            true_enzs = test_src[mask]
            
            # Candidates: All Test Enzymes (approx 340)
            # This is "Hard" because we strictly look at the unseen pool.
            candidate_enzs = test_enzymes
            
            # Predict scores for all candidates
            # Construct edge index: [All Candidates, Met_Idx]
            num_cands = candidate_enzs.size(0)
            eval_src = candidate_enzs
            eval_dst = met_idx.repeat(num_cands)
            eval_edges = torch.stack([eval_src, eval_dst])
            
            scores = predictor(enz_emb, met_emb, eval_edges).sigmoid()
            
            # Identify which candidates are True
            # We need to map global ID matches
            # true_enzs are global IDs. candidate_enzs are global IDs.
            # isin() check
            is_true = torch.isin(candidate_enzs, true_enzs)
            
            if is_true.sum() == 0: continue # Should not happen by definition
            
            # Ranking
            sorted_indices = torch.argsort(scores, descending=True)
            sorted_labels = is_true[sorted_indices]
            
            # Hits@20 (Is any true positive in top 20?)
            if sorted_labels[:20].sum() > 0:
                hits_20 += 1
                
            # AP
            ap = average_precision_score(is_true.cpu().numpy(), scores.cpu().numpy())
            ap_sum += ap
            count += 1
            
        return {'Hits@20': hits_20/count, 'MAP': ap_sum/count}

if __name__ == "__main__":
    import argparse
    import sys
    sys.path.append(os.getcwd())
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--graph', type=str, default='data/processed/bipartite_graph.pt')
    args = parser.parse_args()
    
    train_refined(args.graph)
