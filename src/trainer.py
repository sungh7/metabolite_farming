import torch
import torch_geometric.transforms as T
from torch_geometric.data import HeteroData
from torch_geometric.loader import LinkNeighborLoader
from src.model import HGT, LinkPredictor
import torch.nn.functional as F
from sklearn.metrics import roc_auc_score
from tqdm import tqdm

def train():
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Using device: {device}")
    
    # Load Graph
    print("Loading graph...")
    data = torch.load('data/processed/graph.pt')
    
    # Add reverse edges and self-loops (for HAN)
    # HANConv requires metadata
    data = T.ToUndirected()(data)
    data = T.AddSelfLoops()(data) # Optional, HAN often handles it
    
    # Ensure all nodes have features (GraphBuilder initialized them)
    # Move to device
    data = data.to(device)
    
    # Define Edge Type to Predict
    # Let's predict TF -> Enzyme interactions
    # If the graph is sparse, we might need to check if this edge type exists
    target_edge_type = None
    if ('TF', 'interacts', 'Enzyme') in data.edge_types:
        target_edge_type = ('TF', 'interacts', 'Enzyme')
    else:
        # Fallback to Protein-Protein if types are scarce
        # But our builder creates typed edges.
        target_edge_type = list(data.edge_types)[0]
        
    print(f"Target Edge Type for Prediction: {target_edge_type}")
    
    # Split Data (Simple random split for now, strictly link prediction)
    # Because HeteroData split is complex, we will manually mask edges
    # Or use RandomLinkSplit
    transform = T.RandomLinkSplit(
        num_val=0.1,
        num_test=0.1,
        is_undirected=True,
        edge_types=[target_edge_type],
        add_negative_train_samples=False  # We add them in loader or manually
    )
    
    train_data, val_data, test_data = transform(data)
    
    # Model Init
    # Using HGT instead of HAN for better stability
    model = HGT(
        metadata=data.metadata(),
        in_channels=64,
        hidden_channels=64,
        out_channels=64,
        num_heads=4,
        num_layers=2
    ).to(device)
    
    predictor = LinkPredictor(64).to(device)
    
    optimizer = torch.optim.Adam(
        list(model.parameters()) + list(predictor.parameters()), 
        lr=0.01
    )
    
    # Training Loop
    epochs = 10 # Short run for verification
    
    for epoch in range(1, epochs + 1):
        model.train()
        optimizer.zero_grad()
        
        # Forward pass
        x_dict = model(train_data.x_dict, train_data.edge_index_dict)
        
        # Get embeddings for target edge
        src_type, _, dst_type = target_edge_type
        edge_label_index = train_data[target_edge_type].edge_label_index
        edge_label = train_data[target_edge_type].edge_label
        
        # Negative Sampling (On the fly if not provided)
        # RandomLinkSplit provides edge_label = 1 for supervision
        # We need negative samples
        # For simplicity in this script, we just assume RandomLinkSplit provided positives
        # We need to manually add negatives or use a Loader that supports it.
        # Let's stick to what's in train_data (positives) and sample negatives.
        
        # Simple negative sampling:
        # Permute dst indices
        src_idx = edge_label_index[0]
        dst_idx = edge_label_index[1]
        
        # Positive Scores
        pos_out = predictor(x_dict[src_type], x_dict[dst_type], edge_label_index)
        
        # Negative Scores (Random shuffle)
        neg_dst_idx = torch.randint(0, x_dict[dst_type].size(0), (src_idx.size(0),), device=device)
        neg_out = predictor(x_dict[src_type], x_dict[dst_type], torch.stack([src_idx, neg_dst_idx]))
        
        loss = -torch.log(torch.sigmoid(pos_out) + 1e-15).mean() - torch.log(1 - torch.sigmoid(neg_out) + 1e-15).mean()
        
        loss.backward()
        optimizer.step()
        
        # Validation
        model.eval()
        with torch.no_grad():
            x_dict = model(val_data.x_dict, val_data.edge_index_dict)
            val_edge_label_index = val_data[target_edge_type].edge_label_index
            val_edge_label = val_data[target_edge_type].edge_label # 1s and 0s from Split
            
            out = predictor(x_dict[src_type], x_dict[dst_type], val_edge_label_index)
            pred = torch.sigmoid(out).cpu().numpy()
            label = val_edge_label.cpu().numpy()
            
            auc = roc_auc_score(label, pred)
            
        print(f"Epoch {epoch:02d}, Loss: {loss.item():.4f}, Val AUC: {auc:.4f}")
    
    # Save Model
    os.makedirs('data/models', exist_ok=True)
    torch.save(model.state_dict(), 'data/models/hgt_model.pth')
    torch.save(predictor.state_dict(), 'data/models/predictor.pth')
    print("Model saved to data/models/")

if __name__ == "__main__":
    import os
    # Fix import issue
    import sys
    sys.path.append(os.getcwd())
    
    train()
