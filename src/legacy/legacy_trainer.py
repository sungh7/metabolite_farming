"""
Legacy Trainer (Original Implementation)

This is the original trainer.py preserved for backward compatibility.
Uses RandomLinkSplit and random negative sampling.

DO NOT USE FOR PRODUCTION - Use src/train.py instead.

Original behavior:
- PyG RandomLinkSplit for data splitting
- Random permutation negative sampling
- Basic BCE loss without early stopping
- Fixed 10 epochs
"""

import torch
import torch_geometric.transforms as T
from torch_geometric.data import HeteroData
from torch_geometric.loader import LinkNeighborLoader
import torch.nn.functional as F
from sklearn.metrics import roc_auc_score
import os
import sys

# Add parent directory for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from src.model import HGT, LinkPredictor


def train(graph_path='data/processed/graph.pt'):
    """
    LEGACY training function.

    Uses original random negative sampling and RandomLinkSplit.
    """
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"[LEGACY] Using device: {device}")

    # Load Graph
    print("Loading graph...")
    data = torch.load(graph_path)

    # Add reverse edges and self-loops
    data = T.ToUndirected()(data)
    data = T.AddSelfLoops()(data)

    data = data.to(device)

    # Define Edge Type to Predict
    target_edge_type = None
    if ('TF', 'interacts', 'Enzyme') in data.edge_types:
        target_edge_type = ('TF', 'interacts', 'Enzyme')
    else:
        target_edge_type = list(data.edge_types)[0]

    print(f"Target Edge Type: {target_edge_type}")

    # Split Data using RandomLinkSplit
    transform = T.RandomLinkSplit(
        num_val=0.1,
        num_test=0.1,
        is_undirected=True,
        edge_types=[target_edge_type],
        add_negative_train_samples=False
    )

    train_data, val_data, test_data = transform(data)

    # Model Init
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
    epochs = 10

    for epoch in range(1, epochs + 1):
        model.train()
        optimizer.zero_grad()

        x_dict = model(train_data.x_dict, train_data.edge_index_dict)

        src_type, _, dst_type = target_edge_type
        edge_label_index = train_data[target_edge_type].edge_label_index

        src_idx = edge_label_index[0]
        dst_idx = edge_label_index[1]

        # Positive Scores
        pos_out = predictor(x_dict[src_type], x_dict[dst_type], edge_label_index)

        # LEGACY: Random shuffle negative sampling
        neg_dst_idx = torch.randint(
            0, x_dict[dst_type].size(0),
            (src_idx.size(0),),
            device=device
        )
        neg_out = predictor(
            x_dict[src_type], x_dict[dst_type],
            torch.stack([src_idx, neg_dst_idx])
        )

        loss = (
            -torch.log(torch.sigmoid(pos_out) + 1e-15).mean()
            - torch.log(1 - torch.sigmoid(neg_out) + 1e-15).mean()
        )

        loss.backward()
        optimizer.step()

        # Validation
        model.eval()
        with torch.no_grad():
            x_dict = model(val_data.x_dict, val_data.edge_index_dict)
            val_edge_label_index = val_data[target_edge_type].edge_label_index
            val_edge_label = val_data[target_edge_type].edge_label

            out = predictor(x_dict[src_type], x_dict[dst_type], val_edge_label_index)
            pred = torch.sigmoid(out).cpu().numpy()
            label = val_edge_label.cpu().numpy()

            auc = roc_auc_score(label, pred)

        print(f"Epoch {epoch:02d}, Loss: {loss.item():.4f}, Val AUC: {auc:.4f}")

    # Save Model
    os.makedirs('data/models', exist_ok=True)
    torch.save(model.state_dict(), 'data/models/legacy_hgt_model.pth')
    torch.save(predictor.state_dict(), 'data/models/legacy_predictor.pth')
    print("Model saved to data/models/")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('--graph', type=str, default='data/processed/graph.pt')
    args = parser.parse_args()

    train(args.graph)
