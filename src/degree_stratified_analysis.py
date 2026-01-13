"""
Degree-Stratified Performance Analysis

Compares HGT vs Heuristic baselines across enzyme degree quartiles.
Demonstrates GNN value in sparse-connectivity regions.
"""

import torch
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
sys.path.append(os.getcwd())

from src.model import HGT, LinkPredictor

def set_seed(seed):
    torch.manual_seed(seed)
    np.random.seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)

def compute_enzyme_degrees(data):
    """Compute degree from catalyzes edges."""
    edge_index = data['Enzyme', 'catalyzes', 'Metabolite'].edge_index
    degrees = torch.zeros(data['Enzyme'].num_nodes, device=edge_index.device)
    for src in edge_index[0]:
        degrees[src] += 1
    return degrees

def adamic_adar_score(neighbors, src_idx, dst_idx):
    """Simple Adamic-Adar for enzyme-metabolite."""
    src_neighbors = neighbors.get(src_idx, set())
    dst_neighbors = neighbors.get(dst_idx, set())
    common = src_neighbors & dst_neighbors
    
    if len(common) == 0:
        return 0.0
    
    score = 0.0
    for n in common:
        degree = len(neighbors.get(n, set()))
        if degree > 1:
            score += 1.0 / np.log(degree)
    return score

def train_model(data, train_mask, device, epochs=20, seed=42):
    """Train HGT model."""
    set_seed(seed)
    
    train_data = data.clone()
    edge_index = data['Enzyme', 'catalyzes', 'Metabolite'].edge_index
    mask = train_mask[edge_index[0]]
    train_data['Enzyme', 'catalyzes', 'Metabolite'].edge_index = edge_index[:, mask]
    
    rev_index = data['Metabolite', 'rev_catalyzes', 'Enzyme'].edge_index
    rev_mask = train_mask[rev_index[1]]
    train_data['Metabolite', 'rev_catalyzes', 'Enzyme'].edge_index = rev_index[:, rev_mask]
    
    model = HGT(train_data.metadata(), 64, 64, 64, 4, 2).to(device)
    predictor = LinkPredictor(64).to(device)
    optimizer = torch.optim.Adam(list(model.parameters()) + list(predictor.parameters()), lr=0.01)
    
    for epoch in range(epochs):
        model.train()
        optimizer.zero_grad()
        
        x_dict = model(train_data.x_dict, train_data.edge_index_dict)
        pos_edge = train_data['Enzyme', 'catalyzes', 'Metabolite'].edge_index
        
        if pos_edge.size(1) == 0:
            continue
        
        num_pos = pos_edge.size(1)
        num_met = train_data['Metabolite'].num_nodes
        neg_src = pos_edge[0]
        offset = torch.randint(1, 6, (num_pos,), device=device) * (2 * torch.randint(0, 2, (num_pos,), device=device) - 1)
        neg_dst = (pos_edge[1] + offset) % num_met
        
        pos_out = predictor(x_dict['Enzyme'], x_dict['Metabolite'], pos_edge)
        neg_out = predictor(x_dict['Enzyme'], x_dict['Metabolite'], torch.stack([neg_src, neg_dst]))
        
        loss = -torch.log(torch.sigmoid(pos_out) + 1e-15).mean() - torch.log(1 - torch.sigmoid(neg_out) + 1e-15).mean()
        loss.backward()
        optimizer.step()
    
    return model, predictor, train_data

def evaluate_by_degree_bin(data, train_mask, model, predictor, degrees, degree_bins, device):
    """Evaluate Hits@20 for each degree bin."""
    model.eval()
    
    train_data = data.clone()
    edge_index = data['Enzyme', 'catalyzes', 'Metabolite'].edge_index
    
    # Create mask for edges from training enzymes
    edge_src = edge_index[0].cpu()
    train_mask_cpu = train_mask.cpu()
    mask = train_mask_cpu[edge_src]
    
    train_data['Enzyme', 'catalyzes', 'Metabolite'].edge_index = edge_index[:, mask.to(device)]
    
    test_mask = ~train_mask
    test_enzymes = torch.where(test_mask)[0]
    test_degrees = degrees[test_enzymes].cpu().numpy()
    
    # Get test edges
    test_edge_mask = ~mask
    test_edges = edge_index[:, test_edge_mask.to(device)]
    test_src = test_edges[0]
    test_dst = test_edges[1]
    unique_mets = torch.unique(test_dst)
    
    with torch.no_grad():
        x_dict = model(train_data.x_dict, train_data.edge_index_dict)
        enz_emb = x_dict['Enzyme']
        met_emb = x_dict['Metabolite']
    
    results = {bin_label: {'hits': 0, 'total': 0} for bin_label in degree_bins.keys()}
    test_enz_degrees = {enz.item(): degrees[enz].item() for enz in test_enzymes}
    
    for met_idx in unique_mets:
        met_mask = (test_dst == met_idx)
        true_enzs = set(test_src[met_mask].cpu().numpy())
        
        # Score all test enzymes
        num_cands = test_enzymes.size(0)
        eval_edges = torch.stack([test_enzymes, met_idx.repeat(num_cands)])
        
        with torch.no_grad():
            scores = predictor(enz_emb, met_emb, eval_edges).sigmoid().cpu().numpy()
        
        for enz, score in zip(test_enzymes.cpu().numpy(), scores):
            enz_degree = test_enz_degrees[enz]
            
            # Determine bin
            bin_label = None
            for label, (low, high) in degree_bins.items():
                if low <= enz_degree < high:
                    bin_label = label
                    break
            
            if bin_label is None:
                continue
            
            # Check if this enzyme is in top-20 for this metabolite
            sorted_indices = np.argsort(-scores)
            top_20_enzymes = test_enzymes.cpu().numpy()[sorted_indices[:20]]
            
            if enz in true_enzs:
                results[bin_label]['total'] += 1
                if enz in top_20_enzymes:
                    results[bin_label]['hits'] += 1
    
    # Compute Hits@20 per bin
    bin_hits20 = {}
    for label, data in results.items():
        if data['total'] > 0:
            bin_hits20[label] = data['hits'] / data['total']
        else:
            bin_hits20[label] = 0.0
    
    return bin_hits20

def compute_heuristic_by_degree(data, train_mask, degrees, degree_bins, device):
    """Compute Adamic-Adar performance by degree bin."""
    edge_index = data['Enzyme', 'catalyzes', 'Metabolite'].edge_index
    
    # Create mask for edges from training enzymes
    edge_src = edge_index[0].cpu()
    train_mask_cpu = train_mask.cpu()
    mask = train_mask_cpu[edge_src]
    
    # Build neighbor dict from training edges
    train_edges = edge_index[:, mask.to(device)]
    neighbors = {}
    for i in range(train_edges.size(1)):
        enz = train_edges[0, i].item()
        met = train_edges[1, i].item()
        if enz not in neighbors:
            neighbors[enz] = set()
        neighbors[enz].add(met)
    
    test_mask = ~train_mask
    test_enzymes = torch.where(test_mask)[0]
    
    # Get test edges
    test_edge_mask = ~mask
    test_edges = edge_index[:, test_edge_mask.to(device)]
    test_src = test_edges[0]
    test_dst = test_edges[1]
    unique_mets = torch.unique(test_dst)
    
    test_enz_degrees = {enz.item(): degrees[enz].item() for enz in test_enzymes}
    
    results = {bin_label: {'hits': 0, 'total': 0} for bin_label in degree_bins.keys()}
    
    for met_idx in unique_mets:
        met_mask = (test_dst == met_idx)
        true_enzs = set(test_src[met_mask].cpu().numpy())
        
        # Score all test enzymes with AA
        scores = []
        for enz in test_enzymes.cpu().numpy():
            score = adamic_adar_score(neighbors, enz, met_idx.item())
            scores.append(score)
        
        scores = np.array(scores)
        sorted_indices = np.argsort(-scores)
        top_20_enzymes = test_enzymes.cpu().numpy()[sorted_indices[:20]]
        
        for enz in test_enzymes.cpu().numpy():
            enz_degree = test_enz_degrees[enz]
            
            bin_label = None
            for label, (low, high) in degree_bins.items():
                if low <= enz_degree < high:
                    bin_label = label
                    break
            
            if bin_label is None:
                continue
            
            if enz in true_enzs:
                results[bin_label]['total'] += 1
                if enz in top_20_enzymes:
                    results[bin_label]['hits'] += 1
    
    bin_hits20 = {}
    for label, data in results.items():
        if data['total'] > 0:
            bin_hits20[label] = data['hits'] / data['total']
        else:
            bin_hits20[label] = 0.0
    
    return bin_hits20

def main():
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Using device: {device}")
    
    set_seed(42)
    
    data = torch.load('data/processed/strict_bipartite_v2.pt')
    data = data.to(device)
    
    # Compute degrees
    degrees = compute_enzyme_degrees(data)
    
    # Define degree bins (quartiles)
    degree_values = degrees.cpu().numpy()
    q25, q50, q75 = np.percentile(degree_values[degree_values > 0], [25, 50, 75])
    
    degree_bins = {
        'Q1 (0-25%)': (0, q25),
        'Q2 (25-50%)': (q25, q50),
        'Q3 (50-75%)': (q50, q75),
        'Q4 (75-100%)': (q75, float('inf'))
    }
    
    print(f"Degree quartiles: Q25={q25:.1f}, Q50={q50:.1f}, Q75={q75:.1f}")
    
    # Create train/test split
    num_enzymes = data['Enzyme'].num_nodes
    indices = torch.randperm(num_enzymes)
    split = int(0.9 * num_enzymes)
    train_mask = torch.zeros(num_enzymes, dtype=torch.bool, device=device)
    train_mask[indices[:split]] = True
    
    # Train HGT
    print("Training HGT...")
    model, predictor, train_data = train_model(data, train_mask, device)
    
    # Evaluate by degree bin
    print("Evaluating HGT by degree bin...")
    hgt_results = evaluate_by_degree_bin(data, train_mask, model, predictor, degrees, degree_bins, device)
    
    print("Evaluating Adamic-Adar by degree bin...")
    aa_results = compute_heuristic_by_degree(data, train_mask, degrees, degree_bins, device)
    
    print("\n" + "="*60)
    print("DEGREE-STRATIFIED RESULTS")
    print("="*60)
    
    for bin_label in degree_bins.keys():
        hgt = hgt_results.get(bin_label, 0) * 100
        aa = aa_results.get(bin_label, 0) * 100
        print(f"{bin_label}: HGT={hgt:.1f}%, AA={aa:.1f}%")
    
    # Create figure
    fig, ax = plt.subplots(figsize=(10, 6))
    
    labels = list(degree_bins.keys())
    hgt_vals = [hgt_results.get(l, 0) * 100 for l in labels]
    aa_vals = [aa_results.get(l, 0) * 100 for l in labels]
    
    x = np.arange(len(labels))
    width = 0.35
    
    bars1 = ax.bar(x - width/2, hgt_vals, width, label='HGT', color='#2A9D8F')
    bars2 = ax.bar(x + width/2, aa_vals, width, label='Adamic-Adar', color='#E9C46A')
    
    ax.set_xlabel('Enzyme Degree Percentile', fontsize=12)
    ax.set_ylabel('Hits@20 (%)', fontsize=12)
    ax.set_title('Performance by Enzyme Connectivity: HGT vs. Heuristic', fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=10)
    ax.legend()
    ax.set_ylim(0, max(max(hgt_vals), max(aa_vals)) * 1.3)
    
    # Add value labels
    for bar, val in zip(bars1, hgt_vals):
        ax.annotate(f'{val:.1f}%', xy=(bar.get_x() + bar.get_width()/2, bar.get_height()),
                   ha='center', va='bottom', fontsize=9)
    for bar, val in zip(bars2, aa_vals):
        ax.annotate(f'{val:.1f}%', xy=(bar.get_x() + bar.get_width()/2, bar.get_height()),
                   ha='center', va='bottom', fontsize=9)
    
    plt.tight_layout()
    
    os.makedirs('results/figures', exist_ok=True)
    plt.savefig('results/figures/degree_stratified_performance.png', dpi=150, bbox_inches='tight')
    plt.close()
    
    print("\nFigure saved to results/figures/degree_stratified_performance.png")
    
    # Save data
    with open('results/gnn/degree_stratified.tsv', 'w') as f:
        f.write("Degree_Bin\tHGT_Hits20\tAA_Hits20\n")
        for label in labels:
            hgt = hgt_results.get(label, 0)
            aa = aa_results.get(label, 0)
            f.write(f"{label}\t{hgt:.4f}\t{aa:.4f}\n")
    print("Data saved to results/gnn/degree_stratified.tsv")

if __name__ == "__main__":
    main()
