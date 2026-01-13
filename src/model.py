import torch
import torch.nn as nn
from torch_geometric.nn import HGTConv, Linear

class HGT(nn.Module):
    def __init__(self, metadata, in_channels, hidden_channels, out_channels, num_heads=4, num_layers=2):
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
        # Initial projection to hidden_channels
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

class SimpleMLP(nn.Module):
    def __init__(self, metadata, in_channels, hidden_channels, out_channels):
        super().__init__()
        self.lin_dict = nn.ModuleDict()
        # Separate projection for each node type
        for node_type in metadata[0]:
            self.lin_dict[node_type] = nn.Sequential(
                Linear(in_channels, hidden_channels),
                nn.ReLU(),
                Linear(hidden_channels, hidden_channels),
                nn.ReLU(),
                Linear(hidden_channels, out_channels) # Output embedding dim
            )

    def forward(self, x_dict, edge_index_dict):
        # Ignore edge_index, just transform features
        out_dict = {}
        for node_type, x in x_dict.items():
            if node_type in self.lin_dict:
                out_dict[node_type] = self.lin_dict[node_type](x)
            else:
                out_dict[node_type] = x # Should not happen usually
        return out_dict

# GraphSAGE (Hetero via to_hetero or SAGEConv wrapper)
from torch_geometric.nn import SAGEConv, to_hetero

from torch_geometric.nn import SAGEConv, HeteroConv, HANConv

class HAN(nn.Module):
    def __init__(self, metadata, in_channels, hidden_channels, out_channels, num_heads=2, num_layers=2):
        super().__init__()
        self.lin_dict = nn.ModuleDict()
        for node_type in metadata[0]:
            self.lin_dict[node_type] = Linear(in_channels, hidden_channels)

        self.convs = nn.ModuleList()
        for _ in range(num_layers):
            conv = HANConv(hidden_channels, hidden_channels, metadata, heads=num_heads)
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

class HeteroSAGE(nn.Module):
    def __init__(self, metadata, in_channels, hidden_channels, out_channels, num_layers=2):
        super().__init__()
        # 1. Projections
        self.lin_dict = nn.ModuleDict()
        for node_type in metadata[0]:
            self.lin_dict[node_type] = Linear(in_channels, hidden_channels)

        # 2. HeteroConvs
        self.convs = nn.ModuleList()
        for _ in range(num_layers):
            # Define convolution for each edge type
            conv_dict = {}
            for edge_type in metadata[1]:
                # SAGEConv per edge type
                # ('Enzyme', 'catalyzes', 'Metabolite')
                src_type, _, dst_type = edge_type
                conv_dict[edge_type] = SAGEConv(hidden_channels, hidden_channels)
                
            self.convs.append(HeteroConv(conv_dict, aggr='sum'))

        # 3. Output
        self.out_lin = Linear(hidden_channels, out_channels)

    def forward(self, x_dict, edge_index_dict):
        x_dict = {
            node_type: self.lin_dict[node_type](x).relu_()
            for node_type, x in x_dict.items()
        }
        for conv in self.convs:
            x_dict = conv(x_dict, edge_index_dict)
            x_dict = {k: v.relu_() for k, v in x_dict.items()}
            
        return x_dict
