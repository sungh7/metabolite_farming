import pandas as pd
import numpy as np
from scipy.stats import ttest_ind, ranksums
import os

def simulate_metabolomics_data(n_samples=6, n_metabolites=200):
    """
    Simulates metabolomics data for Control vs Ethylene.
    Features key Soy isoflavonoids with upregulated patterns.
    """
    np.random.seed(42)
    
    # 1. Generate Metadata
    compounds = [f"Metabolite_{i}" for i in range(n_metabolites)]
    
    # Pathway Mapping (Simulated)
    # Phenylpropanoid: 0-19
    # Flavonoid: 20-49
    # Others: 50+
    pathways = {}
    for i in range(20): pathways[compounds[i]] = 'Phenylpropanoid'
    for i in range(20, 50): pathways[compounds[i]] = 'Flavonoid'
    for i in range(50, n_metabolites): pathways[compounds[i]] = 'Other'
    
    # 2. Generate Expression Data (Log2 Intensity)
    # Background: Normal(10, 1)
    control = np.random.normal(10, 1, (n_metabolites, n_samples))
    treatment = np.random.normal(10, 1, (n_metabolites, n_samples))
    
    # 3. Inject Signal (Ethylene Effect)
    # Increase Flavonoids/Phenylpropanoids
    # Effect size: Normal(1.5, 0.5) (Log2 scale -> 2.8x fold change)
    effect_indices = list(range(10)) + list(range(20, 35)) # Selecting specific targets
    
    treatment[effect_indices] += np.random.normal(1.5, 0.3, (len(effect_indices), n_samples))
    
    # DataFrame
    df = pd.DataFrame(data=np.hstack([control, treatment]), index=compounds, 
                      columns=[f"Ctrl_{i}" for i in range(n_samples)] + [f"Trt_{i}" for i in range(n_samples)])
    
    df['Pathway'] = df.index.map(pathways)
    return df

def perform_bootstrap_analysis(df, n_boot=1000):
    """
    Bootstrap resampling to verify pathway enrichment stability.
    """
    print(f"Running Bootstrap Analysis (n={n_boot})...")
    
    results = []
    
    # Original Calculation
    # Exclude 'Pathway' column
    original_fc = df.iloc[:, 6:-1].mean(axis=1) - df.iloc[:, :6].mean(axis=1)
    
    # Bootstrap Loop
    pathway_scores = {'Phenylpropanoid': [], 'Flavonoid': [], 'Other': []}
    
    n_samples = 6
    cols_ctrl = [f"Ctrl_{i}" for i in range(n_samples)]
    cols_trt = [f"Trt_{i}" for i in range(n_samples)]
    
    for _ in range(n_boot):
        # Resample columns with replacement
        boot_ctrl = np.random.choice(cols_ctrl, n_samples, replace=True)
        boot_trt = np.random.choice(cols_trt, n_samples, replace=True)
        
        # Calculate FC for this bootstrap
        boot_fc = df[boot_trt].mean(axis=1) - df[boot_ctrl].mean(axis=1)
        
        # Calculate Pathway Enrichment (Simple Mean FC of members)
        # More rigorous: GSEA or ORA p-value. Using Mean FC for stability check.
        p_fc = boot_fc[df['Pathway'] == 'Phenylpropanoid'].mean()
        f_fc = boot_fc[df['Pathway'] == 'Flavonoid'].mean()
        o_fc = boot_fc[df['Pathway'] == 'Other'].mean()
        
        pathway_scores['Phenylpropanoid'].append(p_fc)
        pathway_scores['Flavonoid'].append(f_fc)
        pathway_scores['Other'].append(o_fc)
        
    # Summarize
    summary = []
    for path, scores in pathway_scores.items():
        scores = np.array(scores)
        mean_score = np.mean(scores)
        ci_lower = np.percentile(scores, 2.5)
        ci_upper = np.percentile(scores, 97.5)
        stability_score = mean_score / np.std(scores) # Signal-to-Noise
        
        summary.append({
            'Pathway': path,
            'Enrichment_Score_Mean': mean_score,
            'CI_Lower': ci_lower,
            'CI_Upper': ci_upper,
            'Stability_Score': stability_score,
            'Significant': ci_lower > 0 # Simple significance check
        })
        
    return pd.DataFrame(summary)

if __name__ == "__main__":
    df = simulate_metabolomics_data()
    stats = perform_bootstrap_analysis(df)
    
    print("\n--- Metabolomics Bootstrap Results ---")
    print(stats)
    
    os.makedirs('results', exist_ok=True)
    stats.to_csv('results/table1_metabolomics.csv', index=False)
    print("Saved to results/table1_metabolomics.csv")
