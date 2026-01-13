# Top-1 Novel TF-Enzyme Case Study Panel (Template)

**Figure 1** (or X): **Multi-Layered Evidence for the Novel Candidate [TF Name] - [Enzyme Name]**

## Panel A: Subgraph Context (Network Logic)
*   **Visual**: A small subgraph showing the TF (Source), the Enzyme (Target), and any shared neighbors or intermediate nodes (e.g., TF -> Gene -> Enzyme).
*   **Highlight**: Edges present in **Strict Graph** (solid lines) vs edges only in Full Graph (dashed, if any).
*   **Caption**: "Predicted Link (HGT Score: 0.98) in the absence of text-mining evidence."

## Panel B: Ranking & Confidence (Statistical Rigor)
*   **Table**:
    | Ranking Metric | Value | Baseline (Random) |
    | :--- | :--- | :--- |
    | Rank (Strict) | **#2** | ~1700 |
    | Score (Prob) | **0.98** | 0.001 |
    | Bootstrap CI | **[0.95-0.99]** | - |
*   **Caption**: "Consistently Top-ranked across 100 bootstrap iterations."

## Panel C: Orthogonal Support (Triangulation)
*   **Left Y-Axis**: Protein Expression (Log2FC) from Proteomics (PXD).
*   **Right Y-Axis**: Metabolite Abundance (Log2FC) of the enzyme's product from Metabolomics (MTBLS).
*   **Bar Chart**: Side-by-side bars for [TF], [Enzyme], [Product].
*   **Motif Inset**: "Promoter Scan: 2 hits of [TF-Family] motif within -2kb TSS (q < 0.05)."
*   **Caption**: "Co-upregulation in Ethylene treatment + Binding site support."

---
**Why is this "Novel"?**
1.  **No direct link** in STRING/KEGG/Araport.
2.  **Strict Graph prediction**: Survived exclusion of text-mining.
3.  **Orthogonal Validation**: Supported by independent Omics (Condition A) + Sequence Motif (Static).
