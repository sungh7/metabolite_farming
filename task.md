# Tasks

- [x] **Phase 0: Initial Setup & Exploration** (Completed)
    - [x] Setup Environment & STRING Data
    - [x] Initial HGT Model Training (AUC 0.97)
    - [x] Preliminary Candidate Inference

- [x] **Phase 1: Metabolomics Evidence (Statistical Rigor)**
    - [x] Implement `src/metabolomics_sim.py`: Simulate Data & Bootstrap Analysis <!-- id: 15 -->
    - [x] Calculate Pathway Enrichment Confidence Intervals (95% CI) <!-- id: 16 -->
    - [x] **Deliverable**: `results/table1_metabolomics.csv` <!-- id: 17 -->

- [x] **Phase 2: Protein Layer Mapping**
    - [x] Map Metabolites to Reactions (EC) -> Genes (Glyma) (Simulated Mapping) <!-- id: 18 -->
    - [x] Construct Bipartite Graph (Metabolite-Enzyme) <!-- id: 19 -->
    - [x] **Deliverable**: `data/processed/bipartite_graph.pt` <!-- id: 20 -->

- [x] **Phase 3: GNN Task Refinement (Rigorous Eval)**
    - [x] Implement `src/refined_trainer.py`: Node-Disjoint Split & Hard Negative Sampling <!-- id: 21 -->
    - [x] Disable `add_negative_train_samples` in Loader (Handled logic manually in loop) <!-- id: 22 -->
    - [x] Train GNN with Ranking Loss <!-- id: 23 -->
    - [x] Evaluate with Precision@20, Hits@20, AUPRC <!-- id: 24 -->

- [x] **Phase 4: Structural Auxiliary Scoring**
    - [x] Implement Structural Similarity Score (`src/structural_ranking.py`) <!-- id: 25 -->
    - [x] Re-rank Candidates <!-- id: 26 -->

- [x] **Phase 5: Integrated Hypothesis & Reporting**
    - [x] Integrate all results into a cohesive narrative (`walkthrough.md`).
    - [x] Finalize the "Experimental-Free" defense strategy.
    - [x] Refine task definitions and baselines for academic rigor.

- [x] **Phase 6: Real Data Integration (MTBLS531)**
    - [x] Process `maf.tsv` to calculate differential abundance (Log2FC, P-value).
    - [x] Map ChEBI IDs to KEGG/PlantCyc IDs.
    - [x] Re-run Pathway Enrichment Analysis with real data.
    - [x] Update `walkthrough.md` to replace "Synthetic" claims with "Real Data" results.

- [x] **Phase 7: Proteomics Integration (PXD006989)**
    - [x] Download and process `MaxQuant_resutls_txt.zip`.
    - [x] Extract `proteinGroups.txt` and identify differentially expressed proteins.
    - [x] Map Protein IDs to STRING/Glyma IDs.
    - [x] Cross-validate: Do significant metabolites (Phase 6) interact with significant proteins (Phase 7)?

- [x] **Phase 8: Final Refinements & Rigor**
    - [x] Implement "Direct-Reaction Baseline" (Stronger than Pathway Baseline).
    - [x] Calculate Hits@20 against "All Enzyme" pool.
    - [x] Document "MaxQuant/LFQ" details in Reproducibility.
    - [x] Create Supplementary Table for IFS/IFR/CHI evidence.
    - [x] Finalize `walkthrough.md` with all defensive edits.

- [x] **Phase 9: Alternative Time-Series Dataset Search** (Completed - No suitable matches)
    - [x] Search SoyMetDB and MetaboLights for soybean time-series datasets.
    - [x] Download and inspect candidate datasets (MTBLS9424 rejected).
    - [x] Feasibility assessment: Decided to proceed with MTBLS531 only.
