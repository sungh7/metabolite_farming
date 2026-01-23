"""
Pathway Mapper Module

Handles KEGG pathway-gene-protein-metabolite mappings for pathway activation analysis.
"""

import json
import pandas as pd
from pathlib import Path
from typing import List, Dict, Set, Optional


# Known compound-pathway mappings for soybean isoflavonoid biosynthesis
# Based on KEGG pathway definitions and actual KEGG compound IDs
COMPOUND_PATHWAY_MAPPING = {
    # Isoflavonoid biosynthesis (gmx00943)
    'C00814': ['gmx00943'],  # Genistein
    'C00858': ['gmx00943'],  # Formononetin (detected in dataset)
    'C05376': ['gmx00943'],  # Biochanin A
    'C02495': ['gmx00943'],  # Daidzein (detected in dataset)
    'C10216': ['gmx00943'],  # Daidzin (detected in dataset)
    'C10205': ['gmx00943'],  # Daidzein (alternative ID)
    'C04895': ['gmx00943'],  # Glycitein
    'C00786': ['gmx00943'],  # 2'-Hydroxyisoflavanone
    'C00891': ['gmx00943'],  # Daidzin (alternative ID)
    'C00899': ['gmx00943'],  # Genistin
    'C10510': ['gmx00943'],  # Genistein 7-O-glucoside
    'C16191': ['gmx00943'],  # 2'-Hydroxyformononetin
    'C16195': ['gmx00943'],  # 2'-Hydroxygenistein
    'C05623': ['gmx00943', 'gmx00941'],  # Liquiritigenin (shared)

    # Flavonoid biosynthesis (gmx00941)
    'C00509': ['gmx00941'],  # Naringenin
    'C01477': ['gmx00941'],  # Apigenin
    'C00389': ['gmx00941'],  # Quercetin
    'C01514': ['gmx00941'],  # Luteolin
    'C00540': ['gmx00941'],  # Catechin
    'C05631': ['gmx00941'],  # Eriodictyol
    'C00774': ['gmx00941'],  # Phloretin
    'C09320': ['gmx00941'],  # Kaempferol
    'C01604': ['gmx00941'],  # Phlorizin
    'C09762': ['gmx00941'],  # Naringenin 7-O-glucoside
    'C10107': ['gmx00941'],  # Naringenin chalcone

    # Phenylpropanoid biosynthesis (gmx00940)
    'C00079': ['gmx00940'],  # Phenylalanine
    'C00423': ['gmx00940'],  # trans-Cinnamic acid
    'C00156': ['gmx00940'],  # 4-Coumaric acid
    'C00323': ['gmx00940'],  # Caffeic acid
    'C00482': ['gmx00940'],  # Sinapic acid
    'C01772': ['gmx00940'],  # trans-Ferulic acid
    'C00406': ['gmx00940'],  # Feruloyl-CoA
    'C00223': ['gmx00940'],  # p-Coumaroyl-CoA
    'C00811': ['gmx00940'],  # 4-Coumarate
    'C12203': ['gmx00940'],  # Coniferyl alcohol
    'C12204': ['gmx00940'],  # Sinapyl alcohol

    # Plant hormone signal transduction (gmx04075)
    'C00805': ['gmx04075'],  # Salicylic acid
    'C00220': ['gmx04075'],  # Auxin (IAA)
    'C00179': ['gmx04075'],  # Zeatin
    'C00147': ['gmx04075'],  # Adenine
    'C01592': ['gmx04075'],  # ACC (ethylene precursor)
    'C06082': ['gmx04075'],  # Jasmonic acid
    'C08468': ['gmx04075'],  # Auxin-related (detected in dataset)
}

# Additional compounds found in isoflavonoid pathway
ISOFLAVONOID_COMPOUNDS = {
    'C00814', 'C00858', 'C05376', 'C10205', 'C04895',
    'C00786', 'C00891', 'C00899', 'C05623'
}


class PathwayMapper:
    """
    Maps KEGG pathways to genes, proteins, and metabolites.
    """

    def __init__(
        self,
        kegg_pathway_path: str,
        kegg_uniprot_mapping_path: str,
        enzyme_string_mapping_path: str
    ):
        """
        Initialize PathwayMapper with data file paths.

        Args:
            kegg_pathway_path: Path to kegg_pathway_data.json
            kegg_uniprot_mapping_path: Path to kegg_uniprot_mapping.csv
            enzyme_string_mapping_path: Path to enzyme_string_mapping.csv
        """
        self.kegg_pathway_path = Path(kegg_pathway_path)
        self.kegg_uniprot_mapping_path = Path(kegg_uniprot_mapping_path)
        self.enzyme_string_mapping_path = Path(enzyme_string_mapping_path)

        # Load data
        self._load_pathway_data()
        self._load_id_mappings()
        self._build_reverse_mappings()

    def _load_pathway_data(self):
        """Load KEGG pathway definitions."""
        with open(self.kegg_pathway_path, 'r') as f:
            data = json.load(f)

        self.organism = data.get('organism', 'gmx')
        self.pathways = data.get('pathways', {})
        self.gene_families = data.get('genes', {})
        self.compounds = data.get('compounds', {})

        # Build pathway -> genes mapping
        self.pathway_genes = {}
        for pathway_id, pathway_info in self.pathways.items():
            self.pathway_genes[pathway_id] = set(pathway_info.get('genes', []))

    def _load_id_mappings(self):
        """Load ID mapping files."""
        # KEGG gene -> UniProt
        self.kegg_to_uniprot = {}
        self.uniprot_to_kegg = {}

        if self.kegg_uniprot_mapping_path.exists():
            df = pd.read_csv(self.kegg_uniprot_mapping_path)
            for _, row in df.iterrows():
                kegg_gene = str(row['kegg_gene'])
                uniprot = str(row['uniprot'])
                self.kegg_to_uniprot[kegg_gene] = uniprot
                self.uniprot_to_kegg[uniprot] = kegg_gene

        # STRING ID -> UniProt
        self.string_to_uniprot = {}
        self.uniprot_to_string = {}

        if self.enzyme_string_mapping_path.exists():
            df = pd.read_csv(self.enzyme_string_mapping_path)
            for _, row in df.iterrows():
                string_id = str(row['string_id'])
                uniprot = str(row['uniprot_id'])
                self.string_to_uniprot[string_id] = uniprot
                self.uniprot_to_string[uniprot] = string_id

    def _build_reverse_mappings(self):
        """Build reverse mappings for efficient lookup."""
        # Pathway -> UniProt IDs
        self.pathway_uniprots = {}
        for pathway_id, genes in self.pathway_genes.items():
            uniprots = set()
            for gene in genes:
                if gene in self.kegg_to_uniprot:
                    uniprots.add(self.kegg_to_uniprot[gene])
            self.pathway_uniprots[pathway_id] = uniprots

        # UniProt -> Pathways
        self.uniprot_to_pathways = {}
        for pathway_id, uniprots in self.pathway_uniprots.items():
            for uniprot in uniprots:
                if uniprot not in self.uniprot_to_pathways:
                    self.uniprot_to_pathways[uniprot] = set()
                self.uniprot_to_pathways[uniprot].add(pathway_id)

        # Compound -> Pathways (from predefined mapping)
        self.compound_to_pathways = COMPOUND_PATHWAY_MAPPING.copy()

        # Pathway -> Compounds
        self.pathway_compounds = {}
        for compound_id, pathways in self.compound_to_pathways.items():
            for pathway_id in pathways:
                if pathway_id not in self.pathway_compounds:
                    self.pathway_compounds[pathway_id] = set()
                self.pathway_compounds[pathway_id].add(compound_id)

    def get_pathway_proteins(self, pathway_id: str) -> List[str]:
        """
        Get UniProt IDs of proteins in a pathway.

        Args:
            pathway_id: KEGG pathway ID (e.g., 'gmx00943')

        Returns:
            List of UniProt IDs
        """
        return list(self.pathway_uniprots.get(pathway_id, set()))

    def get_pathway_string_ids(self, pathway_id: str) -> List[str]:
        """
        Get STRING IDs of proteins in a pathway.

        Args:
            pathway_id: KEGG pathway ID

        Returns:
            List of STRING IDs
        """
        uniprots = self.pathway_uniprots.get(pathway_id, set())
        string_ids = []
        for uniprot in uniprots:
            if uniprot in self.uniprot_to_string:
                string_ids.append(self.uniprot_to_string[uniprot])
        return string_ids

    def get_pathway_metabolites(self, pathway_id: str) -> List[str]:
        """
        Get KEGG compound IDs of metabolites in a pathway.

        Args:
            pathway_id: KEGG pathway ID

        Returns:
            List of KEGG compound IDs
        """
        return list(self.pathway_compounds.get(pathway_id, set()))

    def get_pathway_genes(self, pathway_id: str) -> List[str]:
        """
        Get KEGG gene IDs in a pathway.

        Args:
            pathway_id: KEGG pathway ID

        Returns:
            List of KEGG gene IDs
        """
        return list(self.pathway_genes.get(pathway_id, set()))

    def get_protein_pathways(self, uniprot_id: str) -> List[str]:
        """
        Get pathways containing a protein.

        Args:
            uniprot_id: UniProt ID

        Returns:
            List of pathway IDs
        """
        return list(self.uniprot_to_pathways.get(uniprot_id, set()))

    def get_metabolite_pathways(self, compound_id: str) -> List[str]:
        """
        Get pathways containing a metabolite.

        Args:
            compound_id: KEGG compound ID

        Returns:
            List of pathway IDs
        """
        return list(self.compound_to_pathways.get(compound_id, []))

    def get_pathway_info(self, pathway_id: str) -> Dict:
        """
        Get pathway information.

        Args:
            pathway_id: KEGG pathway ID

        Returns:
            Dictionary with pathway info (name, gene_count, etc.)
        """
        return self.pathways.get(pathway_id, {})

    def get_coverage(
        self,
        pathway_id: str,
        detected_proteins: Set[str],
        detected_metabolites: Set[str]
    ) -> Dict[str, float]:
        """
        Calculate pathway coverage based on detected molecules.

        Args:
            pathway_id: KEGG pathway ID
            detected_proteins: Set of detected UniProt IDs
            detected_metabolites: Set of detected KEGG compound IDs

        Returns:
            Dictionary with protein_coverage, metabolite_coverage, overall_coverage
        """
        pathway_proteins = self.pathway_uniprots.get(pathway_id, set())
        pathway_metabolites = self.pathway_compounds.get(pathway_id, set())

        # Protein coverage
        if len(pathway_proteins) > 0:
            protein_overlap = len(detected_proteins & pathway_proteins)
            protein_coverage = protein_overlap / len(pathway_proteins)
        else:
            protein_coverage = 0.0

        # Metabolite coverage
        if len(pathway_metabolites) > 0:
            metabolite_overlap = len(detected_metabolites & pathway_metabolites)
            metabolite_coverage = metabolite_overlap / len(pathway_metabolites)
        else:
            metabolite_coverage = 0.0

        # Overall coverage (geometric mean)
        if protein_coverage > 0 and metabolite_coverage > 0:
            overall_coverage = (protein_coverage * metabolite_coverage) ** 0.5
        else:
            overall_coverage = max(protein_coverage, metabolite_coverage) * 0.5

        return {
            'protein_coverage': protein_coverage,
            'metabolite_coverage': metabolite_coverage,
            'overall_coverage': overall_coverage,
            'n_proteins_detected': len(detected_proteins & pathway_proteins),
            'n_proteins_total': len(pathway_proteins),
            'n_metabolites_detected': len(detected_metabolites & pathway_metabolites),
            'n_metabolites_total': len(pathway_metabolites)
        }

    def string_to_uniprot_id(self, string_id: str) -> Optional[str]:
        """Convert STRING ID to UniProt ID."""
        return self.string_to_uniprot.get(string_id)

    def uniprot_to_kegg_gene(self, uniprot_id: str) -> Optional[str]:
        """Convert UniProt ID to KEGG gene ID."""
        return self.uniprot_to_kegg.get(uniprot_id)

    def get_all_pathways(self) -> List[str]:
        """Get all pathway IDs."""
        return list(self.pathways.keys())

    def get_pathway_name(self, pathway_id: str) -> str:
        """Get pathway name."""
        info = self.pathways.get(pathway_id, {})
        return info.get('name', pathway_id)

    def summary(self) -> Dict:
        """Get summary statistics of the mapper."""
        return {
            'n_pathways': len(self.pathways),
            'n_kegg_genes': len(self.kegg_to_uniprot),
            'n_string_proteins': len(self.string_to_uniprot),
            'n_compound_mappings': len(self.compound_to_pathways),
            'pathways': {
                pid: {
                    'name': info.get('name', ''),
                    'n_genes': info.get('gene_count', 0),
                    'n_uniprots_mapped': len(self.pathway_uniprots.get(pid, set())),
                    'n_compounds_mapped': len(self.pathway_compounds.get(pid, set()))
                }
                for pid, info in self.pathways.items()
            }
        }
