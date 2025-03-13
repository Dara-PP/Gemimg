from dna_graph.bio.error_correction import add_code_correcteur
from dna_graph.bio.genetic_code import (
    add_bases,
    add_classification,
    add_complementarity,
    add_segments_3mer,
    add_motifs,
    add_processes,
    add_cis_elements,
    add_trans_factors,
    add_degeneracy_edges,
    add_epigenetics_layer,
    add_splicing_layer
)

def add_bases_layer(G):
    """Configure la couche des bases nucléotidiques et leur complémentarité."""
    add_bases(G)
    add_classification(G)
    add_complementarity(G)

def add_motif_layer(G):
    """Configure la couche des motifs et des processus biologiques."""
    add_motifs(G)
    add_processes(G)
    add_code_correcteur(G)

def add_segment_layer(G):
    """Configure la couche des segments (3-mers) et la connexion via la dégénérescence."""
    add_segments_3mer(G)
    add_degeneracy_edges(G)

def add_regulation_layer(G):
    """Configure la couche de régulation avec les éléments cis et trans."""
    add_cis_elements(G)
    add_trans_factors(G)

def add_additional_layers(G): # temp peut etre test
    add_splicing_layer(G)
    add_epigenetics_layer(G)
