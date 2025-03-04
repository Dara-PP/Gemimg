import networkx as nx
from dna_graph.core.graph_layers import (
    add_bases_layer,
    add_motif_layer,
    add_segment_layer,
    add_regulation_layer,
)
def init_graph():
    """
    Initialise et assemble le graphe complet en appelant
    les fonctions par couche.
    Assemble le graphe de connaissance pour l'ADN en 6 couches + Noeuds additionnels :
    - Couche 1 : Bases nucléotidiques
    - Couche 2 : Classification des bases
    - Couche 3 : Complémentarité
    - Couche 4 : Motifs fonctionnels et regroupement en "Gene"
    - Couche 5 : Processus physiologiques
    - Couche 6 : Code correcteur
    - Couche 7 : Segments 3-mers & clusters
    - Couche 8 : Enhancer & Silencer
    - Couche 9 : Trans
    - Edge : Degeneracy
    """
    G = nx.Graph()
    
    add_bases_layer(G)
    add_motif_layer(G)
    add_segment_layer(G)
    add_regulation_layer(G)

    return G
