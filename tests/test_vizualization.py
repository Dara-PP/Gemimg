import networkx as nx
from dna_graph.core.visualization import set_positions_by_layer

def test_set_positions():
    """
    Vérifie que les noeuds des types définis recoivent la bonne coordonnée x,
    et que les autres noeuds obtiennent la position par défaut default_pos.
    """
    G = nx.Graph()
    # Ajout de quelques noeuds avec différents types
    G.add_node("A", type="base")
    G.add_node("B", type="motif")
    G.add_node("C", type="unknown")
    
    layer_config = {"base": (1, 10), "motif": (2, 5)}
    positions = set_positions_by_layer(G, layer_config, default_pos=(0, 0))
    
    assert positions["A"][0] == 1
    assert positions["B"][0] == 2
    assert positions["C"] == (0, 0)
