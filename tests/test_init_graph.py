from dna_graph.core.init_graph import init_graph
from config.config import MANDATORY_NODES

def test_graph_contains_mandatory_nodes():
    G = init_graph()
    for node in MANDATORY_NODES:
        assert node in G.nodes, f"Le noeud obligatoire {node} manque dans le graphe."