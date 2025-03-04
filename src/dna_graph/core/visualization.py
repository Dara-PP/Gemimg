import networkx as nx
import matplotlib.pyplot as plt
from config.config import COLOR_MAP, NODE_SIZE, FONT_SIZE_NODE, FONT_COLOR, FONT_SIZE_INTERACTION, OPTI_PATH, XSIZE, YSIZE

def draw_graph(G, pos, path=None):
    # Fenetre Graphe
    plt.figure(figsize=(XSIZE, YSIZE))
    
    # Creer une liste de couleurs pour chaque noeud basée sur son type
    node_colors = [COLOR_MAP.get(d.get("type"), "gray") for n, d in G.nodes(data=True)]
    
    # Node
    nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=NODE_SIZE)

    # Interaction 
    nx.draw_networkx_edges(G, pos)

    # Nom des node
    labels = {n: d.get("label", n) for n, d in G.nodes(data=True)}
    nx.draw_networkx_labels(G, pos, labels, font_size=FONT_SIZE_NODE)
    
    # Interaction ligne (nom)
    edge_labels = nx.get_edge_attributes(G, "interaction")
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_color=FONT_COLOR, font_size=FONT_SIZE_INTERACTION)

    # Si un chemin est fourni, le dessiner par-dessus en le colorant différemment
    if path:
        # On extrait la liste d'arêtes correspondant au chemin
        path_edges = list(zip(path, path[1:]))
        nx.draw_networkx_edges(G, pos, edgelist=path_edges, width=3, edge_color=OPTI_PATH)

    plt.title("Graphe de connaissance ADN (configuration couche)")
    plt.axis("off")
    plt.show()
    #plt.savefig("graph.png") # pour sav le graph en png


def set_positions_by_layer(G, layer_config, default_pos=(10, 0)):
    """
    Attribue des positions aux noeuds du graphe G en colonnes.

    layer_config (main.py) : 
        Dictionnaire ou chaque clé est un type de noeud (base, motif, process etc...)
        la valeur est un tuple (x, y) x indiquant la colonne et y l'espacement vertical.

    default_pos : tuple, optionnel
        Position par défaut pour les noeuds dont le type n'est pas défini dans layer_config.
    
    Retourne :
        pos : dict
            Un dictionnaire des positions pour chaque noeud.
    """

    pos = {}
    for node_type, (x_coord, y_spacing) in layer_config.items():
        # Extraire tous les noeuds de ce type
        layer_nodes = [n for n, d in G.nodes(data=True) if d.get("type") == node_type]
        # Positionner les noeuds verticalement avec un espacement régulier
        for i, n in enumerate(sorted(layer_nodes)):
            pos[n] = (x_coord, 1 - i * y_spacing)
    
    # Pour les noeuds non assignés, leur attribuer une position par défaut
    for n in G.nodes():
        if n not in pos:
            pos[n] = default_pos
    return pos