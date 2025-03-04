import networkx as nx
from dna_graph.bio.gene_expression import simulate_gene_expression

def multi_criteria_weight(u, v, data, alpha, beta, gamma):
    """
    Calcule un poids unique à partir de : cost, stability, error.
    Minimiser = alpha*cost + beta*(1 - stability) + gamma*error.
    """
    c = data.get("weight_cost", 0.0)
    s = data.get("weight_stability", 0.0)
    e = data.get("weight_error", 0.0)
    return alpha * c + beta * (1 - s) + gamma * e


def compute_path_weight(G, path, alpha, beta, gamma, dna_sequence):
    """
    Calcule la somme des poids multi-critères sur les arêtes d'un chemin donné.
    En cas d'échec de la simulation de l'expression génétique, une pénalité est ajoutée.
    """
    total = 0.0
    for i in range(len(path) - 1):
        edge_data = G[path[i]][path[i + 1]]
        total += alpha * edge_data.get("weight_cost", 0.0)
        total += beta * (1 - edge_data.get("weight_stability", 0.0))
        total += gamma * edge_data.get("weight_error", 0.0)
    try:
        simulate_gene_expression(dna_sequence)
    except ValueError:
        total += 1  # Pénalité en cas d'échec
    return total

# ---- Bellman Ford ---- #
def bellman_ford(G, start_node, end_node, mandatory_nodes, alpha, beta, gamma):
    """
    Calcule un chemin passant par tous les noeuds obligatoires dans l'ordre indiqué en utilisant Bellman-Ford.
    
    Paramètres :
      - G : Le graphe.
      - start_node : Le noeud de départ.
      - end_node : Le noeud d'arrivée.
      - mandatory_nodes : Liste ordonnée des noeuds obligatoires.
      - alpha, beta, gamma : Pondérations pour la fonction de coût.
    
    Retourne :
      - full_path : Liste de noeuds formant le chemin complet.
    """
    full_path = []
    current_node = start_node

    for mandatory in mandatory_nodes:
        segment = nx.bellman_ford_path(
            G,
            source=current_node,
            target=mandatory,
            weight=lambda u, v, d: multi_criteria_weight(u, v, d, alpha, beta, gamma)
        )
        if full_path:
            full_path.extend(segment[1:])
        else:
            full_path.extend(segment)
        current_node = mandatory

    segment = nx.bellman_ford_path(
        G,
        source=current_node,
        target=end_node,
        weight=lambda u, v, d: multi_criteria_weight(u, v, d, alpha, beta, gamma)
    )
    full_path.extend(segment[1:])
    return full_path

# ---- Djikstra ---- #
def dijkstra(G, start_node, end_node, mandatory_nodes, alpha, beta, gamma):
    """
    Calcule un chemin passant par tous les noeuds obligatoires dans l'ordre indiqué.
    Pour les segments vers les noeuds obligatoires, Bellman-Ford est utilisé, et Dijkstra pour le segment final.
    
    Paramètres :
      - G : Le graphe.
      - start_node : Le noeud de départ.
      - end_node : Le noeud d'arrivée.
      - mandatory_nodes : Liste ordonnée des noeuds obligatoires.
      - alpha, beta, gamma : Pondérations pour la fonction de coût.
    
    Retourne :
      - full_path : Liste de noeuds formant le chemin complet.
    """
    full_path = []
    current_node = start_node

    for mandatory in mandatory_nodes:
        segment = nx.bellman_ford_path(
            G,
            source=current_node,
            target=mandatory,
            weight=lambda u, v, d: multi_criteria_weight(u, v, d, alpha, beta, gamma)
        )
        if full_path:
            full_path.extend(segment[1:])
        else:
            full_path.extend(segment)
        current_node = mandatory

    segment = nx.dijkstra_path(
        G,
        source=current_node,
        target=end_node,
        weight=lambda u, v, d: multi_criteria_weight(u, v, d, alpha, beta, gamma)
    )
    full_path.extend(segment[1:])
    return full_path


# ---- Algorithme A* ---- #

def astar(G, start_node, end_node, mandatory_nodes, alpha, beta, gamma, heuristic=lambda u, v: 0):
    """
    Calcule un chemin passant par tous les noeuds obligatoires en utilisant l'algorithme A*.
    Le paramètre 'heuristic' permet de fournir une fonction heuristique.
    Si aucune heuristique n'est fournie, une heuristique nulle est utilisée (équivalent à Dijkstra).
    """
    full_path = []
    current_node = start_node

    for mandatory in mandatory_nodes:
        segment = nx.astar_path(
            G,
            source=current_node,
            target=mandatory,
            heuristic=heuristic,
            weight=lambda u, v, d: multi_criteria_weight(u, v, d, alpha, beta, gamma)
        )
        if full_path:
            full_path.extend(segment[1:])
        else:
            full_path.extend(segment)
        current_node = mandatory

    segment = nx.astar_path(
        G,
        source=current_node,
        target=end_node,
        heuristic=heuristic,
        weight=lambda u, v, d: multi_criteria_weight(u, v, d, alpha, beta, gamma)
    )
    full_path.extend(segment[1:])
    return full_path
