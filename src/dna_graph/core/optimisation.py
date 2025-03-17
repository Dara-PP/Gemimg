import networkx as nx
from dna_graph.bio.gene_expression import simulate_gene_expression
import numpy as np

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

def compute_on_layered_graph(G, global_alpha, global_beta, global_gamma, algorithm,
                             heuristic=lambda u, v: 0, noise_scale=0.01, add_noise=True):
    """
    Calcule le chemin optimal sur un graphe en couches en utilisant une fonction de coût multi-critères.
    
    Paramètres :
      - G : le graphe (NetworkX)
      - global_alpha, global_beta, global_gamma : pondérations globales
      - algorithm : algorithme à utiliser ('dijkstra', 'bellman_ford', 'astar', 'bfs', 'dfs')
      - heuristic : fonction heuristique pour A* (défaut : retourne 0)
      - noise_scale : échelle du bruit aléatoire pour casser l'homogénéité (défaut : 0.01)
      - add_noise : booléen pour activer/désactiver l'ajout de bruit (défaut : True)
    
    Retourne :
      - best_path : chemin optimal trouvé ou None en cas d'erreur.
    """

    # Regrouper les nœuds par couche à partir de leur nom au format "i_..."
    layers = {}
    for node in G.nodes():
        try:
            layer_index = int(node.split('_')[0])
            layers.setdefault(layer_index, []).append(node)
        except (ValueError, IndexError):
            continue  # Ignorer les nœuds dont le nom ne suit pas le format

    if not layers:
        print("Aucune information de couche n'a pu être extraite.")
        return None

    # Définir les nœuds fictifs de départ et d'arrivée et les ajouter s'ils n'existent pas
    start = "start_fictif"
    end = "end_fictif"
    if start not in G:
        G.add_node(start, type="virtual", label="Start Fictif")
    if end not in G:
        G.add_node(end, type="virtual", label="End Fictif")

    # Relier le nœud fictif de départ aux nœuds de la première couche et de la dernière couche au nœud fictif d'arrivée
    start_layer = min(layers.keys())
    end_layer = max(layers.keys())
    for node in layers[start_layer]:
        if not G.has_edge(start, node):
            G.add_edge(start, node, weight_cost=0.01, weight_stability=1.0, weight_error=0.0)
    for node in layers[end_layer]:
        if not G.has_edge(node, end):
            G.add_edge(node, end, weight_cost=0.01, weight_stability=1.0, weight_error=0.0)

    # Définir la fonction de coût multi-critères
    def multi_criteria(u, v, data):
        # Récupérer les poids de l'arête avec des valeurs par défaut
        cost_edge = data.get("weight_cost", 1.0)
        stability_edge = data.get("weight_stability", 1.0)
        error_edge = data.get("weight_error", 0.0)
        # Récupérer les attributs du nœud de destination avec des valeurs par défaut
        node_data = G.nodes[v]
        node_alpha = node_data.get("alpha", 1.0)
        node_beta = node_data.get("beta", 1.0)
        node_gamma = node_data.get("gamma", 1.0)
        # Calcul du coût : on minimise cost et error, et on maximise la stabilité via (1 - stabilité)
        cost = (global_alpha * node_alpha) * cost_edge + \
               (global_beta * node_beta) * (1 - stability_edge) + \
               (global_gamma * node_gamma) * error_edge
        # Ajouter un bruit pour casser l'homogénéité
        noise = np.random.uniform(0, noise_scale) if add_noise else 0
        return cost + noise

    # Sélectionner l'algorithme d'optimisation en fonction du paramètre 'algorithm'
    try:
        algo = algorithm.lower()
        if algo == "dijkstra":
            best_path = nx.dijkstra_path(G, start, end, weight=lambda u, v, d: multi_criteria(u, v, d))
        elif algo == "bellman_ford":
            best_path = nx.bellman_ford_path(G, start, end, weight=lambda u, v, d: multi_criteria(u, v, d))
        elif algo == "astar":
            best_path = nx.astar_path(G, start, end, heuristic=heuristic,
                                      weight=lambda u, v, d: multi_criteria(u, v, d))
        elif algo == "bfs":
            best_path = nx.shortest_path(G, start, end)
        elif algo == "dfs":
            dfs_nodes = list(nx.dfs_tree(G, source=start).nodes())
            if end not in dfs_nodes:
                raise ValueError(f"Le nœud {end} n'est pas accessible en DFS depuis {start}.")
            best_path = nx.shortest_path(G, start, end)
        else:
            raise ValueError(f"Algorithme non supporté : {algorithm}")
        return best_path
    except Exception as e:
        print("Erreur lors de la recherche sur le second graphe :", e)
        return None
