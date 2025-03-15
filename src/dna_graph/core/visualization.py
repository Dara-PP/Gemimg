import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from sklearn.decomposition import PCA
from config.config import COLOR_MAP, NODE_SIZE, FONT_SIZE_NODE, FONT_COLOR, FONT_SIZE_INTERACTION, OPTI_PATH, XSIZE, YSIZE, OPTI_PSIZE

def draw_graph(G, pos, path=None, show_codon_nodes=False):
    plt.figure(figsize=(XSIZE, YSIZE))
    
    # Si on ne veut pas afficher tous les codons, on conserve néanmoins ceux qui font partie du chemin optimal.
    if not show_codon_nodes:
        codon_nodes_in_path = set()
        if path:
            codon_nodes_in_path = {n for n in path if G.nodes[n].get("type") == "segment_3mer"}
        nodes_to_draw = [n for n, d in G.nodes(data=True) 
                         if d.get("type") != "segment_3mer" or n in codon_nodes_in_path]
    else:
        nodes_to_draw = list(G.nodes())
    
    # Couleurs pour les noeuds à afficher
    node_colors = [COLOR_MAP.get(G.nodes[n].get("type"), "gray") for n in nodes_to_draw]
    nx.draw_networkx_nodes(G, pos, nodelist=nodes_to_draw, node_color=node_colors, node_size=NODE_SIZE)
    
    # Filtrer les arêtes pour n'afficher que celles reliant des noeuds affichés
    edges_to_draw = [(u, v) for (u, v, d) in G.edges(data=True)
                     if d.get("display", True) and u in nodes_to_draw and v in nodes_to_draw]
    nx.draw_networkx_edges(G, pos, edgelist=edges_to_draw)
    
    # Labels de noeuds pour ceux affichés
    labels = {n: d.get("label", n) for n, d in G.nodes(data=True) if n in nodes_to_draw}
    nx.draw_networkx_labels(G, pos, labels, font_size=FONT_SIZE_NODE)
    
    # Labels d'arêtes filtrés
    edge_labels = nx.get_edge_attributes(G, "interaction")
    filtered_edge_labels = {(u, v): label for (u, v), label in edge_labels.items()
                            if u in nodes_to_draw and v in nodes_to_draw}
    nx.draw_networkx_edge_labels(G, pos, edge_labels=filtered_edge_labels, font_color=FONT_COLOR, font_size=FONT_SIZE_INTERACTION)
    
    # Si un chemin optimal est passé, on le dessine (en filtrant également)
    if path:
        path_edges = list(zip(path, path[1:]))
        filtered_path_edges = [(u, v) for u, v in path_edges if u in nodes_to_draw and v in nodes_to_draw]
        nx.draw_networkx_edges(G, pos, edgelist=filtered_path_edges, width=OPTI_PSIZE, edge_color=OPTI_PATH)
    
    plt.title("Graphe de connaissance ADN (configuration par couche)")
    plt.axis("off")
    plt.show()

    #plt.savefig("graph.png") # pour sav le graph en png


def set_positions_by_layer(G, layer_config, default_pos=(10, 0), random_seed=None):
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
    for node_type, (x_coord, base_y, y_spacing) in layer_config.items():
        # Sélectionne les noeuds de ce type
        layer_nodes = [n for n, d in G.nodes(data=True) if d.get("type") == node_type]
        # Trie éventuellement par nom pour un placement cohérent
        for i, n in enumerate(sorted(layer_nodes)):
            # On place le 1er noeud à y = base_y, 
            # le 2e à y = base_y - y_spacing, etc.
            y = base_y - i * y_spacing
            pos[n] = (x_coord, y)

    # Position par défaut pour les noeuds non assignés
    for n in G.nodes():
        if n not in pos:
            pos[n] = default_pos

    return pos

def plot_clusters(test_results, dimensions=2):
    """
    Visualise les clusters sur une projection en 2D ou 3D.
    
    Paramètres :
      - test_results (list) : Liste de dictionnaires contenant 'alpha', 'beta', 'gamma', 'mutation_rate' et 'cluster'.
      - dimensions (int) : Dimension de la projection (2 ou 3).
    """
    features = np.array([[res['alpha'], res['beta'], res['gamma'], res['mutation_rate']] 
                         for res in test_results])
    clusters = np.array([res.get('cluster', -1) for res in test_results])
    
    pca = PCA(n_components=dimensions)
    projected = pca.fit_transform(features)
    
    if dimensions == 2:
        plt.figure(figsize=(8,6))
        scatter = plt.scatter(projected[:, 0], projected[:, 1], c=clusters, cmap='viridis', marker='o')
        plt.xlabel("PCA 1")
        plt.ylabel("PCA 2")
        plt.title("Projection des clusters en 2D")
        plt.colorbar(scatter, label="Cluster")
        plt.show()
    elif dimensions == 3:
        fig = plt.figure(figsize=(8,6))
        ax = fig.add_subplot(111, projection='3d')
        scatter = ax.scatter(projected[:, 0], projected[:, 1], projected[:, 2], c=clusters, cmap='viridis', marker='o')
        ax.set_xlabel("PCA 1")
        ax.set_ylabel("PCA 2")
        ax.set_zlabel("PCA 3")
        plt.title("Projection des clusters en 3D")
        fig.colorbar(scatter, ax=ax, label="Cluster")
        plt.show()

def plot_gaussian_distribution(default_value, sigma, num_points=1000):
    """
    Affiche la courbe de la distribution gaussienne pour le taux de mutation.
    
    Paramètres :
      - default_value (float) : La moyenne (mu).
      - sigma (float) : L'écart-type.
      - num_points (int) : Nombre de points pour tracer la courbe.
    """
    # Définir l'intervalle de x en respectant [0, 1]
    x = np.linspace(max(0, default_value - 4 * sigma), min(1, default_value + 4 * sigma), num_points)
    y = (1 / (sigma * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((x - default_value) / sigma) ** 2)
    
    plt.figure(figsize=(8,6))
    plt.plot(x, y, label=f"mu={default_value}, sigma={sigma}")
    plt.xlabel("Taux de mutation")
    plt.ylabel("Densité de probabilité")
    plt.title("Distribution Gaussienne pour le taux de mutation")
    plt.legend()
    plt.show()

def plot_gaussian_with_histogram(test_results, default_value, sigma, num_points=1000, bins=30):
    """
    Affiche la distribution gaussienne et superpose un histogramme des taux de mutation obtenus.
    
    Paramètres :
      - test_results (list) : Liste de dictionnaires issus de gaussian_kernel_test.
      - default_value (float) : Valeur moyenne (mu) pour la courbe.
      - sigma (float) : Écart-type de la distribution.
      - num_points (int) : Nombre de points pour tracer la courbe.
      - bins (int) : Nombre de bins pour l'histogramme.
    """
    rates = [res['mutation_rate'] for res in test_results]
    x = np.linspace(max(0, default_value - 4 * sigma), min(1, default_value + 4 * sigma), num_points)
    y = (1 / (sigma * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((x - default_value) / sigma) ** 2)
    
    plt.figure(figsize=(8,6))
    plt.plot(x, y, label=f"Distribution Gaussienne (mu={default_value}, sigma={sigma})", color="blue")
    plt.hist(rates, bins=bins, density=True, alpha=0.5, label="Histogramme des taux", color="orange")
    plt.xlabel("Taux de mutation")
    plt.ylabel("Densité")
    plt.title("Courbe gaussienne et histogramme")
    plt.legend()
    plt.show()

def draw_layered_sequence_graph(sentence, word_results, n_best=10):
    # Séparer la phrase en mots
    words = sentence.split()
    
    # Créer un graphe orienté pour représenter le chemin
    G = nx.DiGraph()
    layers = []  # Liste de listes de noeuds par couche
    
    for i, word in enumerate(words):
        results = word_results.get(word, [])
        if not results:
            node_id = f"{i}_0"
            G.add_node(node_id, word=word, sequence="—", score=0, alpha=None, beta=None, gamma=None, mutation_rate=None)
            layers.append([node_id])
            continue

        sorted_results = sorted(results, key=lambda r: r.get('score', 0), reverse=True)
        best_candidates = sorted_results[:n_best]
        current_layer = []
        for j, res in enumerate(best_candidates):
            node_id = f"{i}_{j}"
            G.add_node(node_id,
                       word=word,
                       sequence=res.get('sequence', ''),
                       score=res.get('score', 0),
                       alpha=res.get('alpha'),
                       beta=res.get('beta'),
                       gamma=res.get('gamma'),
                       mutation_rate=res.get('mutation_rate'))
            current_layer.append(node_id)
        layers.append(current_layer)
    
    # Connecter chaque couche à la suivante
    for i in range(len(layers) - 1):
        for node_u in layers[i]:
            for node_v in layers[i+1]:
                G.add_edge(node_u, node_v)
    
    # Générer des positions et attribuer une couleur différente par couche
    pos = {}
    layer_spacing = 4
    vertical_spacing = 2
    layer_colors = plt.cm.viridis(np.linspace(0, 1, len(layers)))
    node_colors = {}
    
    for i, layer_nodes in enumerate(layers):
        x = i * layer_spacing
        k = len(layer_nodes)
        ys = np.linspace((k-1)*vertical_spacing/2, -(k-1)*vertical_spacing/2, k) if k > 0 else [0]
        for node, y in zip(layer_nodes, ys):
            pos[node] = (x, y)
            node_colors[node] = layer_colors[i]
    
    plt.figure(figsize=(12, 8))
    nx.draw_networkx_nodes(G, pos, node_color=[node_colors[n] for n in G.nodes()], node_size=1500)
    nx.draw_networkx_edges(G, pos, arrowstyle='->', arrowsize=20)
    node_labels = {}
    for node, data in G.nodes(data=True):
        label = (f"{data.get('word')}\nSeq: {data.get('sequence')}\nScore: {data.get('score'):.2f}\n"
                 f"(α={data.get('alpha'):.2f}, β={data.get('beta'):.2f}, γ={data.get('gamma'):.2f})\n"
                 f"Mut: {data.get('mutation_rate'):.2f}")
        node_labels[node] = label
    nx.draw_networkx_labels(G, pos, labels=node_labels, font_size=8)
    plt.title("Graphe en couches des meilleures séquences par mot")
    plt.axis('off')
    plt.show()
    
    return G  # Retourner le graphe pour l'utiliser dans la recherche
