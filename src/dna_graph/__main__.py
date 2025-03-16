import argparse
import logging
import time
import dna_graph.bio.genetic_code as gen_code

from dna_graph.core.gauss import (
    gaussian_kernel_test,
    gaussian_kernel_test_sentence,
    cluster_results,
    get_word_test_results
)
from dna_graph.core.visualization import (
    draw_graph,
    set_positions_by_layer,
    draw_layered_sequence_graph,
    plot_clusters,
    plot_gaussian_with_histogram
)
from dna_graph.codec.codon_graph import add_codon_subgraph_bio, build_aa_to_codons
from dna_graph.codec.encode_decode import convert_message_to_bases, decode_message_from_path, extract_base_path
from dna_graph.core.optimisation import compute_on_layered_graph, compute_path_weight, dijkstra, bellman_ford, astar
from dna_graph.bio.gene_expression import simulate_gene_expression
from dna_graph.contraintes.gene_contraintes import validate_gene_expression_constraints
from dna_graph.core.init_graph import init_graph
from config.config import (
    LOG_FILE, LOG_LEVEL, LOG_FORMAT, LOG_FILE_MODE,
    ALPHA, BETA, GAMMA, DEFAULT_MESSAGE, MANDATORY_NODES, LAYER_CONFIG,
    PROMOTER, TERMINATION_SIGNAL, ADRN, DEFAULT_MUTATION_RATE, NUMB_TEST, SEED,
    NBR_BEST, NUMBER_TEST
)


def setup_logging():
    """
    Configure le logging pour enregistrer tous les messages dans un fichier unique.
    """
    # Supprimez les handlers existants pour éviter des doublons
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)
        
    logging.basicConfig(
        level=LOG_LEVEL,                   # Capture tous les niveaux de log
        format=LOG_FORMAT,                 # Format du log 
        filename=LOG_FILE,                 # Tous les logs vont dans ce fichier
        filemode=LOG_FILE_MODE             # Réécrit le fichier à chaque démarrage
    )
    
    # Réduit le niveau de log pour certains modules
    logging.getLogger("PIL.PngImagePlugin").setLevel(logging.WARNING)
    logging.getLogger("matplotlib.font_manager").setLevel(logging.WARNING)


def parse_arguments():
    """
    Parse les arguments de la ligne de commande pour GenImg.
    """
    parser = argparse.ArgumentParser(
        description="Exécute GenImg : conversion d'un message texte en séquence ADN."
    )
    parser.add_argument(
        "-m", "--message",
        type=str,
        default=DEFAULT_MESSAGE,
        help="Message à encoder. Par défaut : '%(default)s'."
    )
    parser.add_argument(
        "--alpha",
        type=float,
        default=ALPHA,
        help="Pondération alpha pour l'optimisation (coût). Par défaut : %(default)s."
    )
    parser.add_argument(
        "--beta",
        type=float,
        default=BETA,
        help="Pondération beta pour l'optimisation (stabilité). Par défaut : %(default)s."
    )
    parser.add_argument(
        "--gamma",
        type=float,
        default=GAMMA,
        help="Pondération gamma pour l'optimisation (erreur). Par défaut : %(default)s."
    )
    parser.add_argument(
        "--version",
        action="version",
        version="GenImg 1.0",
        help="Affiche la version du programme."
    )
    return parser.parse_args()


def initialize_graph():
    """
    Initialise le graphe de connaissance via le module.
    """
    logging.info("Initialisation du graphe...")
    G = init_graph()
    return G

def encode_message(message):
    """
    Convertit le message en une séquence de bases.
    """
    logging.info("Conversion du message en bases...")
    return convert_message_to_bases(message)

def add_codon_graph(G, base_list):
    """
    Ajoute le sous-graphe des codons au graphe.
    """
    logging.info("Ajout du sous-graphe codon...")
    aa_to_codons = build_aa_to_codons(gen_code.GENETIC_CODE, include_stop=True) # On inclus les stop plus fidele realité
    start, end = add_codon_subgraph_bio(G, base_list, gen_code.GENETIC_CODE, aa_to_codons)
    return start, end

def simulate_expression(G, base_list):
    """
    Simule l'expression génique à partir de la séquence de bases.
    """
    logging.info("Simulation de l'expression genique...")
    if len(base_list) % 3 != 0:
        padding_needed = 3 - (len(base_list) % 3)
        logging.info(f"Padding : Ajout de {padding_needed} base(s) pour atteindre un multiple de 3.")
        base_list += ["A"] * padding_needed
    # Construction d'une séquence ADN avec promoteur et signal de terminaison ATG pour ADRn
    dna_sequence = PROMOTER + ADRN + ''.join(base_list) + TERMINATION_SIGNAL

    constraints_status = validate_gene_expression_constraints(dna_sequence, G)
    if not constraints_status["is_valid"]:
        for key, error in constraints_status["errors"].items():
            logging.error(error)
        raise ValueError("Les contraintes d'expression genique ne sont pas respectées.")
    
    try:
        protein = simulate_gene_expression(dna_sequence)
        logging.info(f"Proteine synthetisee : {protein}")
    except ValueError as e:
        logging.error(f"Erreur lors de la simulation de l'expression génique: {e}")

def compute(G, start, end, ALPHA, BETA, GAMMA, base_list, message):
    """
    Calcule trois chemins optimisés à l'aide de Bellman-Ford, A* et Dijkstra,
    loggue chacun d'eux avec leur poids, et retourne le chemin ayant le poids minimal.
    """
    logging.info("Recherche du chemin contraint...")
    start_time = time.perf_counter()

    # Calculer les trois chemins et gérer d'éventuelles erreurs
    try:
        bf_path = bellman_ford(G, start, end, MANDATORY_NODES, ALPHA, BETA, GAMMA)
    except Exception as e:
        logging.error("Bellman-Ford a échoué: %s", e)
        bf_path = None

    try:
        astar_path = astar(G, start, end, MANDATORY_NODES, ALPHA, BETA, GAMMA, heuristic=lambda u, v: 0)
    except Exception as e:
        logging.error("A* a échoué: %s", e)
        astar_path = None

    try:
        dj_path = dijkstra(G, start, end, MANDATORY_NODES, ALPHA, BETA, GAMMA)
    except Exception as e:
        logging.error("Dijkstra a échoué: %s", e)
        dj_path = None

    end_time = time.perf_counter()
    execution_time = end_time - start_time
    print(f"Temps d'exécution : {execution_time:.4f} secondes")

    # Construire la séquence ADN complète utilisée pour calculer les poids
    dna_sequence = PROMOTER + ADRN + ''.join(base_list) + TERMINATION_SIGNAL

    # Stocker les candidats avec leur nom, chemin et poids
    candidate_paths = []
    if bf_path is not None:
        weight_bf = compute_path_weight(G, bf_path, ALPHA, BETA, GAMMA, dna_sequence)
        candidate_paths.append(('Bellman-Ford', bf_path, weight_bf))
    if astar_path is not None:
        weight_astar = compute_path_weight(G, astar_path, ALPHA, BETA, GAMMA, dna_sequence)
        candidate_paths.append(('A*', astar_path, weight_astar))
    if dj_path is not None:
        weight_dj = compute_path_weight(G, dj_path, ALPHA, BETA, GAMMA, dna_sequence)
        candidate_paths.append(('Dijkstra', dj_path, weight_dj))

    # Loguer tous les chemins candidats
    for algo_name, path, weight in candidate_paths:
        logging.info(f"{algo_name} path: {path} avec poids {weight}")

    # Sélectionner le chemin optimal (celui avec le poids minimal)
    if candidate_paths:
        best_candidate = min(candidate_paths, key=lambda x: x[2])
        best_path = best_candidate[1]
        logging.info(f"Chemin optimal choisi par {best_candidate[0]}: {best_path} avec poids {best_candidate[2]}")
    else:
        raise ValueError("Aucun chemin n'a pu être calculé.")

    base_path = extract_base_path(best_path)
    logging.info(f"Chemin simple (ATCG) : {base_path}")
    print(f"Chemin base (ATCG) : {base_path}")
    return best_path


def draw(G, best_path, base_list, message, alpha, beta, gamma):    
    """
    Dessine le graph avec le best_path pour le message
    """
    pos = set_positions_by_layer(G, LAYER_CONFIG, default_pos=(10, 0))
    logging.info("Dessin du graphe...")
    for node in best_path:
        if node not in pos:
            pos[node] = (10, 0)
    # Appel à draw_graph avec show_codon_nodes=False pour masquer les codons
    draw_graph(G, pos, best_path, show_codon_nodes=False)
    
    dna_sequence = PROMOTER + ADRN + ''.join(base_list) + TERMINATION_SIGNAL
    total_weight = compute_path_weight(G, best_path, alpha, beta, gamma, dna_sequence)
    logging.info(f"Poids total du chemin: {total_weight}")
    print(f"Poids total du chemin: {total_weight}")
    decoded_message = decode_message_from_path(best_path, original_bases=base_list, original_length=len(message))
    logging.info(f"Message decode: {decoded_message}")
    print(f"Message decode: {decoded_message}")

def gauss_kernel(message):
    """
    Test pour chaque mot de la phrase gauss kernel
    """
    all_results = gaussian_kernel_test_sentence(message, ALPHA, BETA, GAMMA, DEFAULT_MUTATION_RATE, NUMB_TEST, random_seed=SEED)
    for word, results in all_results.items():
        print(f"Résultats pour le mot '{word}':")
        for res in results:
            print(f"  Mutation rate: {res['mutation_rate']:.3f}, Score: {res['score']:.3f}")
    

def main():
    setup_logging()
    args = parse_arguments()
    logging.info("Demarrage de Genimg ...")
    
    try:
        # Initialisation du graphe
        G = initialize_graph()
    except Exception as e:
        logging.error(f"Erreur lors de l'initialisation du graphe : {e}")
        return

    try:
        # Conversion du message en bases
        base_list = encode_message(args.message)
        logging.info(f"Bases generees pour le message '{args.message}': {base_list}")
    except Exception as e:
        logging.error(f"Erreur lors de la conversion du message : {e}")
        return

    try:
        # Ajout du sous-graphe codon
        start, end = add_codon_graph(G, base_list)
    except Exception as e:
        logging.error(f"Erreur lors de l'ajout du sous-graphe codon : {e}")
        return

    try:
        # Simulation de l'expression génétique (transcription, correction et traduction)
        simulate_expression(G, base_list)
    except Exception as e:
        logging.error(f"Erreur lors de la simulation de l'expression genetique : {e}")
        return

    try:
        # Recherche du chemin optimal
        best_path = compute(G, start, end, args.alpha, args.beta, args.gamma, base_list, args.message)
    except Exception as e:
        logging.error(f"Erreur lors du calcul du chemin : {e}")
        return
    
    try:
        # Dessine le graph
        draw(G, best_path, base_list, args.message, args.alpha, args.beta, args.gamma)
    except Exception as e:
        logging.error(f"Erreur lors du dessin du graphe : {e}")
        return

    try:
        # Test avec plusieurs parametre
        gauss_kernel(args.message)
        # Affichage optionnel de la distribution gaussienne avec histogramme pour l'ensemble du message
        some_results = gaussian_kernel_test(args.message, ALPHA, BETA, GAMMA, DEFAULT_MUTATION_RATE, NUMBER_TEST, random_seed=SEED)
        plot_gaussian_with_histogram(some_results,  DEFAULT_MUTATION_RATE, sigma=0.005)

        # Générer des tests pour un mot unique et clusteriser les résultats
        test_results = gaussian_kernel_test(args.message, ALPHA, BETA, GAMMA, DEFAULT_MUTATION_RATE, NUMBER_TEST, random_seed=SEED)
        best_results = cluster_results(test_results, n_clusters=5)
        for idx, res in enumerate(best_results):
            print(f"Cluster {res['cluster']}: alpha={res['alpha']:.3f}, beta={res['beta']:.3f}, gamma={res['gamma']:.3f}, mutation_rate={res['mutation_rate']:.3f}")
            print(f"  Séquence: {res['sequence']}")
            print(f"  Score de similarité: {res['score']:.3f}")
        
        # Visualisation des clusters en 3D
        plot_clusters(test_results, dimensions=3)
        
        # On exécute à nouveau un test gaussien sur le message pour une utilisation ultérieure 
        test_results = gaussian_kernel_test(args.message, ALPHA, BETA, GAMMA, DEFAULT_MUTATION_RATE, NUMBER_TEST, random_seed=SEED)
            
        # Traitement pour le second graphe
        word_results = get_word_test_results(args.message, NUMB_TEST, NBR_BEST, ALPHA, BETA, GAMMA, DEFAULT_MUTATION_RATE, SEED)
        G2 = draw_layered_sequence_graph(args.message, word_results, NBR_BEST) 
    except Exception as e:
        logging.exception("Erreur lors des tests gaussiens : %s", e)

    try:
        best_path_G2 = compute_on_layered_graph(G2, ALPHA, BETA, GAMMA, algorithm="dijkstra")
        #  dijkstra   bellman_ford    astar    bfs    dfs
        if best_path_G2 is not None:
            logging.info(f"Chemin optimal sur le second graphe : {best_path_G2}")
            print(f"Chemin optimal sur le second graphe : {best_path_G2}")
            
        else:
            logging.error("Aucun chemin n'a pu être trouvé sur le second graphe.")
    except Exception as e:
        logging.error(f"Erreur lors du calcul du chemin sur le second graphe : {e}")

if __name__ == '__main__':
    main()
