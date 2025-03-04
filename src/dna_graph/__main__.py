import argparse
import logging
import dna_graph.bio.genetic_code as gen_code
from dna_graph.contraintes.gene_contraintes import validate_gene_expression_constraints
from dna_graph.core.init_graph import init_graph
from dna_graph.core.visualization import draw_graph, set_positions_by_layer
from dna_graph.codec.codon_graph import add_codon_subgraph_bio, build_aa_to_codons
from dna_graph.codec.encode_decode import convert_message_to_bases, decode_message_from_path, extract_base_path
from dna_graph.core.optimisation import compute_path_weight, dijkstra, bellman_ford, astar
from dna_graph.bio.gene_expression import simulate_gene_expression
from config.config import LOG_FILE, LOG_LEVEL,LOG_FORMAT, LOG_FILE_MODE, ALPHA, BETA, GAMMA, DEFAULT_MESSAGE, MANDATORY_NODES, LAYER_CONFIG, PROMOTER, TERMINATION_SIGNAL, ADRN

import logging

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

def compute(G, start, end, ALPHA , BETA, GAMMA, base_list, message):
    """
    Utilise les algorithmes d'optimisation pour recherche le chemin le plus optimisé
    """
    logging.info("Recherche du chemin contraint...")
    #best_path = bellman_ford(G, start, end, MANDATORY_NODES, ALPHA, BETA, GAMMA)
    #best_path = astar(G, start, end, MANDATORY_NODES, ALPHA, BETA, GAMMA, heuristic=lambda u, v: 0)
    best_path = dijkstra(G, start, end, MANDATORY_NODES, ALPHA, BETA, GAMMA)

    logging.info(f"Chemin contraint trouve: {best_path}")
    print(f"Chemin complet trouvé: {best_path}")
    base_path = extract_base_path(best_path)
    logging.info(f"Chemin simple (ATCG) : {base_path}")
    print(f"Chemin base (ATCG) : {base_path}")
    return best_path

def draw(G, best_path, base_list, message, alpha, beta, gamma):    
    """
    Dessine le graph avec le best_path pour le message
    """
    pos = set_positions_by_layer(G,LAYER_CONFIG , default_pos=(10, 0))
    logging.info("Dessin du graphe...")
    for node in best_path:
        if node not in pos:
            pos[node] = (10, 0)
    draw_graph(G, pos, best_path)
    dna_sequence = PROMOTER + ADRN + ''.join(base_list) + TERMINATION_SIGNAL
    total_weight = compute_path_weight(G, best_path, alpha, beta, gamma, dna_sequence)
    logging.info(f"Poids total du chemin: {total_weight}")
    print(f"Poids total du chemin: {total_weight}")
    decoded_message = decode_message_from_path(best_path, original_bases=base_list, original_length=len(message))
    logging.info(f"Message decode: {decoded_message}")

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
    
if __name__ == '__main__':
    main()
