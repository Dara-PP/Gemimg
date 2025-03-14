import networkx as nx
from config.config import PROMOTER, TERMINATION_SIGNAL, MANDATORY_NODES

def validate_restriction_sites(dna_sequence: str, restriction_sites: list = None) -> bool:
    """
    Vérifie que la séquence d'ADN ne contient pas de motifs de restriction enzymatique indésirables.
    Par défaut, on vérifie pour le motif 'GAATTC' (site de EcoRI).
    """
    if restriction_sites is None:
        restriction_sites = ["GAATTC"]  
    for site in restriction_sites:
        if site in dna_sequence:
            return False
    return True


def validate_gc_ratio(dna_sequence: str, lower_bound: float = 0.40, upper_bound: float = 0.60) -> bool:
    """
    Vérifie que le pourcentage de GC de la séquence d'ADN est compris entre lower_bound et upper_bound.
    """
    if len(dna_sequence) == 0:
        return False
    gc_count = dna_sequence.count("G") + dna_sequence.count("C")
    gc_ratio = gc_count / len(dna_sequence)
    return lower_bound <= gc_ratio <= upper_bound

def validate_promoter(dna_sequence: str) -> bool:
    """Vérifie que le promoteur est présent dans la séquence ADN."""
    return PROMOTER in dna_sequence

def validate_termination_signal(dna_sequence: str) -> bool:
    """Vérifie que le signal de terminaison est présent dans la séquence ADN."""
    return TERMINATION_SIGNAL in dna_sequence

def validate_length_for_codons(dna_sequence: str) -> bool:
    """Vérifie que la séquence (hors promoteur et terminaison) est un multiple de 3."""
    # On retire le promoteur et le signal de terminaison, si présents
    seq = dna_sequence
    if dna_sequence.startswith(PROMOTER):
        seq = dna_sequence[len(PROMOTER):]
    if dna_sequence.endswith(TERMINATION_SIGNAL):
        seq = seq[:-len(TERMINATION_SIGNAL)]
    return len(seq) % 3 == 0

def validate_classification_nodes(G: nx.Graph) -> dict:
    """
    Vérifie que les nœuds de classification 'Purines' et 'Pyrimidines' sont présents dans le graphe.
    Retourne un dictionnaire d'erreurs s'ils sont manquants.
    """
    errors = {}
    if "Purines" not in G.nodes:
        errors["purines"] = "Le noeud 'Purines' est manquant dans le graphe."
    if "Pyrimidines" not in G.nodes:
        errors["pyrimidines"] = "Le noeud 'Pyrimidines' est manquant dans le graphe."
    return errors

def validate_complementarity_edges(G: nx.Graph) -> dict:
    """
    Vérifie que les arêtes de complémentarité entre A-T et C-G existent dans le graphe.
    Retourne un dictionnaire d'erreurs si l'une de ces arêtes est manquante.
    """
    errors = {}
    # Vérifier la présence des nœuds de base
    for base in ["A", "T", "C", "G"]:
        if base not in G.nodes:
            errors[f"no_{base}"] = f"Le nœud '{base}' est manquant dans le graphe."
    # Vérifier la complémentarité entre A et T
    if not G.has_edge("A", "T"):
        errors["edge_AT"] = "L'arête de complémentarité entre A et T est manquante."
    # Vérifier la complémentarité entre C et G
    if not G.has_edge("C", "G"):
        errors["edge_CG"] = "L'arête de complémentarité entre C et G est manquante."
    return errors


def validate_mandatory_nodes(G: nx.Graph) -> bool:
    """Vérifie que tous les nœuds obligatoires sont présents dans le graphe."""
    return all(node in G.nodes for node in MANDATORY_NODES)

def validate_gene_expression_constraints(dna_sequence, G) -> dict:
    """
    Valide l'ensemble des contraintes pour l'expression génique.
    Renvoie un dictionnaire avec le statut et les messages d'erreur.
    """
    errors = {}
    
    if not validate_promoter(dna_sequence):
        errors["promoteur"] = "Le promoteur n'est pas present dans la sequence ADN."
    if not validate_termination_signal(dna_sequence):
        errors["termination"] = "Le signal de terminaison n'est pas present dans la sequence ADN."
    if not validate_length_for_codons(dna_sequence):
        errors["longueur"] = "La sequence entre le promoteur et le signal de terminaison n'est pas un multiple de 3."
    if not validate_mandatory_nodes(G):
        errors["noeuds"] = "Tous les noeuds obligatoires ne sont pas presents dans le graphe."

    errors.update(validate_classification_nodes(G))
    errors.update(validate_complementarity_edges(G))

        # Nouvelle contrainte : éviter les sites de restriction enzymatique
    if not validate_restriction_sites(dna_sequence):
        errors["restriction_site"] = "La séquence ADN contient des sites de restriction enzymatique indésirables."
    
    # Nouvelle contrainte : vérifier le ratio GC (40-60%)
    if not validate_gc_ratio(dna_sequence):
        errors["gc_ratio"] = "La séquence ADN ne respecte pas le ratio GC requis (40-60%)."

    
    return {"is_valid": len(errors) == 0, "errors": errors}

