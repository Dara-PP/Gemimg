import dna_graph.bio.constants as const

# Code genetic DEFAULT
GENETIC_CODE = {
    "TTT": "F", "TTC": "F",
    "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I",
    "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y",
    "TAA": "Stop", "TAG": "Stop",
    "CAT": "H", "CAC": "H",
    "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N",
    "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D",
    "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C",
    "TGA": "Stop",
    "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S",
    "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"
}

def add_bases(G):
    """Ajoute les bases nucléotidiques à G (graph), avec code binaire"""

    bases = [("A", "00"), ("C", "01"), ("T", "11"), ("G", "10")]
    for b, code in bases:
        G.add_node(b, type="base", label=b, binaire=code)

def add_classification(G):
    """Ajoute la classification des bases et relie chaque base à sa catégorie"""

    G.add_node("Purines", type="classification", label="Purines")
    G.add_node("Pyrimidines", type="classification", label="Pyrimidines")

    G.add_edge("A", "Purines", interaction="Appartient",
               weight_cost=const.DEFAULT_COST_BUILD_CLASSIFICATION,
               weight_stability=const.DEFAULT_STABILITY_BUILD_CLASSIFICATION,
               weight_error=const.DEFAULT_ERROR_BUILD_CLASSIFICATION)
    G.add_edge("G", "Purines", interaction="Appartient",
               weight_cost=const.DEFAULT_COST_BUILD_CLASSIFICATION,
               weight_stability=const.DEFAULT_STABILITY_BUILD_CLASSIFICATION,
               weight_error=const.DEFAULT_ERROR_BUILD_CLASSIFICATION)
    G.add_edge("C", "Pyrimidines", interaction="Appartient",
               weight_cost=const.DEFAULT_COST_BUILD_CLASSIFICATION,
               weight_stability=const.DEFAULT_STABILITY_BUILD_CLASSIFICATION,
               weight_error=const.DEFAULT_ERROR_BUILD_CLASSIFICATION)
    G.add_edge("T", "Pyrimidines", interaction="Appartient",
               weight_cost=const.DEFAULT_COST_BUILD_CLASSIFICATION,
               weight_stability=const.DEFAULT_STABILITY_BUILD_CLASSIFICATION,
               weight_error=const.DEFAULT_ERROR_BUILD_CLASSIFICATION)
    
    # Connecter la classification au noeud 'Gene' pour intégrer ces informations dans le parcours
    if "Gene" in G:
        G.add_edge("Purines", "Gene", interaction="Classification->Gene",
                   weight_cost=0.1, weight_stability=0.9, weight_error=0.05)
        G.add_edge("Pyrimidines", "Gene", interaction="Classification->Gene",
                   weight_cost=0.1, weight_stability=0.9, weight_error=0.05)


def add_complementarity(G):
    """Ajoute les arêtes de complémentarité avec des poids, cout de build, stabilité, erreur"""

    G.add_edge("A", "T", interaction="Complémentarité", 
               weight_cost=const.DEFAULT_COST_BUILD_AT, 
               weight_stability=const.DEFAULT_STABILITY_BUILD_AT, 
               weight_error=const.DEFAULT_ERROR_BUILD_AT)
    G.add_edge("C", "G", interaction="Complémentarité", 
               weight_cost=const.DEFAULT_COST_BUILD_CG, 
               weight_stability=const.DEFAULT_STABILITY_BUILD_CG , 
               weight_error=const.DEFAULT_ERROR_BUILD_CG)


def add_degeneracy_edges(G):
    """
    Ajoute des aretes entre les segments 3-mers codant pour le même acide aminé (hors codons STOP),
    avec un coût minimal pour modeliser la degeneracy
    """
    genetic_code = GENETIC_CODE

    segments = [n for n, d in G.nodes(data=True) if d.get("type") == "segment_3mer"]
    groups = {}
    for seg_node in segments:
        codon = seg_node.split("(")[1].split(")")[0]
        aa = genetic_code.get(codon)
        if aa and aa != "Stop":
            groups.setdefault(aa, []).append(seg_node)
    for aa, seg_list in groups.items():
        if len(seg_list) > 1:
            for i in range(len(seg_list)):
                for j in range(i+1, len(seg_list)):
                    s1, s2 = seg_list[i], seg_list[j]
                    if not G.has_edge(s1, s2):
                        G.add_edge(s1, s2, interaction="Degeneracy",
                                   weight_cost=const.DEFAULT_COST_BUILD_DEGEN, 
                                   weight_stability=const.DEFAULT_STABILITY_BUILD_DEGEN, 
                                   weight_error=const.DEFAULT_ERROR_BUILD_DEGEN)

def add_segments_3mer(G):
    """Ajoute des noeuds pour chaque segment de 3 bases (3-mers) et les relie aux bases"""
    
    bases = ["A", "C", "T", "G"]
    for b1 in bases:
        for b2 in bases:
            for b3 in bases:
                seg = f"{b1}{b2}{b3}"
                seg_node = f"Seg({seg})"
                G.add_node(seg_node, type="segment_3mer", label=f"{seg}")
                G.add_edge(b1, seg_node, interaction="Build",
                           weight_cost=const.DEFAULT_COST_BUILD_3MER_BUILD, 
                           weight_stability=const.DEFAULT_STABILITY_3MER_BUILD, 
                           weight_error=const.DEFAULT_ERROR_3MER_BUILD)
                G.add_edge(seg_node, b3, interaction="Complete",
                           weight_cost=const.DEFAULT_COST_BUILD_3MER_COM, 
                           weight_stability=const.DEFAULT_STABILITY_3MER_COM, 
                           weight_error=const.DEFAULT_ERROR_3MER_COM)

def add_motifs(G):
    """Ajoute les motifs fonctionnels et crée le noeud 'Gene' regroupant ces motifs"""

    for motif in ["Promoteur", "Exon", "Intron", "Site_de_liaison"]:
        G.add_node(motif, type="motif", label=motif)
    G.add_edge("Promoteur", "A", interaction="Contient", 
               weight_cost=const.DEFAULT_COST_BUILD_MOTIF_A, 
               weight_stability=const.DEFAULT_STABILITY_BUILD_MOTIF_A, 
               weight_error=const.DEFAULT_ERROR_BUILD_MOTIF_A)
    
    G.add_edge("Promoteur", "T", interaction="Contient", 
               weight_cost=const.DEFAULT_COST_BUILD_MOTIF_T, 
               weight_stability=const.DEFAULT_STABILITY_BUILD_MOTIF_T, 
               weight_error=const.DEFAULT_ERROR_BUILD_MOTIF_T)
    G.add_node("Gene", type="gene", label="Gène")


    for motif in ["Promoteur", "Exon", "Intron", "Site_de_liaison"]:
        G.add_edge(motif, "Gene", interaction="Fait partie de", 
                   weight_cost=const.DEFAULT_COST_BUILD_MOTIF_EXON_INTRO, 
                   weight_stability=const.DEFAULT_STABILITY_BUILD_MOTIF_EXON_INTRO,
                     weight_error=const.DEFAULT_ERROR_BUILD_MOTIF_EXON_INTRO)
        
    for b in ["A", "C", "T", "G"]:
        G.add_edge("Gene", b, interaction="Contient", 
                   weight_cost=const.DEFAULT_COST_BUILD_MOTIF_ACTG, 
                   weight_stability=const.DEFAULT_STABILITY_BUILD_MOTIF_ACTG, 
                   weight_error=const.DEFAULT_ERROR_BUILD_MOTIF_ACTG)

def add_processes(G):
    """Ajoute les processus biologiques et leurs interactions avec les motifs."""
    # Création du noeud de méthylation
    G.add_node("Methy", type="process", label="Méthylation")
    # Lien initial vers Promoteur
    G.add_edge("Methy", "Promoteur", interaction="Modifie", 
               weight_cost=const.DEFAULT_COST_BUILD_PRO_MOD, 
               weight_stability=const.DEFAULT_STABILITY_BUILD_PRO_MOD, 
               weight_error=const.DEFAULT_ERROR_BUILD_PRO_MOD)
    # Pour éviter une impasse, on ajoute plusieurs liaisons sortantes depuis Methy :
    if "Gene" in G and not G.has_edge("Methy", "Gene"):
        G.add_edge("Methy", "Gene", interaction="Transmet",
                   weight_cost=0.25, weight_stability=0.85, weight_error=0.1)
    if "Purines" in G and not G.has_edge("Methy", "Purines"):
        G.add_edge("Methy", "Purines", interaction="Transmet",
                   weight_cost=0.3, weight_stability=0.80, weight_error=0.1)
    # Création du noeud de réparation
    G.add_node("Reparation", type="process", label="Réparation de l'ADN")
    # Lien initial vers Gene
    G.add_edge("Reparation", "Gene", interaction="Corrige", 
               weight_cost=const.DEFAULT_COST_BUILD_PRO_COR, 
               weight_stability=const.DEFAULT_STABILITY_BUILD_PRO_COR, 
               weight_error=const.DEFAULT_ERROR_BUILD_PRO_COR)
    # Ajouter une liaison de Reparation vers Code_Correcteur pour plus d'options
    if "Code_Correcteur" in G and not G.has_edge("Reparation", "Code_Correcteur"):
        G.add_edge("Reparation", "Code_Correcteur", interaction="Transmet",
                   weight_cost=0.15, weight_stability=0.90, weight_error=0.05)
    # Ajouter une liaison de Reparation vers Promoteur pour diversifier le parcours
    if "Promoteur" in G and not G.has_edge("Reparation", "Promoteur"):
        G.add_edge("Reparation", "Promoteur", interaction="Corrige->Promoteur",
                   weight_cost=0.2, weight_stability=0.85, weight_error=0.1)


def add_cis_elements(G):
    """Ajoute des éléments régulateurs cis (Enhancer, Silencer) et les relie aux motifs"""
    
    G.add_node("Enhancer", type="cis", label="Enhancer")
    G.add_node("Silencer", type="cis", label="Silencer")
    if "Promoteur" in G:
        G.add_edge("Enhancer", "Promoteur", interaction="Cis-regulation",
                   weight_cost=const.DEFAULT_COST_BUILD_ELEM_ENHANCER_PRO, 
                   weight_stability=const.DEFAULT_STABILITY_BUILD_ELEM_ENHANCER_PRO, 
                   weight_error=const.DEFAULT_ERROR_BUILD_ELEM_ENHANCER_PRO)
        G.add_edge("Silencer", "Promoteur", interaction="Cis-regulation",
                   weight_cost=const.DEFAULT_COST_BUILD_ELEM_SILENCER_PRO, 
                   weight_stability=const.DEFAULT_STABILITY_BUILD_ELEM_SILENCER_PRO, 
                   weight_error=const.DEFAULT_ERROR_BUILD_ELEM_SILENCER_PRO)
    if "Gene" in G:
        G.add_edge("Enhancer", "Gene", interaction="Cis-regulation",
                   weight_cost=const.DEFAULT_COST_BUILD_ELEM_ENHANCER_GEN, 
                   weight_stability=const.DEFAULT_STABILITY_BUILD_ELEM_ENHANCER_GEN, 
                   weight_error=const.DEFAULT_ERROR_BUILD_ELEM_ENHANCER_GEN)
        G.add_edge("Silencer", "Gene", interaction="Cis-regulation",
                   weight_cost=const.DEFAULT_COST_BUILD_ELEM_SILENCER_GEN, 
                   weight_stability=const.DEFAULT_STABILITY_BUILD_ELEM_SILENCER_GEN, 
                   weight_error=const.DEFAULT_ERROR_BUILD_ELEM_SILENCER_GEN)

def add_trans_factors(G):
    """ Ajoute des facteurs régulateurs trans (TF1, TF2) et les relie au noeud Gene"""

    G.add_node("TF1", type="trans", label="TF1")
    G.add_node("TF2", type="trans", label="TF2")
    if "Gene" in G:
        G.add_edge("TF1", "Gene", interaction="Trans-regulation",
                   weight_cost=const.DEFAULT_COST_BUILD_TF1, 
                   weight_stability=const.DEFAULT_STABILITY_TF1, 
                   weight_error=const.DEFAULT_ERROR_TF1)
        G.add_edge("TF2", "Gene", interaction="Trans-regulation",
                   weight_cost=const.DEFAULT_COST_BUILD_TF2, 
                   weight_stability=const.DEFAULT_STABILITY_TF2, 
                   weight_error=const.DEFAULT_ERROR_TF2)
