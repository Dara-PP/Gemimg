def build_aa_to_codons(genetic_code, include_stop=False):
    """
    Construit un dictionnaire associant chaque acide aminé à une liste de codons.
    
    Paramètres :
      - genetic_code : dictionnaire codon -> acide aminé.
      - include_stop : si False, on ignore les codons de stop.
      
    Retourne :
      - aa_to_codons : dictionnaire de la forme {acide_aminé: [codon1, codon2, ...]}.
    """
    aa_to_codons = {}
    for codon, aa in genetic_code.items():
        if not include_stop and aa == "Stop":
            continue
        aa_to_codons.setdefault(aa, []).append(codon)
    return aa_to_codons

def add_codon_subgraph_bio(G, base_list, genetic_code, aa_to_codons):
    """
    Ajoute un sous-graphe codon au graphe G en utilisant la séquence de bases.
    
    La séquence est découpée en groupes de 3 bases (codons). Pour chaque codon,
    un nœud est créé avec un label identifiant le codon.
    
    En plus, cette fonction crée explicitement deux nœuds "start" et "end" de type "virtual"
    qui servent de points de raccord pour le chemin retour vers le stop de l'ADN.
    
    Le chemin final est constitué de la suite : [start] + liste des nœuds codon + [end].
    
    Paramètres :
      - G : graphe NetworkX.
      - base_list : liste des bases générées.
      - genetic_code : dictionnaire codon -> acide aminé.
      - aa_to_codons : dictionnaire acide aminé -> liste de codons.
      
    Retourne :
      - (start, end) : tuple contenant le nœud "start" et le nœud "end".
    """
    # Création explicite des nœuds "start" et "end" avec type "virtual"
    start = "start"
    end = "end"
    G.add_node(start, type="virtual", label="Start")
    G.add_node(end, type="virtual", label="End")
    
    # Découper la séquence en codons (groupes de 3 bases)
    codon_nodes = []
    for i in range(0, len(base_list), 3):
        codon = ''.join(base_list[i:i+3])
        if len(codon) < 3:
            break
        node = f"Seg({codon})"
        G.add_node(node, type="segment_3mer", label=codon)
        codon_nodes.append(node)
    
    # Construire le chemin complet : start -> codon_nodes -> end
    full_path_nodes = [start] + codon_nodes + [end]
    
    # Ajouter les arêtes entre chaque noeud successif du chemin
    for i in range(len(full_path_nodes) - 1):
        G.add_edge(full_path_nodes[i], full_path_nodes[i+1],
                   interaction="Codon_Path",
                   weight_cost=0.5,
                   weight_stability=0.9,
                   weight_error=0.1)
    
    return start, end

