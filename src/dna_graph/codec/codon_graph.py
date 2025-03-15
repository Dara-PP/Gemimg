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

def add_edge_if_not_exists(G, u, v, **attrs):
    if not G.has_edge(u, v):
        G.add_edge(u, v, **attrs)

def add_codon_subgraph_bio(G, base_list, genetic_code, aa_to_codons):
    """
    Ajoute un sous-graphe codon au graphe G en utilisant la séquence de bases.
    Crée un graphe en couches permettant des chemins alternatifs (via la dégénérescence)
    tout en évitant les doublons.
    
    Retourne :
      - (start, end) : tuple contenant le nœud "start" et le nœud "end".
    """
    start, end = "start", "end"
    for node, label in [(start, "Start"), (end, "End")]:
        if node not in G:
            G.add_node(node, type="virtual", label=label)

    num_codons = len(base_list) // 3
    layers = []
    for pos in range(num_codons):
        # Extraire le codon du message pour cette position
        msg_codon = "".join(base_list[pos*3: pos*3+3])
        aa = genetic_code.get(msg_codon)
        # Récupère les options pour l'acide aminé, ou garde le codon du message s'il n'est pas reconnu
        codon_options = aa_to_codons.get(aa, [msg_codon]) if aa is not None else [msg_codon]
        # Enlever les doublons tout en conservant l'ordre
        seen = set()
        codon_options = [c for c in codon_options if c not in seen and not seen.add(c)]
        
        layer_nodes = []
        for codon in codon_options:
            node_name = f"Seg({codon})_pos{pos}"
            # Marquer comme "préféré" si c'est le codon issu du message
            preferred = (codon == msg_codon)
            if node_name not in G:
                G.add_node(node_name, type="segment_3mer", label=codon, preferred=preferred)
            layer_nodes.append(node_name)
        layers.append(layer_nodes)

    # Connexion du nœud start à la première couche
    if layers:
        for node in layers[0]:
            add_edge_if_not_exists(G, start, node,
                interaction="Codon_Path",
                weight_cost=0.5,
                weight_stability=0.9,
                weight_error=0.1,
                display=False)

    # Connexion entre chaque couche
    for i in range(len(layers) - 1):
        for u in layers[i]:
            for v in layers[i+1]:
                penalty = 0.1 if not G.nodes[v].get("preferred", False) else 0.0
                add_edge_if_not_exists(G, u, v,
                    interaction="Codon_Path",
                    weight_cost=0.5 + penalty,
                    weight_stability=0.9,
                    weight_error=0.1)

    # Connexion de la dernière couche au nœud end
    if layers:
        for node in layers[-1]:
            add_edge_if_not_exists(G, node, end,
                interaction="Codon_Path",
                weight_cost=0.5,
                weight_stability=0.9,
                weight_error=0.1,
                display=False)

    # Connecter 'Promoteur' aux nœuds préférés de la première et dernière couche
    if "Promoteur" in G and layers:
        first_pref = next((node for node in layers[0] if G.nodes[node].get("preferred", False)), layers[0][0])
        add_edge_if_not_exists(G, "Promoteur", first_pref,
            interaction="Codon_Path",
            weight_cost=0.5,
            weight_stability=0.9,
            weight_error=0.1,
            display=True)
        last_pref = next((node for node in layers[-1] if G.nodes[node].get("preferred", False)), layers[-1][0])
        add_edge_if_not_exists(G, last_pref, "Promoteur",
            interaction="Codon_Path",
            weight_cost=0.5,
            weight_stability=0.9,
            weight_error=0.1,
            display=True)

    return start, end
