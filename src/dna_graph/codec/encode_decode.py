def convert_message_to_bases(message):
    """
    Convertit un message en une liste de bases.
    Chaque caractère est encodé en 4 bases en convertissant sa valeur ASCII en base 4.
    
    Exemple :
      'h' (ASCII 104) en base 4 → [1,2,2,0] → ["C", "G", "G", "A"]
    """
    bases = []
    mapping = {0: "A", 1: "C", 2: "G", 3: "T"}
    for char in message:
        ascii_val = ord(char)
        # Conversion en base 4 sur 4 chiffres
        digits = []
        for _ in range(4):
            digits.append(ascii_val % 4)
            ascii_val //= 4
        digits = digits[::-1]  # pour remettre dans l'ordre
        for d in digits:
            bases.append(mapping[d])
    return bases

def decode_message_from_path(path, original_bases, original_length):
    """
    Reconstruit un message à partir des bases originales.
    
    Ici, on considère que 'original_bases' contient la séquence complète des bases générées initialement.
    On découpe cette séquence en groupes de 4 bases pour reconstituer chaque caractère.
    
    Le paramètre original_length permet de limiter le décodage au nombre de caractères originaux.
    """
    # On prend les n*4 premières bases
    relevant_bases = original_bases[:original_length * 4]
    mapping_inv = {"A": 0, "C": 1, "G": 2, "T": 3}
    message = ""
    for i in range(0, len(relevant_bases), 4):
        digits = relevant_bases[i:i+4]
        num = 0
        for d in digits:
            num = num * 4 + mapping_inv[d]
        message += chr(num)
    return message

def extract_base_path(full_path):
    """Extrait du chemin complet en décomposant les noeuds de bases et codons.
       Cette version gère les nœuds au format 'Seg(codon)_posX'."""
    base_path = []
    for node in full_path:
        if node in {"A", "T", "C", "G"}:
            base_path.append(node)
        elif node.startswith("Seg("):
            # Trouver la parenthèse fermante pour extraire le codon
            end_index = node.find(")")
            if end_index != -1:
                codon = node[4:end_index]
                base_path.extend(list(codon))
    return base_path

