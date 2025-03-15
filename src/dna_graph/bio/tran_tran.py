import logging
from config.config import PROMOTER, TERMINATION_SIGNAL
import dna_graph.bio.genetic_code as gen_code
import random

def codon_to_amino_acid(codon: str) -> str:
    """
    Convertit un codon (ARN, 3 nucléotides) en son acide aminé correspondant (lettre).
    Si le codon n'est pas reconnu, retourne "?".
    La conversion U->T permet d'utiliser le dictionnaire en notation ADN.
    """
    codon_dna = codon.replace("U", "T")
    return gen_code.GENETIC_CODE.get(codon_dna, "?")

def transcribe(dna_sequence: str) -> str:
    """
    Transcrit une séquence d'ADN en ARNm en respectant le promoteur et le signal de terminaison.
    
    Contraintes :
      - Le promoteur ("TATAATG") doit être présent dans la séquence.
      - La transcription démarre juste après le promoteur.
      - La transcription se termine au signal de terminaison ("ATT") s'il est présent, sinon à la fin de la séquence.
    
    Retourne :
      mRNA (str) : La chaîne d'ARN messager obtenue (où T est remplacé par U).
    """
    
    if PROMOTER not in dna_sequence:
        raise ValueError("La sequence d'ADN doit commencer par le promoteur 'TATAATG'.")
    
    promoter_index = dna_sequence.index(PROMOTER)
    if promoter_index != 0:
        logging.info("Passage par le Promoteur detecte a l'index %d. Transcription initiee apres le promoteur.", promoter_index)
    else:
        logging.info("La sequence commence par le Promoteur.")
    
    # Début de la transcription juste après le promoteur
    start_transcription = promoter_index + len(PROMOTER)
    termination_index = dna_sequence.find(TERMINATION_SIGNAL, start_transcription)
    if termination_index == -1:
        termination_index = len(dna_sequence)
    
    # Extraction de la portion à transcrire
    transcript_dna = dna_sequence[start_transcription:termination_index]
    
    # Transcription : conversion de T en U pour simuler l'ARNm
    mRNA = transcript_dna.replace("T", "U")
    
    return mRNA

def translate(mRNA: str) -> str:
    """
    Traduit une chaîne d'ARNm en protéine.
    
    Contraintes :
      - La traduction doit commencer par le codon START (AUG) au tout début de l'ARNm.
      - La traduction se poursuit codon par codon jusqu'à la rencontre d'un codon stop (UAA, UAG, UGA).
    
    Retourne :
      protein (str) : La protéine synthétisée.
    """
    # Vérifier que l'ARNm commence bien par le codon START (AUG)
    if not mRNA.startswith("AUG"):
        raise ValueError("L'ARNm ne commence pas par le codon START (AUG).")
    
    stop_codons = {"UAA", "UAG", "UGA"}
    protein = []
    
    # Puisque l'ARNm commence par AUG, on démarre la traduction dès le début
    for i in range(0, len(mRNA), 3):
        codon = mRNA[i:i+3]
        if len(codon) < 3:
            break  # Fin de la séquence ou padding insuffisant
        if codon in stop_codons:
            break  # Arrêt dès le codon stop
        protein.append(codon_to_amino_acid(codon))
    
    return ''.join(protein)

def introduce_mutations(dna_sequence: str, mutation_rate) -> str:
    """
    Introduit des mutations dans la séquence d'ADN en substituant aléatoirement des nucléotides,
    avec un taux de mutation donné (par exemple, 1%).
    
    Paramètres :
      - dna_sequence (str) : La séquence d'ADN originale.
      - mutation_rate (float) : Probabilité qu'une position soit mutée.
      
    Retourne :
      - mutated_sequence (str) : La séquence d'ADN après mutation.
    """
    nucleotides = ["A", "T", "C", "G"]
    mutated_sequence = list(dna_sequence)
    for i in range(len(mutated_sequence)):
        if random.random() < mutation_rate:
            current = mutated_sequence[i]
            alternatives = [n for n in nucleotides if n != current]
            mutated_sequence[i] = random.choice(alternatives)
    return "".join(mutated_sequence)

def introduce_insertion(dna_sequence: str, insertion_rate) -> str:
    """
    Introduit des insertions aléatoires dans la séquence d'ADN.
    
    Paramètres :
      - dna_sequence (str) : La séquence d'ADN originale.
      - insertion_rate (float) : Probabilité d'insertion après chaque nucléotide.
      
    Retourne :
      - new_sequence (str) : La séquence d'ADN avec des insertions.
    """
    nucleotides = ["A", "T", "C", "G"]
    result = []
    for nucleotide in dna_sequence:
        result.append(nucleotide)
        if random.random() < insertion_rate:
            result.append(random.choice(nucleotides))
    return "".join(result)

def introduce_deletion(dna_sequence: str, deletion_rate) -> str:
    """
    Introduit des délétions aléatoires dans la séquence d'ADN.
    
    Paramètres :
      - dna_sequence (str) : La séquence d'ADN originale.
      - deletion_rate (float) : Probabilité de supprimer un nucléotide.
      
    Retourne :
      - new_sequence (str) : La séquence d'ADN après suppression de certains nucléotides.
    """
    result = []
    for nucleotide in dna_sequence:
        if random.random() >= deletion_rate:
            result.append(nucleotide)
    return "".join(result)

def modify_dna_sequence(dna_sequence: str, mutation_rate: float = 0.01,
                        insertion_rate: float = 0.005, deletion_rate: float = 0.005) -> str:
    #dna_sequence: str, mutation_rate: float = 0.01, insertion_rate: float = 0.005, deletion_rate: float = 0.005
    """
    Applique successivement des mutations par substitution, insertion et délétion sur la séquence d'ADN.
    Cela modifie directement la séquence et peut altérer des régions critiques (promoteur, codons, etc.).
    """
    if dna_sequence.startswith(PROMOTER):
        promoter_part = dna_sequence[:len(PROMOTER)]
        # On suppose que le start codon (ATG) occupe les 3 bases suivant le promoteur
        start_codon = dna_sequence[len(PROMOTER):len(PROMOTER)+3]
        rest = dna_sequence[len(PROMOTER)+3:]
        # Appliquer les modifications uniquement sur la partie après le start codon
        rest = introduce_mutations(rest, mutation_rate)
        rest = introduce_insertion(rest, insertion_rate)
        rest = introduce_deletion(rest, deletion_rate)
        return promoter_part + start_codon + rest
    else:
        # Si le promoteur n'est pas détecté, on applique les modifications sur toute la séquence
        seq = introduce_mutations(dna_sequence, mutation_rate)
        seq = introduce_insertion(seq, insertion_rate)
        seq = introduce_deletion(seq, deletion_rate)
        return seq