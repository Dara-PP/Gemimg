import logging
from config.config import PROMOTER, TERMINATION_SIGNAL
import dna_graph.bio.genetic_code as gen_code

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

