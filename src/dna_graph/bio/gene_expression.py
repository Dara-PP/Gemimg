# fichier: dna_graph/contraintes/gene_expression.py
from dna_graph.bio.tran_tran import transcribe, translate

def simulate_gene_expression(dna_sequence: str) -> str:
    """
    Simule l'expression génique en appliquant d'abord la transcription, et la traduction.
    
    Étapes :
      1. Transcription : conversion de l'ADN en ARNm après passage par le promoteur.
      2. Traduction : conversion de l'ARNm (corrigé le cas échéant) en protéine.
    
    Paramètres :
      - dna_sequence (str) : La séquence d'ADN à traiter.
      
    Retourne :
      - protein (str) : La protéine synthétisée.
    """
    # Étape 1 : Transcription
    mRNA = transcribe(dna_sequence)
    
    # Étape 2 : Traduction
    protein = translate(mRNA)

    return protein
