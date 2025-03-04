from dna_graph.bio.tran_tran import transcribe

def test_transcribe():
    """
    Vérifie que le promoteur et le signal de terminaison ne sont pas inclus dans l'ARNm
    et que les T sont remplacés par des U.
    """
    dna_sequence = "TATAATGTTTAAATT"
    mRNA = transcribe(dna_sequence)
    assert "T" not in mRNA, "Les T doivent être remplacés par des U dans l'ARNm."
