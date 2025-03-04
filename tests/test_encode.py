from dna_graph.codec.encode_decode import convert_message_to_bases

def test_convert_message_to_bases():
    """
    Vérifie que le résultat est une liste et qu'elle a la bonne longueur (4 bases par caractère)
    """
    message = "ABC"
    bases = convert_message_to_bases(message)
    assert isinstance(bases, list)
    assert len(bases) == len(message) * 4