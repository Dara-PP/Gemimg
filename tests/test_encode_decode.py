from dna_graph.codec.encode_decode import convert_message_to_bases, decode_message_from_path

def test_encode_decode():
    """
    La fonction decode_message_from_path n'utilise pas le paramètre "path"
    et se base uniquement sur la liste des bases et la longueur du message.
    """
    message = "Hello"
    bases = convert_message_to_bases(message)
    decoded = decode_message_from_path([], bases, len(message))
    assert decoded == message
