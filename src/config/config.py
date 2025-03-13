import os
import logging

# ----- Parametres Draw -----#
# Taille fenetre
XSIZE = 12
YSIZE = 7

# Palette couleur du graphe
COLOR_MAP = {
        "base": "lightblue",
        "classification": "lightgreen",
        "motif": "orange",
        "process": "yellow",
        "correcteur": "plum",
        "segment_3mer": "pink",
        "cis": "violet",
        "trans": "grey",
        "virtual": "red"
    }

# Organisation des couches x, y du graphe
LAYER_CONFIG = {
        "virtual": (-8, 5),
        "segment_3mer": (-4, 1),
        "base": (1, 15),
        "classification": (3, 7),
        "splicing": (5, 8),             
        "motif": (7, 8),
        "gene": (9, 8),
        "splicing": (11, 7), 
        "process": (13, 6),
        "correcteur": (15, 6),
        "epigenetics": (17, 4),         
        "post_transcription": (19, 4),  
        "chromatin": (21, 4),
        "cis": (23, 6),
        "trans": (25, 10),       
    }
"""
LAYER_CONFIG = {
    "virtual": (-10, 5),             # Nœuds virtuels en marge
    "segment_3mer": (-8, 15),         # Couche codons
    "base": (-6, 25),                # Bases nucléotidiques
    "classification": (-4, 20),      # Purines / Pyrimidines
    "motif": (-2, 18),               # Motifs (Promoteur, Exon, Intron, etc.)
    "gene": (0, 16),                 # Nœud Gene
    "splicing": (2, 14),             # Épissage
    "process": (4, 12),              # Processus (Méthylation, Réparation)
    "correcteur": (6, 12),           # Code correcteur
    "epigenetics": (8, 10),          # Modifications épigénétiques
    "post_transcription": (10, 10),  # Régulation post-transcriptionnelle
    "chromatin": (12, 10),           # Structure chromatinienne
    "cis": (14, 12),                 # Régulation cis
    "trans": (16, 14)                # Régulation trans
}
"""


# Parametres des nodes et textes
NODE_SIZE = 400
FONT_SIZE_NODE = 8
FONT_SIZE_INTERACTION = 5
FONT_COLOR = "red"
OPTI_PATH = "blue"

# ----- Chemins et répertoires -----
# Définir le dossier racine pour les logs (ici, un dossier 'logs' au même niveau que le projet)
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
LOG_DIR = os.path.join(BASE_DIR, '..', 'logs')

# Créer le dossier des logs s'il n'existe pas déjà
if not os.path.exists(LOG_DIR):
    os.makedirs(LOG_DIR)

# Chemin complet du fichier de log
LOG_FILE = os.path.join(LOG_DIR, 'genimg.log')

# ----- Configuration du Logging -----
LOG_LEVEL = logging.DEBUG
LOG_FORMAT = '%(asctime)s - %(levelname)s - %(name)s - %(filename)s -%(funcName)s - %(lineno)d - %(message)s'
LOG_FILE_MODE = 'w'  # 'w' pour écraser à chaque démarrage, 'a' pour ajouter

# ----- Paramètres d'Optimisation -----
# Coefficients pour la fonction de coût dans l'optimisation du graphe (moyenne), build, stabilité relation, taux erreur
ALPHA = 0.2
#BETA = 0.8
BETA = 100
GAMMA = 0.5

# Message à traduire par défaut
DEFAULT_MESSAGE = "world"

# Couche par lasquel le chemin doit obligatoirement passer afin de respecter les contraintes biologique
MANDATORY_NODES = ["Promoteur", "Code_Correcteur", "Gene", "Purines", "Pyrimidines","Enhancer", "Silencer", "TF1", "TF2"]

# ----- Paramètres Biologiques Contraintes -----
# Constantes pour la simulation d'expression génétique
PROMOTER = "TATAATG"
ADRN = "ATG"
TERMINATION_SIGNAL = "ATT"


