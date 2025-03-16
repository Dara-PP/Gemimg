import os
import logging

# ----- Parametres Draw -----#
# Taille fenetre
XSIZE = 12
YSIZE = 7

# Palette couleur du graphe
COLOR_MAP = {
    "virtual": "red",
    "segment_3mer": "pink",
    "base": "lightblue",
    "classification": "lightgreen",
    "splicing": "lime",
    "motif": "orange",
    "gene": "teal",
    "process": "yellow",
    "correcteur": "plum",
    "epigenetics": "gold",
    "post_transcription": "darkcyan",
    "chromatin": "brown",
    "cis": "violet",
    "trans": "grey",
}


# Organisation des couches(x_coord, base_y, y_spacing) du graphe
LAYER_CONFIG = {
    "virtual":           (-12,  0, 100),   # start/end
    "segment_3mer":      ( -8,  0, 10.5),  # codons (3-mers)
    "base":              ( -4,  6, 200),   # A, C, G, T
    "classification":    (  0,  5, 200),   # Purines, Pyrimidines
    "motif":             (  4,  8, 200),   # Promoteur, Exon, Intron
    "splicing":          (  8,  -100, 200),   # Spliceosome
    "gene":              ( 12,  -100, 200),   # Gène
    "process":           ( 16,  -50, 200),    # Méthylation, Réparation
    "correcteur":        ( 20,  -105, 200),   # Code_Correcteur
    "epigenetics":       ( 24,  3, 200),      # Histone_Acetylation
    "post_transcription":( 28,  3, 200),
    "chromatin":         ( 32,  7, 200),
    "cis":               ( 33,  5, 200),   # Enhancer, Silencer
    "trans":             ( 38,  5, 200),   # TF1, TF2
}

# Parametres des nodes et textes 
NODE_SIZE = 400
FONT_SIZE_NODE = 8
FONT_SIZE_INTERACTION = 5
FONT_COLOR = "red"
OPTI_PATH = "blue"
OPTI_PSIZE = 1.5

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
ALPHA = 0.1
BETA = 0.1
GAMMA = 0.5

# ----- Paramètres Gaussien -----
DEFAULT_MUTATION_RATE = 0.01
NUMB_TEST = 50
SEED = 42
NUMBER_TEST = 100

# ----- Paramètres Graph 2-----
NBR_BEST = 10

# Message à traduire par défaut
DEFAULT_MESSAGE = "world"

# Couche par lasquel le chemin doit obligatoirement passer afin de respecter les contraintes biologique
MANDATORY_NODES = ["Promoteur", "Code_Correcteur", "Gene", "Purines", "Pyrimidines","Enhancer", "Silencer", "TF1", "TF2"]

# ----- Paramètres Biologiques Contraintes -----
# Constantes pour la simulation d'expression génétique
PROMOTER = "TATAATG"
ADRN = "ATG"
TERMINATION_SIGNAL = "ATT"


