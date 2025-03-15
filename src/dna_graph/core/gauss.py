import numpy as np
import difflib
from dna_graph.codec.encode_decode import convert_message_to_bases
from dna_graph.bio.tran_tran import modify_dna_sequence
from sklearn.cluster import KMeans

def gaussian_kernel_test_sentence(sentence, default_alpha, default_beta, default_gamma, default_mutation_rate,
                                  num_tests=10, sigma_alpha=0.1, sigma_beta=0.1, sigma_gamma=0.1, sigma_mutation=0.005,
                                  random_seed=None):
    """
    Découpe une phrase en mots et effectue le test pour chacun d'eux.
    Seed pour des résultats reproductibles pour les tests et la clustering.
    Retourne un dictionnaire où chaque mot (sans ponctuation) est associé à sa liste de résultats.
    """
    if random_seed is not None:
        np.random.seed(random_seed)
    words = sentence.split()
    sentence_results = {}
    for word in words:
        # Pour chaque mot, on génère les tests sans réinitialiser le seed
        results = gaussian_kernel_test(word, default_alpha, default_beta, default_gamma,
                                    default_mutation_rate, num_tests,
                                    sigma_alpha, sigma_beta, sigma_gamma, sigma_mutation,
                                    random_seed=None)
        sentence_results[word] = results
    return sentence_results


def gaussian_kernel_test(word, default_alpha, default_beta, default_gamma, default_mutation_rate,
                         num_tests=10, sigma_alpha=0.1, sigma_beta=0.1, sigma_gamma=0.1, sigma_mutation=0.005,
                         random_seed=None):
    """
    Génère plusieurs variantes de la séquence correspondant à un mot en perturbant les paramètres
    alpha, beta, gamma et le taux de mutation avec une distribution gaussienne.
    
    Paramètres :
      - word (str) : Le mot à tester.
      - default_alpha, default_beta, default_gamma (float) : Valeurs par défaut pour les paramètres.
      - default_mutation_rate (float) : Taux de mutation par défaut.
      - num_tests (int) : Nombre d'échantillons à générer.
      - sigma_alpha, sigma_beta, sigma_gamma (float) : Écart-type des perturbations pour chaque paramètre.
      - sigma_mutation (float) : Écart-type de la perturbation pour le taux de mutation.
      - random_seed (int, optionnel) : Seed pour la reproductibilité.
      
    Retourne :
      - results (list) : Liste de dictionnaires pour chaque test contenant :
          * 'alpha', 'beta', 'gamma' : Paramètres utilisés.
          * 'mutation_rate' : Taux de mutation (compris entre 0 et 1).
          * 'sequence' : Séquence modifiée.
          * 'score' : Similarité (0 à 1) entre la séquence originale et la séquence modifiée.
    """
    #if random_seed is not None:
        #np.random.seed(random_seed)
        
    results = []
    # Conversion du mot en liste de bases puis en chaîne
    base_list = convert_message_to_bases(word)
    original_sequence = ''.join(base_list)
    
    for _ in range(num_tests):
        # Perturber les paramètres par bruit gaussien
        alpha = np.random.normal(default_alpha, sigma_alpha)
        beta = np.random.normal(default_beta, sigma_beta)
        gamma = np.random.normal(default_gamma, sigma_gamma)
        mutation_rate = np.random.normal(default_mutation_rate, sigma_mutation)
        # S'assurer que le taux de mutation reste dans [0, 1]
        mutation_rate = np.clip(mutation_rate, 0, 1)
        
        # Appliquer les modifications sur la séquence
        mutated_sequence = modify_dna_sequence(original_sequence, mutation_rate=mutation_rate)
        
        # Calcul du score de similarité entre séquence originale et modifiée
        score = difflib.SequenceMatcher(None, original_sequence, mutated_sequence).ratio()
        
        results.append({
            'alpha': alpha,
            'beta': beta,
            'gamma': gamma,
            'mutation_rate': mutation_rate,
            'sequence': mutated_sequence,
            'score': score
        })
        
    return results

def cluster_results(test_results, n_clusters=5):
    """
    Regroupe les résultats en n_clusters et sélectionne, pour chaque cluster, le test avec le meilleur score.
    
    Paramètres :
      - test_results (list) : Liste de dictionnaires issus de gaussian_kernel_test.
      - n_clusters (int) : Nombre de clusters souhaités.
    
    Retourne :
      - best_representatives (list) : Liste des meilleurs résultats (un par cluster).
    """
    # Ajuster le nombre de clusters si le nombre de tests est inférieur
    n_clusters = min(n_clusters, len(test_results))
    
    features = np.array([
        [res['alpha'], res['beta'], res['gamma'], res['mutation_rate']]
        for res in test_results
    ])
    
    kmeans = KMeans(n_clusters=n_clusters, random_state=42).fit(features)
    labels = kmeans.labels_
    
    for i, res in enumerate(test_results):
        res['cluster'] = labels[i]
    
    best_in_cluster = {}
    for res in test_results:
        cluster = res['cluster']
        if cluster not in best_in_cluster or res['score'] > best_in_cluster[cluster]['score']:
            best_in_cluster[cluster] = res
            
    best_representatives = list(best_in_cluster.values())
    return best_representatives

def get_word_test_results(sentence, num_tests, n_best, alpha, beta, gamma, mutation_rate, seed):
    """
    Pour chaque mot de la phrase, effectue les tests gaussiens et retourne un dictionnaire
    où chaque clé est un mot et la valeur est une liste (de longueur n_best) des meilleurs
    résultats sous forme de dictionnaires (contenant 'sequence' et 'score').
    
    Paramètres:
      - sentence (str) : La phrase à traiter.
      - num_tests (int) : Nombre d'échantillons générés par mot.
      - n_best (int) : Nombre de meilleurs résultats à conserver par mot.
      - alpha, beta, gamma (float) : Paramètres de l'algorithme gaussien.
      - mutation_rate (float) : Taux de mutation par défaut.
      - seed (int) : Seed pour la reproductibilité.
    
    Retourne:
      - filtered_results (dict) : Dictionnaire avec pour chaque mot une liste de dictionnaires.
    """
    # Utilise la fonction gaussienne déjà définie qui retourne un dictionnaire {mot: [résultats]}
    sentence_results = gaussian_kernel_test_sentence(
        sentence, alpha, beta, gamma, mutation_rate,
        num_tests=num_tests, random_seed=seed
    )
    
    # Pour chaque mot, on trie les résultats par score décroissant et on garde les n_best
    filtered_results = {}
    for word, results in sentence_results.items():
        sorted_results = sorted(results, key=lambda r: r.get('score', 0), reverse=True)
        filtered_results[word] = sorted_results[:n_best]
    return filtered_results

