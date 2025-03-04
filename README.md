# Genimg

Genimg est un projet qui permet de convertir un message texte en une séquence de bases nucléotidiques (A, T, C, G), tout en respectant de manière basique les réactions génétiques et les interactions biologiques. 

## Objectifs du Projet

- **Conversion du Message**  
  Transformer un message en une séquence de bases nucléotidiques en utilisant une représentation en base 4. Afin de garantir que le message converti respecte assez fidèlement les réactions génétiques.

- **Simulation basique des Réactions Génétiques**  
  Respecter les étapes clés de la transcription (conversion ADN → ARNm) et de la traduction (ARNm → Protéine) en intégrant les contraintes biologiques telles que :
  - Passage par un **Promoteur** pour initier la transcription.
  - Utilisation d’un **signal de terminaison** pour clore la transcription.
  - Démarrage de la traduction par le **codon START** (ATG) et arrêt au codon stop.

- **Représentation par Graphe**  
  Construire un graphe de connaissances qui modélise :
  - Les bases nucléotidiques et leurs interactions.
  - Les motifs, segments et processus biologiques.
  - Les chemins de conversion optimisés à travers d'algorithmes

## Prérequis

- **Dépendances**  
  - `pytest`
  - `networkx`
  - `matplotlib`

- **Installation**  
  ```bash
  pip install -r requirements.txt
  pip install -e .

## Utilisation
```bash
  dna_graph [-h] [-m MESSAGE] [--alpha ALPHA] [--beta BETA] [--gamma GAMMA] [--version]

- **Usage**  
  ```bash
  dna_graph -m "Votre message ici" --alpha 0.2 --beta 0.8 --gamma 0.5    

- **Désinstallation**
  ```bash
  pip uninstall dna_graph

  
