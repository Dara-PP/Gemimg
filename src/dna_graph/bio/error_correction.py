import dna_graph.bio.constants as const

def add_code_correcteur(G):
    """Ajoute le noeud Code_Correcteur et le relie au Promoteur et aux processus."""
    G.add_node("Code_Correcteur", type="correcteur", label="Code Correcteur")

    G.add_edge("Code_Correcteur", "Promoteur", interaction="Corrige", 
               weight_cost=const.DEFAULT_COST_BUILD_COR_COR, 
               weight_stability=const.DEFAULT_STABILITY_BUILD_COR_COR, 
               weight_error=const.DEFAULT_ERROR_BUILD_COR_COR)
    
    if "Methy" in G:
        G.add_edge("Code_Correcteur", "Methy", interaction="Corrige",
                   weight_cost=const.DEFAULT_COST_BUILD_CODE, 
                   weight_stability=const.DEFAULT_STABILITY_BUILD_CODE, 
                   weight_error=const.DEFAULT_ERROR_BUILD_CODE)
    if "Reparation" in G:
        G.add_edge("Code_Correcteur", "Reparation", interaction="Corrige",
                   weight_cost=const.DEFAULT_STABILITY_BUILD_CODE, 
                   weight_stability=const.DEFAULT_STABILITY_BUILD_CODE, 
                   weight_error=const.DEFAULT_ERROR_BUILD_CODE)