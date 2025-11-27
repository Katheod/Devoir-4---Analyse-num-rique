import numpy as np
from tridiagonal import tridiagonal

def problimite(h, P, Q, R, a, b, alpha, beta):
    """
    Résout un problème à valeurs aux limites par différences finies
    
    Arguments:
    h : pas entre deux points consecutifs
    P : vecteur des valeurs de p 
    Q : vecteur des valeurs de q 
    R : vecteur des valeurs de r
    a, b : bornes de l'intervalle
    alpha, beta : conditions aux limites
    
    Sortie:
    y : vecteur solution
    """

    N = int((b - a) / h) - 1  # nombre de noeuds intérieurs

    I = -1.0 - P[1:] * (h / 2.0)       # taille N-1 (a_2..a_N)
    D = 2.0 + Q * (h ** 2)            # taille N
    S = -1.0 + P[:-1] * (h / 2.0)     # taille N-1 (c_1..c_{N-1})
    b_vec = -R * (h ** 2)             # taille N

    # Conditions aux limites incorporées au second membre
    b_vec[0] = b_vec[0] + (1.0 + P[0] * (h / 2.0)) * alpha
    b_vec[N-1] = b_vec[N-1] + (1.0 - P[N-1] * (h / 2.0)) * beta

    # Résoudre
    y = tridiagonal(I, D, S, b_vec)
    return y