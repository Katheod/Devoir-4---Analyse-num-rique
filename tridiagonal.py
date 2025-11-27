import numpy as np

def tridiagonal(I, D, S, b):
    """
    Résout un système linéaire tridiagonal Ax = b avec l'agorithme de Thomas.
    Où l'agortihme a été trouvé -> Url: https://www.thevisualroom.com/01_barba_theory/tri_diagonal_matrix.html

    arguments:
        I : vecteur diagonale inférieure 
        D : vecteur diagonale principale 
        S : vecteur diagonale supérieure 
        b : vecteur 

    returns:
        x : solution du système linéaire
    """
    n = len(D) # taille du système
    Ic, Dc, Sc, bc = map(np.array, (I, D, S, b)) # copie des vecteurs
    for i in range(1, n):
        # Modification de Ic car le message suivant se lève (IndexError: index 1 is out of bounds for axis 0 with size 1) si Ic est découpé par i (comme dans l'algortihme de Thomas) au lieu de Ic[i-1]
        mc = Ic[i-1]/Dc[i-1]
        Dc[i] = Dc[i] - mc*Sc[i-1]
        bc[i] = bc[i] - mc*bc[i-1]

    x = np.zeros(n)
    x[n-1] = bc[n-1] / Dc[n-1]

    for j in range(n-2, -1, -1):
        x[j] = (bc[j]-Sc[j]*x[j+1])/Dc[j]

    del Dc, Sc, bc # effacer les variables temporaires

    return x
