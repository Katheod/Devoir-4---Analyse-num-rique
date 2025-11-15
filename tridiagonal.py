import numpy as np

def tridiagonal(D, I, S, b):
    """
    Résoudre la matrice tridiagonale Ax = b avec l'algorithme de Thomas.
    url = https://gist.github.com/vuddameri/75212bfab7d98a9c75861243a9f8f272
    
    Arguments:
    D : Vecteur
        La diagonale principal.
    I : Vecteur
        La diagonale inférieure.
    S : Vecteur
        La diagonal supérieur.
    b : Vecteur
        Le vecteur du coté droit.
        
    Return:
    X : Vecteur
        La solution du x.
    """
    N = len(I)
    # Create copies of the arrays to avoid modifying the original ones
    S_prime = np.zeros(N, dtype=float)
    b_prime = np.zeros(N, dtype=float)
    X = np.zeros(N, dtype=float)
    
    # Balayage avant
    S_prime[0] = S[0] / D[0]
    b_prime[0] = b[0] / D[0]
    
    for i in np.arange(1, N, 1):
        denom = D[i] - I[i] * S_prime[i - 1]
        S_prime[i] = S[i] / denom
        b_prime[i] = (b[i] - I[i] * b_prime[i - 1]) / denom
    
    # Substitution arrière
    X[N - 1] = b_prime[N - 1] # Pour obtenir le dernier xn
    
    for i in np.arange(N - 2, -1, -1):
        X[i] = (b_prime[i] - S_prime[i] * X[i + 1]) # Utiliser x[i + 1] pour trouver x[i]

    return X

I = [0, -1, -1, -1]
D = [2, 2, 2, 1]
S = [-1, -1, -1, 0]
b = [0, 0, 1, 0]

#x = tridiagonal(D, I, S, b)
#print(x)






    
