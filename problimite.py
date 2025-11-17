from tridiagonal import tridiagonal

def problimite(h, P, Q, R, a, b, alpha, beta):
    """Résoudre les problèmes à valeurs aux limites avec la méthode des différences finies.
    Arguements:
    h : float
    P : Vecteur
        contient les évaluations de la fonctions p aux noeuds xi.
    Q : Vecteur
        contient les évaluations de la fonctions q aux noeuds xi.
    R : Vecteur
        contient les évaluations de la fonctions r aux noeuds xi.
    a : scalaire
        limite de l'invertalle gauche.
    b : scalaire
        limite de l'intervalle droit.
    alpha : scalaire
        Condition aux limites en x = a (y(a) = alpha).
    beta : scalaire
        Condition aux limites en x = b (y(b) = beta).
    Returns:
    y : vecteur
        La solution approchée du problème aux limites.
    """
    # n = n+1
    N = int((b - a) / h) - 1
    I = [0] * N
    D = [0] * N
    S = [0] * N
    B = [0] * N

    # n équivaut à n-1 dans ce cas car range(1,n) donne n éléments de 1 à n-1
    for i in range(1, N):
        x_i = a + i * h
        if i == 1:
            I[i] = 0
            D[i] = 2 + Q(x_i)*(h ** 2)
            S[i] = -1 - P(x_i)*(h / 2)
            B[i] =  -R(x_i)*(h**2) + (1 + P(x_i)*(h/2))* alpha
        elif i == N-1:
            S[i] = 0
            I[i] = -1 - P(x_i)*(h / 2)
            D[i] = 2 + Q(x_i)*(h ** 2)
            B[i] = -R(x_i)*(h**2) + (1 - P(x_i)*(h/2)) * beta
        else:
            I[i] = -1 - P(x_i)*(h / 2)
            D[i] = 2 + Q(x_i)*(h ** 2)
            S[i] = -1 - P(x_i)*(h / 2)
            B[i] = -R(x_i)*(h ** 2)

    Y = alpha + tridiagonal(D[1:N], I[1:N], S[1:N], B[1:N]) + beta
    # Trouver le vecteur y en résolvant le système tridiagonal
    return Y