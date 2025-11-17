from tridiagonal import tridiagonal

def problimite(h, P, Q, R, a, b, alpha, beta):
    """Solve the boundary value problem using the finite difference method.
    Parameters:
    h : float
        Step size for the finite difference method.
    P : function
        Coefficient function P(x).
    Q : function
        Coefficient function Q(x).
    R : function
        Coefficient function R(x).
    a : scalaire
        Left boundary of the interval.
    b : scalaire
        Right boundary of the interval.
    alpha : scalaire
        Boundary condition at x = a (y(a) = alpha).
    beta : scalaire
        Boundary condition at x = b (y(b) = beta).
    Returns:
    y : vecteur
        La solution approchée du problème aux limites.
    """

    n = int((b - a) / h) - 1
    I = [0] * (n - 1)
    D = [0] * (n - 1)
    S = [0] * (n - 1)
    b = [0] * (n - 1)

    for i in range(1, n):
        x_i = a + i * h
        I[i - 1] = -1 - P(x_i)*(h / 2)
        D[i - 1] = 2 + Q(x_i)*(h ** 2)
        S[i - 1] = -1 - P(x_i)*(h / 2)
        b[i - 1] = -R(x_i)*(h ** 2)

    
    # Adjust the first equation for boundary condition at a
    x_i = a + 0*h
    b[0] =  -R(x_i) + (1 + P(x_i)*(h/2))* alpha
    I[0] = 0

    # Adjust the last equation for boundary condition at b
    x_i = a + n*h
    b[-1] = -R(x_i) + (1 + P(x_i)*(h/2)) * beta
    S[-1] = 0

    return tridiagonal(D, I, S, b)