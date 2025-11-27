import numpy as np
import matplotlib.pyplot as plt
from problimite import problimite

a = 0.9
b = 1.0
alpha = 0.0
beta = 0.0

def y(x):
    """y(x)"""
    return (0.4 - 0.4 / x**2) - (0.4 - 0.4 / 0.81) * np.log(x) / np.log(0.9)

def p(x):
    """Fonction p(x)"""
    return -1.0 / x

def q(x):
    """Fonction q(x)"""
    return 0.0

def r(x):
    """Fonction r(x)"""
    return -1.6 / x**4


# h = 1/30
h1 = 1.0 / 30
N1 = int((b - a) / h1)
points1 = np.array([a + i * h1 for i in range(1, N1+1)])
P1 = np.array([p(x) for x in points1])
Q1 = np.array([q(x) for x in points1])
R1 = np.array([r(x) for x in points1])
y1 = problimite(h1, P1, Q1, R1, a, b, alpha, beta)

# h = 1/100
h2 = 1.0 / 100
N2 = int((b - a) / h2)
points2 = np.array([a + i * h2 for i in range(1, N2+1)])
P2 = np.array([p(x) for x in points2])
Q2 = np.array([q(x) for x in points2])
R2 = np.array([r(x) for x in points2])
y2 = problimite(h2, P2, Q2, R2, a, b, alpha, beta)

# y(x)
x_sol = np.linspace(a, b, 100)
y_sol = y(x_sol)

# Question a)
plt.plot(x_sol, y_sol, 'k-', label='y(x)')
plt.plot(points1, y1, 'ro-', label='h = 1/30')
plt.plot(points2, y2, 'bs-', label='h = 1/100')
plt.xlabel('x (distance à l\'axe)')
plt.ylabel('y (température)')
plt.title('Distribution de température entre cylindres coaxiaux')
plt.legend()
plt.grid(True)
plt.figure()

# ordre de convergence
h_valeurs = [1e-2, 1e-3, 1e-4, 1e-5]
erreurs = []

for h in h_valeurs:
    N = int((b - a) / h)
    points = np.array([a + i * h for i in range(1, N + 1)])
    
    P = np.array([p(x) for x in points])
    Q = np.array([q(x) for x in points])
    R = np.array([r(x) for x in points])
    
    y_num = problimite(h, P, Q, R, a, b, alpha, beta)
    y_sol = y(points)
    
    # Calcul de l'erreur
    erreur = np.max(np.abs(y_num - y_sol))
    erreurs.append(erreur)

erreurs = np.array(erreurs)
h_valeurs = np.array(h_valeurs)

plt.loglog(h_valeurs, erreurs, 'bo-', label="E(h)")
plt.xlabel('h')
plt.ylabel('Erreur')
plt.title('Convergence de la méthode des différences finies')
plt.legend()
plt.grid(True)

plt.show()
