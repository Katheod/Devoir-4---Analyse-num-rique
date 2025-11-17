from tridiagonal import tridiagonal
from problimite import problimite

P = lambda x: 1/x
Q = lambda x: 0
R = lambda x: -1.6/(x**4)
a = 0.9
b = 1.0
alpha = 0
beta = 0
h = 1/30

y = problimite(h, P, Q, R, a, b, alpha, beta)
print(y)