from data import *
import scipy as sp
from scipy.optimize import linprog


def section(A, a):
    B = []
    for i in range(0, len(A)):
        B.append(list(range(0, len(A))))
        for j in range(0, len(A)):
            if A[i][j] >= 1:
                B[i][j] = [fs[A[i][j]][1] - (fs[A[i][j]][1] - fs[A[i][j]][0]) * (1 - a),
                           fs[A[i][j]][1] + (fs[A[i][j]][2] - fs[A[i][j]][1]) * (1 - a)]
            else:
                B[i][j] = [1 / fs[A[j][i]][1] - (1 / fs[A[j][i]][1] - 1 / fs[A[j][i]][2]) * (1 - a),
                           1 / fs[A[j][i]][1] + (1 / fs[A[j][i]][0] - 1 / fs[A[j][i]][1]) * (1 - a)]
    return B


def interval_weights(A):
    n = len(A)
    c = [1] * n * (n - 1)
    c += [0] * n
    B = [0] * n * (n - 1)
    b = [0] * n * (n - 1)
    for i in range(0, n * (n - 1)):
        B[i] = [0] * n * n
    l = 0
    for i in range(0, n * (n - 1)):
        B[i][i] = -1
    for i in range(0, n-1):
        for j in range(i+1, n):
            B[l][n * (n - 1) + i] = -1
            B[l][n * (n - 1) + j] = 1
            B[round(l + n * (n - 1) / 2)][n * (n - 1) + i] = 1
            B[round(l + n * (n - 1) / 2)][n * (n - 1) + j] = -1
            b[l] = - sp.log(A[i][j][0])
            b[round(l + n * (n - 1) / 2)] = sp.log(A[i][j][1])
            l += 1
    l = [0] * n * (n - 1)
    l += [None] * n
    u = [None] * n * (n - 1)
    u += [0] * n
    delta = linprog(c=c, A_ub=B, b_ub=b, bounds=tuple([(l[i], u[i]) for i in range(0, n*n)]))
    l = 0
    for i in range(0, n - 1):
        for j in range(i + 1, n):
            A[i][j][0] *= sp.e ** (- delta.x[l])
            A[i][j][1] *= sp.e ** delta.x[round(l + n * (n - 1) / 2)]
            l += 1
    B = [0] * n * (n - 1)
    for i in range(0, n * (n - 1)):
        B[i] = [0] * n
    l = 0
    for i in range(0, n - 1):
        for j in range(i + 1, n):
            B[l][i] = -1
            B[l][j] = 1
            B[round(l + n * (n - 1) / 2)][i] = 1
            B[round(l + n * (n - 1) / 2)][j] = -1
            b[l] = - sp.log(A[i][j][0])
            b[round(l + n * (n - 1) / 2)] = sp.log(A[i][j][1])
            l += 1
    l = [None] * n
    u = [0] * n
    w = [0] * n
    for i in range(0, n):
        w[i] = [0] * 2
        c = [0] * n
        c[i] = 1
        w[i][0] = linprog(c=c, A_ub=B, b_ub=b, bounds=tuple([(l[i], u[i]) for i in range(0, n)]))
        w[i][0] = sp.e ** w[i][0].x[i]
        c[i] = -1
        w[i][1] = linprog(c=c, A_ub=B, b_ub=b, bounds=tuple([(l[i], u[i]) for i in range(0, n)]))
        w[i][1] = sp.e ** w[i][1].x[i]
    return w

