from data import *


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


