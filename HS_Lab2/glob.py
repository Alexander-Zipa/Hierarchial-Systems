from scipy import log
from scipy.optimize import linprog


def generate_matr(A):
    B = [0] * len(A)
    for k in range(0, len(A)):
        B[k] = [[0] * len(A) for i in [0] * len(A)]
        for i in range(0, len(A)):
            for j in range(0, len(A)):
                if (i == k) or (j == k):
                    B[k][i][j] = [A[i][j][0], A[i][j][1]]
                else:
                    B[k][i][j] = [A[i][k][0] * A[k][j][0], A[i][k][1] * A[k][j][1]]
    return B


def spectrum(v):
    scale = [0] * 11
    for k in range(0, 11):
        scale[k] = k * 0.1
    s = [[0] * len(scale) for i in [0] * len(v)]
    delta = [0] * len(v)
    for i in range(len(v)):
        delta[i] = [0] * len(v[i])
        for j in range(len(v[i])):
            delta[i][j] = (((v[i][j][0] + v[i][j][1]) / 2) ** 2 + (1 / 3) * ((v[i][j][1] - v[i][j][0]) / 2) ** 2) ** 0.5
            tmp = [abs(k - delta[i][j]) for k in scale]
            k = tmp.index(min(tmp))
            s[j][k] += 1
    return s


def ky(s):
    m = sum(s)
    N = len(s)
    a = (1 / m) * sum([(k + 1) * s[k] for k in range(0, N)])
    G = m / (N * log(m) * log(N))
    entropy = 0
    for k in s:
        if k > 0:
            entropy -= (k / m) * log(k / m)
    K = 1 - ((1 / m) * sum([abs(k + 1 - a) * s[k] for k in range(0, N)]) +
             entropy) / (G * sum([abs(k + 1 - (N + 1) / 2) for k in range(0, N)]) + log(N))
    return K


def global_weights(w, wc):
    w_glob = [0] * len(w[0])
    for i in range(0, len(w[0])):
        w_glob[i] = [0] * 2
        cl = [w[k][i][0] for k in range(0, len(wc))]
        cu = [-w[k][i][1] for k in range(0, len(wc))]
        wcl = [wc[k][0] for k in range(0, len(wc))]
        wcu = [wc[k][1] for k in range(0, len(wc))]
        x = linprog(c=cl, A_eq=[[1] * len(wc)], b_eq=[1], bounds=tuple([(wcl[i], wcu[i]) for i in range(0, len(wc))]))
        y = linprog(c=cu, A_eq=[[1] * len(wc)], b_eq=[1], bounds=tuple([(wcl[i], wcu[i]) for i in range(0, len(wc))]))
        w_glob[i][0] = x.fun
        w_glob[i][1] = -y.fun
    return w_glob
