from loc import *
from glob import *

w = [0] * len(D)
for i in range(len(D)):
    w[i] = [0] * 2
    print('Interval matrix', i, 'section', 0,)
    w[i][0] = interval_weights(section(D[i], 0))
    print('Interval matrix', i, 'section', 0.5)
    w[i][1] = interval_weights(section(D[i], 0.5))
    print('w[0]=')
    disp_inter_matr([w[i][0]])
    print('w[0.5]=')
    disp_inter_matr([w[i][1]])
    print()

s = [0] + [1] * 10
s[6] = 2
print('T0 = %0.4f' % ky(s))
s = [1] * 2 + [0] * 9
print('Tu = %0.4f' % ky(s))
print()
for i in range(len(D)):
    print('Matrix', i)
    B = generate_matr(section(D[i], 0.5))
    v = [0] * len(B)
    for j in range(len(B)):
        v[j] = interval_weights(B[j], out=False)
    s = spectrum(v)
    for k in s:
        print(k)
        print('Ky = %0.4f' % ky(k))
    print()
wc = interval_weights(section(DC, 0.5), out=False)
w = [i[1] for i in w]
print('wc =')
disp_inter_matr([wc])
print('w =')
disp_inter_matr(w)
print()
w_glob = global_weights(w, wc)
print('w_glob = ')
disp_inter_matr([w_glob])
