import numpy as np

L = 3
S = 2
for l in range(L):
    for s in np.arange(0, S, 0.5):
        jmin = min(abs(l-s), abs(l+s))
        jmax = max(abs(l-s), abs(l+s))
        for j in np.arange(jmin, jmax+1, 1):
            print("{},{},{}: {}".format(j,l,s, (j*(j+1)-l*(l+1)-s*(s+1))/2))
