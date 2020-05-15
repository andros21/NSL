import numpy as np
import pandas as pd
import os


def error(AV, AV2, n):
    if n == 0:
        return 0
    else:
        return np.sqrt((AV2[n] - AV[n]**2)/n)


def blocking(vct, L):
    M = 5e5
    N = int(M/L)
    ave = np.zeros(N)

    for i in range(N):
        sum = 0
        for j in range(L):
            k = j+i*L
            sum += vct[k]
        ave[i] = sum/L  # r_i

    ave = np.asarray(ave)

    avg = np.sum(ave)/N
    err = np.sqrt((np.sum(ave**2)/N - avg**2)/N)

    return (avg, err)


Lrange = [10, 25, 50, 100, 250, 500, 1000, 2500, 5000, 10000]
states = ['solid', 'liquid', 'gas']
vars = "PV"
result = dict()

os.chdir('/home/rox/bktemp')

for state in states:
    print(state)
    result[state] = dict()
    for var in vars:
        print(var)
        result[state][var] = dict()
        result[state][var]['avg'] = list()
        result[state][var]['err'] = list()
        vct = np.asarray(pd.read_csv(f"{var}-{state}.out", header=None))
        for L in Lrange:
            print(L)
            avg = blocking(vct, L)
            result[state][var]['avg'].append(avg[0])
            result[state][var]['err'].append(avg[1])

for state in states:
    for var in vars:
        df = pd.DataFrame(
            {'avg': result[state][var]['avg'],
             'err': result[state][var]['err']})
        df.to_csv(f"smart-{var}-{state}-blk-evo.csv", header=None, index=False)
