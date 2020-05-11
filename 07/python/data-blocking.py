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
    av2 = np.zeros(N)
    sum_prog = np.zeros(N)
    su2_prog = np.zeros(N)
    err_prog = np.zeros(N)

    for i in range(N):
        sum = 0
        for j in range(L):
            k = j+i*L
            sum += vct[k]
        ave[i] = sum/L  # r_i
        av2[i] = (ave[i])**2  # (r_i)^2

    for i in range(N):
        for j in range(i+1):
            sum_prog[i] += ave[j]  # SUM_{j=0,i} r_j
            su2_prog[i] += av2[j]  # SUM_{j=0,i} (r_j)^2
        sum_prog[i] /= (i+1)  # Cumulative average
        su2_prog[i] /= (i+1)  # Cumulative square average
        err_prog[i] = error(sum_prog, su2_prog, i)  # Statistical uncertainty

    return (sum_prog[-1], err_prog[-1])


Lrange = [10, 25, 50, 100, 250, 500, 1000, 2500, 5000]
states = ['solid', 'liquid', 'gas']
vars = "PV"
result = dict()

os.chdir('../out/')

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
        df.to_csv(f"{var}-{state}-blk-evo.csv", header=None, index=False)
