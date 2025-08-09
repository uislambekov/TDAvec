from tdavec.tdavec_core import computePersistenceSilhouette
from tdavec import TDAvectorizer

import numpy as np
import timeit
import statistics as stats

import sys

def benchmark(func, repeats, v):
    time = timeit.repeat(lambda: func(v), number = 1, repeat = repeats)
    return {
        'min': min(time),
        'max': max(time),
        'mean': stats.mean(time),
        'median': stats.median(time),
        'stddev': stats.stdev(time),
        'n': n}
    return time

def myFunc(v):
    v.transform(homDim = 0)

vect_name = "ps"
if len(sys.argv) > 1:
    vect_name = sys.argv[1]

v = TDAvectorizer()
v.setParams({"scale":np.linspace(0, 2, 101), "output":vect_name})

out_str = "vect_name, n, min, max, mean, median, stddev\n"
for n in [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]:
    b = np.random.random(size = (n, 1))
    d = b + np.random.random(size = (n, 1))
    pd = np.array([b, d]).transpose()
    v.diags = [pd]
    

    res = benchmark(myFunc, repeats = 50, v = v)
    out_str += f"{vect_name}, {n}, {res['min']:.6f}, {res['max']:.6f}, {res['mean']:.6f}, {res['median']:.6f}, {res['stddev']:.6f}\n"

with open(f"benchmark_{vect_name}.csv", "w") as f:
    f.write(out_str)
print(f"Benchmark results saved to benchmark_{vect_name}.csv")    
# print(out_str)



