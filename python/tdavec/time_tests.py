from tdavec.tdavec_core import computePersistenceSilhouette
from ripser import ripser
import numpy as np
import timeit
import statistics as stats


def benchmark(func, repeats, pd, scaleSeq):
    time = timeit.repeat(lambda: func(pd, scaleSeq), number = 1, repeat = repeats)
    return {
        'min': min(time),
        'max': max(time),
        'mean': stats.mean(time),
        'median': stats.median(time),
        'stddev': stats.stdev(time),
        'n': n}
    return time

def myFunc(pd, scaleSeq):
    computePersistenceSilhouette(pd, 0, scaleSeq)


print("n, min, max, mean, median, stddev")
for n in [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]:
    b = np.random.random(size = (n, 1))
    d = b + np.random.random(size = (n, 1))
    pd = np.array([b, d]).transpose()
    scaleSeq = np.linspace(0, 2, 101)

    res = benchmark(myFunc, repeats = 10, pd = pd, scaleSeq = scaleSeq)
    print(f"{n}, {res['min']:.6f}, {res['max']:.6f}, {res['mean']:.6f}, {res['median']:.6f}, {res['stddev']:.6f}")



