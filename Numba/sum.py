# Compares sum vs pure sum
from numba import d
from numba.decorators import jit as jit

def sum(arr):
    M, N = arr.shape
    result = 0.0
    for i in range(M):
        for j in range(N):
            result += arr[i,j]
    return result

csum = jit(ret_type=d, arg_types=[d[:,:]])(sum)

from numpy import random
arr = random.randn(1000,1000)

import time
start = time.time()

## Bad python
res=0
start = time.time()
res = sum(arr)
duration = time.time() - start
print "Result from 'bad' python is %s in %s (msec)" % (res, duration*1000)
raw_input("press enter")

## Normal C based sum
start = time.time()
res = arr.sum()
duration2 = time.time() - start
print "Result from numpy (Cbased) is %s in %s (msec)" % (res, duration2*1000)
print "Speed up is %s" % (duration/duration2)
raw_input("press enter")


## Numba
start = time.time()
res = csum(arr)
duration3 = time.time() - start
print "Result from compiled is %s in %s (msec)" % (res, duration3*1000)

print "Speed up is %s (vs bad python) %s (vs Cbased)" % (duration / duration3, duration2 / duration3)
