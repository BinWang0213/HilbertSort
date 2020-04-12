import pyHilbertSort as m
import numpy as np

assert m.__version__ == '0.0.1'
assert m.subtract(1, 2) == -1

pts_sort,idx_sort=m.hilbertSort(3, np.random.rand(5,3),True)
print(pts_sort,idx_sort,True)
