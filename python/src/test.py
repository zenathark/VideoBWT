'''
Created on 03/06/2012

@author: zenathar
'''
import numpy as np
from WavePy.lwt import *

l = np.zeros((32,32))
for i in range(0,32):
   # l[:,i] = 5+i+0.4*i*i-0.02*i*i*i
    l[:,i] = i
m = l.copy()

l = cdf53(l)
l = icdf53(l)

print m - l

