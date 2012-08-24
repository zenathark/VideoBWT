'''
Created on 03/06/2012

@author: zenathar
'''
import numpy as np
#from WavePy.lwt import *
#import WavePy.wavelet as wv
#reload(wv)

l = np.zeros((8,8))
for i in range(8):
    for j in range(8):
        l[i,j] = i * 8 + j

#a = wv.wavelet2D(l,1)
#print a.getLH(l)

'''
m = l.copy()

l = cdf53(l)
l = icdf53(l)
'''
#print(m - l)

