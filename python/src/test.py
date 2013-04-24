'''
Created on 03/06/2012

@author: zenathar
'''
from scipy import *
from cv2 import *
from WavePy.lwt import *
from WavePy.tools import *
from WavePy.speck import *

im = imread('d:/standard_test_images/cameraman.tif',IMREAD_GRAYSCALE)
w = cdf97(im,4,False)
q = quant(w,0.001)
sp = speck()
sp.compress(q, 8)

