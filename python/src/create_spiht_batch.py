import WavePy.wavelet as wv
import WavePy.spiht as spiht
from WavePy.lwt import *
import pickle
from pylab import *

lim = [ "cameraman", 
    "house", 
    "jetplane", 
    "lake", 
    "lena_gray_512", 
    "livingroom", 
    "mandril_gray", 
    "peppers_gray", 
    "pirate", 
    "walkbridge", 
    "woman_blonde", 
    "woman_darkhair"]
path = "/home/zenathar/Pictures/test/standard_test_images/{0}.tif"
output = "/home/zenathar/Pictures/test/spiht/{0}.spt"

fil = open("/home/zenathar/Documents/src/VideoBWT/python/512x512_5.dol","r+")
dct = pickle.load(fil)
first = True
for bpp in linspace(0,8,100):
    if bpp >= 1.85858585859:
        first = False
    if not first:
        for i in lim:
            r,g,b = wv.LoadImageRGB(path.format(i))
            k = spiht.spiht_image_pack([r],"cdf97",5,[bpp],"bi.orth",0.01, True,i,dct,True)
            outputfile = open(output.format(i) + str(bpp),"w+")
            pickle.dump(k["test_data"],outputfile)
            outputfile.close()
            print path.format(i) + str(bpp)
