import WavePy.wavelet as wv
import WavePy.spiht as spiht
from WavePy.lwt import *
import pickle
from pylab import *

lim = [
        #"cameraman", 
        #"house", 
        #"jetplane", 
        #"lake", 
        #"lena_gray_512", 
        #"livingroom", 
        #"mandril_gray", 
    "peppers_gray"#, 
    #"pirate", 
    #"walkbridge", 
    #"woman_blonde", 
    #"woman_darkhair"
    ]
path = "/home/zenathar/Pictures/test/standard_test_images/{0}.tif"
output = "/home/zenathar/Pictures/test/spiht/fvht_{0}.spti"

#fil = open("/home/zenathar/Documents/src/VideoBWT/python/512x512_5.dol","r+")
#dct = pickle.load(fil)
first = False 
for bpp in linspace(0,8,100):
    #if bpp >2.9090909090:
    #    first = False
    if not first:
        for i in lim:
            r,g,b = wv.LoadImageRGB(path.format(i))
            k = spiht.fvht_image_pack([r],"cdf97",5,(256,256),bpp,0.06,0.3,1,1,"bi.orth",0.0001, True,i,{},True)
            outputfile = open(output.format(i) + "_" + str(bpp),"w+")
            pickle.dump(k["test_data"],outputfile)
            outputfile.close()
            print path.format(i) + str(bpp)
