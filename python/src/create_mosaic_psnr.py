import pickle
from pylab import *
from WavePy.lwt import *
import math
import cv


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
fvhtoutput = "/home/zenathar/Pictures/test/spiht/fvht_{0}.spti_{1}"
spihtoutput = "/home/zenathar/Pictures/test/spiht/{0}.spti{1}"
path = "/home/zenathar/Pictures/test/standard_test_images/{0}.tif"
psnrc = 20 * math.log(255,10)
fvlm = {}
for i in lim:
    fvlm[i] = []
#import pdb; pdb.set_trace()
center = 256
siz = 76
for i in lim:
    spfil = open(spihtoutput.format(i,"4.0404040404"),"r")#"1.05050505051"),"r")
    sp = pickle.load(spfil)
    fvfil = open(fvhtoutput.format(i,"4.0404040404"),"r")#"1.05050505051"),"r")
    fv = pickle.load(fvfil)
    spfil.close()
    fvfil.close()
    im = cv.LoadImageM(path.format(i))
    im = asarray(im)[:,:,0]
    fv.data *= 0.0001
    sp.data *= 0.00001
    fv.level = 5
    sp.level = 5
    ifv = icdf97(fv,False)
    isp = icdf97(sp,False)
    for size in linspace(76,255):
        MSEifv = sum((im[center-size:center+size,center-size:center+size]- ifv[center-size:center+size,center-size:center+size]) ** 2) / ((size*2)**2)
        MSEisp = sum((im[center-size:center+size,center-size:center+size]- isp[center-size:center+size,center-size:center+size]) ** 2) / ((size*2)**2)
        PSNRifv = psnrc - 10 * math.log(MSEifv,10)
        PSNRisp = psnrc - 10 * math.log(MSEisp,10)
        fvlm[i] += [(size,PSNRisp,PSNRifv)]
        print i + str(size)

