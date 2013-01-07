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
fvl = {}
banned = []
for i in lim:
    fvl[i] = []
#import pdb; pdb.set_trace()
for bpp in linspace(0,8,100):
    if bpp != 0:
        for i in lim:
            try:
                spfil = open(spihtoutput.format(i,str(bpp)),"r")
                sp = pickle.load(spfil)
                fvfil = open(fvhtoutput.format(i,str(bpp)),"r")
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
                MSEifv = sum((im - ifv) ** 2) / (512. * 512.)
                MSEisp = sum((im - isp) ** 2) / (512. * 512.)
                PSNRifv = psnrc - 10 * math.log(MSEifv,10)
                PSNRisp = psnrc - 10 * math.log(MSEisp,10)
                fvl[i] += [(bpp,PSNRisp,PSNRifv)]
                print i + str(bpp)
            except:
                banned += [bpp]
                print("BANNED!!" + str(bpp))
                break





