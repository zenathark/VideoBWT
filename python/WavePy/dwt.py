import pywt
import wavelet as wv
import numpy as np

def dwt2(image, wavelet, mode, level):
    signal = np.asarray(image)
    if signal.ndim == 2:
        r = pack_wave_coeff(pywt.wavedec2(signal,wavelet,mode,level))
        return r 
    elif signal.ndim == 3:
        r,g,b = wv.splitRGB(image)
        rw = pack_wave_coeff(pywt.wavedec2(r,wavelet,mode,level))
        gw = pack_wave_coeff(pywt.wavedec2(g,wavelet,mode,level))
        bw = pack_wave_coeff(pywt.wavedec2(b,wavelet,mode,level))
        return (rw,gw,bw)

def pack_wave_coeff(coeff):
    mtx = np.zeros((2*len(coeff[-1][0]),2*len(coeff[-1][0][0])))
    rw = len(coeff[0])
    cl = len(coeff[0][0])
    mtx[:rw,:cl] = coeff[0].copy()
    y = rw
    x = cl
    for i in range(1,len(coeff)):
        rw = len(coeff[i][0])
        cl = len(coeff[i][0][0])
        mtx[:rw,x:x+cl] = coeff[i][0].copy()
        mtx[y:y+rw,:cl] = coeff[i][1].copy()
        mtx[y:y+rw,x:x+cl] = coeff[i][2].copy()
        y += rw
        x += cl
    wavelet = wv.wavelet2D(mtx,len(coeff) - 1)
    return wavelet
