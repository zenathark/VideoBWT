'''
@author: zenathar
'''

import numpy as np
import math
import wavelet as wv

if __name__ == '__main__':
    pass

def cdf97(signal, level = 1,  in_place = True):
    if not isinstance(signal, np.ndarray) or signal.ndim != 2:
        raise TypeError, "Signal expected as 2D ndarray (numpy)"
    if not in_place:
        signal = signal.copy()
    if signal.dtype == np.uint8:
        sig_i8 = np.zeros((len(signal),2*len(signal[0])),np.uint8)
        sig_i8[:,0:-1:2] = signal
        sig_i16 = sig_i8.view(np.uint16)
        sig_f16 = sig_i16.view(np.float16)
        sig_f16[:] = sig_i16
        signal = sig_f16
    signal = normal_forward(signal,level,1/1.149604398,(-1.586134342,-0.05298011854,0.8829110762,0.4435068522))
    return wv.wavelet2D(signal,level,"cdf97")

def cdf53(signal, level = 1, in_place = True):
    if not isinstance(signal, np.ndarray) or signal.ndim != 2:
        raise TypeError, "Signal expected as 2D ndarray (numpy)"
    if not in_place:
        signal = signal.copy()
    if signal.dtype == np.uint8:
        sig_i8 = np.zeros((len(signal),2*len(signal[0])),np.uint8)
        sig_i8[:,0:-1:2] = signal
        sig_i16 = sig_i8.view(np.uint16)
        sig_f16 = sig_i16.view(np.float16)
        sig_f16[:] = sig_i16
        signal = sig_f16
    signal = normal_forward(signal,level,math.sqrt(2),(-0.5,0.25))
    return wv.wavelet2D(signal,level,"cdf53")

def icdf97(wave, in_place = True):
    if not isinstance(wave, wv.wavelet2D):
        raise TypeError, "Signal expected as wavelet2D"
    signal = wave.data.copy()
    signal.dtype = np.float32
    signal[:] = wave.data
    signal = normal_inverse(signal,wave.level,1.149604398,(-0.4435068522,-0.8829110762,0.05298011854,1.586134342))
    return signal

def icdf53(wave, in_place = True):
    if not isinstance(wave, wv.wavelet2D):
        raise TypeError, "Signal expected as wavelet2D"
    signal = wave.data.copy()
    signal.dtype = np.float32
    signal[:] = wave.data
    signal = normal_inverse(signal,wave.level,1/math.sqrt(2),(-0.25,0.5))
    return signal

def normal_forward(signal, level,  scale_coeff, coeff):
    decomposed_signal = signal
    for x in range(level):
        forward(decomposed_signal, scale_coeff, coeff)
        decomposed_signal = forward(decomposed_signal.T, scale_coeff, coeff).T
        updated_rows = int(len(decomposed_signal) / 2)
        updated_cols = int(len(decomposed_signal[0]) / 2)
        decomposed_signal = decomposed_signal[:updated_rows,:updated_cols]
    return signal

def normal_inverse(signal, level, scale_coeff, coeff):
    updated_rows = len(signal) / 2 **(level-1)
    updated_cols = len(signal[0]) / 2 **(level-1)
    for x in range(level):
        recomposed_signal = signal[:updated_rows,:updated_cols]
        recomposed_signal = inverse(recomposed_signal.T, scale_coeff, coeff).T
        recomposed_signal = inverse(recomposed_signal, scale_coeff, coeff)
        updated_rows = updated_rows * 2
        updated_cols = updated_cols * 2
    return signal

def forward(signal, scale_coeff, coeff):
    if not isinstance(signal, np.ndarray) or signal.ndim != 2:
        raise TypeError, "Signal expected as 2D ndarray (numpy)"
    if len(coeff) <= 0:
        raise TypeError, "Filter bank empty"
    if len(coeff) % 2 != 0:
        raise TypeError, "Filter bank expected to be Predict-Update pairs"
    for i in np.arange(0,len(coeff),2):
        #predict
        signal[:,1:-1:2] += coeff[i] * (signal[:,:-2:2] + signal[:,2::2])
        signal[:,-1] += 2 * coeff[i] * signal[:,-2]
        #update
        signal[:,2::2] += coeff[i+1] * (signal[:,1:-1:2] + signal[:,3::2])
        signal[:,0] += 2 * coeff[i+1] * signal[:,1]
    #scale
    signal[:,:2] *= scale_coeff
    signal[:,1::2] /= scale_coeff
    #sort lazy wavelet
    signal = sort(signal)
    return signal

def inverse(signal, scale_coeff, coeff):
    if not isinstance(signal, np.ndarray) or signal.ndim != 2:
        raise TypeError, "Signal expected as 2D ndarray (numpy)"
    if len(coeff) <= 0:
        raise TypeError, "Filter bank empty"
    if len(coeff) % 2 != 0:
        raise TypeError, "Filter bank expected to be Predict-Update pairs"
    #undo scale
    signal = unsort(signal)
    signal[:,:2] *= scale_coeff
    signal[:,1::2] /= scale_coeff
    for i in np.arange(0,len(coeff),2):
        #Undo update
        signal[:,2::2] += coeff[i] * (signal[:,1:-1:2] + signal[:,3::2])
        signal[:,0] += 2 * coeff[i] * signal[:,1] 
        #Undo predict
        signal[:,1:-1:2] += coeff[i+1] * (signal[:,:-2:2] + signal[:,2::2])
        signal[:,-1] += 2 * coeff[i+1] * signal[:,-2]
    return signal

def sort(signal):
    to = signal.shape[1]
    for i in range(1, to/2+1, 1):
        temp = signal[:,i].copy()
        signal[:,i:to-1] = signal[:,i+1:to]
        signal[:,-1] = temp
    return signal
        
def unsort(signal):
    to = signal.shape[1]
    for i in range(to/2, 0, -1):
        temp = signal[:,-1].copy()
        signal[:,i+1:] = signal[:,i:-1]
        signal[:,i] = temp
    return signal
