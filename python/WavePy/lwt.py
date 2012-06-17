'''
Created on 03/06/2012

@author: zenathar
'''

import numpy as np

if __name__ == '__main__':
    pass

def cdf97(signal, in_place = True):
    if not isinstance(signal, np.ndarray) or signal.ndim != 2:
        raise TypeError, "Signal expected as 2D ndarray (numpy)"
    if not in_place:
        signal = signal.copy()
    signal = forward(signal,1/1.149604398,-1.586134342,-0.05298011854,0.8829110762,0.4435068522)
    return forward(signal.T,1/1.149604398,-1.586134342,-0.05298011854,0.8829110762,0.4435068522).T

def icdf97(signal, in_place = True):
    if not isinstance(signal, np.ndarray) or signal.ndim != 2:
        raise TypeError, "Signal expected as 2D ndarray (numpy)"
    if not in_place:
        signal = signal.copy()
    signal = inverse(signal.T,1.149604398,-0.4435068522,-0.8829110762,0.05298011854,1.586134342).T
    return inverse(signal,1.149604398,-0.4435068522,-0.8829110762,0.05298011854,1.586134342)

def forward(signal, scale_coeff, *coeff):
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

def inverse(signal, scale_coeff, *coeff):
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