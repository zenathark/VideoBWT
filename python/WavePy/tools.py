'''
Created on 18/04/2013

@author: JuanCarlos
'''
import scipy as np


class CircularStack(object):
    '''
    classdocs
    '''
    data = []
    size = 0
    index = 0
    tail = 0

    def __init__(self, size=0):
        '''
        Constructor
        '''
        if size > 0:
            self.data = [0] * size
            self.size = size
        self.index = 0
        self.tail = 0

    def push(self, data):
        next_tail = self._get_next_tail()
        if next_tail == self.index:
            raise OverflowError
        else:
            self.data[self.tail] = data
            self.tail = next_tail

    def pop(self):
        if self.index == self.tail:
            raise IndexError
        else:
            result = self.data[self.index]
            self.index = (self.index + 1) % self.size
            return result

    def _get_next_tail(self):
        return (self.tail + 1) % self.size


def quant(wavelet, delta):
    #iw = zeros((len(wavelet.data),len(wavelet.data[0])),int64)
    iw = np.array(np.trunc(wavelet.data / delta), np.int64)
    wavelet.data = iw
    return wavelet


def aquant(wavelet, delta):
    #iw = zeros((len(wavelet.data),len(wavelet.data[0])),int64)
    iw = np.array(np.trunc(wavelet / delta), dtype=np.uint8)
    return iw


def zero_padding(signal, squared=True):
    rows, cols = signal.shape
    pow_rows = int(np.ceil(np.log2(rows)))
    pow_cols = int(np.ceil(np.log2(cols)))
    if squared:
        if pow_cols > pow_rows:
            pow_rows = pow_cols
        else:
            pow_cols = pow_rows
    padded_signal = np.zeros((2 ** pow_rows, 2 ** pow_cols),
                             dtype=signal.dtype)
    y_0 = np.trunc((2 ** pow_rows - rows) / 2)
    y_t = y_0 + signal.shape[0]
    x_0 = np.trunc((2 ** pow_cols - cols) / 2)
    x_t = x_0 + signal.shape[1]
    padded_signal[y_0:y_t, x_0:x_t] = signal
    return padded_signal

def unpadding(signal, dim):
    y, x = dim
    sy, sx = signal.shape
    y0 = np.trunc((sy-y) / 2)
    x0 = np.trunc((sx-x) / 2)
    yt = y0 + y
    xt = x0 + x

    return signal[y0:yt,x0:xt]

