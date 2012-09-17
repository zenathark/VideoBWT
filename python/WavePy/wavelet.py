'''
Created on 02/06/2012

@author: zenathar
'''
import numpy as np
from numpy import array 
import math
import cv
import threading

def get_z_order(dim):
    mtx = []
    n = int(math.log(dim,2))
    pows = range(int(n/2))
    for i in range(dim):
        x = 0
        y = 0
        for j in pows:
            x |= ((i >> 2*j) & 1) << j 
            y |= ((i >> 2*j+1) & 1) << j 
        mtx += [vector((y,x))]
    return mtx


class wavelet(object):
    def __init__(self):
        print("echo")

class wavelet2D(object):
    '''
    classdocs
    '''
    
    level = 1
    w_type = float


    def __init__(self, data, level, w_type = float):
        '''
        Constructor
        '''
        if not isinstance(data, np.ndarray) or data.ndim != 2:
            raise TypeError, "data must be a 2D numpy array matrix"
        self.level = level
        self.w_type = w_type
        self.data = data
        
    def getHH(self):
        '''
        Returns the HH band of the wavelet
        '''
        data1 = self.data[0:len(self.data) / 2**self.level,0:len(self.data[0]) / 2 ** self.level]
        return data1

    def getLH(self,subdata):
        '''
        Returns the subband LH of a wavelet assuming one level of decomposition
        '''
        data1 = subdata[len(subdata) / 2:,:len(subdata[0]) / 2]
        return data1

    def getHL(self,subdata):
        '''
        Returns the subband HL of a wavelet assuming one level of decomposition
        '''
        data1 = subdata[:len(subdata) / 2,len(subdata[0]) / 2:]
        return data1

    def getLL(self,subdata):
        '''
        Returns the subband LH of a wavelet assuming one level of decomposition
        '''
        data1 = subdata[int(len(subdata) / 2):,int(len(subdata[0]) / 2):]
        return data1

    def getLH_n(self, level):
        '''
        Returns the LH subband of level "level" of the wavelet decomposition
        '''
        level = level - 1
        #calculate the number of rows and columns that the decomposition must have at level - 1
        rows = len(self.data) / 2 ** level
        cols = len(self.data[0]) / 2 ** level
        _subband = self.getLH(self.data[:rows, :cols])
        return _subband
    
    def getHL_n(self, level):
        '''
        Returns the HL subband of level "level" of the wavelet decomposition
        '''
        level = level - 1
        #calculate the number of rows and columns that the decomposition must have at level - 1
        rows = len(self.data) / 2 ** level
        cols = len(self.data[0]) / 2 ** level
        _subband = self.getHL(self.data[:rows, :cols])
        return _subband

    def getLL_n(self, level):
        '''
        Returns the LL subband of level "level" of the wavelet decomposition
        '''
        level = level - 1
        #calculate the number of rows and columns that the decomposition must have at level - 1
        rows = len(self.data) / 2 ** level
        cols = len(self.data[0]) / 2 ** level
        _subband = self.getLL(self.data[:rows, :cols])
        return _subband

    def show(self):
        temp = self.data
        #temporal wavelet copy
        w = self.data.copy()
        rows = len(w)
        cols = len(w[0])
        #scaling coefficients for display
        #HH first
        #hh = self.getHH()
        #m = w.min() 
        #M = w.max() - m
        #hh[:] = (hh - m) / M * 255
        m = w.min()
        M = w.max() - m
        w = (w - m) / M * 255
        #Now all bands
        for i in reversed(range(1,self.level+1)):
            #HL
            hl = w[rows/2**i:rows/2**(i-1),:cols/2**i]
            m = hl.min()
            M = hl.max() - m
            hl[:] = (hl - m) / M * 255
            #LH
            lh = w[:rows/2**i,cols/2**i:cols/2**(i-1)]
            m = lh.min()
            M = lh.max() - m
            lh[:] = (lh - m) / M * 255
            #HL
            hh = w[rows/2**i:rows/2**(i-1),cols/2**i:cols/2**(i-1)]
            m = hh.min()
            M = hh.max() - m
            hh[:] = (hh - m) / M * 255

        #change type to uint8
        w_ui8 = w.view(np.uint8)
        w_ui8_d = w_ui8[:,::2] 
        w_ui8_d[:] = w
        w = w_ui8_d
        self.data = temp
        #show with a thread
        wm = WindowManager(1)
        wm.img = cv.fromarray(w.copy())
        wm.name = "ASDADS"
        wm.start()
        wm.join()

class WindowManager(threading.Thread):

    img = 0
    name = "Empty"

    def __init__(self, num):
        threading.Thread.__init__(self)
        self.num = num

    def run(self):
        cv.ShowImage(self.name,self.img)
        while True:
            key = cv.WaitKey(0)
            if key == -1:
                cv.DestroyWindow(self.name)
                print self.name + " destroyed..."
                break

class vector(object):
    def __init__(self, data, entry_type = "-"):
        if isinstance(data,np.ndarray):
            self.data = data
        else:
            self.data = np.array(data)
        self.entry_type = type

    def __hash__(self):
        return hash(tuple(self.data))

    def __add__(self, other):
        if isinstance(other,np.ndarray):
            return vector(self.data + other.data)
        else:
            return vector(self.data + np.array(other))

    def __str__(self):
        return self.data.__str__()

    def __repr__(self):
        return self.data.__repr__()

    def __mul__(self,other):
        return vector(self.data * other)

    def __rmul__(self,other):
        return vector(self.data * other)

    def __lt__(self, other):
        if (isinstance(other,np.ndarray)):
            return np.all(self.data < other)
        else:
            return np.all(self.data < np.array(other))

    def tolist(self):
        return self.data.tolist()

    def __getitem__(self,index):
        return self.data[index]
