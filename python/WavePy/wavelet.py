'''
Created on 02/06/2012

@author: zenathar
'''
import numpy as np
from numpy import array 
import math

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
        data1 = subdata[len(subdata) / 2:,len(subdata[0]) / 2:]
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
