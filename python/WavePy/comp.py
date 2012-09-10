'''
Created on 07/09/2012

@author: zenathar
'''

import numpy as np
from wavelet import *
import itertools as it
import math

if __name__ == '__main__':
    pass

def bitplane_encoding(wavelet, LSP, power, output):
    for p in LSP:
        output.append((wavelet[p[0]][p[1]] & 2 ** power) >> power)
    return output

def inv_bitplane_encoding(stream, LSP, wavelet, power):
    for p in LSP:
        bt = stream.pop()
        wavelet[p[0]][p[1]] = wavelet[p[0]][p[1]] | bt << power
        if not stream:
            break
    return wavelet

def bitplane_encoding_n(wavelet, varsize):
    output = []
    LSP = list(it.product(range(len(wavelet)) , range(len(wavelet[0])) ))
    for x in range(varsize):
        output = bitplane_encoding(wavelet, LSP, x, output)
    return output

def inv_bitplane_encoding_n(stream,rows,cols,varsize):
    wavelet = np.zeros((rows,cols),dtype = np.int16)
    LSP = list(it.product(range(rows),range(cols)))
    import pdb; pdb.set_trace();
    stream.reverse()
    for i in range(varsize):
        wavelet = inv_bitplane_encoding(stream, LSP, wavelet, i)
    return wavelet

def has_sons(ij,wavelet, typ = bool):
    if scalar_tuple(2,ij) < (len(wavelet.data),len(wavelet.data[0])):
        if typ == bool:
            return True
        else:
            return 1
    else:
        if typ == bool:
            return False
        else:
            return 0

def tuple_to_Array(a):
    return np.array([a[0],a[1]])

def add_tuple(a,b):
    return list(tuple_to_array(a),tuple_to_array(b))

def scalar_tuple(k,a):
    return list(k * tuple_to_array(a))

def fill_zerotree(root,wavelet,threshold):
    l = [root]
    c = 1
    while True:
        if has_sons(l[-c]):
            for i in range(c).reverse():
                son = get_sons(l[i])
                l += son[1]
                l += son[2]
                l += son[3]
                l += son[4]
            c *= 2
        else:
            break
    tree = {}
    l.reverse()
    for i in range(c):
        p = l.pop()
        tree[p] =  (((wavelet.data[p[0],p[1]] >> threshold) & 1) == 0)
    if not not l:
        for c in l:
            son = get_sons(c)
            tree[c] = (((wavelet.data[c[0],c[1]] >> threshold) & 1) == 0) && \\
                    tree[son[1]] && tree[son[2]] && tree[son[3]] && tree[son[4]]
    return tree

def is_zerotree(tree,ij):
    return tree[ij[0:2]] 

def get_sons(c):
    son1 += scalar_tuple(2,c)
    son2 += add_tuple(scalar_tuple(2,c), (1,0))
    son3 += add_tuple(scalar_tuple(2,c), (0,1))
    son4 += add_tuple(scalar_tuple(2,c), (1,0))
    return [son1,son2,son3,son4]
 
class SPIHT(object):
    '''
    '''

    output_stream = []
    def __init__(self, wavelet):
        self.wavelet = wavelet

    def init(self):
        maxs = abs(self.wavelet.data)
        self.n = math.log(maxs.max(),2)
        self.LSP = np.array([], np.int16)
        rows = len(self.wavelet.data) / 2 ** (self.wavelet.level-1 )
        cols = len(self.wavelet.data[0]) / 2 ** (self.wavelet.level-1)
        self.LIP = list(it.product(range(rows),range(cols)))
        self.LIS = list(it.product(range(rows/2),range(cols/2,cols),'A'))
#Append merges all sublists into one big vector, remember using LIS[1::2] for columns
#and LIS[0::2] for rows
        self.LIS += list(it.product(range(rows/2,rows),range(cols/2),'A'))
        self.LIS += list(it.product(range(rows/2,rows),range(cols/2,cols),'A'))
        return

    def S_n(self, Tau):
        T = np.append(Tau,np.array([],np.int16))
        return (abs(self.wavelet.data[T[0::2],T[1::2]]).max() >> int(self.n)) & 1

    def sorting(self):
        nextLSP = []
        for ij in self.LIP:
            out = self.S_n(ij)
            self.output_stream += out
            if out == 1:
                nextLSP += [ij]
        for ij in nextLSP:
            self.LIP.remove(ij)
        zero_tree = {}
        for ij in self.LIS:
            zero_tree = dict(zero_tree.items() + is_zerotree(ij,self.wavelet,self.n).items())
        for ij in self.LIS:
            if ij[3] == 'A':
                out = self.S_n(ij[0:2])
                self.output_stream += out


