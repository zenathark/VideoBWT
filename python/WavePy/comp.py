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
    stream.reverse()
    for i in range(varsize):
        wavelet = inv_bitplane_encoding(stream, LSP, wavelet, i)
    return wavelet

#Get descendants, Pearlman named this set D
def get_D(ij,wavelet):
    D = []
    O = get_O(ij,wavelet)
    D += O
    window_size = 4
    while O:
    #while np.all((4 * root) < max_root):
        reverse_D = D[::-1]
        for i in range(window_size):
            O = get_O(reverse_D[i],wavelet)
            D += O
        window_size *= 4 #Window size increases as rows x 2 and cols x 2
    return D

#Get offsprings of ij, Pearlman named this set O
def get_O(ij,wavelet):
    O = []
    if has_offspring(ij,wavelet):
        O += [2*ij, 2*ij+(1,0), 2*ij+(0,1), 2*ij+(1,1)]
    return O

#Get all subsets of ij decendants O,D and L
def get_DOL(ij,wavelet):
    O = get_O(ij,wavelet)
    D = get_D(ij,wavelet)
    L = D[4:]
    return D, O, L

def has_offspring(ij,wavelet, typ = bool):
    if np.all(2*ij < (len(wavelet.data),len(wavelet.data[0]))):
        if typ == bool:
            return True
        else:
            return 1
    else:
        if typ == bool:
            return False
        else:
            return 1

def fill_zerotree(root,wavelet,threshold):
    l = [root]
    c = 1
    while True:
        if has_offspring(l[-c]):
            for i in range(c).reverse():
                son = get_sons(l[i])
                l += son[1]
                l += son[2]
                l += son[3]
                l += son[4]
            c *= 2
        else:
            break
    tree = tuple(root)
    l.reverse()
    for i in range(c):
        p = l.pop()
        tree[p] =  (((wavelet.data[p[0],p[1]] >> threshold) & 1) == 0)
    if not not l:
        for c in l:
            son = get_offspring(c)
            tree[c] = ((((wavelet.data[c[0],c[1]] >> threshold) & 1) == 0) & tree[son[1]] & tree[son[2]] & tree[son[3]] & tree[son[4]])
    return tree

def is_zerotree(tree,ij):
    return tree[ij[0:2]] 
 
class SPIHT(object):
    '''
    '''

    output_stream = []
    def __init__(self, wavelet):
        self.wavelet = wavelet

    def init(self):
        maxs = abs(self.wavelet.data)
        self.n = int(math.log(maxs.max(),2))
        self.LSP = []
        rows = len(self.wavelet.data) / 2 ** (self.wavelet.level-1 )
        cols = len(self.wavelet.data[0]) / 2 ** (self.wavelet.level-1)
        self.LIP = []
        for i in  list(it.product(range(rows),range(cols))):
            self.LIP += [np.array(i)]
        LIS = []
        LIS = list(it.product(range(rows/2),range(cols/2,cols)))
        LIS += list(it.product(range(rows/2,rows),range(cols/2)))
        LIS += list(it.product(range(rows/2,rows),range(cols/2,cols)))
        self.LIS = []
        for i in LIS:
            self.LIS += [np.array(i)]
        self.LIS_type = ["A" for i in range(len(LIS))]
        return

    def S_n(self, Tau):
        T = np.array([i.tolist() for i in Tau])
        return (abs(self.wavelet.data[T[:,0],T[:,1]]).max() >> int(self.n)) & 1

#outputs coefficient sign
    def out_sign(self, coeff):
        if self.wavelet.data[coeff[0],coeff[1]] > 0:
            return 0
        else:
            return 1

    def sorting(self):
        nextLSP = []                                                    
        #Check for each significant pixel on the LIP
        for ij in self.LIP:
            out = self.S_n([ij])
            self.output_stream += [out]
            if out == 1:
                nextLSP += [ij]
                self.ouput_stream += [output_sign(ij)]
        #Remove new Significant pixels from LIP list
        #This is done after the output so the for loop doesnt have any problem
        for ij in nextLSP:
            self.LIP.remove(ij)
        remove_from_LIS = []
        c = -1
        for ij in self.LIS:
            c+=1
        #Check for zerotree roots (2.2.1)
            D, O, L = get_DOL(ij,self.wavelet)
            if self.LIS_type == 'A':
                out = self.S_n(D)
                self.output_stream += out
                if out == 1:
                    for kl in O:
                        out = self.S_n([kl])
                        self.output_stream += [out]
                        if out == 1:
                            nextLSP += [kl]
                            self.output_stream += [output_sign(kl)]
                        else:
                            self.LIP += [o]
                #Check if ij has grandsons, if not remove from LIS
                    if L:
                        self.LIS += [ij]
                        self.LIS_type += "B"
                    else:
                        pass
                    remove_from_LIS += [c]
            else: #Entry is type B
                out = self.S_n(L)
                self.output_stream += [out]
                if out == 1:
                    self.LIS += O
                    self.LIS_type += ["A" for i in range(4)]
                    remove_from_LIS += [c]
        return remove_from_LIS

    def compress(self):
        self.init()
        r = range(self.n + 1)
        r.reverse()
        for i in r:
            erase = self.sorting()
            self.output_stream += bitplane_encoding_n(self.wavelet,i)
            for c in erase:
                self.LIS.pop(c)
                self.LIS_type.pop(c)

