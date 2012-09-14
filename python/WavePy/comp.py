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
        output.append((wavelet.data[p[0]][p[1]] & 2 ** power) >> power)
    return output

def inv_bitplane_encoding(stream, LSP, wavelet, power):
    for p in LSP:
        bt = stream.pop()
        wavelet.data[p[0]][p[1]] = wavelet.data[p[0]][p[1]] | bt << power
        if not stream:
            break
    return wavelet

def bitplane_encoding_n(wavelet, varsize):
    output = []
    LSP = list(it.product(range(len(wavelet.data)) , range(len(wavelet.data[0])) ))
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
    D = get_D(ij,wavelet)
    O = D[:4]
    L = D[4:]
    return D, O, L

def has_offspring(ij,wavelet, typ = bool):
    if ij*2 < (len(wavelet.data),len(wavelet.data[0])):
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

    bpp = 0
    def __init__(self, wavelet):
        self.wavelet = wavelet

    def init(self):
        self.LSP = []
        rows = len(self.wavelet.data) / 2 ** (self.wavelet.level-1 )
        cols = len(self.wavelet.data[0]) / 2 ** (self.wavelet.level-1)
        self.LIP = get_z_order(rows*cols)
        self.LIS = get_z_order(rows*cols)[rows*cols/4:]
        for i in self.LIS:
            i.entry_type = "A"
        return

    def S_n(self, Tau, n):
        T = np.array([i.tolist() for i in Tau])
        return (abs(self.wavelet.data[T[:,0],T[:,1]]).max() >> int(n)) & 1

#outputs coefficient sign
    def out_sign(self, coeff):
        if self.wavelet.data[coeff[0],coeff[1]] > 0:
            return 0
        else:
            return 1

    def sorting(self, n):
        nextLSP = []                                                    
        removeLIP = []
        #Check for each significant pixel on the LIP
        c = -1
        for ij in self.LIP:
            c += 1
            out = self.S_n([ij], n)
            self.output_stream += [out]
            if out == 1:
                nextLSP += [ij]
                self.output_stream += [self.out_sign(ij)]
                removeLIP += [ij]
        #Remove new Significant pixels from LIP list
        #This is done after the output so the for loop doesnt have any problem
        for i in removeLIP:
            self.LIP.remove(i)
        remove_from_LIS = []
        for ij in self.LIS:
        #Check for zerotree roots (2.2.1)
            D, O, L = get_DOL(ij,self.wavelet)
            if ij.entry_type == 'A':
                out = self.S_n(D,n)
                self.output_stream += [out]
                if out == 1:
                    for kl in O:
                        out = self.S_n([kl],n)
                        self.output_stream += [out]
                        if out == 1:
                            nextLSP += [kl]
                            self.output_stream += [self.out_sign(kl)]
                        else:
                            self.LIP += [kl]
                #Check if ij has grandsons, if not remove from LIS
                    if L:
                        ij.entry_type = "B"
                        self.LIS += [ij]
                    else:
                        pass
                    remove_from_LIS += [ij]
            else: #Entry is type B
                out = self.S_n(L,n)
                self.output_stream += [out]
                if out == 1:
                    for k in O:
                        k.entry_type = "A"
                    self.LIS += O
                    remove_from_LIS += [ij]
        for i in remove_from_LIS:
            self.LIS.remove(i)
        remove_from_LIS = [] 
        return nextLSP 

    def compress(self):
        maxs = abs(self.wavelet.data)
        self.n = int(math.log(maxs.max(),2))
        n = self.n
        self.init()
        bit_bucket = self.bpp * len(self.wavelet.data) * len(self.wavelet.data[0])
        self.output_stream=buffer([],bit_bucket)
        #self.output_stream = []
        while n >= 0:
            try:
                newLSP = self.sorting(n)
                bitplane_encoding(self.wavelet,self.LSP,n,self.output_stream)
                self.LSP += newLSP
                n -= 1
            except NameError:
                return

    def inv_sorting(self,n):
        nextLSP = []
        #Fill each significant pixel
        removeLIP = []
        for ij in self.LIP:
            out = self.output_stream.pop()
            if out == 1:
                nextLSP += [ij]
                self.wavelet.data[ij[0],ij[1]] |= (1 << n)
                sign = self.output_stream.pop()
                if sign:
                    self.wavelet.data[ij[0],ij[1]] *= -1
                removeLIP += [ij]
        for i in removeLIP:
            self.LIP.remove(i)
        remove_from_LIS = []
        for ij in self.LIS:
            D, O, L = get_DOL(ij,self.wavelet)
            if ij.entry_type == "A":
                out = self.output_stream.pop()
                if out == 1:
                    for kl in O:
                        out = self.output_stream.pop()
                        if out == 1:
                            nextLSP += [kl]
                            sign = self.output_stream.pop()
                            self.wavelet.data[kl[0],kl[1]] |= (1 << n)
                            if sign:
                                self.wavelet.data[kl[0],kl[1]] *= -1
                        else:
                            self.LIP += [kl]
                    if L:
                        ij.entry_type = "B"
                        self.LIS += [ij]
                    else:
                        pass
                    remove_from_LIS += [ij]
            else:
                out = self.output_stream.pop() 
                if out == 1:
                    for i in O:
                        i.entry_type = "A"
                    self.LIS += O
                    remove_from_LIS += [ij]
        remove_from_LIS.reverse()
        for i in remove_from_LIS:
            self.LIS.remove(i)
        remove_from_LIS = []
        return nextLSP

    def uncompress(self):
        self.init()
        self.output_stream.reverse()
        n = self.n
        while n >= 0:
            try:
                newLSP = self.inv_sorting(n)
                inv_bitplane_encoding(self.output_stream,self.LSP,self.wavelet, n)
                self.LSP += newLSP
                n -= 1
            except IndexError:
                break
                    
class buffer(list):
    size = 1024
    
    def __init__(self, arg = [], size = 1024):
        super(buffer, self).__init__(arg)
        self.size = size

    def __iadd__(self, other):
        if len(self) <= self.size:
            return buffer(self + other,self.size)
        else:
            raise NameError("Stream full")

    def append(self,other):
        if len(self) <= self.size:
            super(buffer,self).append(other)
        else:
            raise NameError("Stream full")
