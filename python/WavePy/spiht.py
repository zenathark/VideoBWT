'''
Created on Sep 28, 2012

@author: jcg112
'''

from collections import deque
import numpy as np
import wavelet as wv
import scipy.weave as weave
from scipy.weave import converters
import pickle
import lwt
import dwt
import math
import cv

def spiht_image_pack(img, wavename, level, bpp, mode = "bi.orth", delta = 0.01, display_progress = True, str_pr = "", d = {}, handle = True):
    if not isinstance(img,tuple):
        img = (img)
        bpp = (bpp)
    ch_stream = {
            "size": 0,
            "channels":len(img),
            "bpp":bpp,
            "quant_delta":delta,
            "wave_type": wavename,
            "mode":mode,
            "decomp_level":level,
            "wise_bit":[],
            "payload":[]
            }
    for i in range(len(img)):
        m = np.asarray(img[i])
        if mode == "bi.orth":
            if wavename == "cdf97": 
                w = lwt.cdf97(m,level,False)
            elif wavename == "cdf53":
                w = lwt.cdf53(m,level,False)
        else:
            w = dwt.dwt2(m,wavename,mode, level)
        stream = spiht_pack(w,bpp[i],delta,display_progress,str_pr + "[channel "+str(i)+"]",d,handle)
        ch_stream["payload"] += stream["payload"]
        ch_stream["wave_type"] = w.name
        ch_stream["size"] += stream["size"]
        ch_stream["wise_bit"] += [stream["wise_bit"]]
        ch_stream["rows"] = stream["rows"]
        ch_stream["cols"] = stream["cols"]
        ch_stream["test_data"] = stream["test_data"]
    return ch_stream

def spiht_image_unpack(frame, display_progress = True, str_pr = "", d = {}, handle = True):
    rgb = []
    for i in range(frame["channels"]):
        pack = {
            "rows":frame["rows"],
            "cols":frame["cols"],
            "channels": 1,
            "wave_type":frame["wave_type"],
            "bpp":frame["bpp"],
            "quant_delta":frame["quant_delta"],
            "decomp_level":frame["decomp_level"],
            "wise_bit":frame["wise_bit"][i],
            "payload":frame["payload"][i]
            }
        wave = spiht_unpack(pack,True,str_pr + "[channel "+str(i)+"]", d, handle)
        if frame["wave_type"] == 'cdf97':
            ch = lwt.icdf97(wave,False)
        elif frame["wave_type"] == 'cdf53':
            ch = lwt.icdf53(wave,False)
        else:
            ch = dwt.dwt2(wave,frame["wave_type"],frame["decomp_level"],frame["mode"])
        ch = ch - ch.min()
        ch = ch / ch.max() * 255
        ch_i = np.zeros((len(ch),len(ch[0])),np.uint8)
        ch_i[:] = ch
        rgb += [cv.fromarray(ch_i)]
    return tuple(rgb)

def spiht_pack(wave,bpp,delta = 0.01, display_progress=True, str_pr = "", d = {}, handle = True):
    """Compresses a wavelet with SPIHT.

    Runs the original SPIHT algorithm from Dr. Pearlsman paper over the 
    given wavelet.

    Args:
        wavelet: A wavelet to be compressed, must be wavelet2D data type
        bpp: Bits per pixel compression ratio
        delta: Quantization delta if wavelet data is on floating point. 
            This leads to lossy compression always if coefficients are not 
            int
        display_progress: True if you want to print on screen algorithm 
                    progress

    Returns:
        A dictionary with a structure with the compressed data and some 
            decoding info. For example:
        {
            "size":0,
            "wave_type":"db7",
            "bpp":0.5,
            "quant_delta":0.001,
            "payload":[]
        }
    """
    #filename = str(len(wave.data)) + "x" + str(len(wave.data[0])) + "_" + str(wave.level) +".dol"
    #if not d:
    #    try:
    #        f = open(filename,"r")
    #        d = pickle.load(f)
    #        f.close()
    #        handle = True
    #    except IOError as e:
    #        d = {}
    #        handle = True
    #update = len(d)
    codec = SPIHT()
    codec.wavelet = wave
    codec.bpp = bpp
    codec.delta = delta
    codec.check_floating_point()
    codec.str_pr = str_pr
    codec.display_progress = display_progress
    codec.d_memory = d
    codec.compress()
    stream = list(codec.output_stream)[:int(bpp*wave.rows*wave.cols)]
    pack = {
            "size":len(stream),
            "rows":len(wave.data),
            "cols":len(wave.data[0]),
            "channels": 1,
            "wave_type":wave.name,
            "bpp":(bpp),
            "quant_delta":delta,
            "decomp_level":wave.level,
            "wise_bit":codec.n,
            "payload":[stream],
            "test_data":codec.test_data
            }
    #if len(codec.d_memory) > update and handle:
    #    f = open(filename,"w+")
    #    pickle.dump(codec.d_memory,f)
    #    f.close()
    return pack

def spiht_unpack(frame, display_progress=True, str_pr = "",d = {}, handle = True):
    """De-compresses a wavelet with SPIHT.

    Runs the original SPIHT algorithm from Dr. Pearlsman paper over the 
    given wavelet.

    Args:
        frame: A previously compressed frame
        display_progress: True if you want to print on screen algorithm 
                    progress

    Returns:
        An uncompressed wavelet2D instance 
    """
    filename = str(frame["rows"]) + "x" + str(frame["cols"]) + "_" + str(frame["decomp_level"]) +".dol" 
    if not d:
        try:
            f = open(filename,"r")
            d = pickle.load(f)
            f.close()
            handle = True
        except IOError as e:
            d = {}
            handle = True
    update = len(d)
    codec = SPIHT()
    data = np.zeros((frame["rows"],frame["cols"]),np.int32)
    codec.wavelet = wv.wavelet2D(data,frame["decomp_level"],frame["wave_type"])
    codec.delta = frame["quant_delta"]
    codec.str_pr = str_pr
    codec.n = frame["wise_bit"]
    codec.d_memory = d
    codec.output_stream = deque(frame["payload"])
    codec.uncompress()
    if len(codec.d_memory) > update and handle:
        f = open(filename,"w+")
       # pickle.dump(codec.d_memory,f)
        f.close()
    return codec.wavelet

class SPIHT(object):
    '''
    This class represents a SPIHT codec. It can compress and decompress
    data from Wavelet2D data type.

    This class is based on the original paper of Dr. Pearlsman
    '''

    #bpp resolution for compression
    bpp = 0
    #inner wavelet data type structure
    wavelet = 0
    #Quantization delta for floating point wavelets
    delta = 0.01
    #outputdata
    output_stream = 0
    
    zero_tree_data = 0

    #Extra info for progress display
    str_pr = ""
    display_progress = True
    d_memory = {}
    
    def __init__(self):
        pass
    
    def init(self):
        N = self.wavelet.rows * self.wavelet.cols
        self.LSP = deque([])
        rows = self.wavelet.rows / 2 ** (self.wavelet.level-1)
        cols = self.wavelet.cols / 2 ** (self.wavelet.level-1)
        self.LIP = wv.get_morton_order(rows*cols)
        self.LIS = wv.get_morton_order(rows*cols/4,rows*cols)
        for i in self.LIS:
            i.entry_type = "A"
        self.test_data = wv.wavelet2D(np.zeros((512,512)),self.wavelet.level)
        return
    
    def check_floating_point(self):
        ints = [np.int, np.int16, np.int32, np.int64]
        if not self.wavelet.data.dtype in ints:
            self.wavelet = quant(self.wavelet,self.delta)

    def S_n(self, Tau, n):
        T = np.array([i.tolist() for i in Tau])
        return int((abs(self.zero_tree_data[T[:,0],T[:,1]]).max() >= n))

    def bitplane_encoding(self, power, output):
#        if self.LSP:
#            T = np.array([i.tolist() for i in self.LSP])
#            stream = self.wavelet.data[T[:,0], T[:,1]] & 2 ** power >> power
#            self.output_stream.extend(list(stream))
        data = self.wavelet.data
        accum = self.test_data.data
        for p in self.LSP:
            if (accum[p[0],p[1]] + 2 ** power) <=  data[p[0],p[1]]:
                out = 1
                accum[p[0],p[1]] += 2 ** power
                exit = False
            else:
                out = 0
            #out = (self.wavelet.data[p[0],p[1]] & (2 ** power))
            #self.test_data.data[p[0],p[1]] |= out
            #out = out >> power
            self.output_stream.append(out)
            
            
        
    def inv_bitplane_encoding(self, power):
        for p in self.LSP:
            bt = self.output_stream.popleft()
            self.wavelet.data[p[0],p[1]] = self.wavelet.data[p[0],p[1]] | (bt << power)
    
    def get_D(self, ij):
        try:
            return self.d_memory[ij.__hash__()]
        except KeyError:
            Di = []
            O = self.get_O(ij)
            Di += O
            window_size = 4
            while O:
                reverse_D = Di[::-1]
                for i in range(window_size):
                    O = self.get_O(reverse_D[i])
                    Di += O
                window_size *= 4 #Window size increases as rows x 2 and cols x 2
            self.d_memory[ij.__hash__()] = Di
            return Di
    #Get offsprings of ij, Pearlman named this set O
    def get_O(self, ij):
        if self.has_offspring(ij):
            return [2*ij, 2*ij+(1,0), 2*ij+(0,1), 2*ij+(1,1)]
        else:
            return []

    #Get all subsets of ij decendants O,D and L
    def get_DOL(self, ij):
        D = self.get_D(ij)
        O = D[:4]
        L = D[4:]
        return D, O, L

    def has_offspring(self, ij, typ = bool):
        if ij*2 < (len(self.wavelet.data),len(self.wavelet.data[0])):
            if typ == bool:
                return True
            else:
                return 1
        else:
            if typ == bool:
                return False
            else:
                return 1

    def has_grandchilds(self, ij, typ = bool):
        if ij*4 < (len(self.wavelet.data),len(self.wavelet.data[0])):
            if typ == bool:
                return True
            else:
                return 1
        else:
            if typ == bool:
                return False
            else:
                return 1

    #outputs coefficient sign
    def out_sign(self, coeff):
        if self.wavelet.data[coeff[0],coeff[1]] > 0:
            return 0
        else:
            return 1
        
    def sorting(self, n):
        nextLSP = deque()
        newLIP = deque()
        data = self.wavelet.data
        zt_matrix = self.zero_tree_data
        #Check for each significant pixel on the LIP
        for ij in self.LIP:
            #out = self.S_n([ij], n)
            out = int((abs(data[ij[0],ij[1]]) & 2**n) > 0)
            self.output_stream.append(out)
            if out == 1:
                self.test_data.data[ij[0],ij[1]] = 2**n
                nextLSP.append(ij)
                self.output_stream.append(self.out_sign(ij))
                if self.out_sign(ij) == 1:
                    self.test_data.data[ij[0],ij[1]] *= -1
            else:
                newLIP.append(ij)
        #Update LIP
        self.LIP = newLIP
        newLIS = deque([],self.wavelet.rows*self.wavelet.cols)
        while self.LIS:
            ij = self.LIS.popleft()
        #Check for zerotree roots (2.2.1)
            #D, O, L = self.get_DOL(ij)
            O = self.get_O(ij)
            if ij.entry_type == 'A':
                #out = self.S_n(D,n)
                out = int((zt_matrix[ij[0],ij[1]] >= n))
                self.output_stream.append(out)
                if out == 1:
                    for kl in O:
                        #out = self.S_n([kl],n)
                        out = int((abs(data[kl[0],kl[1]]) & 2**n) > 0)
                        self.output_stream.append(out)
                        if out == 1:
                            nextLSP.append(kl)
                            self.test_data.data[kl[0],kl[1]] = 2**n
                            self.output_stream.append(self.out_sign(kl))
                            if self.out_sign(kl) == 1:
                                self.test_data.data[kl[0],kl[1]] *= -1
                        else:
                            self.LIP.append(kl)
                #Check if ij has grandsons, if not remove from LIS
                    if self.has_grandchilds(ij):
                        ij.entry_type = "B"
                        self.LIS.append(ij)
                else:
                    newLIS.append(ij)
            else: #Entry is type B
                out = self.S_n(O,n)
                self.output_stream += [out]
                if out == 1:
                    for k in O:
                        k.entry_type = "A"
                        self.LIS.append(k)
                else:
                    newLIS.append(ij) 
        self.LIS = newLIS
        return nextLSP
    
    def resolve_zero_tree(self):
        data = self.wavelet.data
        rows = self.wavelet.rows
        cols = self.wavelet.cols
        maxr = rows
        maxc = cols
        levels = self.wavelet.level
        zt_matrix = np.zeros((rows,cols),np.uint32)
        self.zero_tree_data = zt_matrix
        code = """
                #line 393 "spiht.py"
                int ix,iy;
                for(int i = 0; i < levels; ++i) {
                    for(int c = (int)(cols / 2); c < cols; ++c) {
                        for(int r = (int)(rows / 2); r < rows; ++r) {
               zt_matrix(r,c) = floor(log(abs(data(r,c)))/log(2));
                            ix = c * 2;
                            iy = r * 2;
                            if (ix < maxc && iy < maxr) {
                                if(zt_matrix(iy,ix)>zt_matrix(r,c))
                                    zt_matrix(r,c) = zt_matrix(iy,ix);
                                ++iy;
                                if(zt_matrix(iy,ix)>zt_matrix(r,c))
                                    zt_matrix(r,c) = zt_matrix(iy,ix);
                                ++ix;
                                if(zt_matrix(iy,ix)>zt_matrix(r,c))
                                    zt_matrix(r,c) = zt_matrix(iy,ix);
                                --iy;
                                if(zt_matrix(iy,ix)>zt_matrix(r,c))
                                    zt_matrix(r,c) = zt_matrix(iy,ix);
                            }
                        }
                    }
                    for(int c = 0; c < (int)(cols / 2); ++c) {
                        for(int r = (int)(rows / 2); r < rows; ++r) {
               zt_matrix(r,c) = floor(log(abs(data(r,c)))/log(2));
                            ix = c * 2;
                            iy = r * 2;
                            if (ix < maxc && iy < maxr) {
                                if(zt_matrix(iy,ix)>zt_matrix(r,c))
                                    zt_matrix(r,c) = zt_matrix(iy,ix);
                                ++iy;
                                if(zt_matrix(iy,ix)>zt_matrix(r,c))
                                    zt_matrix(r,c) = zt_matrix(iy,ix);
                                ++ix;
                                if(zt_matrix(iy,ix)>zt_matrix(r,c))
                                    zt_matrix(r,c) = zt_matrix(iy,ix);
                                --iy;
                                if(zt_matrix(iy,ix)>zt_matrix(r,c))
                                    zt_matrix(r,c) = zt_matrix(iy,ix);
                            }
                        }
                    }
                    for(int c = (int)(cols / 2); c < cols; ++c) {
                        for(int r = 0; r < (int)(rows / 2); ++r) {
               zt_matrix(r,c) = floor(log(abs(data(r,c)))/log(2));
                            ix = c * 2;
                            iy = r * 2;
                            if (ix < maxc && iy < maxr) {
                                if(zt_matrix(iy,ix)>zt_matrix(r,c))
                                    zt_matrix(r,c) = zt_matrix(iy,ix);
                                ++iy;
                                if(zt_matrix(iy,ix)>zt_matrix(r,c))
                                    zt_matrix(r,c) = zt_matrix(iy,ix);
                                ++ix;
                                if(zt_matrix(iy,ix)>zt_matrix(r,c))
                                    zt_matrix(r,c) = zt_matrix(iy,ix);
                                --iy;
                                if(zt_matrix(iy,ix)>zt_matrix(r,c))
                                    zt_matrix(r,c) = zt_matrix(iy,ix);
                            }
                        }
                    }
                    rows /= 2;
                    cols /= 2;
                }
                for(int c = 0; c < rows; ++c) {
                    for(int r = 0; r < cols; ++r) {
               zt_matrix(r,c) = floor(log(abs(data(r,c)))/log(2));
                    }
                }

        """
        weave.inline(code, 
                     ['rows','cols','maxr','maxc','zt_matrix','data',
                     'levels'],
                     type_converters = converters.blitz,
                     compiler = 'gcc',
                     headers = ["<math.h>"])
        #for i in range(level):
        #    for c in range(cols / 2,cols):
        #        for r in range(rows / 2,rows):
        #            zt_matrix[r,c] = 2 ** int(math.log(abs(data[r,c]),2))
        #            O = self.get_O((r,c))
        #            for o in O:
        #                if zt_matrix[o[0],o[1]] > zt_matrix[r,c]:
        #                    zt_matrix[r,c] = zt_matrix[o[0],o[1]]
        #    for c in range(cols / 2):
        #        for r in range(rows / 2,rows):
        #            zt_matrix[r,c] = 2 ** int(math.log(max(abs(data)),2))
        #            O = self.get_O((r,c))
        #            for o in O:
        #                if zt_matrix[o[0],o[1]] > zt_matrix[r,c]:
        #                    zt_matrix[r,c] = zt_matrix[o[0],o[1]]
        #    for c in range(cols / 2,cols):
        #        for r in range(rows / 2):
        #            zt_matrix[r,c] = 2 ** int(math.log(max(abs(data)),2))
        #            O = self.get_O((r,c))
        #            for o in O:
        #                if zt_matrix[o[0],o[1]] > zt_matrix[r,c]:
        #                    zt_matrix[r,c] = zt_matrix[o[0],o[1]]
        #    rows /= 2
        #    cols /= 2

    def compress(self):
        maxs = abs(self.wavelet.data)
        self.n = int(math.log(maxs.max(),2))
        n = self.n
        self.init()
        self.resolve_zero_tree()
        bit_bucket = self.bpp * self.wavelet.rows * self.wavelet.cols
        self.output_stream=deque([],2*bit_bucket)
        while n >= 0 and len(self.output_stream) < bit_bucket:
            print "Sorting...{0}".format(len(self.output_stream)/float(bit_bucket))
            newLSP = self.sorting(n)
            print "Bitplane Encoding...{0}".format(len(self.output_stream)/float(bit_bucket))
            self.bitplane_encoding(n,self.output_stream)
            self.LSP.extend(newLSP)
            n -= 1
    
    def inv_sorting(self,n):
        nextLSP = deque()                                                    
        newLIP = deque()
        for ij in self.LIP:
            out = self.output_stream.popleft()
            if out == 1:
                nextLSP.append(ij)
                self.wavelet.data[ij[0],ij[1]] = (1 << n)
                sign = self.output_stream.popleft()
                if sign == 1:
                    self.wavelet.data[ij[0],ij[1]] *= -1
            else:
                newLIP.append(ij)
        self.LIP = newLIP
        newLIS = deque([],self.wavelet.rows*self.wavelet.cols)
        while self.LIS:
            ij = self.LIS.popleft()
            D, O, L = self.get_DOL(ij)
            if ij.entry_type == "A":
                out = self.output_stream.popleft()
                if out == 1:
                    for kl in O:
                        out = self.output_stream.popleft()
                        if out == 1:
                            nextLSP.append(kl)
                            sign = self.output_stream.popleft()
                            self.wavelet.data[kl[0],kl[1]] = (1 << n)
                            if sign == 1:
                                self.wavelet.data[kl[0],kl[1]] *= -1
                        else:
                            self.LIP.append(kl)
                    if L:
                        ij.entry_type = "B"
                        self.LIS.append(ij)
                else:
                    newLIS.append(ij)
            else:
                out = self.output_stream.popleft() 
                if out == 1:
                    for k in O:
                        k.entry_type = "A"
                        self.LIS.append(k)
                else:
                    newLIS.append(ij)
        self.LIS = newLIS
        return nextLSP

    def uncompress(self):
        self.init()
        n = self.n
        while n >= 0:
            try:
                newLSP = self.inv_sorting(n)
                self.inv_bitplane_encoding(n)
                self.LSP.extend(newLSP)
                n -= 1
            except IndexError:
                break

def fvht_image_pack(img, wavename, level, f_center, Lbpp, lbpp, alpha, c, gamma, mode = "bi.orth", delta = 0.01, display_progress = True, str_pr = "", d = {}, handle = True):
    if not isinstance(img,tuple):
        img = (img)
    ch_stream = {
            "size": 0,
            "channels":len(img),
            "bpp":Lbpp,
            "quant_delta":delta,
            "wave_type": wavename,
            "mode":mode,
            "decomp_level":level,
            "wise_bit":[],
            "Lbpp": Lbpp,
            "lbpp": lbpp,
            "alpha": alpha,
            "c": c,
            "gamma": gamma,
            "fovea_center": f_center,
            "payload":[]
            }
    for i in range(len(img)):
        m = np.asarray(img[i])
        if mode == "bi.orth":
            if wavename == "cdf97":
                w = lwt.cdf97(m,level,False)
            elif wavename == "cdf53":
                w = lwt.cdf53(m,level,False)
        else:
            w = dwt.dwt2(m,wavename,level,mode)
        stream = fvht_pack(w, f_center, Lbpp, lbpp, alpha, c, gamma, delta, display_progress,str_pr + "[channel "+str(i)+"]",d,handle)
        ch_stream["payload"] += stream["payload"]
        ch_stream["wave_type"] = w.name
        ch_stream["size"] += stream["size"]
        ch_stream["wise_bit"] += [stream["wise_bit"]]
        ch_stream["rows"] = stream["rows"]
        ch_stream["cols"] = stream["cols"]
        ch_stream["bits"] = stream["bits"]
        ch_stream["test_data"] = stream["test_data"]
    return ch_stream

def fvht_image_unpack(frame, display_progress = True, str_pr = "", d = {}, handle = True):
    rgb = []
    for i in range(frame["channels"]):
        pack = {
            "rows":frame["rows"],
            "cols":frame["cols"],
            "channels": 1,
            "wave_type":frame["wave_type"],
            "bpp":frame["bpp"],
            "quant_delta":frame["quant_delta"],
            "decomp_level":frame["decomp_level"],
            "wise_bit":frame["wise_bit"][i],
            "payload":frame["payload"][i],
            "Lbpp": frame["Lbpp"],
            "lbpp": frame["lbpp"],
            "alpha": frame["alpha"],
            "c": frame["c"],
            "gamma": frame["gamma"],
            "fovea_center": frame["fovea_center"]
            }
        wave = fvht_unpack(pack,True,str_pr + "[channel "+str(i)+"]", d, handle)
        if frame["wave_type"] == 'cdf97':
            ch = lwt.icdf97(wave,False)
        elif frame["wave_type"] == 'cdf53':
            ch = lwt.icdf53(wave,False)
        else:
            ch = dwt.dwt2(wave,frame["wave_type"],frame["decomp_level"],frame["mode"])
        ch = ch - ch.min()
        ch = ch / ch.max() * 255
        ch_i = np.zeros((len(ch),len(ch[0])),np.uint8)
        ch_i[:] = ch
        rgb += [cv.fromarray(ch_i)]
    return tuple(rgb)

def fvht_pack(wave, f_center, Lbpp, lbpp, alpha, c, gamma, delta = 0.01, display_progress=True, str_pr = "", d = {}, handle = True):
    """Compresses a wavelet with SPIHT.

    Runs the original SPIHT algorithm from Dr. Pearlsman paper over the 
    given wavelet.

    Args:
        wavelet: A wavelet to be compressed, must be wavelet2D data type
        bpp: Bits per pixel compression ratio
        delta: Quantization delta if wavelet data is on floating point. 
            This leads to lossy compression always if coefficients are not 
            int
        Lbpp: minimum bit rate and final bit rate compression
        lbpp: maximum bit rate compression at the fovea edges
        alpha: fovea area in percent
        c: constant of the power law function
        gamma: pow of the powerlaw function
        display_progress: True if you want to print on screen algorithm 
                    progress

    Returns:
        A dictionary with a structure with the compressed data and some 
            decoding info. For example:
        {
            "size":0,
            "wave_type":"db7",
            "bpp":0.5,
            "quant_delta":0.001,
            "payload":[]
        }
    """
    #filename = str(len(wave.data)) + "x" + str(len(wave.data[0])) + "_" + str(wave.level) +".dol"
    #if not d:
    #    try:
    #        f = open(filename,"r")
    #        d = pickle.load(f)
    #        f.close()
    #        handle = True
    #    except IOError as e:
    #        d = {}
    #        handle = True
    #update = len(d)
    codec = FVHT()
    codec.wavelet = wave
    codec.bpp = Lbpp
    codec.delta = delta
    codec.check_floating_point()
    codec.str_pr = str_pr
    codec.display_progress = display_progress
    codec.d_memory = d
    codec.Lbpp = Lbpp
    codec.lbpp = lbpp
    codec.alpha = alpha
    codec.c = c
    codec.gamma = gamma
    codec.P = f_center
    codec.compress()
    stream = list(codec.output_stream)[:int(Lbpp*wave.rows*wave.cols)]
    pack = {
            "size":len(stream),
            "rows":len(wave.data),
            "cols":len(wave.data[0]),
            "channels": 1,
            "wave_type":wave.name,
            "quant_delta":delta,
            "decomp_level":wave.level,
            "wise_bit":codec.n,
            "Lbpp": Lbpp,
            "lbpp": lbpp,
            "alpha": alpha,
            "c": c,
            "gamma": gamma,
            "fovea_center" : f_center,
            "bits": codec.amount_of_bits,
            "payload":[stream],
            "test_data":codec.test_data
            }
    #if len(codec.d_memory) > update and handle:
    #    f = open(filename,"w+")
    #    pickle.dump(codec.d_memory,f)
    #    f.close()
    return pack

def fvht_unpack(frame, display_progress=True, str_pr = "",d = {}, handle = True):
    """De-compresses a wavelet with SPIHT.

    Runs the original SPIHT algorithm from Dr. Pearlsman paper over the 
    given wavelet.

    Args:
        frame: A previously compressed frame
        display_progress: True if you want to print on screen algorithm 
                    progress

    Returns:
        An uncompressed wavelet2D instance 
    """
    #filename = str(frame["rows"]) + "x" + str(frame["cols"]) + "_" + str(frame["decomp_level"]) +".dol" 
    #if not d:
    #    try:
    #        f = open(filename,"r")
    #        d = pickle.load(f)
    #        f.close()
    #        handle = True
    #    except IOError as e:
    #        d = {}
    #        handle = True
    #update = len(d)
    codec = FVHT()
    data = np.zeros((frame["rows"],frame["cols"]),np.int32)
    codec.wavelet = wv.wavelet2D(data,frame["decomp_level"],frame["wave_type"])
    codec.bpp = frame["Lbpp"]
    codec.delta = frame["quant_delta"]
    codec.str_pr = str_pr
    codec.n = frame["wise_bit"]
    codec.d_memory = d
    codec.output_stream = deque(frame["payload"])
    codec.Lbpp = frame["Lbpp"]
    codec.lbpp = frame["lbpp"]
    codec.alpha = frame["alpha"]
    codec.c = frame["c"]
    codec.gamma = frame["gamma"]
    codec.P = frame["fovea_center"]
    codec.uncompress()
    #if len(codec.d_memory) > update and handle:
    #    f = open(filename,"w+")
    #    pickle.dump(codec.d_memory,f)
    #    f.close()
    return codec.wavelet

class FVHT(SPIHT):
    #Highest bpp resolution for compression
    Lbpp = 0
    #Lowest bpp resolution for compression
    lbpp = 0 
    #Area of the fovea to be compressed at highest bpp
    fovea_tap = 0
    #Distance of the observer
    dist = 0
    #inner wavelet data type structure
    wavelet = 0
    #Quantization delta for floating point wavelets
    delta = 0.01
    #Extra info for progress display
    str_pr = ""
    display_progress = True
    #DOL sets memory cache
    d_memory = {}
    #Fovea center
    P = (0,0)
    #fovea area
    alpha = 0
    #calculated fovea length
    fovea_length = 0
    
    #power law parameters
    c = 1
    gamma = 1

    #TESTING Variables
    amount_of_bits = 0
    
    def bitplane_encoding(self, power):
        newLSP = deque()
        data = self.wavelet.data
        accum = self.test_data.data
        for p in self.LSP:
            if self.calculate_fovea_w(p) >= self.get_current_bpp():
                #out = (self.wavelet.data[p[0],p[1]] & (2 ** power))
                #self.test_data.data[p[0],p[1]] |= out
                #out = out >> power
                if (accum[p[0],p[1]] + 2 ** power) <=  data[p[0],p[1]]:
                    out = 1
                    accum[p[0],p[1]] += 2 ** power
                    exit = False
                else:
                    out = 0
                self.output_stream.append(out)
                self.amount_of_bits[tuple(p)] += 1
                newLSP.append(p)
        self.LSP = newLSP
        
    def inv_bitplane_encoding(self, power):
        newLSP = deque()
        for p in self.LSP:
            max_limit = self.bit_bucket/float(self.wavelet.cols*self.wavelet.rows)
            if self.calculate_fovea_w(p) >= (max_limit - self.get_current_bpp()):
                bt = self.output_stream.popleft()
                self.wavelet.data[p[0],p[1]] = self.wavelet.data[p[0],p[1]] | (bt << power)
                newLSP.append(p)
        self.LSP = newLSP
    
    def sorting(self, n):
        nextLSP = deque()                                                    
        newLIP = deque()
        data = self.wavelet.data
        zt_matrix = self.zero_tree_data

        #Check for each significant pixel on the LIP
        dif = deque()
        for ij in self.LIP:
            if self.calculate_fovea_w(ij) >= self.get_current_bpp():
                dif.append(ij)
                out = int((abs(data[ij[0],ij[1]]) & 2**n) > 0)
                #out = self.S_n([ij], n)
                self.output_stream.append(out)
                if out == 1:
                    self.test_data.data[ij[0],ij[1]] = (2 ** n)
                    nextLSP.append(ij)
                    self.output_stream.append(self.out_sign(ij))
                    if self.out_sign(ij) == 1:
                        self.test_data.data[ij[0],ij[1]] *= -1
                    self.amount_of_bits[tuple(ij)] = 1
                else:
                    newLIP.append(ij)
            #self.check.append((ij,len(self.output_stream)))
        #Update LIP
        self.LIP = newLIP
        newLIS = deque([],self.wavelet.rows*self.wavelet.cols)
        dif = deque()
        while self.LIS:
            ij = self.LIS.popleft()
        #Check for zerotree roots (2.2.1)
            #D, O, L = self.get_DOL(ij)
            O = self.get_O(ij)
            if ij.entry_type == 'A':
                #out = self.S_n(D,n)
                out = int((zt_matrix[ij[0],ij[1]] >= n))
                self.output_stream.append(out)
                if out == 1:
                    for kl in O:
                        #out = self.S_n([kl],n)
                        out = int((abs(data[kl[0],kl[1]]) & 2**n) > 0)
                        self.output_stream.append(out)
                        if out == 1:
                            if self.calculate_fovea_w(kl) >= self.get_current_bpp():
                                nextLSP.append(kl)
                                dif.append(kl)
                                self.test_data.data[kl[0],kl[1]] = 2 ** n
                                self.output_stream.append(self.out_sign(kl))
                                if self.out_sign(kl) == 1:
                                    self.test_data.data[kl[0],kl[1]] *= -1
                        else:
                            self.LIP.append(kl)
                #Check if ij has grandsons, if not remove from LIS
                    if self.has_grandchilds(ij):
                        ij.entry_type = "B"
                        self.LIS.append(ij)
                else:
                    newLIS.append(ij)
            else: #Entry is type B
                out = self.S_n(O,n)
                self.output_stream += [out]
                if out == 1:
                    for k in O:
                        k.entry_type = "A"
                        self.LIS.append(k)
                        dif.append(k)
                else:
                    newLIS.append(ij) 
        self.LIS = newLIS
        return nextLSP
    
    def compress(self):
        maxs = abs(self.wavelet.data)
        self.n = int(math.log(maxs.max(),2))
        n = self.n
        self.init()
        bit_bucket = self.Lbpp * self.wavelet.rows * self.wavelet.cols
        self.output_stream=deque([],2*bit_bucket)
        self.calculate_fovea_length()
        self.resolve_zero_tree()
        self.amount_of_bits = np.zeros((len(self.wavelet.data), len(self.wavelet.data[0])))
        self.check = deque()
        test = deque()
        #k = self.printFoveaWindow()
        while n >= 0 and len(self.output_stream) < bit_bucket:
            print "Sorting...{0}".format(len(self.output_stream)/float(bit_bucket))
            newLSP = self.sorting(n)
            test.append((len(self.output_stream),len(self.LIP),len(self.LSP)))
            print "Bitplane Encoding...{0}".format(len(self.output_stream)/float(bit_bucket))
            self.bitplane_encoding(n)
            test.append((len(self.output_stream),len(self.LIP),len(self.LSP)))
            self.LSP.extend(newLSP)
            n -= 1
        print "Done"
    
    def inv_sorting(self,n):
        nextLSP = deque()                                                    
        newLIP = deque()
        dif = deque()
        for ij in self.LIP: 
            max_limit = self.bit_bucket/float(self.wavelet.cols*self.wavelet.rows)
            if self.calculate_fovea_w(ij) >= (max_limit - self.get_current_bpp()):
                try:
                    kk = dif.popleft()
                    if not kk == ij:
                        print "error"
                except:
                    pass
                out = self.output_stream.popleft()
                if out == 1:
                    nextLSP.append(ij)
                    self.wavelet.data[ij[0],ij[1]] = (1 << n)
                    sign = self.output_stream.popleft() 
                    if sign == 1:
                        self.wavelet.data[ij[0],ij[1]] *= -1
                else:
                    newLIP.append(ij)
        self.LIP = newLIP
        newLIS = deque([],self.wavelet.rows*self.wavelet.cols)
        dif = deque()
        while self.LIS:
            ij = self.LIS.popleft()
            D, O, L = self.get_DOL(ij)
            if ij.entry_type == "A":
                out = self.output_stream.popleft()
                if out == 1:
                    for kl in O:
                        out = self.output_stream.popleft()
                        if out == 1:
                            max_limit = self.bit_bucket/float(self.wavelet.cols*self.wavelet.rows)
                            if self.calculate_fovea_w(kl) >= (max_limit - self.get_current_bpp()):
                                nextLSP.append(kl)
                                try:
                                    kk = dif.popleft()
                                    if not kk == ij:
                                        print "error"
                                except:
                                    pass
                                sign = self.output_stream.popleft()
                                self.wavelet.data[kl[0],kl[1]] = (1 << n)
                                if sign:
                                    self.wavelet.data[kl[0],kl[1]] *= -1
                        else:
                            self.LIP.append(kl)
                    if L:
                        ij.entry_type = "B"
                        self.LIS.append(ij)
                else:
                    newLIS.append(ij)
            else:
                out = self.output_stream.popleft() 
                if out == 1:
                    for k in O:
                        k.entry_type = "A"
                        self.LIS.append(k)
                        try:
                            kk = dif.popleft()
                            if not kk == k:
                                print "error"
                        except:
                            pass
                else:
                    newLIS.append(ij)
        self.LIS = newLIS
        return nextLSP

    def uncompress(self):
        self.init()
        n = self.n
        self.calculate_fovea_length()
        self.bit_bucket = len(self.output_stream)
        fil = open("check.t","r+")
        test = deque()
        while n >= 0:
            try:
                print "Sorting...{0}".format(1-len(self.output_stream)/float(self.bit_bucket))
                newLSP = self.inv_sorting(n)
                test.append((self.bit_bucket-len(self.output_stream),len(self.LIP),len(self.LSP)))
                print test[-1]
                print "Bitplane Encoding...{0}".format(1-len(self.output_stream)/float(self.bit_bucket))
                self.inv_bitplane_encoding(n)
                test.append((self.bit_bucket-len(self.output_stream),len(self.LIP),len(self.LSP)))
                print test[1]
                self.LSP += newLSP
                n -= 1
            except IndexError:
                break
        print "Done"
    
    def get_current_bpp(self):
        bpp = len(self.output_stream)
        bpp /= float(self.wavelet.rows*self.wavelet.cols)
        return bpp
    
    def calculate_fovea_w(self, ij):
        try:
            P = self.get_center(ij)
        except NameError:
            return self.Lbpp
        H = len(self.wavelet.data)
        W = len(self.wavelet.data[0])
        d = self.norm(P[1] - ij[1],P[0] - ij[0]) * 2**P[2] / self.fovea_length
        if d<self.alpha:
            return self.Lbpp
        elif d>=1:
            return self.lbpp
        else:
            return self.powerlaw(d) * (self.Lbpp - self.lbpp) + self.lbpp
    
    def norm(self,x,y):
        mx = abs(x)
        if mx<abs(y):
            mx = abs(y)
        return mx#math.sqrt(float(x**2 + y ** 2))
    
    def powerlaw(self,n):
        return self.c * (1 - ((n-self.alpha) / (1-self.alpha))) ** self.gamma
        
    def get_center(self,ij):
        if (ij[0] == 0 and ij[1] == 0):
            raise NameError("ij on HH")
        else:
            if ij[0] == 0:
                aprx_level_r = self.wavelet.level + 1
            else:
                aprx_level_r = math.ceil(math.log(self.wavelet.rows/float(ij[0]),2))
                if aprx_level_r > self.wavelet.level:
                    aprx_level_r = self.wavelet.level + 1
            if ij[1] == 0:
                aprx_level_c = self.wavelet.level + 1
            else:
                aprx_level_c = math.ceil(math.log(self.wavelet.rows/float(ij[1]),2))
                if aprx_level_c > self.wavelet.level:
                    aprx_level_c = self.wavelet.level + 1
            if (aprx_level_r > self.wavelet.level) and (aprx_level_c > self.wavelet.level):
#                raise NameError("ij on HH")
                y = float(self.P[0]) / 2**(aprx_level_r-1)
                x = float(self.P[1]) / 2**(aprx_level_r-1)
                return (y,x,aprx_level_r-1) 
            if aprx_level_r <= aprx_level_c:
                aprx_level = aprx_level_r
            else:
                aprx_level = aprx_level_c
            y = float(self.P[0]) / 2**aprx_level
            x = float(self.P[1]) / 2**aprx_level 
            if aprx_level_r == aprx_level:
                y += float(self.wavelet.rows) / 2**aprx_level
            if aprx_level_c == aprx_level:
                x += float(self.wavelet.cols) / 2**aprx_level
            return (y, x, aprx_level)
    
    def calculate_fovea_length(self):
        H = len(self.wavelet.data)
        W = len(self.wavelet.data[0])
        k = np.zeros(4)
        k[0] = self.norm(self.P[0],H-self.P[1])
        k[1] = self.norm(W-self.P[0],self.P[1])
        k[2] = self.norm(W-self.P[0],H-self.P[1])
        k[3] = self.norm(self.P[0],H-self.P[1])
        self.fovea_length = k.max()
        
    def printFoveaWindow(self):
        window = np.zeros((self.wavelet.rows, self.wavelet.cols))
        points = wv.get_z_order(self.wavelet.rows * self.wavelet.cols)
        for i in points:
            window[tuple(i)] = self.calculate_fovea_w(i)
        return window
            
def quant(wavelet,delta):
    iw = np.zeros((len(wavelet.data),len(wavelet.data[0])),np.int64)
    iw= np.array(np.trunc(wavelet.data / delta), np.int64)
    wavelet.data = iw
    return wavelet
