import copy
from PIL import Image
import numpy as np
from numpy.core.multiarray import ndarray
from scipy.fftpack import dctn
from scipy.fftpack import idctn
import math
import matplotlib.pyplot as plt


class Dct:


# bs_n is blocksize n.
#imagepath is the input of n.
    def __init__(self,imagepath,bs_n):
        self.bs_n = bs_n
        image_int8 = Image.open(imagepath)
        # image is the input of n.
        self.image = np.asarray(image_int8, dtype="int32")
        self.Row = len(self.image)
        self.Col = len(self.image[0])
        self.procImage = np.asarray(self.image, dtype="float")

    def save(self, path):
        im = Image.fromarray(np.uint8(self.procImage))
        im.save(path)

    def CompressionRate(self):
        n = self.bs_n
        return 8*n*n/(3*int((n*n-1)/10)+2*int((n*n-1)/10)+4)

    def Snr(self):
        self.noise = abs(self.signal-self.image)
        pnoise = np.mean(self.noise**2)
        psignal = np.mean(self.image**2)
        return 10*math.log10(psignal/pnoise)



    def blockImage(self):
        DC_term=[]
        n = self.bs_n
        for row in range(0, self.Row, n):
            for col in range(0, self.Col, n):
                tempblock = self.procImage[row: row+n, col: col+n]
                dctblock = dctn(tempblock, norm='ortho')
                DC_term.append(dctblock[0][0])
                self.getACterm(dctblock)
                self.procImage[row: row + n, col: col + n] = dctblock
        self.QdeQ_array(DC_term, 16)
        for row in range(0, self.Row, n):
            for col in range(0, self.Col, n):
                dcColLen = int(self.Col/n)
                tempblock = self.procImage[row: row + n, col: col + n]
                index=int(row / n * dcColLen + col / n)
                tempblock[0][0] = DC_term[index]
                idctblock = idctn(tempblock, norm='ortho')
                self.procImage[row: row + n, col: col + n] = idctblock
        self.signal =  np.uint8(self.procImage)
        pass



    def getACterm(self,mat):
        n = self.bs_n
        termNum = int((n * n - 1) / 10)
        firstTerm = []
        secondTerm = []
        zzIndex = 0
        for i in range(1, 2 * n - 1):
            if i % 2 == 1:
                # down left
                x = 0 if i < n else i - n + 1
                y = i if i < n else n - 1
                while x < n and y >= 0:
                    if zzIndex < termNum:
                        firstTerm.append(mat[x][y])
                    elif termNum <= zzIndex < 2 * termNum:
                        secondTerm.append(mat[x][y])
                    else:
                        mat[x][y] = 0.0
                    x += 1
                    y -= 1
                    zzIndex += 1
            else:
                x = i if i < n else n - 1
                y = 0 if i < n else i - n + 1
                while x >= 0 and y < n:
                    if zzIndex < termNum:
                        firstTerm.append(mat[x][y])
                    elif termNum <= zzIndex < 2 * termNum:
                        secondTerm.append(mat[x][y])
                    else:
                        mat[x][y] = 0.0
                    x -= 1
                    y += 1
                    zzIndex += 1
        if firstTerm:
            self.QdeQ_array(firstTerm, 8)
        if secondTerm:
            self.QdeQ_array(secondTerm, 4)
        zzIndex = 0
        fi = 0
        si = 0
        for i in range(1, 2 * n - 1):
            if zzIndex >= 2 * termNum:
                break
            if i % 2 == 1:
                # down left
                x = 0 if i < n else i - n + 1
                y = i if i < n else n - 1
                while x < n and y >= 0:
                    if zzIndex < termNum:
                        mat[x][y] = firstTerm[fi]
                        fi += 1
                    elif termNum <= zzIndex < 2 * termNum:
                        mat[x][y] = secondTerm[si]
                        si += 1
                    else:
                        break
                    x += 1
                    y -= 1
                    zzIndex += 1
            else:
                x = i if i < n else n - 1
                y = 0 if i < n else i - n + 1
                while x >= 0 and y < n:
                    if zzIndex < termNum:
                        mat[x][y] = firstTerm[fi]
                        fi += 1
                    elif termNum <= zzIndex < 2 * termNum:
                        mat[x][y] = secondTerm[si]
                        si += 1
                    else:
                        break
                    x -= 1
                    y += 1
                    zzIndex += 1

    def QdeQ_array(self, array, level):
        H = np.max(array) + 0.000001
        L = np.min(array)
        for i in range(len(array)):
            array[i] = self.QdeQ(H, L, array[i], level)

    def QdeQ(self,H,L,value,level):
        gap = (H-L)/level
        return (int((value - L)/gap)+0.5)*gap+L





#quantized then dequantized





if __name__ == '__main__':
    path = '/Users/shibowen/Documents/data compression/newpicture.png'
    pathlena = '/Users/shibowen/Documents/data compression/lena.png'
    n =[2,4,8,16,32,64]
    #bird
    snrplot=[]
    compplot = []
    for element in n:
        path2 = f'/Users/shibowen/Documents/data compression/newimage{element}.png'
        imagea = Dct(path,element)
        imagea.blockImage()
        imagea.save(path2)
        snr = imagea.Snr()
        comp = imagea.CompressionRate()
        snrplot.append(snr)
        compplot.append(comp)
    print('sky and bird result')
    print(n)
    print(snrplot)
    print(compplot)
    plt.plot(n,snrplot)
    plt.title("the sky and bird snr to n")
    plt.show()

    #lena
    snrplotlena=[]
    compplotlena = []
    for element in n:
        path2 = f'/Users/shibowen/Documents/data compression/lena{element}.png'
        imagea = Dct(pathlena,element)
        imagea.blockImage()
        imagea.save(path2)
        snr = imagea.Snr()
        comp = imagea.CompressionRate()
        snrplotlena.append(snr)
        compplotlena.append(comp)
    print('lena result')
    print(n)
    print(snrplotlena)
    print(compplotlena)
    plt.plot(n, snrplotlena)
    plt.title("the lena snr to n")
    plt.show()




