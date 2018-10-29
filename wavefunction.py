# -*- coding:utf-8 -*-
import numpy as np
import sys
sys.path.append(r"/Users/koda/PyModules/SimpleQmap")
sys.path.append(r"/Users/koda/PyModules")
sys.path.append(r"C:\PyModules\SimpleQmap")
sys.path.append(r"C:\PyModules")
import matplotlib.pyplot as plt
import SimpleQmap as sq
from maps import StandardMap as st
import numpy as np
import pylab
import matplotlib.pyplot as plt
import scipy as sp
import pandas as pd


Seeds = [[],[],[]]
psi = [[],[],[]]#[orbitnumber,value,p]
density = [[],[]]#[value,p]
twopi = 2.0*np.pi
grid=500
regrid = 500
omega = 1
initp = 0
imsize = 3
resize = 4
realscale = 1
imscale = 1 #(realscale,imscale)=(1,1) is used as usual.
step = 6#これ任意のステップでなんとかなるように.実際に出力したマップとの整合性に気をつけて.
radius = 0
diskgrid = 200
orbit = [[],[],[]]
orbitnumber = 0
originseed = []
originseed2 =[]
origin = np.array([])
arg2 = 0
k,d = 1.2,5
radius = 0.01
class ModifiedKMap(st):
    def __init__(self, k,d):
        st.__init__(self,k)
        self.d = d
    def func0(self,x,k):
        return k*np.sin(x)
    def func1(self, x):
        return (pow(x,7)/(pow(x,6)+pow(d,6)))*(1+(3*pow(d,6)/(pow(x,6)+pow(d,6))))
    def func2(self,x):
        return (x**2/2)*( x**6/(x**6 + d**6)  )
    def func3(self,x,k):
        return k * np.cos(x)

def probdensity(S,p,cmap):# S,p are got from "Action".
     #orbit has information of inittheta
     #cmap here means tMap
     densit = psi[1]*np.conj(psi[1])
     #density = [density, p]
     return densit,p

#点の情報を代入して作用を求める.
def Action(omega,theta,p,cmap,step,orbitnumber):#theta,p is initial condition.
    S = 0
    ppt = initp
    for i in range(step):
         theta, p = Map(theta,ppt,cmap)
         S = S + cmap.func2(p) + omega * p + cmap.func3(p,k) - theta*(ppt-p)
         ppt = p
    #print(S)
    return  S,p#is that all we need evaluating that

def Map(theta,p,cmap):
    thetaa = theta+(omega+cmap.func1(p))
    pp = p-cmap.func0(thetaa,cmap.k)
    return [thetaa,pp]

def tMap(theta,p,cmap,step):#まとめてマップ
    for i in range(step):
        theta,p  = Map(theta,p,cmap)
    return [theta,p]


def WaveFunction(S,p):
    psi = [[],[],[]]
    psi[0].extend([orbitnumber])
    psi[1].extend( [np.exp(1j * S)] )
    #print(psi[1])
    psi[2].extend([p])
    return psi
def appendorbit(orbit):
    global key,Seeds
    print("Which orbit do you append?")
    key = int(input())#番号を入力
    orbit = np.loadtxt("orbit{0}.txt".format(key))
    for i in range(3):
        Seeds[i].extend( [orbit[i][k] for k in range(len(orbit[0])) if orbit[0][k] == key])#orbitnumber(number) is that
    return Seeds

cmap = ModifiedKMap(k,d)
#初期条件の塊を持ってくる.
Seeds = appendorbit(orbit)
#波動関数の描画
fig =  plt.figure(figsize = (6,6))
ax = fig.add_subplot(111)
for i in range(len(Seeds[0])):
    theta = Seeds[1][i] + 1j * Seeds[2][i]
    S,p = Action(omega,theta,initp,cmap,step,orbitnumber)
    #theta,p= tMap(theta,initp,cmap,step)  #for check
    psi = WaveFunction(S,p)
    density[0].extend([probdensity(S,p,cmap)[0]])
    density[1].extend([(probdensity(S,p,cmap)[1]).real/twopi])
print(density[0])
plt.plot(density[1],density[0],',k')
plt.xlim(-2,2)
plt.yscale("log")
plt.show()
