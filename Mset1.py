# cording : utf-8
import numpy as np
import sys
import os.path
sys.path.append(r"/Users/koda/PyModules/SimpleQmap")
sys.path.append(r"/Users/koda/PyModules")
sys.path.append(r"C:\PyModules\SimpleQmap")
sys.path.append(r"C:\PyModules")
import matplotlib.pyplot as plt
import SimpleQmap as sq
from maps import StandardMap as st
from mpl_toolkits.mplot3d import Axes3D
import sympy as sym
from sympy import Array
from numba import jit,autojit
import math
#The initial letter 'p' in this code means pre- .
twopi = 2.0*np.pi
grid=100
regrid = 100
omega = 1
initp = 0
imsize = 3
resize = 5
reorigin = 0.0
imorigin = 0.0 #Origin on the complex surface is (reorigin,imorigin)
realscale =1
imscale = 1 #(realscale,imscale)=(1,1) is used as usual.
step = 6
radius = 0
diskgrid = 200
orbit = [[],[],[]]
diskcount = True
originseed = []
originseed2 =[]
origin = np.array([])
arg2 = 0
switch = False
switch3 = False
trycount = 0
radius = 0.01
cradius = 0.01
orbitnumber = 17
#print("Number?")
#orbitnumber = int(input())
#while os.path.exists('orbit{0}.txt'.format(orbitnumber)) == True:
#    print("Reserved.Number again please?")
#    orbitnumber = int(input())


#################
#もし点がうまく取れなければ、
#半径(radius)の値を小さくしてみてください。
#スピードは少々遅くなってしまいますが、小さいほど正確に取れます。
#また、自然分岐に関してはうまくやってください.
#今のところステップ数に応じて変更しなければなりませんが、
#ヘッダーでステップ数を読み込ませる機能を後ですぐに実装します。
#################
class ModifiedKMap(st):
    def __init__(self, k,d):
        st.__init__(self,k)
        self.d = d
    def func0(self,x,k):
        return k*np.sin(x)
    def func1(self, x):
        return (pow(x,7)/(pow(x,6)+pow(d,6)))*(1+(3*pow(d,6)/(pow(x,6)+pow(d,6))))


def Map(theta,p,cmap):
    thetaa = theta+(omega+cmap.func1(p))
    pp = p-cmap.func0(thetaa,cmap.k)
    return [thetaa,pp]
######################################
def draw(p,cmap):#等高線
    seed = getseed(p,cmap,step,0)
    fig = plt.figure(figsize=(6,6))
    ax = fig.add_subplot(1,1,1)
    ax.plot(seed[0],seed[1],',k')
    np.savetxt('Mset_step1.txt',seed)
    #plot xi and eta
    plt.xlim(reorigin,reorigin+(1/realscale))
    plt.ylim(imorigin,imorigin+(1/imscale))
    plt.show()
    return


def tMap(theta,p,step,cmap):
    for i in range(step):
        theta,p = Map(theta,p,cmap)
        theta = theta.real-np.floor(theta.real/twopi)*twopi+1j*theta.imag
    return theta,p

def getseed(p,cmap,step,x):
    a = [[],[]]
    for i in range(imsize * grid):
        if i%10 == 0:
            print(i)
        for k in range(resize*regrid): #k is real axis.
            ppt = p
            p = initp
            theta = reorigin+k*(1/(realscale*regrid)) + 1j*((1/(imscale*grid))*i+imorigin)
            theta = theta.real-np.floor(theta.real/twopi)*twopi+1j*theta.imag
            inittheta = theta
            theta, p = tMap(theta,p,step,cmap)
            if np.sign((ppt.imag-x)*(p.imag-x)) < 0 :
                a[0].extend([inittheta.real-1/(2*realscale*regrid)])
                a[1].extend([inittheta.imag])
        for k in range(resize*regrid): #k is real axis.
            ppt = p
            p = initp
            theta = -reorigin-k*(1/(realscale*regrid)) - 1j*((1/(imscale*grid))*i-imorigin)
            theta = theta.real-np.floor(theta.real/twopi)*twopi+1j*theta.imag
            inittheta = theta
            theta, p = tMap(theta,p,step,cmap)
            if np.sign((ppt.imag-x)*(p.imag-x)) < 0 :
                a[0].extend([inittheta.real-1/(2*realscale*regrid)])
                a[1].extend([inittheta.imag])
    return a

def v(self,p,cmap):#等高線
    b = [[],[]]
    for x in range(5):
       constant = x * 0.4
       b = getseed(p,cmap,step,x)
       seed[0].extend(b[0])
       seed[1].extend(b[1])
    np.savetxt('Mset_step1.txt',seed)#出力するファイルの名前
    fig = plt.figure(figsize=(6,6))
    ax = fig.add_subplot(1,1,1)
    ax.plot(seed[0],seed[1],',k')
    #plot xi and eta
    plt.xlim(reorigin,reorigin+(1/realscale))
    plt.ylim(imorigin,imorigin+(1/imscale))
    plt.show()
##############################################################
def Press(event):
    global x1,y1,DragFlag,realscale,regrid

    if(event.xdata is None) or (event.ydata is None):
        return
    cx = event.xdata
    cy = event.ydata

    x1 = cx-1/(2*realscale*regrid)
    y1 = cy
    if event.button == 2:
        DragFlag = True
    else:
        return

def Drag(event):
    global x1,y1,x2,y2,DragFlag

    if DragFlag == False:
        return

    if (event.xdata is None) or (event.ydata is None):
        return

    cx = event.xdata
    cy = event.ydata

    x2 = cx-1/(2*realscale*regrid)
    y2 = cy

    ix1,ix2 = sorted([x1,x2])
    iy1,iy2 = sorted([y1,y2])
    DrawRect(x1,x2,y1,y2)
    plt.draw()

def Release(event):#Get x1,y1(initialoriginx,initialoriginy)
    global x1,y1,x2,y2
    global seed
    global orbit
    global DragFlag
    global trycount,diskcount
    orbitnumber = 0
    switch  = False
    switch2 = False
    switch3 = False
    radius  = cradius 
    if event.button == 2:
        DragFlag = False
    #print(x1,y1)
    #plt.plot(x1,y1,'.')
    #pointexistjudge(seed)
    originx = x1
    originy = y1
    alpha =0
    beta = 0
    preal = 0
    pradius = 0
    originseed,arg2 = initialdisksize(alpha,beta,originx,originy,radius,diskgrid,cmap,switch)
    pradius = cradius#仮の値
    radius = originseed[0]
    orbit = getpoint(alpha,beta,originx,originy,0,0,radius,pradius,cmap,switch)#originseed = [radius,arg]
    diskcount = not(diskcount)
    radius = originseed[0]
    arg = originseed[1]
    #print("originx = %f,originy = %f " %(originx,originy))
    poriginx,poriginy  = originx,originy
    originx,originy = bisectionOndisk(originx,originy,radius,arg,alpha,beta,cmap)#bisectionOndisk(originx,originy,radius,arg,alpha,beta,cmap)    return
    #print("originx = %f,originy = %f " %(originx,originy))
    origin = initialdisksize(alpha,beta,originx,originy,radius,diskgrid,cmap,switch)[0]#initialdisksize(alpha,beta,originx,originy,radius,diskgrid,cmap)
    #print(origin)
    pradius = radius
    radius = origin[0]
    #print("radius = %f,arg = %f"%(radius,arg))
    alpha,beta = disksections(poriginx,poriginy,originx,originy,radius,pradius)
    origin = []
    origin = disksize(alpha,beta,originx,originy,radius,diskgrid,cmap,switch)
    arg = origin[1]
    orbit = getpoint(alpha,beta,originx,originy,poriginx,poriginy,radius,pradius,cmap,switch)
    print(arg)
    poriginx = originx
    poriginy = originy
    originx,originy = bisectionOndisk(originx,originy,radius,arg,alpha,beta,cmap)#二個目の原点を確定させる.
    plt.plot(originx,originy,'.')
    inittheta = originx + 1j * originy
    preal = tMap(inittheta,initp,step,cmap)[1].real
    print('Re(p) = %f' %preal)
    while abs(tMap(inittheta,initp,step,cmap)[1].real) < 5 :
        #plt.plot(originx,originy,'.')
        #前の円の半径の情報
        trycount += 1
        origin =[]
        #decidenext disksize
        origin = initialdisksize(0,0,originx,originy,radius,diskgrid,cmap,switch)[0] #次の円の大きさを決める（交点が求められる.）
        #print(origin)
        #print('origin = {0}'.format(origin))
        pradius = radius
        radius  = origin[0]
        #print("radius = %f,arg = %f"%(radius,arg))
        alpha,beta = disksections(poriginx,poriginy,originx,originy,radius,pradius)
        orbit = getpoint(alpha,beta,originx,originy,poriginx,poriginy,radius,pradius,cmap,switch)
        poriginx = originx
        poriginy = originy
        origin = []
        origin = disksize(alpha,beta,originx,originy,radius,diskgrid,cmap,switch)
        arg = origin[1]
        originx,originy = bisectionOndisk(originx,originy,radius,arg,alpha,beta,cmap)#二個目の原点を確定させる.
        #print('radius = %f ,arg = %f'%(radius,arg))
        #a,b is provisional
        #originx,originy.
        #次の円との交点.
        #前の円の情報を
        preal = tMap(inittheta,initp,step,cmap)[1].real
        print('Re(p) = %f' %preal)
        inittheta = originx + 1j * originy
            #if trycount >  10:
            #trycount = 0
            #break
    #disksection evaluate alpha(minarg),beta(maxarg)
    trycount = 0
    print('逆の方向')
    alpha = 0
    beta = 0
    radius = cradius 
    pradius = 0
    originx = x1
    originy = y1
    switch = True
    origin =  initialdisksize(0,0,originx,originy,radius,diskgrid,cmap,switch)[0] #次の円の大きさを決める（交点が求められる.）
    radius = origin[0]
    arg = origin[1]
    print("radius = %f,arg = %f"%(radius,arg))
    poriginx = originx
    poriginy  = originy
    originx,originy = bisectionOndisk(originx,originy,radius,arg,alpha,beta,cmap)
    print("originx =%f,originy = %f" %(originx,originy))
    origin = []
    origin =  initialdisksize(0,0,originx,originy,radius,diskgrid,cmap,switch)[0]#2つ目
    pradius = radius
    radius = origin[0]
    origin = []
    alpha,beta = disksections(poriginx,poriginy,originx,originy,radius,pradius)
    origin = disksize(alpha,beta,originx,originy,radius,diskgrid,cmap,switch)
    arg = origin[1]
    orbit = getpoint(alpha,beta,originx,originy,poriginx,poriginy,radius,pradius,cmap,switch)#originseed = [radius,arg]
    poriginx = originx
    poriginy= originy
    originx,originy = bisectionOndisk(originx,originy,radius,arg,alpha,beta,cmap)
    inittheta = originx + 1j * originy
    plt.plot(originx,originy,'.')
    origin = []
    #decidenext disksize
    while abs(tMap(inittheta,initp,step,cmap)[1].real) < 5 :
            trycount += 1
            origin =  initialdisksize(0,0,originx,originy,radius,diskgrid,cmap,switch)[0]#2つ目
            pradius = radius
            radius = origin[0]
            alpha,beta = disksections(poriginx,poriginy,originx,originy,radius,pradius)#disksize(alpha,beta,originx,originy,poriginx,poriginy,radius,pradius,diskgrid,cmap):            origin = []
            origin  = disksize(alpha,beta,originx,originy,radius,diskgrid,cmap,switch)
            arg = origin[1]
            orbit = getpoint(alpha,beta,originx,originy,poriginx,poriginy,radius,pradius,cmap,switch)
            poriginx,poriginy = originx,originy
            originx,originy = bisectionOndisk(originx,originy,radius,arg,alpha,beta,cmap)
            #print("poriginx = %f,poriginy = %f" %(poriginx,poriginy))
            inittheta = originx + 1j * originy
            plt.plot(originx,originy,'.')
            preal = tMap(inittheta,initp,step,cmap)[1].real
            print('Re(p) = %f' %preal)
            #print(arg,arg2)
            origin = []
                #   if trycount > 1:
                #trycount = 0
                #break
        #disksection evaluate alpha(minarg),beta(maxarg)
        #alpha,beta = disksections(poriginx,poriginy,originx,originy,radius,pradius) #disksections(poriginx,poriginy,originx,originy,radius,pradius):
    print(orbit)
    #np.savetxt('orbit.txt',orbit)
    np.savetxt("orbit{0}.txt".format(orbitnumber),orbit)
        #with open('orbit.txt',mode = 'wb') as f:
        #f.write('\n'.join(orbit[][i] for i in len(orbit[0])))
    orbitnumber += 1
     #########################################################
#########################################################
#def pointexistjudge(seed):#中にちゃんと点があるのかどうか.判定する場所
#    global x1,y1
#    check = True
#    for i in range(len(seed[0])):

#        if  (seed[0][i]-x1)**2 + (seed[1][i] - y1)**2 < radius**2:
#            check = False
#            break
#
#    if check == True:
#        print("please click the place where orbit is there")
#        return 1
#
#        return
########################################################
def initialdisksize(alpha,beta,originx,originy,radius,diskgrid,cmap,switch): #this method is used to one of the directions of the two.
        global orbit
        initp = 0
        count = 0
        i = 0
        disk = [[],[]]
        originseed =np.array([])
        originseed2 = np.array([])
        originseedx = 0
        originseedy = 0
        arg2 = 0
        startarg = beta - twopi
        endarg = alpha
        startarg,endarg =sorted([startarg,endarg])
        ppointx = originx + radius * np.exp(1j*startarg).real
        ppointy = originy + radius * np.exp(1j*startarg).imag
        while i < diskgrid :
            i += 1
            arg = (endarg-startarg) * i/diskgrid + startarg
            #arg = arg - np.floor(arg/twopi)*twopi
            pointx = originx + radius * np.exp(1j*arg).real
            pointy =  originy + radius * np.exp(1j*arg).imag # point on inittheta plane
            disk[0].extend([pointx])
            disk[1].extend([pointy])
            inittheta = pointx + 1j * pointy
            pinittheta = ppointx + 1j * ppointy
            #print('Im(p)=%f,Im(pp)=%f'%(tMap(inittheta,initp,step,cmap)[1].imag,tMap(pinittheta,initp,step,cmap)[1].imag))
            if  tMap(inittheta,initp,step,cmap)[1].imag * tMap(pinittheta,initp,step,cmap)[1].imag < 0 : #p is in tMap[1]
                if count == 0:
                    if switch == False:
                        originseed = np.array([radius,arg])
                        originseedx =  originx + radius * np.exp(1j * arg).real
                        originseedy = originy + radius * np.exp(1j * arg).imag
                    else:
                        arg2 = arg
                    count += 1
                elif count == 1:
                    if switch == True:
                        originseed = np.array([radius,arg])
                        originseedx =  originx + radius * np.exp(1j * arg).real
                        originseedy =  originy + radius * np.exp(1j * arg).imag
                    else:
                        arg2 = arg
                    count += 1
                else: count += 1
            ppointx = pointx
            ppointy = pointy
            if count >= 3  and originy > 0.1 ** 3 :
                   radius = radius/2
                   i = 0
                   count = 0
                   disk = [[],[]]
           # elif count >= 3 and originy < 0.1 ** 3 :
           #         print('自然分岐')
           #         radius  = radius  * 2
           #         i = 0
           #         count = 0
           #         disk = [[],[]]
            if arg >= endarg:
                if count ==1 or count == 0:
                    print('aaaa')
                    np.savetxt('orbit{0}.txt'.format(orbitnumber),orbit)
                    return
                    i = 0
                break
        plt.plot(disk[0],disk[1],',k')
        return originseed,arg2
        #originseed = (radius,arg), 2 be that of another side point

def disksize(alpha,beta,originx,originy,radius,diskgrid,cmap,switch): #this method is used to one of the directions of the two.
        initp = 0
        count = 0
        i = 0
        disk = [[],[]]
        originseed =np.array([])
        originseed2 = np.array([])
        originseedx = 0
        originseedy = 0
        originseed2x = 0
        originseed2y = 0
        startarg = alpha
        endarg = beta
        ppointx = originx + radius * np.exp(1j*startarg).real
        ppointy = originy + radius * np.exp(1j*startarg).imag
        #print("originx = %f,originy = %f" %(originx,originy))
        #print("startarg(disksize) = %f, endarg(disksize) = %f" %(startarg,endarg))#
        while i < diskgrid :
            i +=1
            arg = (endarg-startarg)*i/diskgrid+startarg
            pointx = originx + radius * np.exp(1j*arg).real
            pointy =  originy + radius * np.exp(1j*arg).imag # point on inittheta plane
            disk[0].extend([pointx])
            disk[1].extend([pointy])
            inittheta = pointx + 1j * pointy
            pinittheta = ppointx + 1j * ppointy
            #print('Im(p)=%f,Im(pp)=%f'%(tMap(inittheta,initp,step,cmap)[1].imag,tMap(pinittheta,initp,step,cmap)[1].imag))
            if  tMap(inittheta,initp,step,cmap)[1].imag * tMap(pinittheta,initp,step,cmap)[1].imag < 0 : #p is in tMap[1]
                #print(originx,originy)
                if count == 0  :
                    #plt.plot(pointx,pointy,'.')
                    originseed = np.array([radius,arg])
                    #print('originseed = %f' %originseed[1])
                    originseedx =  originx + radius * np.exp(1j * arg).real
                    originseedy = originy + radius * np.exp(1j * arg).imag
                    count += 1
                elif count == 1 :
                     count += 1
            ppointx = pointx
            ppointy = pointy
            if (count >= 2 and alpha != 0 and beta != twopi ) or ( alpha == 0 and beta == twopi and count >= 3) and originy > 0.1 ** 3:
                radius =radius/2
                count = 0
                i = 0
                disk = [[],[]]
                print(radius)
           # elif (count >= 2 and alpha != 0 and beta != twopi ) or ( alpha == 0 and beta == twopi and count >= 3) and originy > 0.1 ** 3 :
           #      print('自然分岐')
           #      radius  = radius * 2
           #      count = 0
            if arg >= endarg:
                if count == 0 :
                 break
        plt.plot(disk[0],disk[1],',k')
        return originseed
        #originseed = (radius,arg), 2 be that of another side point


def disksections(poriginx,poriginy,originx,originy,radius,pradius):
    global switch
    a = 0
    b = 0
    arg  = 0
    alpha = 0
    beta = 0
    xa = 0
    ya = 0
    i=0
    initx = originx + radius * math.cos(arg)
    inity = originy + radius * math.sin(arg)
    itheta = initx + 1j * inity
    a = math.atan2(poriginy-originy,poriginx-originx)
    #print("originx = %f ,originy = %f" %(originx,originy))
    #print("poriginx = %f ,poriginy = %f" %(poriginx,poriginy))
    print("pradius = %f , radius = %f " %(pradius ,radius ))
    d = math.sqrt((originx - poriginx) ** 2 + (originy - poriginy) ** 2 )
    cos = (d**2 + radius **2 - pradius **2)/(2*d*radius)
    if (abs(cos)>1):
        alpha = 0
        beta = twopi
        return alpha,beta
    print((d**2 + radius **2 - pradius **2)/(2*d*radius))
    arg = math.acos((d**2 + radius **2 - pradius **2)/(2*d*radius))
    print("arg = %f"%arg)
    alpha =  a + arg
    beta  = a - arg
    alpha = alpha - np.floor(alpha/twopi)*twopi
    beta = beta - np.floor(beta/twopi)*twopi
    alpha,beta = sorted([alpha,beta])
    print('alpha = %f,beta = %f' %(alpha,beta))
    xa = originx + radius  * np.cos(alpha)
    ya = originy + radius * np.sin(alpha)
    if abs(alpha - beta) < np.pi :
        beta = beta - twopi
        tmp = alpha
        alpha = beta
        beta = tmp
    plt.plot(xa,ya,'.')
    return alpha,beta
##################################################bisections
def bisection(x,y,x3,y3,cmap):
    global originx,originy
    epsilon = 0.001
    i = 0
    if x > x3:
        tmpx = x
        tmpy = y
        x = x3
        x3 = tmpx
        y =y3
        y3 = tmpy
    #print("x= %f, y= %f"%(x,y))
    dinittheta = x + 1j *y
    uinittheta = x3 + 1j * y3
    while i < 100:
        midx = (x+x3)/2
        midy = (y+y3)/2
        inittheta = x + 1j * y
        midinittheta = midx + 1j * midy
        #print(x,x3)
        if tMap(dinittheta,initp,step,cmap)[1].imag * tMap(uinittheta,initp,step,cmap)[1].imag > 0:#tMap(theta,p,step,cmap)
            return [np.nan,np.nan]
        if tMap(inittheta,initp,step,cmap)[1].imag * tMap(midinittheta,initp,step,cmap)[1].imag < 0 :
            x3 = midx
            y3 = midy
        else:
             x = midx
             y = midy
        if math.sqrt((x-x3)**2+(y-y3)**2) < 0.5 ** 20:
            #print('break')
            break
            i += 1
    return [midx,midy]

def bisectionOndisk(originx,originy,radius,arg,alpha,beta,cmap):
    #print(arg)
    epsilon = twopi/diskgrid #epsilon is  angular variable
    if alpha > 0 and  arg  < 0 :
        beta = beta - twopi
        tmp = alpha
        alpha = beta
        beta = tmp
    startarg = min(alpha,beta)
    endarg = max(alpha,beta)
    #print(beta,alpha)
    #print("startarg(bisectionOndisk) = %f, endarg(bisectionOndisk) = %f" %(startarg,endarg))
    upperarg = (arg  + epsilon)- np.floor((arg  + epsilon)/ twopi)*twopi
    downerarg =(arg  - epsilon) - np.floor((arg  - epsilon) / twopi)*twopi
    #plt.plot(originx + radius * math.cos(upperarg),originy + radius * math.sin(upperarg),'.')
    i = 0
    x = 0
    y = 0
    while i < 100:
       # print("startarg(bisectionOndisk) = %f, endarg(bisectionOndisk) = %f" %(startarg,endarg))
        arg = (upperarg + downerarg)/2
        pointx =  originx +radius * np.exp(1j * arg).real
        pointy =  originy + radius * np.exp(1j * arg).imag #upointx,dpointx..etc means the poit on the disk near Im0p point
        upointx = originx + radius * np.exp(1j * upperarg).real
        dpointx = originx +radius *np.exp(1j * downerarg).real
        upointy = originy + radius *np.exp(1j * upperarg).imag
        dpointy = originy + radius * np.exp(1j * (downerarg)).imag
        plt.plot(pointx, pointy,'.')
        dinittheta = dpointx + 1j * dpointy
        uinittheta = upointx + 1j * upointy
        midinittheta = pointx + 1j * pointy
        if tMap(dinittheta,initp,step,cmap)[1].imag * tMap(uinittheta,initp,step,cmap)[1].imag >= 0 :
                upperarg += epsilon
                downerarg -= epsilon
                i = 0
            #print('a')
        elif tMap(midinittheta,initp,step,cmap)[1].imag * tMap(uinittheta,initp,step,cmap)[1].imag < 0:#tMap(theta,p,step,cmap):
            downerarg = arg
            #print('b')
        else:
            upperarg = arg
            #print('c')
        if abs(upperarg-downerarg) < 0.5**10 :
            #print('break')
            x = originx + radius * np.cos(arg)
            y = originy + radius * np.sin(arg)
            return originx + radius * np.cos(arg) , originy + radius * np.sin(arg)
        #print(arg)
        i += 1
        plt.plot(x,y,'.')
    return x,y

###################################################
def imp0seed(x,y):#Evaluate  Im(p)=0 accurately by using bisection
    epsilon = 2.5
    xmin =x-epsilon
    xmax =x+epsilon
    initreal,initimag = biseciton(xmin,y,xmax,ymax,cmap)
    return initreal,initimag
@jit
def getpoint(alpha,beta,originx,originy,poriginx,poriginy,radius,pradius,cmap,switch):#少なくともこれは最適化しないとまずい.
    global orbit,orbitnumber,diskcount,x1,y1
    a = 0
    b = 0
    count = 0
    seeds = int(radius * 5000)
    #print(switch)
    atan = math.atan2(originy-poriginy,originx-poriginx)
    print("atan = %f " %atan)
    #if switch == True :
    #    atan = -atan
    startarg = atan + np.pi/2
    endarg = startarg + np.pi
    print("alpha = %f,beta = %f" %(alpha,beta))
    print("startarg = %f,endarg = %f"%(startarg , endarg))
    if poriginx == 0 and poriginy == 0  :
        startarg = 0
        seeds = int(radius * 7000)
        endarg = twopi
        alpha = -0.01
        beta = twopi+0.01
        #print("alpha = %f,beta = %f" %(alpha,beta))
    elif endarg > beta  :
        endarg = endarg -twopi
        tmp = startarg
        startarg = endarg
        endarg = tmp
    if alpha > 0 and startarg < 0:
        #print("cc")
        endarg = endarg -twopi
        tmp = startarg
        startarg = endarg
        endarg = tmp
        if switch == True and abs(atan) > np.pi/2:
            #print("breeze")
            startarg = startarg - np.pi/2
            endarg = startarg + np.pi
            startarg = startarg - atan/2
            endarg = startarg+ np.pi
    #print("startarg = %f,endarg = %f"%(startarg , endarg))
    if poriginx != 0 and poriginy != 0:
        for i in range(seeds):
            slice = (i+1) *  radius * (1/seeds)
            d = np.sqrt(radius ** 2  - slice ** 2 )
            theta = np.pi/2 -math.acos(slice/radius)
            x = originx + radius * math.cos(startarg  + theta)
            y = originy + radius  * np.sin(startarg + theta)
            x3 = originx + radius * np.cos(endarg -theta)
            y3 = originy + radius * np.sin(endarg - theta)
            #print("x= %f,y = %f" %(x,y))
            #plt.plot(x,y,'.')
            #plt.plot(x3,y3,'.')
            a,b = bisection(x,y,x3,y3,cmap)
            plt.plot(a,b,'.')
            if a != np.nan and a != 0 and b != 0 and b != np.nan:
                orbit[0].extend([orbitnumber])
                orbit[1].extend([a])
                orbit[2].extend([b])
            count += 1
    else:
        seeds = int(radius  * 7000)
        for i in range(seeds):
            slice = (i+1) *  radius * (1/seeds)
            theta = np.pi/2 -math.acos(slice/radius)
            x = originx - radius * np.cos(theta)
            y = originy - radius  * np.sin(theta)
            x3 = originx - radius * np.cos(theta + np.pi)
            #print("x= %f,y = %f" %(x,y))
            #plt.plot(x,y,'.')
            #plt.plot(x3,y,'.')
            a,b = bisection(x,y,x3,y,cmap)
            plt.plot(a,b,'.')
            if a != np.nan and a != 0 and b != 0 and b != np.nan:
                orbit[0].extend([orbitnumber])
                orbit[1].extend([a])
                orbit[2].extend([b])
            count += 1
        for i in range(seeds):
            slice = (i+1) *  radius * (1/seeds)
            theta = np.pi/2 - math.acos(slice/radius)
            x = originx + radius * np.cos(theta)
            y = originy + radius  * np.sin(theta)
            x3 = originx + radius * np.cos(theta + np.pi)
            #print("x= %f,y = %f" %(x,y))
            #plt.plot(x,y,'.')
            #plt.plot(x3,y,'.')
            a,b = bisection(x,y,x3,y,cmap)
            plt.plot(a,b,'.')
            if a != np.nan and a != 0 and b != 0 and b != np.nan:
                orbit[0].extend([orbitnumber])
                orbit[1].extend([a])
                orbit[2].extend([b])
            count += 1

    #if tMap(a + 1j * b,initp,step,cmap)[1].real > 1000 : #p,theta = tMap(theta,p,step,cmap):
        #    epsilon = a
        #    originx = x1
        #    originy = y1
        #print("atan= %f" %atan)
    count += 1
    return orbit
####################################
k,d=1.2,5
cmap =ModifiedKMap(k,d)
seed = [[],[]]
p=initp

#draw(p,cmap) #Mset描きたい時は使ってください.その時は下のloadtxtはoff
#v(p,cmap)

seed = np.loadtxt('Msetimp0.txt')#ここに読み込みたいMsetのファイルを持ってくる。
fig = plt.figure(figsize = (6,6))
ax = fig.add_subplot(1,1,1)
ax.plot(seed[0],seed[1],',k')

x1=0
y1=0

plt.connect('button_press_event',Press)
plt.connect('button_release_event',Release)

plt.show()
