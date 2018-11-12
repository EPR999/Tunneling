#cording : utf-8

import numpy as np

def dotV(x):  #kicked potential の時間微分
     return k * x * np.exp(-8 * x**2) - e2 * (np.exp(-8 * pow(x - xb, 2)) - np.exp(-8 * pow(x + xb, 2)))
def U(qp): #正の時間方向の写像
    return [qp[0] + qp[1] - 0.5 * dotV(qp[0]), qp[1] - 0.5 * dotV(qp[0]) - 0.5 * dotV(qp[0] + qp[1] - 0.5 * dotV(qp[0]))]
def Ui(qp): #inverse of
     return qp[0] - qp[1] - 0.5 * dotV(qp[0]), qp[1] + 0.5 * dotV(qp[0]) + 0.5 * dotV(qp[0] - qp[1] - 0.5 * dotV(qp[0]))
def f(x):#動作確認用
    return x[0]+x[1]  , x[0]**2

def distance(q,p,q0,p0):#距離を求める関数.
    distance =np.array([(p0[i] - p[i]) * np.conj(p0[i] - p[i]) + (q0[i] - q[i]) * np.conj(q0[i] - q[i]) for i in range(len(p)) ])
    return distance

def periodicpoint(q0,p0,step):
    periodicpoint =  [[],[]]
    q = np.array([])
    p = np.array([])
    periodic_q0 = np.array([])
    periodic_p0 = np.array([])
    q,p = tMap(q0,p0,step)
    d = distance(q,p,q0,p0)
    a = np.array([])
    a = [ [q0[i],p0[i]] for i in range(len(q))  if d[i] < 10 ** (-5) ]
    periodic_q0 = [a[i][0] for i in range(len(a))]
    periodic_p0 = [a[i][1] for i in range(len(a))]
    print(len(a))
    return periodic_q0,periodic_p0   
    
def Map(q,p): #写像
    q,p = U([q,p])
    return q,p

def tMap(q0,p0,step):#step数だけ写像
    q = q0
    p = p0
    for i in range(dimension):
        q,p  = Map(q,p)
    return q,p    
        
initp = 0 
step = 1
L = 1
N = 10 ** 6
dimension = 2
k = 20
xf = 1.2
xb = 1.0
e2 = k * xf * (np.exp(-8 * xb * (2 * xf - xb)) / (1 - np.exp(-32 * xf * xb)))

def main():
    seeds = 1000
    q0 = np.random.random(seeds)
    p0 = np.full(seeds,0)
    periodic_q,periodic_p = periodicpoint(q0,p0,step)
    print(periodic_q,periodic_p)

if __name__ == '__main__':
    main() 
