import numpy as np
import matplotlib.pyplot as plt

class ScatteringMap: #made by R.Kogawa
    def __init__(self,xb,k,e2):
        self.xb = xb
        self.k = k
        self.e2 = e2
    def U(self,qp): #正の時間方向の写像
        return [qp[0] + qp[1] - 0.5 * dotV(qp[0]), qp[1] - 0.5 * dotV(qp[0]) - 0.5 * dotV(qp[0] + qp[1] - 0.5 * dotV(qp[0]))]
    def Ui(self,qp): #inverse of
        return [qp[0] - qp[1] - 0.5 * dotV(qp[0]), qp[1] + 0.5 * dotV(qp[0]) + 0.5 * dotV(qp[0] - qp[1] - 0.5 * dotV(qp[0]))]

def dotV(x):  #kicked potential の時間微分
        return k * x * np.exp(-8 * x**2) - e2 * (np.exp(-8 * pow(x - xb, 2)) - np.exp(-8 * pow(x + xb, 2)))

k = 20
xf = 1.2
xb = 1.0
e2 = k * xf * (np.exp(-8 * xb * (2 * xf - xb)) / (1 - np.exp(-32 * xf * xb)))
twopi = 2 * np.pi
step  = 6
seed = 1000
cmap = ScatteringMap(xb,k,e2)

def main():
    twopi = 2 * np.pi
    xb = 1
    step  = 6
    seed = 1000
    q = [   np.random.random() for i in range(1000)   ]
    p = [   np.random.random() for i in range(1000)   ]
    a = [[],[]] #一時的な格納用
    a = [cmap.U([q[i],p[i]]) for i in range(len(q)) ]
    q = [a[i][0] for i in range(len(a))]
    p = [a[i][1] for i in range(len(a))]
    fig = plt.figure(figsize = (6,6))
    ax = fig.add_subplot(111)
    ax.plot(q,p,',k')
    plt.plot()
    plt.show()

if __name__ == '__main__':
        main()
