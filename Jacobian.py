import numpy as np

#線形化を行うために写像についてヤコビアンを求める.
#This program evaluate Jacobian of the map at fixed and periodic orbits in order to do linearlization around that point.
k = 20
xf = 1.2
xb = 1.0
e2 = k * xf * (np.exp(-8 * xb * (2 * xf - xb)) / (1 - np.exp(-32 * xf * xb)))
dimension = 2

def dotV(x):  #kicked potential の時間微分
     return k * x * np.exp(-8 * x**2) - e2 * (np.exp(-8 * pow(x - xb, 2)) - np.exp(-8 * pow(x + xb, 2)))
def U(qp): #正の時間方向の写像
    return [qp[0] + qp[1] - 0.5 * dotV(qp[0]), qp[1] - 0.5 * dotV(qp[0]) - 0.5 * dotV(qp[0] + qp[1] - 0.5 * dotV(qp[0]))]
def Ui(qp): #inverse of
     return qp[0] - qp[1] - 0.5 * dotV(qp[0]), qp[1] + 0.5 * dotV(qp[0]) + 0.5 * dotV(qp[0] - qp[1] - 0.5 * dotV(qp[0]))

def derivative(func,x,i):#多変数関数の微分(i番目が微分され戻る.) 参考:blog.amedama.jp
    h = 0.1 ** 4
    h_vec = np.zeros_like(x)
    h_vec[i] = h
    h_vec = h_vec
    xp = x + h_vec
    xm = x - h_vec
    print((np.array(func(xp))-np.array(func(xm)))/(2*h))
    return (np.array(func(xp))-np.array(func(xm)))/(2 * h)

def main():
    x = np.array([0.0,0.0]) #x = [Re(q),Re(p)] 不動点・周期点を受け取る.#This get
    jacobian = [np.array([]),np.array([])]
    for i in range(dimension):#theta is one of the misterious
       jacobian =  [np.append(jacobian[j],derivative(U,x,i)[j]) for j in range(dimension)]
       print(jacobian)
    print(jacobian)#ヤコビアンが求まる.次は固有値を
    w,v = np.linalg.eig(jacobian)#wは固有値,vは固有ベクトル.
    print(v)
    

if __name__ == '__main__':
          main()

