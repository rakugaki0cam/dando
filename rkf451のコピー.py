# ルンゲ・クッタ=フェールベルグ法
#   刻み幅自動制御
#
#   連立1階微分方程式
#
# sikinote(https://slpr.sakura.ne/jp/qp) 2020/04/22によります
#

import numpy as np
import math
import matplotlib.pyplot as plt

def grk(n, x, y):
    #計算する関数
    # n:    方程式階数
    # x:    変数
    # y:    解
    f = np.zeros(n)
    ganma = 0.15

    f[0] = y[1]                         #dy/dx = v
    f[1] = -2 * ganma * y[1] - y[0]     #dv/dx = -2 * ganma * v - y
    return f


def dkf45(N, x, y, h, tol):
    # ルンゲ・クッタ=フェールベルグ法
    # 刻み幅自動制御
    # N:    方程式階数
    # x:    変数
    # y:    解
    # h:    刻み幅
    # tol:  精度
    #帰り値
    #info = -2  異常　計算範囲をオーバー
    #     = -1  異常　刻み値が小さくなりすぎ　tolの再検討が必要
    #     =  0  正常計算進行中
    #     =  1  x終了値に到達し計算完了

    #hのチェック
    hmin = 1e-14
    hmax = 1

    info = 0
    flag = 0

    if h > hmax:
        h = hmax
    if h < -hmax:
        h = -hmax    
    
    if abs(x - xe) <= h:
        h = xe - x

    if abs(x - xe) <= hmin:
        info = 1
        flag = 1
        #return x, y, h, info

    #精度値の計算
    Sy = 0
    for n in range(N):
        Sy += y[n] ** 2
    Sy = math.sqrt(Sy)
    if Sy >= 1:
        err = tol * Sy
    else:
        err = tol
    #print("x ={:13.10f} tol ={:13.10f}".format(x, err))

    #ルンゲクッタ法のブッチャーテーブル準備
    #5次6段
    S = 6    #段数
    b = np.zeros((2, S))
    b[0] = (25/216, 0, 1408/2565, 2197/4104, -1/5, 0)           #4次
    b[1] = (16/135, 0, 6656/12825, 28561/56430, -9/50, 2/55)    #5次
    c = (0, 1/4, 3/8, 12/13, 1, 1/2)
    #Rb = b[1] - b[0]
    Rb = (1/360, 0, -128/4275, -2197/75240, 1/50, 2/55)         # 検算用

    ty = np.zeros(N)
    fn = np.zeros(N)
    K = np.zeros((N, S))
    R = np.zeros(N)
    #y5 = np.zeros(N)    #.astype(np.float64)
    #Rcnt = 0

    #計算ループ
    while flag == 0:
        #係数Kの計算
        for s in range(S):
            tx = x + c[s] * h
            for n in range(N):
                if s == 0:
                    ty[n] = y[n]
                else:
                    ty[n] = y[n] + c[s] * K[n][s - 1]
            fn = grk(N, tx, ty)
            for n in range(N):
                K[n][s] = h * fn[n]
                
        #step 4
        #4次での解と5次での解の差を求める
        # R = |x(4次)-x(5次)| / hN
        R = 0
        r = np.zeros(N)
        for n in range(N):
            for s in range(S):
                r[n] += Rb[s] * K[n][s]
            R += r[n] ** 2
        R = abs(math.sqrt(R) / h / N)
        #Rcnt += 1
        #print('{:3d}回目 刻み幅:h ={:10.7f}  差:R ={:13.10f}'.format(Rcnt, h, R))

        #step 5
        #4次での解を計算
        #for n in range(N):
            #y5[n] = y[n]
        if R <= err:
            #print('R計算回数', Rcnt)
            x += h  
            for n in range(N):
                for s in range(S):
                    y[n] += b[0][s] * K[n][s]   #4次
                    #y5[n] += b[1][s] * K[n][s] #5次　比較用
                #print('y4[{:1d}] ={:20.16f}  y5[{:1d}] ={:20.16f} 差={:20.16f}'.format(n, y[n], n, y5[n], y[n] - y5[n]))
            flag = 1

        #step 6
        if R >= 1e-20:
            delta = (err / (2 * R)) ** (1 / 4)
        else:
            #Rがほぼ0の時
            delta = 4
        #print('delta =', delta)

        #step 7
        if delta <= 0.1:
            #現在の刻み幅hは大きすぎるので1/10にする
            h = 0.1 * h
        elif delta >= 4:
            #現在の刻み幅は小さすぎるので4倍する
            h = 4 * h
        else:
            #刻み幅を適正値に修正する
            h = delta * h

        #step 8
        #刻み幅がmaxを超えないように
        if h > hmax:
            h = hmax
        if h < -hmax:
            h = -hmax

        #step 9
        if abs(xe - x) <= abs(h):
            #hをx終端値に合わせる
            h = xe - x
            if abs(h) <= hmin:
                info = 1
                flag = 1

        #計算が範囲外になっている場合
        if h <= 0 and (xe - x) >= 0:
            info = -2
            flag = 1
        elif h > 0 and (xe - x) <= 0:
            info = -2
            flag = 1
        
    
    return x, y, h, info



# MAIN ------------------------------------------------------------

x = 0
xe = 20     #x終了値
#初期値
N = 2
y = np.zeros(N)
y[0] = 1
y[1] = -0.15

gx = []
gy = []

tol = 1e-8  #計算精度
h = xe - x  #刻み
i = 0
info = 0
step = 200

while info <= 0:
    x, y, h, info = dkf45(N, x, y, h, tol)
    i += 1
    if i % step == 0 or info == 1:     #最終データも表示
        print('{:8d} x:{:9.5f} y:{:13.9f} v:{:13.9f}   h:{:9.6f}'.format(i, x, y[0], y[1], h))
        gx.append(x)
        gy.append(y[0])
        
print('計算回数 {:8d}回'.format(i))    


#グラフ表示
ganma = 0.15
v0 = -0.15
y0 = 1
A = 1       #####
alfa = 0    #####

gx2 = []
gy2 = []
for i in range(0, 20000):
    x = float(i / 1000)
    y =  A * math.exp(-ganma * x) * math.cos(math.sqrt(1 - ganma ** 2) * x - alfa)
    gx2.append(x)
    gy2.append(y)

plt.plot(gx, gy, linestyle = 'None', marker = 'x', markeredgecolor = 'r')
plt.plot(gx2, gy2, color = 'g')
plt.xlabel('x')
plt.ylabel('y')
plt.show()


#  stop
#end program main

