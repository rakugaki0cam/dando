# 弾道計算プログラム
#  
#     BB弾の弾道を求める
#     ホップ回転による揚力計算
#     空気抵抗係数は補正して使用
#
#     追加　着弾までの時間を求める
#
#        2022/08/29
#           カム
#
#  pythonへ移植
#  
#
# original by
# sikinote http://slpr.sakura.ne.jp/qp/
# Author : sikino
#        : 2016/08/06
#

import numpy as np
import math
import matplotlib.pyplot as plt


### 諸元 #####
g = 9.80665             # 重力[m/sec^2]
mBb = 0.28 * 1e-3       # BB弾質量[kg] (g)
rBb = 5.95 / 2 * 1e-3   # BB弾半径[m]  (φmm)
v0 = 80.0               # 初速[m/sec]
hop = 200               # ホップ回転数[rps]
## 気象
temp = 25.0             # 気温[°C]
pres = 1013.25 * 100    # 気圧[Pa]  (hPa)
humi = 60.00 / 100      # 湿度[%]   (%RH)
## マトまでの距離
xTarget = 7.5           # 水平距離[m] 
## 射出時の条件
theta = 0.0             # 射出角度[°]
## 風速
vWind = 0.0             # 風速[m/s]
dWind = 3               # 風向き[時計の短針]

## 風速
sWind = np.deg2rad(dWind * 30)   # 3時=90° 右から左、9時=270° 左から右、6時=180° 追い風、12時=360° 向かい風
ux = -vWind * math.cos(sWind)    # 追い風-、向かい風+[m/sec]
uy = -vWind * math.sin(sWind)    # 左からの風-、右からの風+[m/sec]
uz = 0                           # 下から上への風+、上から下への風-[m/sec]
## 初期位置
rx = 0.0                # 距離[m]
ry = 0.0                # 左右[m]
rz = 1.0                # 高さ[m]
## 初期速度
vy = 0                     # 左右方向のブレ[m/sec]
vx = math.sqrt(v0 ** 2 - vy ** 2) * math.cos(math.pi * theta / 180)
                           # 水平方向の初速[m/sec]
vz = math.sqrt(v0 ** 2 - vy ** 2) * math.sin(math.pi * theta / 180)
                           # 鉛直方向の初速[m/sec]
## 回転による周速度/r(半径)
omgx = 0.0 * 2 * math.pi   # カーブ
omgy = -hop * 2 * math.pi  # ホップ
omgz = 0.0 * 2 * math.pi   # ライフリング

## 空気抵抗係数のフィッティング式
fCd = "Morrison"
#fCd = "Clift&Gauvin"
## 空気抵抗係数の実験補正値
kCd = 1.0
## 回転数減衰の計算方法
ik = 0                     #0:積分計算　1:近似計算
## 時間
h = 0.001                  #計算時間刻み初期値[sec]
te = 5.0                   #終了時間[sec]
##ルンゲクッタの計算精度
# 1e-2  < tol           : 粗い　計算は速い
# 1e-10 < tol < 1e-2    : 普通
# 1e-12 < tol < 1e-10   : 精密　計算は遅い
#         tol < 1e-12   : 非推奨
tol = 1e-8



def rho_humid(T, P, H):
   # ρ : 空気密度 [kg/m^3]
   # T : 気温 [°C]
   # P : 気圧 [Pa] (1気圧 = 1013.25 hPa)
   # H : 湿度 [％RH]
   psT = 6.1078 * 10 ** (7.5 * T / (T + 237.3))              # 飽和水蒸気圧 psT[hPa]
   psT *= 100 * H
   #print("et ", et, "Pa")
   rhoair = 0.0034856447 * P / ( T + 273.15 - 0.670)        ##あとで調べる
   rho_humid = rhoair * (1 - 0.378 * psT / P)
   #print(rho_humid)
   return rho_humid 


def eta_air(T):
   #サザーランドの式により粘性係数を求める
   # η : 空気粘度　[Pa s]
   # T : 気温 [°C]
   # 式も2つあり、値もいろいろある   20°C 1.82e-6 [Pa s]程度
   eta_air = 1.487e-6 * ((T + 273.15) ** 1.5) / ((T + 273.15) + 117)
   #print(eta_air)
   #eta_air = 1.458e-6 * ((T + 273.15) ** 1.5) / ((T + 273.15) + 110.4)
   #print(eta_air)
   #eta_air = 1.7894e-5 * ((T + 273.15)/288.15) ** 1.5 * ((288.15+110.4)/(T+273.15+110.4))
   #print(eta_air)
   return eta_air


def rkf45(t, h, x, xe):
   # ルンゲ・クッタ=フェールベルグ法
   # 刻み幅自動制御
   # t:    時刻(+)
   # x[]:  微分方程式の解
   # h:    刻み幅(+)
   # xe:   距離終点値
   # 戻り値
   # info = -2  異常　計算範囲をオーバー
   #      = -1  異常　刻み値が小さくなりすぎ　tolの再検討が必要
   #      =  0  正常計算進行中
   #      =  1  xe終点に到達し、精度調整中
   #      =  2  xe終点値に到達し計算完了

   hMin = 1e-14
   hMax = 1

   info = 0
   loopFlag = 1
   
   #精度値の計算
   #Sy = scalar(x, N)   ## すべての階について評価
   Sy = scalar(x, 3)    ##距離だけを評価　{x,y,z} = {x[0],x[1],x[3]}
   if Sy > 1:
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

   x4 = np.zeros(N)
   y4 = np.zeros(N)
   y5 = np.zeros(N)
   K = np.zeros((N, S))
   R = np.zeros(N)
   Rcnt = 0       #刻み幅をきめる時の回数
   endCnt = 0     #終点位置繰り返しの回数

   #計算ループ
   while loopFlag == 1:
      x4 = x[0:N]    #スタート値をコピー
      #係数Kの計算
      for s in range(S):
         tx = t + c[s] * h
         for n in range(N):
            if s == 0:
               y4[n] = x[n]      #c[0]==0ではないときには成り立たない
            else:
               y4[n] = x[n] + c[s] * K[n][s - 1]
         for n in range(N):      
            K[n][s] = h * rkfd(y4, n)
      #print("t=",t, "h=",h, end=' ')

      #4次での解と5次での解の差Rを求める
      # R = |x(4次)-x(5次)| / hN
      R = 0
      r = np.zeros(N)
      #for n in range(N):        ##全ての階を対象にする
      for n in range(0,3):       ###刻み幅修正計算の対象を位置だけにする
         for s in range(S):
            r[n] += Rb[s] * K[n][s]
         R += r[n] ** 2
      #R = abs(math.sqrt(R) / N / h)      ##全ての階を対象にする
      R = abs(math.sqrt(R) / 3 / h)       ###刻み幅修正計算の対象を位置だけにする
      #print("R:",r)

      if Rcnt >= 1:
         print('刻み幅変更{:3d}回目  差:R ={:13.10f}  h ={:10.7f}'.format(Rcnt, R, h))
         dummy = ""
      Rcnt += 1

      if R <= err:
         #4次での微分方程式の解を計算
         #for n in range(N):  ##5次の計算
         #  y5[n] = y[n]
         Rcnt = 0
         t2 = t + h  
         for n in range(N):
            for s in range(S):
               x4[n] += b[0][s] * K[n][s]   #4次
               #y5[n] += b[1][s] * K[n][s] #5次　比較用
            #print('y4[{:1d}] ={:20.16f}  y5[{:1d}] ={:20.16f} 差={:20.16f}'.format(n, y[n], n, y5[n], y[n] - y5[n]))

         #終点を過ぎたかの判定
         dx = xe - x4[0]
         if dx <= 0:
            #終点を過ぎた時
            endCnt += 1
            print("終点到達", endCnt, "回目     x=", x4[0],"  Δx=",dx, "   h:",h )
            if abs(dx) < 0.0001:
               # 距離誤差が0.1mm以下になったら終了
               info = 2
               loopFlag = 0
               break
            else:
               # まだ誤差が大きい時は再計算
               h *= 0.5    #刻み幅を狭くして再計算
               loopFlag = 1
               continue
         elif endCnt > 0:
            #刻み幅を変えて再計算した後、終点まで届かなくなった場合
            #一度メインループへ戻り、データを表示させる
            endCnt += 1
            print("終点到達", endCnt, "回目     x=",x4[0],"  Δx=",dx, "   h:",h )
            print("終点に届かなくなった")
            h *= 1.5    #刻み幅を広くする
            info = 1
            loopFlag = 0
            continue
         else:
            #終点に届いていない時
            info = 0
            loopFlag = 0
            #刻み幅修正へすすむ
      
      #刻み幅hを修正  
      if R >= 1e-20:
         delta = (err / (2 * R)) ** (1 / 4)
         if delta > 0.1:
            #刻み幅を適正値に修正する
            h *= delta
         else:
            #現在の刻み幅hは大きすぎるので1/10に
            h *= 0.1
      else:
         #Rがほぼ0の時は刻み幅は小さすぎるので4倍に
         h *= 4
         
      #刻み幅のチェック
      if h > hMax:
         #刻み幅がmaxを超えないように
         h = hMax
      if h < hMin:
         #刻み幅が小さくなり過ぎた時は中断
         #コンピュータの計算誤差に埋もれて正常に計算できなくなるため
         info = -1
         loopFlag = 0

   return t2, x4, h, info


def rkfd(x, n):
   # 玉の運動方程式
   #  微分方程式
   #  dx/dt = fn(t,x)
   #
   # t:  時刻
   # x[]:  値
   # n:  計算する階
   # 戻り値
   # fn: 解

   r   = np.zeros(3)    #位置 rx,ry,rz
   v   = np.zeros(3)    #速度 vx,vy,vz
   omg = np.zeros(3)    #回転数 ωx,ωy,ωz
   u   = np.zeros(3)    #風速 ux,uy,uz
   relv = np.zeros(3)   #相対速度　風を考慮した速度
   
   for i in range(3):
      r[i]   = x[i]        #0,1,2
      v[i]   = x[3 + i]    #3,4,5
      omg[i] = x[6 +i]     #6,7,8
      u[i]   = x[9 + i]    #9,10,11
      relv[i] = v[i] - u[i]   # 風がある時の相対速度
 
   if n < 3:
      # 0,1,2
      # 位置の変化=速度（等速運動）
      # d {x,y,z}/dt = v{x,y,z}
      fn = v[n]

   elif n < 6:
      # 3,4,5
      # 速度の変化=加速度(抵抗力による減速度)
      # d v_{x,y,z}/dt = F{x,y,z}
      vn = n - 3  #v[n] = x[vn]
      fn = (gravity(vn) + vis1(vn, relv) + vis2(vn, relv) + mag(vn, omg, relv)) / mBb

   elif n < 9:
      #6,7,8
      # 回転角の変化（回転の減衰）
      # d omega{x,y,z}/dt = N_{x,y,z}/I
      I = 0.4 * mBb * rBb ** 2
      nv = scalar(relv, 3)
      nw = scalar(omg, 3)
      if nw <= 1e-13:
         fn = 0
      else:
         if ik == 0:
            #積分計算
            fn = (Nz(nv, nw) / I) * x[n] / nw
         else:
            #近似計算
            fn = (Nze(nv, nw) / I) * x[n] / nw  

   else:
      #9,10,11
      # 風の計算 = 相対速度
      # d u{x,y,z}/dt = Fu_{x,y,z}
      fn = 0
 
   return fn


#-----------------------------
def scalar(x, n = 3):
   #スカラーを求める
   # x: ベクトル
   # n: 次元 (デフォルト: 3)
   sum = 0
   for i in range(n):
      sum += x[i] ** 2
   scalar = math.sqrt(sum)   
   return scalar   


def energy(m, v):
   #初速エネルギ
   # v[x,y,z]: vx,vy,vz
   energy =  m * (v[0] ** 2+ v[1] ** 2 + v[2] ** 2) / 2
   return energy


def gravity(dir):
   # 重力 -mg z軸方向のみ
   # dir: 0,1,2 = x,y,z 
   if dir == 2:
      gravity = -mBb * g
   else:
      gravity = 0
   return gravity


def vis1(dir, v):
   # 粘性抵抗分
   vis1 = -6 * math.pi * eta * rBb * v[dir]
   return vis1


def vis2(dir, v):
   # 空気抵抗分
   nv = scalar(v)#####################
   vis2 = -1 / 2 * Cd(Reynolds(nv, 2 * rBb)) * rho * math.pi * rBb ** 2 * nv * v[dir]
   return vis2


def mag(dir, omg, v):
   # 揚力分
   mag = 0
   nomg = scalar(omg, 3)
   if nomg < 1e-14:
      return mag
  
   nv = scalar(v, 3)
   L = np.zeros(3)
   L[0] = v[1] * omg[2] - v[2] * omg[1]
   L[1] = v[2] * omg[0] - v[0] * omg[2]
   L[2] = v[0] * omg[1] - v[1] * omg[0]
   nL = scalar(L, 3)
   if nL < 1e-14:
      return mag
  
   Cl = 0.12   # 揚力係数
   mag = - 4 / 3 * Cl * math.pi * rBb ** 3 * 2 * rho * nomg * nv * L[dir] / nL
   return mag


def Reynolds(nv, d):
   #レイノルズ数
   # nv : norm of velocity of object
   # d  : diameter

   #keta means Kinetic viscosity
   keta = eta / rho
   Reynolds = nv * d / keta
   return Reynolds


def Cd(Re):
   #空気抵抗係数
   # Re : レイノルズ数
   # 空気抵抗力: F = 1/2 Cd ρπ R^2 |V|^2

   if fCd == 'Morrison':
      #空気抵抗係数　Morrisonの式
      #乱流域を含めた広範囲な係数を求められる
      #BB弾での領域Re=20000〜40000程度では少し小さめの値となっているよう
      # From
      #http://www.chem.mtu.edu/~fmorriso/DataCorrelationForSphereDrag2013.pdf  
      c1 = 24 / Re
      c2 = Re / 5
      c2 = 2.6 * c2 / (1 + c2 **1.52)
      c3 = Re / 263000
      c3 = 0.411 * c3 ** (-7.94) / (1 + c3 ** (-8))
      c4 = Re ** 0.8 / 461000
      Cd = c1 + c2 + c3 + c4
      #print("Morrison Cd :", Cd)

   elif fCd == 'Clift&Gauvin':
      #空気抵抗係数　　Clift and Gauvinの式
      # 層流域に限定　Re<300000
      # マッハ数  Ma = v/c 音速に対する比率　BB弾では0.3以下
      # なので、マッハ数による補正無し
      c1 = 24 / Re * (1 + 0.15 * Re ** 0.687)
      c2 = 0.42 / (1 + 42500 * Re ** -1.16)
      Cd = c1 + c2
      #print("C&G Cd =",Cd)

   # 補正値をかける
   Cd = Cd * kCd
   #print("Cd =", Cd, "k =",kCd)

   return Cd


def Nz(nv, omg):
   #積分計算
   #Moment of omg direction
   Nz = 1 / 2 * rho * Cf(nv) * rBb ** 3 * Fintegral(nv, omg)
   return Nz


def Nze(nv, omg):
   #近似計算
   pc = math.pi / 5.32065     # magic phi
   tc = math.pi / 3.60475     # magic theta
  
   vu = nv * math.sin(pc) - rBb * omg * math.sin(tc)
   vd = -nv * math.sin(pc) - rBb * omg * math.sin(tc)
   Nze = -0.5 * rho * Cf(nv)
   Nze *= (4 * math.pi * rBb ** 2) * rBb * 0.5
   Nze *= -(abs(vu) * vu + abs(vd) * vd)
   return Nze


def Cf(nv):
   #摩擦抗力係数
   # nv : norm of velocity of object
   
   # 層流の時　Re < 5e+6
   Cf = 1.328 / math.sqrt(Reynolds(nv, 2 * rBb))
   # 乱流の時 Re > 10^7
   #Cf = 0.455 / (log10(Reynolds(nv))**(2.58))
   return Cf


### 積分計算　#####

def Fintegral(u, omega):
   #積分計算
   # 2016/07/23
   #   2π      π
   #  /       /
   #  | dφ    |  |u sinφ - Rωsinθ| (u sinφ - Rωsinθ) sin^2 θ dθ
   #  /       /
   #  0      0

   xx = np.zeros(15)
   ww = np.zeros(15)
   ss = 0
   xx, ww = GaussKronrod15ab(0, math.pi)
   for i in range(15):
      ss += ww[i] * Fphi(u, omega, xx[i])
  
   xx, ww = GaussKronrod15ab(math.pi, 2 * math.pi)
   for i in range(15):
      ss += ww[i] * Fphi(u, omega, xx[i])
       
   return ss


def GaussKronrod15ab(a, b):
   #Gauss-Kronrod Quadrature Nodes and Weights
   #http://www.advanpix.com/2011/11/07/gauss-kronrod-quadrature-nodes-weights/
   xx = np.zeros(15)
   ww = np.zeros(15)
 
   xx[7] = 0
   xx[8] = 2.077849550078984676006894037732449e-1
   xx[9] = 4.058451513773971669066064120769615e-1
   xx[10] = 5.860872354676911302941448382587296e-1
   xx[11] = 7.415311855993944398638647732807884e-1
   xx[12] = 8.648644233597690727897127886409262e-1
   xx[13] = 9.491079123427585245261896840478513e-1
   xx[14] = 9.914553711208126392068546975263285e-1  
   
   ww[7] = 2.094821410847278280129991748917143e-1
   ww[8] = 2.044329400752988924141619992346491e-1
   ww[9] = 1.903505780647854099132564024210137e-1
   ww[10] = 1.690047266392679028265834265985503e-1
   ww[11] = 1.406532597155259187451895905102379e-1
   ww[12] = 1.047900103222501838398763225415180e-1
   ww[13] = 6.309209262997855329070066318920429e-2
   ww[14] = 2.293532201052922496373200805896959e-2
   for i in range(7):
      xx[i] = -xx[14 - i]
      ww[i] = ww[14 - i]

   xx = 0.5 * ((b - a) * xx + (a + b))
   ww = 0.5 * (b - a) * ww
   #print(xx)
   #print(ww)

   return xx, ww


#----------------------------------------
def Fphi(u, omega, phi):
   # using analysis solution.
   #    π
   #   /
   #   |  |u sinφ - Rωsinθ| (u sinφ - Rωsinθ) sin^2 θ dθ
   #   /
   #   0
 
   a = (u / (rBb * omega)) * math.sin(phi)
 
   if a <= 0:
      Fphi = -(0.5 * a * a * math.pi - 8 * a / 3 + 3 * math.pi / 8)
   elif a >= 1:
      Fphi = 0.5 * a * a * math.pi - 8 * a / 3 + 3 * math.pi / 8
   elif a > 0 and a < 1:
      tp = math.asin(a)
      b = -a * a * math.pi * 0.5 - 8 * a / 3 - 3 * math.pi / 8 + (2 * a * a + 3 / 2) * tp
      c = -(a * a + 1) * math.sin(2 * tp ) + math.sin(4 * tp) / 8 + 6 * a * math.cos(tp) - 2 * a * math.cos(3 * tp) / 3
      Fphi = b + c
   else:
      print("undefined at Fphi, program stop", a)
      while 1:
         exit

   Rw = rBb * omega
   Fphi *= Rw * abs(Rw)

   return Fphi






#####  main ###################################



N = 12   #微分方程式の階数
x = [rx, ry, rz, vx, vy, vz, omgx, omgy, omgz, ux, uy, uz]  #v[0]~v[11]
v = [vx, vy, vz]
Ene = energy(mBb, v)
eta = eta_air(temp)
rho = rho_humid(temp, pres, humi)

print("###########################################################################")
print("# 弾道計算")
#print("# BB弾直径:          {:5.3f} mm".format(2 * rBb * 1000))
print("# BB弾質量:          {:5.3f} g".format(mBb * 1000))
print("# 気温:              {:5.2f} °C".format(temp))
print("# 湿度:              {:5.2f} %RH".format(humi * 100))
print("# 気圧:            {:5.2f} hPa".format(pres / 100))
#print("# 重力:              {:5.3f} m/sec^2".format(g))
#print("# 空気密度:          {:5.3f} kg/m^3".format(rho))
#print("# 空気粘性率:    {:.3e} kg/ms".format(eta))
print("# 初速:             {:6.2f} m/sec".format(v0, "[]"))
print("# エネルギ:          {:5.3f} J".format(Ene))
print("# ホップ回転数:      {:5.1f} rps".format(abs(omgy / 2 / math.pi)))
print("# マト距離:         {:5.3f} m".format(xTarget))
print("# 計算精度:       {:.2e} ".format(tol))
print("# 空気抵抗係数:   {} のフィッティング式による".format(fCd))
print("# 空気抵抗補正:      {:5.3f} ".format(kCd))
print()
print("計算回数     時刻     水平距離    着弾高さ      玉速度   ホップ回転数   エネルギ   計算刻み幅")
print("          t[msec]         x[m]      Δz[mm]     vx[m/s]      ωy[rps]        E[J]     Δt[msec]")


##### 表示サブルーチン
def flightData(t, x):
   #飛翔中のデータ表示
   #風のベクトル　まだ未反映
   v = [x[3], x[4], x[5]]         
   Ene = energy(mBb, v)
   Hop = -x[7] / 2 / math.pi
   print("{:6d}  ".format(i), end = '')
   print("{:9.3f}    {:9.4f}    {:+8.2f}      ".format(t * 1000, x[0], (x[2] - rz) * 1000), end = '')
   print("{:6.2f}       {:6.1f}      {:6.3f}     ".format(x[3], Hop, Ene), end = '')
   print("{:8.4f}".format(h * 1000))
   return


def impactData(text, x):
   #着弾、落下時の角度表示
   ke = x[3] / x[5]
   si = np.rad2deg(math.atan(1 / ke))         
   print("           ^^----- {}  角度{:6.1f}° = 1/{:5.1f} (z/x) -----------------------------".format(text, si, abs(ke)))
   return


#####
time = []
cal = []
gx = []
gz = []

t = 0#####
step = 100  #表示周期  100###############
impFlag = 0

for i in range(999999):
   info = 0
   t, x, h, info = rkf45(t, h, x, xTarget)
   time.append(t)
   cal.append(x[0:11])
   #cal.append([x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10],x[11]]) 
   gx.append(x[0])
   gz.append(x[2])
   
   #着弾した時
   #if impFlag == 0 and x[0] >= xTarget:   # x = x[0]
   if impFlag == 0 and info == 2:   # x = x[0]
      #着弾距離に達した時に一度だけ表示
      text = "マトへ着弾"
      flightData(t, cal[i])
      impactData(text, x)
      impFlag = 1
      break             ########## 着地まで見るときはコメントアウトする

   #着地した時
   if x[2] <= 0:   # z = x[2]
      text = "地面に落下"
      flightData(t, x)
      impactData(text, x)
      break
   
   #データを表示
   if i % step == 0 or info == 1:
      #表示周期毎と終点精度調整中に表示
      flightData(t, x)

   if t >= te:
      #時間切れ
      print("タイムオーバー")
      break


#グラフの表示
fig, ax = plt.subplots(figsize = (10, 4))

ax.plot(gx, gz)
#ax.set_xlim(0, 50)
#ax.set_ylim(0, 1.4)
plt.grid()
plt.xlabel('x  [m]')
plt.ylabel('z  [m]')
plt.show()

print("stop")
while 1:
   exit
 



  