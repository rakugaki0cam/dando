# 弾道計算
#
#        2022/08/29
#
#  pythonへ移植
#
# original by
# sikinote http://slpr.sakura.ne.jp/qp/
# Author : sikino
#        : 2016/08/06
#

import numpy as np
import math
import matplotlib.pyplot as plt


#パラメータ設定例
# 1, tol=1e-10  => とても精密　計算は遅い                     
# 2, tol=1e-6   => 十分に実用的な計算結果
#
#ルンゲクッタの計算精度　(tol は固定値)
# 1e-2  < tol           : 粗い　計算は速い
# 1e-10 < tol < 1e-2    : 普通
# 1e-12 < tol < 1e-10   : 精密　計算は遅い
#         tol < 1e-12   : 非推奨
#tol = 1e-10       #桁精度



def rho_humid(T, P, H):
   # ρ : 空気密度 [kg/m^3]
   # T : 気温 [°C]
   # P : 気圧 [Pa] (1気圧 = 101325 Pa)
   # H : 湿度 [％RH]

   et = 6.1078 - 10 ** (7.5 * T / (T + 237.3))              # et[hPa]
   et *= 100 * H
   #print("et ", et, "Pa")
   rhoair = 0.0034856447 * P / ( T + 273.15 - 0.670)        ##あとで調べる
   rho_humid = rhoair * (1 - 0.378 * et / P)

   return rho_humid 


def eta_air(T):
   # η : 空気粘度　[kg/ms]
   # T : 気温 [°C]

   eta_air = 1.487e-6 * ((T + 273.15) ** 1.5) / ((T + 273.15) + 117)       ##あとで調べる
   #eta_air = 1.458e-6 * ((T + 273.15) ** 1.5) / ((T + 273.15) + 110.4)

   return eta_air

#TEST-----------------------------------------------------------------
#rho_humid(20,1013.25,60)
#eta_air(20)



def rkf451_e(x, y, h, xe):
   # ルンゲ・クッタ=フェールベルグ法
   # 刻み幅自動制御
   # x:    変数
   # y:    解
   # h:    刻み幅
   # xe:   変数終値
   #帰り値
   #info = -2  異常　計算範囲をオーバー
   #     = -1  異常　刻み値が小さくなりすぎ　tolの再検討が必要
   #     =  0  正常計算進行中
   #     =  1  x終了値に到達し計算完了

   #h刻み値のチェック
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

      K1 = np.zeros(N)
      K2 = np.zeros(N)
      K3 = np.zeros(N)
      K4 = np.zeros(N)
      K5 = np.zeros(N)
      K6 = np.zeros(N)
      ty0 = np.zeros(N)
      
      for n in range(N):
         tx0 = x
         ty0[n] = y[n]
         K1[n] = h * rkfd(tx0        , ty0, n)                #c[0]=0
         ty0[n] = y[n]+1/4*  K1[n]
         K2[n] = h * rkfd(tx0+1/4*  h, ty0, n )   #c[1]=1/4
         ty0[n] = y[n]+3/8*  K2[n]
         K3[n] = h * rkfd(tx0+3/8*  h, ty0, n )   #c[2]=3/8
         ty0[n] = y[n]+12/13*K3[n]
         K4[n] = h * rkfd(tx0+12/13*h, ty0, n )   #c[3]=12/13
         ty0[n] = y[n]+      K4[n]
         K5[n] = h * rkfd(tx0+      h, ty0, n )   #c[4]=1
         ty0[n] = y[n]+1/2*  K5[n]
         K6[n] = h * rkfd(tx0+1/2*  h, ty0, n )   #c[5]=1/2




      for s in range(S):
         tx = x + c[s] * h
         for n in range(N):
            if s == 0:
               ty[n] = y[n]      #c[s]がゼロでないときは成り立たない
            else:
               ty[n] = y[n] + c[s] * K[n][s - 1]
         for n in range(N):      
            K[n][s] = h * rkfd(tx, ty, n)

      """
      #係数Kの計算 ~9/4
      for s in range(S):
         tx = x + c[s] * h
         for n in range(N):
            if s == 0:
               ty[n] = y[n]
            else:
               ty[n] = y[n] + c[s] * K[n][s - 1]
         for n in range(N):
            K[n][s] = h * rkfd(tx, ty, n)
      """             







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


def rkfd(t, x, n):
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
      relv[i] = v[i] - u[i]
 
   #fn = 0#########
   if n < 3:
      # 0,1,2
      # 位置の変化=速度（慣性力）
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
      # 風の計算
      # d u{x,y,z}/dt = Fu_{x,y,z}
      fn = 0
 
   return fn


#-----------------------------
def scalar(x, n):
   #スカラーを求める
   # x: ベクトル
   # n: 次元=3
   sum = 0
   for i in range(n):
      sum += x[i] ** 2
   scalar = math.sqrt(sum)   
   return scalar   


def energy(x):
   #初速エネルギ
   # x[3,4,5] : vx,vy,vz
   energy =  mBb * (x[3] ** 2 + x[4] ** 2 + x[5] ** 2) / 2
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
   nv = scalar(v, 3)
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
   # From
   #http://www.chem.mtu.edu/~fmorriso/DataCorrelationForSphereDrag2013.pdf
   #Drag coefficient Cd,
   #Cd depend on Reynolds number,Re.
   #Fource of Drug,D is written by
   # D = 1/2 Cd ρπ R^2 |V|^2
 
   c1 = 24 / Re
   c2 = Re / 5
   c2 = 2.6 * c2 / (1 + c2 **1.52)
   c3 = Re / 263000
   c3 = 0.411 * c3 ** (-7.94) / (1 + c3 ** (-8))
   c4 = Re ** 0.8 / 461000
 
   Cd = c1 + c2 + c3 + c4
 
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
      
   Fintegral = ss  
 
   return Fintegral


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
   Fphi = Fphi * Rw * abs(Rw)

   return Fphi




#####  main ###################################

# 諸元
mBb = 0.28 * 1e-3       # BB弾質量[kg] (g)
rBb = 5.95 / 2 * 1e-3   # BB弾半径[m]  (φmm)
v0 = 83.0               # 初速[m/sec]
g = 9.80665             # 重力[m/sec^2]
# 気象
temp = 20.0             # 気温[°C]
pres = 1013.25 * 100    # 気圧[Pa]  (hPa)
humi = 60.0 / 100       # 湿度[%]   (%RH)
## 風速
ux = 0                  # 追い風+、向かい風-[m/sec]
uy = 0                  # 左からの風+、右からの風-[m/sec]
uz = 0                  # 下から上への風+、上から下への風-[m/sec]
#射出時の条件
theta = 0.0             # 射出角度[°]
## 初期位置
rx = 0.0                # 距離[m]
ry = 0.0                # 左右[m]
rz = 1.0                # 高さ[m]
## 初期速度[m/sec]
vy = 0                     # 左右方向のブレ[m/sec]
vx = math.sqrt(v0 ** 2 - vy ** 2) * math.cos(math.pi * theta / 180)
                           # 水平方向の初速[m/sec]
vz = math.sqrt(v0 ** 2 - vy ** 2) * math.sin(math.pi * theta / 180)
                           # 鉛直方向の初速[m/sec]
# 回転速度
omgx = 0.0 * 2 * math.pi   # 回転数[rad/sec]　カーブ
omgy = -210 * 2 * math.pi  # ホップ回転数[rad/sec] (rps)
omgz = 0.0 * 2 * math.pi   # ライフリング[rad/sec]
# 時間
stept = 0.001            #計算時間刻み[sec]
te = 3                   #終了時間[sec]
# 回転数減衰の計算方法
ik = 0                  #0:積分計算　1:近似計算
# 計算精度
tol = 1e-6              #桁精度


N = 12 #微分方程式の階数

eta = eta_air(temp)
rho = rho_humid(temp, pres, humi)
 
x = [rx, ry, rz, vx, vy, vz, omgx, omgy, omgz, ux, uy, uz]  #v[0]~v[11]
Ene = energy(x)

Nt = int(te / stept)
time = np.zeros(Nt + 2)
for i in range(Nt + 2):
   time[i] = i * stept

step = 100  #表示周期 

print("弾道計算")
print("# m:        ", mBb * 1000, "[g]")
print("# g:        ", g, "[m/sec]")
print("# dia:      ", 2 * rBb * 1000, "[mm]")
print("# temp:     ", temp, "[°C]")
print("# humidity: ", humi * 100, "[%RH]")
print("# press:    ", pres / 100, "[hPa]")
print("# eta:      ", eta, "[kg/ms]")
print("# rho:      ", rho, "[kg/m^3]")
print("# tolerance:", tol)
print("# v0:       ", v0,"[m/sec]")
print("# Energy:   ", Ene,"[J]")


print("  t[s]       x[m]         z[m]       vx[m/s]       wy[rot/s]      Energy[J]")

info = 0
gx = []
gz = []

for j in range(Nt + 1):
   t = time[j]
   te = time[j + 1]
   h = te - t
   t, x, _, info = rkf451_e(t, x, h, te)

   gx.append(x[0])
   gz.append(x[2])

   if j % step == 0:          
      Ene = energy(x)
      Hop = -x[7] / 2 / math.pi
      print("{:5.3f}[sec] x:{:7.3f}[m] z:{:6.3f}[m] vx:{:5.1f}[m/s] ωy:{:6.1f}[rps] E:{:6.3f}[J]".format(t, x[0], x[2], x[3], Hop, Ene))

   #Stop calculation when bullet reach ground.
   if x[2] < 0:   # z = x[2]
      print("Reach ground at ", t, "sec")

      
      fig, ax = plt.subplots(figsize = (10, 4))

      ax.plot(gx, gz)

      ax.set_xlim(0, 50)
      ax.set_ylim(0, 1.4)
      plt.grid()
      plt.xlabel('x  [m]')
      plt.ylabel('z  [m]')

      plt.show()




      while 1:
         exit
 



  