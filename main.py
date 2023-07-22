# 弾道計算プログラム
#
#     BB弾の弾道を求める
#     ホップ回転による揚力計算
#     空気抵抗係数は補正して使用
#
#     追加　着弾までの時間を求める
#
#        2022/08/29
#           カムのらくがき帳
#
#  pythonへ移植
#
#   2023/07/08 変数名修正
#
#
# original by
# sikinote http://slpr.sakura.ne.jp/qp/
# Author : sikino
#        : 2016/08/06
#

import math
import numpy as np
import matplotlib.pyplot as plt


### 諸元 #####
# 気象
temp = 29.3             # 気温[°C]
humi = 58.5 / 100      # 湿度[%]   (%RH)
pres = 1016.2 * 100    # 気圧[Pa]  (hPa)
# BB
g = 9.80665             # 重力[m/sec^2]
mBb = 0.280 * 1e-3       # BB弾質量[kg] (g)
rBb = 5.95 / 2 * 1e-3   # BB弾半径[m]  (φmm)
v0 = 20.93               # 初速[m/sec]
hop = 194               # ホップ回転数[rps]
# マトまでの距離
v0meas = 0.031          # 初速測定装置中心[m]　センサ1-2の中間位置
v0Correct = 1           # 初速補正あり:1 　センサ1位置での初速を推測補正する
xTarget1 = 2.5        # 水平距離[m]
xTarget2 = 30          # 次が999でストップ
xTarget3 = 50
# 射出時の角度
elevAngle = 0.0         # 射出時の上下角度[°] 仰角＋、俯角ー
lrAngle = 0.0         # 射出時の左右角度[°] 右＋、左ー （必ずゼロ）
tiltAngle = 0.0         # ホップの傾き 0:正常 右＋、左ー
# 風速
vWind = 0.0             # 風速[m/s]
dWind = 3               # 風向き[時] 3時:右からの風
# 初期位置
x0 = 0.0                # 距離[m]
y0 = 0.0                # 左右[m]
z0 = 1.0                # 高さ[m]
# ブレ（初期速度）
vy = 0                  # 左右方向のブレ[m/sec]
# 空気抵抗係数のフィッティング式
#CdMethod = "Morrison"
#CdMethod = "Clift&Gauvin"
CdMethod = "fixed0.43"
# 空気抵抗係数の実験補正値
kCd = 0.979
# 回転数減衰の計算方法
ik = 0  # 0:積分計算　1:近似計算
# 時間
hh = 0.001  # 計算時間刻み初期値[sec]
te = 3.0    # 終了時間[sec]
# ルンゲクッタの計算精度
# 1e-2  < tol           : 粗い　計算は速い
# 1e-10 < tol < 1e-2    : 普通
# 1e-12 < tol < 1e-10   : 精密　計算は遅い
#         tol < 1e-12   : 非推奨
tol = 1e-7


# 初速の位置(水平距離を求める)
xV0meas = v0meas * np.cos(-np.deg2rad(elevAngle)) * \
    np.cos(np.deg2rad(lrAngle))  # 左右角度は必ずゼロ
# 初期速度


def muzzleVelocity(v1, v2):
    # 左右方向のブレ[m/sec]
    v0y = v2
    v0x = math.sqrt(v1 ** 2 - v0y ** 2) * math.cos(math.pi *
                                                   elevAngle / 180)  # 水平方向の初速[m/sec]
    v0z = math.sqrt(v1 ** 2 - v0y ** 2) * math.sin(math.pi *
                                                   elevAngle / 180)  # 鉛直方向の初速[m/sec]
    return v0x, v0y, v0z


vx, vy, vz = muzzleVelocity(v0, vy)

# ホップ回転の角速度  ーー 軸の向きに注意
omgx = 0.0 * 2 * math.pi                                      # ライフリング
omgy = -hop * math.cos(np.deg2rad(tiltAngle)) * 2 * \
    math.pi    # ホップ       y軸まわりの回転
omgz = hop * math.sin(np.deg2rad(tiltAngle)) * 2 * math.pi    # カーブ
# 風速
# 風向き　3時=90° 右から左、9時=270° 左から右、6時=180° 追い風、12時=360° 向かい風
sWind = np.deg2rad(dWind * 30)
ux = -vWind * math.cos(sWind)    # 追い風-、向かい風+[m/sec]
uy = -vWind * math.sin(sWind)    # 左からの風-、右からの風+[m/sec]
uz = 0                           # 下から上への風+、上から下への風-[m/sec]


def rho_humid(T, P, H):
    # ρ : 空気密度 [kg/m^3]
    # T : 気温 [°C]
    # P : 気圧 [Pa] (1気圧 = 1013.25 hPa)
    # H : 湿度 [％RH]
    psT = 6.1078 * 10 ** (7.5 * T / (T + 237.3))              # 飽和水蒸気圧 psT[hPa]
    psT *= 100 * H
    #print("et ", et, "Pa")
    rhoair = 0.0034856447 * P / (T + 273.15 - 0.670)  # あとで調べる
    rho_h = rhoair * (1 - 0.378 * psT / P)
    # print(rho_h)
    return rho_h


def eta_air(T):
    # サザーランドの式により粘性係数を求める
    # η : 空気粘度　[Pa s]
    # T : 気温 [°C]
    # 式も2つあり、値もいろいろある   20°C 1.82e-6 [Pa s]程度
    etaA = 1.487e-6 * ((T + 273.15) ** 1.5) / ((T + 273.15) + 117)
    # print(etaA)
    # etaA = 1.458e-6 * ((T + 273.15) ** 1.5) / ((T + 273.15) + 110.4)
    # print(etaA)
    # etaA = 1.7894e-5 * ((T + 273.15)/288.15) ** 1.5 * ((288.15+110.4)/(T+273.15+110.4))
    # print(etaA)
    return etaA


def rkfd(x, n):
    # 玉の運動方程式
    #  微分方程式
    #  dx/dt = fn(t,x)
    #
    # t:  時刻　（未使用）
    # x[]:  値
    # n:  計算する階
    # 戻り値
    # fn: 解

    p = np.zeros(3)  # 位置 rx,ry,rz
    v = np.zeros(3)  # 速度 vx,vy,vz
    omg = np.zeros(3)  # 回転数 ωx,ωy,ωz
    u = np.zeros(3)  # 風速 ux,uy,uz
    relV = np.zeros(3)  # 相対速度　風を考慮した速度

    for i in range(3):
        p[i] = x[i]  # 0,1,2
        v[i] = x[3 + i]  # 3,4,5
        omg[i] = x[6 + i]  # 6,7,8
        u[i] = x[9 + i]  # 9,10,11
        relV[i] = v[i] - u[i]   # 風がある時の相対速度

    # 呼び出し番号nによる処理式の区分け
    if n < 3:
        # 0,1,2
        # 位置の変化=速度（等速運動）
        # d {x,y,z}/dt = v{x,y,z}
        fn = v[n]
    elif n < 6:
        # 3,4,5
        # 速度の変化=加速度(抵抗力による減速度)
        # d v_{x,y,z}/dt = F{x,y,z}
        nn = n - 3  # v[n] = x[n-3]
        Fg = fGr(nn)
        Fa = fAirReg(nn, relV)
        Fl = fMagnus(nn, omg, relV)
        fn = (Fg + Fa + Fl) / mBb
    elif n < 9:
        # 6,7,8
        # 回転速度の減衰
        # d omega{x,y,z}/dt = N_{x,y,z}/I
        mOmega = magVector(omg, 3)
        if mOmega < 1e-13:
            fn = 0
        else:
            I = 0.4 * mBb * rBb ** 2
            mRelV = magVector(relV, 3)
            if ik == 0:
                # 積分計算
                fn = (Nz(mRelV, mOmega) / I) * x[n] / mOmega
            else:
                #(ik == 1)
                # 近似計算
                fn = (Nze(mRelV, mOmega) / I) * x[n] / mOmega
    elif n < 12:
        # 9,10,11
        # 風の計算 = 相対速度
        # d u{x,y,z}/dt = Fu_{x,y,z}
        fn = 0
    else:
        # error
        fn = 0

    return fn


def rkf45(t, h, x, xe):
    # ルンゲ・クッタ=フェールベルグ法
    # 刻み幅自動制御
    # t:    時刻(+)
    # h:    刻み幅(+)
    # x[]:  微分方程式の解
    # xe:   距離終点値
    #
    # 戻り値
    # info = -2  異常　刻み幅が大きくなり過ぎ
    #      = -1  異常　刻み値が小さくなりすぎ　tolの再検討が必要
    #      =  0  正常計算進行中
    #      =  1  xe終点に到達し、精度調整中
    #      =  2  xe終点値に到達し計算完了

    kai = 3  # 刻み幅修正計算の対象を位置だけにする　{x,y,z} = {x[0],x[1],x[3]}
    # kai = N       ## 全ての階を対象にする
    deBug001 = 0  # 刻み幅変更と終点検出のデバッグ値　 1:表示する、0:表示しない

    hMin = 1e-14
    hMax = 1

    info = 0
    loopFlag = 1

    # 精度値の計算
    Sy = magVector(x, kai)
    if Sy > 1:
        err = tol * Sy
    else:
        err = tol
    #print("x ={:13.10f} tol ={:13.10f}".format(x, err))

    # ルンゲクッタ法のブッチャーテーブル準備
    # 5次6段
    S = 6  # 段数
    b = np.zeros((2, S))
    b[0] = (25/216, 0, 1408/2565, 2197/4104, -1/5, 0)  # 4次
    b[1] = (16/135, 0, 6656/12825, 28561/56430, -9/50, 2/55)  # 5次
    c = (0, 1/4, 3/8, 12/13, 1, 1/2)
    # Rbは5次-4次　Rb = b[1] - b[0]
    Rb = (1/360, 0, -128/4275, -2197/75240, 1/50, 2/55)  # 差

    f4 = np.zeros(N)  # Nはグローバル　階数
    K = np.zeros((N, S))
    R = np.zeros(N)
    Rcnt = 0  # 刻み幅をきめる時の回数
    endCnt = 0  # 終点位置繰り返しの回数

    # 計算ループ
    while loopFlag == 1:
        x4 = np.copy(x)  # スタート値をコピー
        # 係数Kの計算
        for s in range(S):
            # tx = t + c[s] * h  # いらないよう
            for n in range(N):
                if s == 0:
                    f4[n] = x[n]  # c[0]==0ではないときには成り立たない
                else:
                    f4[n] = x[n] + c[s] * K[n][s - 1]
            for n in range(N):
                K[n][s] = h * rkfd(f4, n)       # rkfd(tx, f4, n) tx不要なよう

        # 4次での解と5次での解の差Rを求める
        # R = |x(4次)-x(5次)| / hN
        R = 0
        r = np.zeros(N)
        for n in range(kai):
            for s in range(S):
                r[n] += Rb[s] * K[n][s]
            R += r[n] ** 2
        R = abs(math.sqrt(R) / kai / h)

        def kizami():
            if deBug001 == 1:
                print(f'刻み幅変更{Rcnt:3d}回目  差:R ={R:13.10f}  h ={h:10.7f} sec')
            return

        def shuten(txt=''):
            if deBug001 == 1:
                print(f'終点に到達{endCnt:2d}回目  x ={x4[0]:9.4f}    Δx = {(dx*1000):5.2f} mm     h ={(h*1000):7.4f} ms  ', end='')
                if txt != '':
                    print(txt)
                else:
                    print()
            return

        if Rcnt >= 1:
            kizami()
        Rcnt += 1

        if R <= err:
            # 4次での微分方程式の解を計算
            Rcnt = 0
            t2 = t + h
            for n in range(N):
                for s in range(S):
                    x4[n] += b[0][s] * K[n][s]  # 4次

            # 水平距離が終点を過ぎたかの判定
            dx = xe - x4[0]
            if dx <= 0:
                # 終点を過ぎた時
                endCnt += 1
                shuten()
                if abs(dx) < 0.0001:
                    # 距離誤差が0.1mm以下になったら終了
                    info = 2
                    break
                else:
                    # まだ誤差が大きい時は再計算
                    h *= 0.5  # 刻み幅を狭くして再計算
                    if h < hMin:
                        info = -1  # エラー終了
                        break
                    loopFlag = 1
                    continue
            elif endCnt == 0:
                # 解の計算をしたけれど終点に届いていない時
                loopFlag = 0  # 続く刻み幅修正処理してからリターン
            else:
                # 終点調整後、終点まで届かなくなった場合、一度メインループへ戻り、データを表示させる
                endCnt += 1
                shuten("終点に届かなくなった")
                h *= 1.5  # 刻み幅を広くする（2倍するとひとつ前の刻み幅と同じになる）
                info = 1
                break
        # 刻み幅hを修正
        if R >= 1e-20:  # ゼロ割り算防止
            delta = (err / (2 * R)) ** (1 / 4)
            if delta > 0.1:
                # 刻み幅を適正値に修正する
                h *= delta
            else:
                # 現在の刻み幅hは大きすぎるので1/10に
                h *= 0.1
        else:
            # Rがほぼ0の時は刻み幅は小さすぎるので4倍に
            h *= 4
        # 刻み幅のチェック
        if h > hMax:
            info = -2
            break
        if h < hMin:
            info = -1
            break

    return t2, h, x4, info


# -----------------------------

def magVector(vec, n=3):
    # 三次元ベクトルの大きさを求める
    # vec: ベクトル
    # n:   次元 (デフォルト: 3)
    v = np.array(vec[0:n])
    mag = np.sqrt(np.sum(v ** 2))
    return mag


def energy(m, v):
    # 初速エネルギ
    # v[x,y,z]: vx,vy,vz
    e = m * (v[0] ** 2 + v[1] ** 2 + v[2] ** 2) / 2
    return e


def fGr(direction):
    # 重力 -mg z軸方向のみ
    # direction: 0,1,2 = x,y,z
    if direction == 2:
        Fgravity = -mBb * g
    else:
        Fgravity = 0
    return Fgravity


def fAirReg(direction, v):
    # 空気抵抗力
    # 粘性抵抗分も含まれる
    # dir: 0,1,2 = x,y,z
    # v[dir] = vx,vy,vz
    mV = magVector(v)
    Fairregistance = -1 / 2 * \
        Cd(reynoldsNum(mV, 2 * rBb)) * rho * math.pi * rBb ** 2 * mV * v[direction]
    return Fairregistance


def fMagnus(direction, omg, v):
    # 揚力分
    Fmagnus = 0
    omgS = magVector(omg, 3)
    if omgS < 1e-14:
        return Fmagnus

    vS = magVector(v, 3)
    l = np.zeros(3)
    l[0] = v[1] * omg[2] - v[2] * omg[1]
    l[1] = v[2] * omg[0] - v[0] * omg[2]
    l[2] = v[0] * omg[1] - v[1] * omg[0]
    lS = magVector(l, 3)
    if lS < 1e-14:
        return Fmagnus

    Cl = 0.12   # 揚力係数
    Fmagnus = - 4 / 3 * Cl * math.pi * rBb ** 3 * \
        2 * rho * omgS * vS * l[direction] / lS
    return Fmagnus


def reynoldsNum(vS, d):
    # レイノルズ数
    # nv : norm of velocity of object
    # d  : diameter

    # keta means Kinetic viscosity
    keta = eta / rho
    reynoldsNumber = vS * d / keta
    return reynoldsNumber


def Cd(Re):
    # 空気抵抗係数
    # Re : レイノルズ数
    # 空気抵抗力: F = 1/2 Cd ρπ R^2 |V|^2

    cd = 0

    if CdMethod == 'Morrison':
        # 空気抵抗係数　Morrisonの式
        # 乱流域を含めた広範囲な係数を求められる
        # BB弾での領域Re=20000〜40000程度では少し小さめの値となっているよう
        # From
        # http://www.chem.mtu.edu/~fmorriso/DataCorrelationForSphereDrag2013.pdf
        c1 = 24 / Re
        c2 = Re / 5
        c2 = 2.6 * c2 / (1 + c2 ** 1.52)
        c3 = Re / 263000
        c3 = 0.411 * c3 ** (-7.94) / (1 + c3 ** (-8))
        c4 = Re ** 0.8 / 461000
        cd = c1 + c2 + c3 + c4
        #print("Morrison Cd :", cd)

    elif CdMethod == 'Clift&Gauvin':
        # 空気抵抗係数　　Clift and Gauvinの式
        # 層流域に限定　Re<300000
        # マッハ数  Ma = v/c 音速に対する比率　BB弾では0.3以下
        # なので、マッハ数による補正無し
        c1 = 24 / Re * (1 + 0.15 * Re ** 0.687)
        c2 = 0.42 / (1 + 42500 * Re ** -1.16)
        cd = c1 + c2
        #print("C&G Cd =", cd)

    elif CdMethod == "fixed0.43":
        cd = 0.43
        #print("Fixed Cd =", cd)
    else:
        print("Cd エラー")
        while 1:
            pass

    # 補正値をかける
    cd *= kCd
    #print("Cd =", cd, "k =",kCd)

    return cd


def Nz(vS, omg):
    # 積分計算
    # Moment of omg direction
    nz = 1 / 2 * rho * Cf(vS) * rBb ** 3 * Fintegral(vS, omg)
    return nz


def Nze(vS, omg):
    # 近似計算
    pc = math.pi / 5.32065     # magic phi
    tc = math.pi / 3.60475     # magic theta

    vu = vS * math.sin(pc) - rBb * omg * math.sin(tc)
    vd = -vS * math.sin(pc) - rBb * omg * math.sin(tc)
    nze = -0.5 * rho * Cf(vS)
    nze *= (4 * math.pi * rBb ** 2) * rBb * 0.5
    nze *= -(abs(vu) * vu + abs(vd) * vd)
    return nze


def Cf(vS):
    # 摩擦抗力係数
    # nv : norm of velocity of object

    # 層流の時　Re < 5e+6
    cf = 1.328 / math.sqrt(reynoldsNum(vS, 2 * rBb))
    # 乱流の時 Re > 10^7
    #cf = 0.455 / (log10(Reynolds(nv))**(2.58))
    return cf


### 積分計算　#####

def Fintegral(u, omega):
    # 積分計算
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
    # Gauss-Kronrod Quadrature Nodes and Weights
    # http://www.advanpix.com/2011/11/07/gauss-kronrod-quadrature-nodes-weights/
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
    # print(xx)
    # print(ww)

    return xx, ww


# ----------------------------------------
def Fphi(u, omega, phi):
    # using analysis solution.
    #    π
    #   /
    #   |  |u sinφ - Rωsinθ| (u sinφ - Rωsinθ) sin^2 θ dθ
    #   /
    #   0

    a = (u / (rBb * omega)) * math.sin(phi)

    if a <= 0:
        fphi = -(0.5 * a * a * math.pi - 8 * a / 3 + 3 * math.pi / 8)
    elif a >= 1:
        fphi = 0.5 * a * a * math.pi - 8 * a / 3 + 3 * math.pi / 8
    elif a > 0 and a < 1:
        tp = math.asin(a)
        b = -a * a * math.pi * 0.5 - 8 * a / 3 - 3 * \
            math.pi / 8 + (2 * a * a + 3 / 2) * tp
        c = -(a * a + 1) * math.sin(2 * tp) + math.sin(4 * tp) / \
            8 + 6 * a * math.cos(tp) - 2 * a * math.cos(3 * tp) / 3
        fphi = b + c
    else:
        print("undefined at Fphi, program stop", a)
        while 1:
            pass

    Rw = rBb * omega
    fphi *= Rw * abs(Rw)

    return fphi


#####  main ###################################

N = 12  # 微分方程式の階数
xArray = np.array([x0, y0, z0, vx, vy, vz, omgx, omgy, omgz, ux, uy, uz])  # 微分方程式の解（初期値）を代入
Ene = energy(mBb, xArray[3:6])
eta = eta_air(temp)
rho = rho_humid(temp, pres, humi)

print("###########################################################################")
print("# 弾道計算")
# print(f"# BB弾直径:          {(2 * rBb * 1000):5.3f} mm")
print(f"# BB弾質量:          {(mBb * 1000):5.3f} g")
print(f"# 気温:              {temp:5.2f} °C")
print(f"# 湿度:              {(humi * 100):5.2f} %RH")
print(f"# 気圧:            {(pres / 100):5.2f} hPa")
# print(f"# 重力:              {g:5.3f} m/sec^2")
# print(f"# 空気密度:          {rho:5.3f} kg/m^3")
# print(f"# 空気粘性率:    {eta:.3e} kg/ms")
print(f"# 初速:             {v0:6.2f} m/sec")
# print(f"# エネルギ:          {Ene:5.3f} J")
print(f"# ホップ回転数:      {(abs(omgy / 2 / math.pi)):5.1f} rps")
print(f"# 射出仰俯角:      {elevAngle:+7.2f} °")
print(f"# ホップ傾斜角:    {tiltAngle:+7.2f} °")
# print(f"# 初速位置:          {xV0meas:5.3f} m".)
print(f"# マト1距離:         {xTarget1:5.3f} m")
print(f"# マト2距離:         {xTarget2:5.3f} m")
print(f"# 計算精度:       {tol:.2e} ")
print(f"# 空気抵抗係数:   {CdMethod} の式による")
print(f"# 空気抵抗補正:      {kCd:5.3f} ")
print()
print("計算回数    時刻  水平距離    玉高さ  玉速度 HOP回転 エネルギ 時間刻み", end='')
print("  LR位置  LR速度  LR回転")
print(
    "         t[msec]      x[m]    Δz[mm] vx[m/s] ωy[rps]    E[J]   Δt[ms]", end='')
print("    y[mm] vy[m/s] ωz[rps]")  # 　y　横方向


# 表示サブルーチン
def flightData(i, t, h, xf):
    # 飛翔中のデータ表示
    dz = (xf[2] - z0) * 1000  # 着弾高さ Δz[mm]
    ymm = xf[1] * 1000  # 横へのブレ y[mm]
    En = energy(mBb, xf[3:6])
    HopY = xf[7] / 2 / math.pi  # 通常ホップ軸
    HopZ = xf[8] / 2 / math.pi  # 傾きがある時
    print(f"{i:6d} ", end='')
    print(f"{(t * 1000):9.3f}  {xf[0]:8.4f} {dz:+9.2f}  ", end='')
    print(f"{xf[3]:6.2f}  {HopY:6.1f}  {En:6.3f}  {(h * 1000):7.4f}", end='')
    print(f" {ymm:+8.2f}  {xf[4]:6.3f}  {HopZ:6.1f}")  # 　y　横方向
    return


def impactData(tmpText):
    # 着弾、落下時の角度
    ke = xArray[3] / xArray[5]
    si = np.rad2deg(math.atan(1 / ke))
    # 重力落下量
    # tg = - 1 / 2 * g * t ** 2 * 1000  # [mm]空気抵抗の考慮がないため値が大きくなりすぎる
    # 左右
    if xArray[1] > 0.000005:  # [m]
        lr = '右'
    elif xArray[1] < -0.000005:
        lr = '左'
    else:
        lr = ''
    print(f"---------^^^^^^^ {tmpText}  角度{si:6.1f}° = 1/{abs(ke):5.1f} (z/x) -----------------------------------")
    # print(f"     重力落下量{tg:+8.2f} mm    ホップアップ量{(x[2] - z0) * 1000 - tg):+8.2f} mm  空気抵抗を考慮していないのでボツ")
    print(f"                 左右の着弾ズレ  {(xArray[1] * 1000):+8.2f} mm {lr}")
    return

#####


tt = 0
hh = 0.001
step = 100  # 表示周期
ii = 0
if v0Correct == 1:  # 初速補正する場合
    xVzero = np.copy(xArray)  # 計算データのコピー
    flightData(ii, tt, hh, xVzero)  # t=0のデータ表示
    for ii in range(1, 999):
        tt, hh, xVzero, ans = rkf45(tt, hh, xVzero, xV0meas)
        if ans < 0:
            print("error stop")
            while 1:
                pass

        # 初速測定位置
        if ans == 2:
            flightData(ii, tt, hh, xVzero)
            dV0 = v0 - magVector(xVzero[3:6], 3)  # ベクトルの大きさで計算
            v00 = v0 + dV0
            xArray[3], xArray[4], xArray[5] = muzzleVelocity(v00, vy)  # 初速値を補正変更
            print(f"-----------^^^^^---- 初速測定位置 -----^^^^^--  {dV0:+6.3f}m/s 修正し再計算 -----------------------")
            break
        # 時間切れ
        if tt >= te:
            print("タイムオーバー")
            break
        # データを表示
        if ii % step == 0 or ans == 1:
            # 表示周期毎と終点精度調整中に表示
            flightData(ii, tt, hh, xVzero)

# 弾道計算メイン
time = np.array([])
gx = np.array([])
gy = np.array([])
gz = np.array([])

xTarget = np.array([xTarget1, xTarget2, xTarget3])
targetNum = 0

tt = 0
hh = 0.001
ii = 0
flightData(ii, tt, hh, xArray)  # tt=0のデータ表示
for ii in range(1, 999999):
    ans = 0
    tt, hh, xArray, ans = rkf45(tt, hh, xArray, xTarget[targetNum])
    if ans < 0:
        print("error stop")
        while 1:
            pass

    time = np.append(time, tt)
    gx = np.append(gx, xArray[0])
    gy = np.append(gy, xArray[1])
    gz = np.append(gz, xArray[2])

    # 着弾した時
    if ans == 2:
        flightData(ii, tt, hh, xArray)
        text = "マトへ着弾"
        impactData(text)
        targetNum += 1
        if xTarget[targetNum] >= 999:
            break
        hh = 0.001
        continue

    # 着地した時
    if xArray[2] <= 0:   # z = x[2]
        flightData(ii, tt, hh, xArray)
        text = "地面に落下"
        impactData(text)
        break

    # 時間切れ
    if tt >= te:
        print("タイムオーバー")
        break

    # データを表示
    if ii % step == 0 or ans == 1:
        # 表示周期毎と終点精度調整中に表示
        flightData(ii, tt, hh, xArray)


# グラフの表示
# フォント設定
plt.rcParams['font.family'] = 'Hiragino Sans'
# 太さ設定
#plt.rcParams['font.weight'] = 'bold'

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(
    8, 6), tight_layout=True, sharex="all")

ax1.plot(gx, gz)
ax1.set_title("BB弾の弾道")
#ax.set_xlim(0, 50)
#ax.set_ylim(0, 1.4)
ax1.grid()
#ax1.set_xlabel('距離  [m]')
ax1.set_ylabel('高さ  [m]')

gy = gy * 1000
ax2.plot(gx, gy)
ax2.sharex(ax1)
#ax.set_xlim(0, 50)
#ax.set_ylim(0, 1.0)
ax2.invert_yaxis()
ax2.grid()
ax2.set_xlabel('距離  [m]')
ax2.set_ylabel('左右  [mm]')
plt.show()

print("stop")
while 1:
    pass
