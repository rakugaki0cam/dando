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


from logging import raiseExceptions
import numpy as np
import math
import matplotlib.pyplot as plt


#パラメータ設定例
# 1, tol=1e-10, eps=1e-9   => とても精密　計算は遅い                     
# 2, tol=1e-6, eps=1e-6    => 十分に実用的な計算結果
#
#ルンゲクッタの計算精度　(tol は固定値)
# 1e-2  < tol           : 粗い　計算は速い
# 1e-10 < tol < 1e-2    : 普通
# 1e-12 < tol < 1e-10   : 精密　計算は遅い
#         tol < 1e-12   : 非推奨
tol = 1e-10       #桁精度
#
#Conversence epsilon at rootfinding method
#     1 < eps           : 粗い　計算は速い
# 1e-10 < eps < 1       : 普通
#         eps < 1e-10   : 精密　計算は遅い
#                          or cannot converge less than 1e-10.
eps = 1e-7



def rho_humid(T, P, H):
   # ρ : 空気密度 [kg/m^3]
   # T : 気温 [°C]
   # P : 気圧 [Pa] (1気圧 = 101325 Pa)
   # H : 湿度 [％RH]

   et = 6.1078 - 10 ** (7.5 * T / (T + 237.3))              # et[hPa]
   et = 100 * et
   et = H * et
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


def rk_preparation(method):
   #ルンゲクッタ法の計算テーブル準備
   #set Butcher table
   # method:  "rk4" or "rkf45"

   global a, b, c, Rc, s

   if(method == "rk4"):
      #4次4段
      s = 4    #次数
      a = np.zeros((s, s))
      a[0] = (0,0,0,0)
      a[1] = (0.5,0,0,0)
      a[2] = (0,0.5,0,0)
      a[3] = (0,0,1,0)

      b = (1/6, 1/3, 1/3, 1/6)
      c = (0, 0.5, 0.5, 1)
      Rc = 0

   elif(method == "rkf45"):
      #5次6段
      s = 6    #次数
      a = np.zeros((s, s))
      a[0] = (0, 0, 0, 0, 0, 0)
      a[1] = (0.25, 0, 0, 0, 0, 0)
      a[2] = (0.09375, 0.28125, 0, 0, 0, 0)
      a[3] = (1932/2197, -7200/2197, 7296/2197, 0, 0, 0)
      a[4] = (439/216, -8, 3680/513, -845/4104, 0, 0)
      a[5] = (-8/27, 2, -3544/2565, 1859/4104, -11/40, 0)
      b = np.zeros((2, s))
      b[0] = (25/216, 0, 1408/2565, 2197/4104, -0.2, 0)  
      b[1] = (16/135, 0, 6656/12825, 28561/56430, -9/50, 2/55)
      c = (0, 0.25, 3/8, 12/13, 1, 0.5)
      Rc = (1/360, 0, -128/4275, -2197/75240, 1/50, 2/55)
    
   else:
      #メソッドエラー
      raise Exception("program stop at rk_preparation")
      
   
   return a, b, c, Rc

##TEST ---------------------------------
#a, b, c, Rc = rk_preparation('rk4')
#print(a, b, c, Rc)


def rkf451(func, N, x, h, y, xbound, info, tol, ik):
   #------------
   #info = -2  (Failed. calclation range has already over.)
   #     = -1  (Failed.
   #             h becomes too small. change tol or check condition of func.)
   #     =  0  (Success. running now)
   #     =  1  (Success. x reach xbound normally)
   #------------

      #implicit none
       #interface
      #   def func(iN,ix,iy,is,iik)
      #     implicit none
      #     integer,intent(in)::iN,is,iik
      #      double precision,intent(in)::ix,iy(1:iN)
      #      double precision::func
      #    end def func
      # end interface
    #integer,intent(in)::N,ik
    #double precision,intent(in)::xbound
    #double precision,intent(in)::tol
    #double precision,intent(inout)::x,h,y(1:N)
    #integer,intent(inout)::info
    #double precision,parameter::
   hmin = 1e-14
   hmax = 1
    #integer::i,j,FLAG
    #double precision::R,delta,tx,tmp(1:N),K(1:s,1:N),Sy,err
   K = np.zeros((s,N))

   if(abs(h) >= hmax):
      if(h <= 0):
         h = -hmax
      else:
         h = hmax
  
   FLAG=1
   if(abs(x - xbound) <= 1e-15):
      info = 1
      FLAG = 0
   else:
      if(abs(h) <= 1e-15):
         print("maybe overflow or underflow, please change tol.")
         print("====Err info====")
         print("x    --> ", x)
         print("h    --> ", h)

         for i in range(N):
            print("y(",i,") --> ",y(i))
            
         print("================")
         info = -1
         FLAG = 0
         raise Exception("stop")

      if((h <= 0) and (xbound - x >= 0)):
         info = -2
         FLAG = 0
      elif((h > 0) and (xbound - x <= 0)):
         info = -2
         FLAG = 0


   while(FLAG == 1):
      tx = x
      for j in range(s):
         tx = x + c[j] * h #############
         tmp = y
         for i in range(j - 1):
            #tmp(1:N) = tmp(1:N) + K(i, 1:N) * a(j, i) ##############################
            temp
         for i in range(N):
            K[j, i] = h * func(N, tx, tmp, i, ik) ###############################FUNC()
            

      #step 4
      R = 0
      for i in range(N):
         R += (Rc[1] * K[1][i] + Rc[3] * K[3][i] + Rc[4] * K[4][i] + Rc[5] * K[5][i] + Rc [6] * K[6][i]) ** 2
      
      R = abs(dsqrt(R) / h)

      Sy = 0
      for i in range(N):
         Sy = Sy + (y[i] * y[i])


      Sy = dsqrt(Sy)
      if(Sy >= 1):
         err = tol * Sy
      else:
         err = tol
      
      #step 5
      if(R <= err):
         x = x + h  
         for i in range(s):
            y(1:N) = y(1:N) + b[0][i] * K(i, 1:N) ##########################
         FLAG = 0
      
      #step 6
      #  Avoid zero deviding.
      if(R >= 1e-20):
         delta = (err/(2 * R)) ** 0.25
      else:
         delta = 4

      #step 7
      if(delta <= 0.1):
         #def changes dramatically.
         h = 0.1 * h
      elif(delta >= 4):
         #def changes loosely.
         h = 4 * h
      else:
         #def changes moderately.
         h = delta * h

      #step 8
      if(abs(h) <= hmax):
         if(h <= 0):
            h = -hmax
         else:
            h = hmax

      #step 9
      if(abs(xbound-x) <= abs(h)):
         h = xbound - x
         if(abs(h) <= hmin):
            info = 1
            FLAG = 0
      
   return #################################

'''
def rkf451_e(func,x,y,xbound,info,tol,ik):
   #sikinote
   #  propagate from y(x) to y(xbound) without interval
   #
   #  info = -1 : h < hmin. Maybe path the singular point.
   #       =  1 : x reach xbound.
   #
   #implicit none
    #interface
   #    def func(iN,ix,iy,is,iik)
   #      implicit none
   #      integer,intent(in)::iN,is,iik
   #      double precision,intent(in)::ix,iy(1:iN)
   #      double precision::func
   #    end def func
   # end interface
   # integer::N,ik
   # double precision,intent(in)::xbound,tol
   # double precision,intent(inout)::x,y(:)
   # integer,intent(inout)::info
   # double precision,parameter::hmin=1d-14,hmax=1.

   # integer::i,j,FLAG,key,disc
   # double precision::R,delta,tx,Sy,err,h,h0
   # double precision,allocatable::tmp(:),K(:,:)

   disc=0
   key=0
   h0=999
   N=size(y,1)
   #allocate(tmp(1:N),K(1:s,1:N))
   tmp=0; K=0  
   
   h=xbound-x   
   if(abs(h) > hmax):
      if(h < 0):
         h = -hmax
      else:
         h = hmax


   FLAG=1
   if(abs(x-xbound) < hmin*0.1):
      info=1
      FLAG=0

   while(FLAG == 1):
      tx=x
      for j in range(s):
         tx=x+c(j)*h
         tmp(1:N)=y(1:N)
         for i in range(j-1):
            tmp(1:N)=tmp(1:N)+K(i,1:N)*a(j,i)
         
         for i in range(N):
            K(j,i)=h*func(N,tx,tmp,i,ik)
         
      

      #step 4
      R=0
      for i in range(N):
         R=R+(Rc(1)*K(1,i)+Rc(3)*K(3,i)+Rc(4)*K(4,i)+Rc(5)*K(5,i)+Rc(6)*K(6,i))**2
      
      R=abs(dsqrt(R)/h)

      Sy=0
      for i in range(N):
         Sy=Sy+(y(i)*y(i))
      
      Sy=dsqrt(Sy)
      if(Sy.ge.1):
         err=tol*Sy
      else
         err=tol
      

      #step 5
      if(R.le.err.or.key == 1):
         x=x+h
         do i=1,s
            y(1:N)=y(1:N)+b1(i)*K(i,1:N)
         
         key=0
      

      #step 6
      #  Avoid zero deviding.
      if(R.ge.1d-20):
         delta=(err/(2*R))**0.25
      else
         delta=4
      

      #step 7
      if(delta.le.0.1):
         #def changes dramatically.
         h=0.1*h
      elseif(delta.ge.4):
         #def changes loosely.
         h=4*h
      else
         #def changes moderately.
         h=delta*h
      

      #step 8
      if(abs(h).ge.hmax):
         h=sign(1,h)*hmax
      elseif(abs(h).lt.hmin):
         h=sign(1,h)*hmin
         key=1
         disc=1
      
       
      #step 9
      if(abs(xbound-x).le.abs(h)):
         h=xbound-x
         if(abs(h).le.hmin):
            info=1
            FLAG=0
         
      
   
   
   if(disc == 1):
      info=-1
   
       
   #deallocate(tmp,K)
   return





#--------------------
def seeklm(xa,xb,za,zb,nodes,poles,r0,v0,omg0,u0,ik)
  #nodes --> number of nodes
  #poles --> number of local maximam/minimam
  #
  #nodes=0
  #              
  #              
  #    x--------------------------------
  #      \_______  
  #               \______
  #                      \
  #                      orbit of bullet
  #
  #nodes=1
  #   za_____________xa
  #               __/\_                
  #      ________/     \_  
  #     /                \  
  #    x------------------\------------
  #                        \
  #                         \
  #                          \
  #                             \orbit of bullet
  #
  #nodes=2
  #
  #  zb________________ xb                
  #                    /\
  #          xa      /    \  
  #    x-----------/--------\------------
  #      \        /          \
  #        \    /             \
  #  za______\/                \
  #                             \orbit of bullet
  #
  use RKmod
  use GBL
  implicit none
  integer,intent(in)::ik
  double precision,intent(in)::omg0(1:3),r0(1:3),v0(1:3),u0(1:3)
  double precision,intent(out)::xa,xb,za,zb
  integer,intent(out)::nodes,poles

  integer::i,j,info
  integer,parameter::N=12
  double precision::vz0,vz1,vz2,t0,t1,t2,ta,tb,tbef,rz0,rz1
  double precision::x(1:N),tt,th,tx(1:N),bt,bh,bx(1:N)
  double precision::t,h,xbound,rkfd
  external::rkfd
 
  x(1:3)=r0(1:3)
  x(4:6)=v0(1:3)
  x(7:9)=omg0(1:3)
  x(10:12)=u0(1:3)
 
  xbound=10
  i=1; t=0; h=1d-1; info=0
  info=0; nodes=0; poles=0
 
  xa=0; xb=0; za=0; zb=0
 
  call rk_preparation("rkf45")
  do while(info == 0)
     tt=t; th=h; tx=x

     rz0=x(3); vz0=x(6); t0=t  
     call rkf451(rkfd,size(x),t,h,x,xbound,info,tol,ik)
     rz1=x(3); vz1=x(6); t1=t

     if(vz0*vz1.lt.0.):
        bt=tt; bh=th; bx=tx

        t0=bt; t2=0.5*(2*bt+th); t1=bt+th
        do j=1,60
           tt=bt; th=bh; tx=bx

           th=t2-bt
           call rkf451(rkfd,size(tx),tt,th,tx,xbound,info,tol,ik)
           vz2=tx(6)
           if(vz2*vz0.ge.0.):
              t2=t2+0.5*(t1-t2)
           else
              tbef=t2
              t2=t2-0.5*(t1-t2)
              t1=tbef
           

           if(abs(t1-t2)/t1.lt.1d-10):
              if(poles == 0):
                 ta=t2; xa=tx(1); za=tx(3)
              else
                 tb=t2; xb=tx(1); zb=tx(3)
              
              exit
           
        
        poles=poles+1
     
     if((rz0-r0(3))*(rz1-r0(3)).lt.0):
        nodes=nodes+1
     
     
     if(x(3).lt.0):
        info=1
     
     i=i+1
  
  call rk_deallocation("rkf45")
 
  return
end def seeklm

def seekzeroin(xa,xb,nodes,r0,v0,omg0,u0,ik)
  #nodes --> number of nodes
  #
  #nodes=0
  #              
  #              
  #    x--------------------------------
  #      \_______  
  #               \______
  #                      \
  #                      orbit of bullet
  #
  #nodes=1
  #
  #               __/\_                
  #      ________/     \_  
  #     /                \ xa
  #    x------------------\------------
  #                        \
  #                         \
  #                          \
  #                             \orbit of bullet
  #
  #nodes=2
  #
  #
  #                    /\
  #              xa  /    \  xb
  #    x-----------/--------\------------
  #      \        /          \
  #        \    /             \
  #          \/                \
  #                             \orbit of bullet
  #
  use RKmod
  use GBL
  implicit none
  integer,intent(in)::ik
  double precision,intent(in)::omg0(1:3),r0(1:3),v0(1:3),u0(1:3)
  double precision,intent(out)::xa,xb
  integer,intent(out)::nodes
 
  integer::i,j,info
  integer,parameter::N=12
  double precision::rz0,rz1,rz2,t0,t1,t2,ta,tb,tbef,rzini
  double precision::x(1:N),tt,th,tx(1:N),bt,bh,bx(1:N)
  double precision::t,h,xbound,rkfd
  external::rkfd

  rzini=r0(3)
  x(1:3)=r0(1:3)
  x(4:6)=v0(1:3)
  x(7:9)=omg0(1:3)
  x(10:12)=u0(1:3)

  xbound=10
  i=1; t=0; h=1d-1; info=0
  info=0; nodes=0
 
  xa=0; xb=0
 
  call rk_preparation("rkf45")

  #-----1step-----
  call rkf451(rkfd,size(x),t,h,x,xbound,info,tol,ik)
  #----------------
  do while(info == 0)
     tt=t; th=h; tx=x

     rz0=x(3); t0=t  
     call rkf451(rkfd,size(x),t,h,x,xbound,info,tol,ik)
     rz1=x(3); t1=t
     
     if(rz0-rzini.lt.0..neqv.rz1-rzini.lt.0.):
        bt=tt; bh=th; bx=tx

        t0=bt; t2=0.5*(2*bt+th); t1=bt+th
        do j=1,60
           tt=bt; th=bh; tx=bx

           th=t2-bt
           call rkf451(rkfd,size(tx),tt,th,tx,xbound,info,tol,ik)
           rz2=tx(3)
           if(rz2-rzini.lt.0..eqv.rz0-rzini.lt.0.):
              t2=t2+0.5*(t1-t2)
           else
              tbef=t2
              t2=t2-0.5*(t1-t2)
              t1=tbef
           

           if(abs(t1-t2)/t1.lt.1d-10):
              if(nodes == 0):
                 ta=t2; xa=tx(1)
              else
                 tb=t2; xb=tx(1)
              
              nodes=nodes+1
              exit
           
        
     

     if(x(3).lt.0):
        info=1
     
     i=i+1
  
  call rk_deallocation("rkf45")
 
  return
end def seekzeroin

def detwy(energy,r,vy,omgx,omgz,u,zeroinx,exwy,exvz,ik)
  implicit none
  integer,intent(in)::ik
  double precision,intent(in)::energy,zeroinx
  double precision,intent(in)::r(1:3),omgx,omgz,u(1:3),vy
  double precision,intent(out)::exwy,exvz 

  integer::Nr,count,Nsplit
  double precision,parameter::pi=dacos(-1)
  double precision,allocatable::rd(:)
  double precision::ad,bd,wycond
  external::wycond
 
  ad=-800*2*pi; bd=2*pi; Nsplit=25; Nr=1
  allocate(rd(1:Nr))
  rd(1:Nr)=0
  count=0
  call rootfindingw(wycond,energy,r,vy,omgx,omgz,u,zeroinx,ad,bd,Nsplit,Nr,rd,exvz,count,"bisection",1)
  exwy=rd(1)

  if(ik == 0):  
     ad=exwy-5*2*pi; bd=exwy+5*2*pi; Nsplit=1; Nr=1
     rd(1:Nr)=0
     count=0
     call rootfindingw(wycond,energy,r,vy,omgx,omgz,u,zeroinx,ad,bd,Nsplit,Nr,rd,exvz,count,"false_position",0)
     exwy=rd(1)
  

  deallocate(rd)
  return
end def detwy
#--------------------------
def wycond(energy,r,vy,omgx,omgy,omgz,u,zeroinx,exvz,ik)
  use GBL, only:m
  implicit none
  integer,intent(in)::ik
  double precision,intent(in)::energy,r(1:3),u(1:3),vy,omgx,omgy,omgz,zeroinx
  double precision,intent(out)::exvz

  integer::nodes
  double precision::zeroinx1,zeroinx2,wycond,v(1:3),omg(1:3)
 
  omg(1)=omgx; omg(2)=omgy; omg(3)=omgz

  call detvz(energy,r,vy,omg,u,exvz,ik)
  v(1)=dsqrt((2*energy/m)-exvz**2-vy**2)
  v(2)=vy
  v(3)=exvz
  call seekzeroin(zeroinx1,zeroinx2,nodes,r,v,omg,u,ik)
  if(zeroinx2.ge.1d-10):
     wycond=zeroinx2-zeroinx
  else
     wycond=2*zeroinx
  
     
  return
end def wycond
#------------------------------
def detwy2(energy,r,vy,omgx,omgz,u,udh,exwy,exvz,ik)  
  implicit none
  integer,intent(in)::ik
  double precision,intent(in)::energy,udh,r(1:3),omgx,omgz,u(1:3),vy
  double precision,intent(out)::exwy,exvz

  integer::Nr,count,Nsplit
  double precision,parameter::pi=dacos(-1)
  double precision,allocatable::rd(:)
  double precision::ad,bd,wycond2
  external::wycond2
 
  ad=-1000*2*pi; bd=2*pi; Nsplit=25; Nr=1
  allocate(rd(1:Nr))
  rd(1:Nr)=0
  count=0 
  call rootfindingw(wycond2,energy,r,vy,omgx,omgz,u,udh,ad,bd,Nsplit,Nr,rd,exvz,count,"bisection",1)
  exwy=rd(1)

  if(ik == 0):
     ad=exwy-5*2*pi; bd=exwy+5*2*pi; Nsplit=1; Nr=1
     count=0; rd(1:Nr)=0
     call rootfindingw(wycond2,energy,r,vy,omgx,omgz,u,udh,ad,bd,Nsplit,Nr,rd,exvz,count,"false_position",0)
     exwy=rd(1)
  

  deallocate(rd) 
  return
end def detwy2
#--------------------------
def wycond2(energy,r,vy,omgx,omgy,omgz,u,udh,exvz,ik)
  use GBL, only:m
  implicit none
  integer,intent(in)::ik
  double precision,intent(in)::energy,r(1:3),u(1:3),vy,omgx,omgy,omgz,udh
  double precision,intent(out)::exvz

  integer::nodes,poles
  double precision::xa,xb,za,zb,z0,wycond2,v(1:3),omg(1:3)
 
  omg(1)=omgx; omg(2)=omgy; omg(3)=omgz

  call detvz(energy,r,vy,omg,u,exvz,ik)
  v(1)=dsqrt((2*energy/m)-exvz**2-vy**2)
  v(2)=vy
  v(3)=exvz
 
  z0=r(3)
  call seeklm(xa,xb,za,zb,nodes,poles,r,v,omg,u,ik)
  if(poles == 0.or.za.lt.0):
     wycond2=10
  elseif(poles == 1):
     wycond2=za-z0
  elseif(poles == 2):
     wycond2=zb-za
  else
     print(6,*)"unknown,poles-->",poles
     stop
  
  wycond2=wycond2-udh
 
  return
end def wycond2
#-------------------------
def detvz(energy,r,vy,omg,u,exvz,ik)  
  implicit none
  integer,intent(in)::ik
  double precision,intent(in)::energy
  double precision,intent(in)::r(1:3),omg(1:3),u(1:3),vy
  integer::Nr,count,Nsplit
  double precision,intent(out)::exvz
 
  double precision,allocatable::rd(:)
  double precision::ad,bd,vzcond
  external::vzcond
 
  ad=-30; bd=1; Nsplit=80; Nr=1
  allocate(rd(1:Nr))
  rd(1:Nr)=0
  count=0
  call vzroot(vzcond,energy,r,vy,omg,u,ad,bd,Nsplit,Nr,rd,count,ik)
  exvz=rd(1)
 
  return
end def detvz
#------------------------
def vzcond(energy,r,vy,vz,omg,u,ik)
  use GBL, only:m
  implicit none
  integer,intent(in)::ik
  double precision,intent(in)::energy,r(1:3),vy,vz,omg(1:3),u(1:3)
  double precision::vzcond,v(1:3)

  integer::nodes,poles
  double precision::xa,xb,za,zb,z0
 
  v(1)=dsqrt((2*energy/m)-vz**2-vy**2)
  v(2)=vy
  v(3)=vz

  z0=r(3)
  call seeklm(xa,xb,za,zb,nodes,poles,r,v,omg,u,ik)
  if(poles == 0):
     vzcond=-z0-1
  elseif(poles == 1):
     vzcond=za-z0
  elseif(poles == 2):
     vzcond=za+zb-2*z0
  else
     print(6,*)"unknown,poles-->",poles
     stop
  
 
  return
end def vzcond
#--------------------
def vzroot(func,energy,r,vy,omg,u,ad,bd,Nsplit,Nr,rd,count,ik)
  #Developer : sikino
  use GBL
  implicit none
  interface
     def func(ienergy,ir,ivy,ivz,iomg,iu,iik)
       implicit none
       integer,intent(in)::iik
       double precision,intent(in)::ienergy,ir(1:3),ivy,ivz,iomg(1:3),iu(1:3)
       double precision::func
     end def func
  end interface
  integer,intent(in)::Nr,Nsplit,ik
  double precision,intent(in)::ad,bd,energy,r(1:3),vy,omg(1:3),u(1:3)
  integer,intent(out)::count
  double precision,intent(out)::rd(1:Nr)
 
  integer::N,i,j 
  double precision::h,hd,x0,x1,x2,y0,y1,y2,Dy0,tx1,ty1

  #ad < bd
  if(ad.ge.bd):
     print(6,'(A,2e15.6e3)')"must be ad < bd, your ad,bd --> ",ad,bd
     print(6,'(A)')"program stop at rootfinding"
     stop
  end if
 
  N=Nsplit
  hd=abs(ad-bd)/dble(N)
  count=0

  x0=ad
  y0=func(energy,r,vy,x0,omg,u,ik)
  for i in range(N):
     x1=x0+hd
     y1=func(energy,r,vy,x1,omg,u,ik)
     if(y0.lt.0..neqv.y1.lt.0.):
        tx1=x1
        ty1=y1
        
        #bisection rule
        do j=1,60
           x2=0.5*(x0+x1)
           y2=func(energy,r,vy,x2,omg,u,ik)
           if(y2.lt.0..eqv.y0.lt.0.):
              x0=x2
           else
              x1=x2
           
           if(abs((x1-x0)/x2).lt.eps):
              if(y2.le.1):
                 count=count+1
                 rd(count)=x2
                 exit
              
           
        
        if(j == 61)print(6,*)"+--cannot convergence at root-finding method bs--+"
        
        if(count.ge.Nr)exit
        x1=tx1
        y1=ty1
     
     x0=x1
     y0=y1
  

  return
end def vzroot

def rootfindingw(func,energy,r,vy, &
     omgx,omgz,u,zeroinx,ad,bd,Nsplit,Nr,rd,exvz,count,method,ik)
  #Date : 2015/07/28
  #Developer : sikino
  use GBL
  implicit none
  interface
     def func(ienergy,ir,ivy,iomgx,iomgy,iomgz,iu,izin,iexvz,iik)
       implicit none
       integer,intent(in)::iik
       double precision,intent(in)::ienergy,ir(1:3),ivy,iomgx,iomgy,iomgz,izin,iu(1:3)
       double precision,intent(out)::iexvz
       double precision::func
     end def func
  end interface
  integer,intent(in)::Nr,Nsplit,ik
  double precision,intent(in)::ad,bd,energy,r(1:3),vy,u(1:3)
  double precision,intent(in)::omgx,omgz,zeroinx
  integer,intent(out)::count
  double precision,intent(out)::rd(1:Nr),exvz
  double precision,parameter::pi=dacos(-1)
  character(*),intent(in)::method
 
  integer::N,i,j
               
  double precision::h,hd,x0,x1,x2,y0,y1,y2,Dy0,tx1,ty1

  #ad < bd
  if(ad.ge.bd):
     print(6,'(A,2e15.6e3)')"must be ad < bd, your ad,bd --> ",ad,bd
     print(6,'(A)')"program stop at rootfinding"
     stop
  end if

  #Announsment
  print(6,'(A,A)')" ---- ",trim(method)
 
  N=Nsplit
  hd=abs(ad-bd)/dble(N)
  count=0

  x0=ad
  y0=func(energy,r,vy,omgx,x0,omgz,u,zeroinx,exvz,ik)
  for i in range(N):
     x1=x0+hd
     y1=func(energy,r,vy,omgx,x1,omgz,u,zeroinx,exvz,ik)
     print(6,'(A,i0,A,i0,A,f14.8,A,f14.8)')" Sequence ",i," of ",N," : omgy/2pi ",x1/2/pi, " : zeroin ",y1
     if(y0.lt.0..neqv.y1.lt.0.):
        print(6,'(A)')"  +--- Root find, go to conversion phase --+  "
        tx1=x1
        ty1=y1
        if(trim(method) == "bisection"):
           #bisection rule
           do j=1,60              
              x2=0.5*(x0+x1)
              y2=func(energy,r,vy,omgx,x2,omgz,u,zeroinx,exvz,ik)
              if(y2.lt.0..eqv.y0.lt.0.):
                 x0=x2
              else
                 x1=x2
              
              print(6,'(A,f16.8,A)')"Conversion => ",eps*abs(x2/(x1-x0))*100,"%"
              if(abs((x1-x0)/x2).lt.eps):
                 #If y2 is large, y2 is singular point.  
                 if(y2.le.1):
                    count=count+1
                    rd(count)=x2
                 
                 exit
              
           
           if(j == 60)print(6,*)"+--cannot convergence at root-finding method bs--+"
        elseif(trim(method) == "newton_raphson"):
           #Newton-Raphson method
           do j=1,10
              h=1.d-4
              y0=func(energy,r,vy,omgx,x0,omgz,u,zeroinx,exvz,ik)
              Dy0=(func(energy,r,vy,omgx,x0+h,omgz,u,zeroinx,exvz,ik)-y0)/h
              x2=x0-y0/Dy0
              if(abs((x0-x2)/x2).lt.eps):
                 if(y0.le.1):
                    count=count+1
                    rd(count)=x2
                 
                 exit
              
              x0=x2
           
           if(j == 10)print(6,*)"+--cannot convergence at root-finding method nr--+"
        elseif(trim(method) == "false_position"):
           #false position method
           if(abs(y1).gt.abs(y0)):
              do j=1,30
                 x2=x0-y0*(x1-x0)/(y1-y0)
                 print(6,'(A,f16.8,A)')"Conversion => ",eps*abs(x2/(x2-x0))*100,"%"
                 if(abs((x2-x0)/x2).lt.eps):
                    if(y0.le.1):
                       count=count+1
                       rd(count)=x2
                    end if
                    exit
                 
                 x0=x2
                 y0=func(energy,r,vy,omgx,x0,omgz,u,zeroinx,exvz,ik)
              
           else
              do j=1,30
                 x2=x0-y0*(x1-x0)/(y1-y0)
                 print(6,'(A,f15.8,A)')"Conversion => ",eps*abs(x2/(x1-x2))*100,"%"
                 if(abs((x2-x1)/x2).lt.eps):
                    if(y1.le.1):
                       count=count+1
                       rd(count)=x2
                    
                    exit
                 
                 x1=x2
                 y1=func(energy,r,vy,omgx,x1,omgz,u,zeroinx,exvz,ik)
              
           
           if(j == 30)print(6,*)"+--cannot convergence at root-finding method fp--+"
        else
           print(6,'(A,A)')"unknown type of method, your method--> ",trim(method)
           print(6,'(A,A)')"program stop"
           stop
        end if
        if(count.ge.Nr)exit
        x1=tx1
        y1=ty1
     
     x0=x1
     y0=y1
  

  return
end def rootfindingw
#---------------------------
def rkfd(N,t,x,s,ik)
  use GBL
  implicit none
  integer,intent(in)::N,s,ik
  double precision,intent(in)::t,x(1:N)
  double precision::r(1:3),v(1:3),relv(1:3),omg(1:3),u(1:3),nw,nv,I
  double precision::rkfd,gravity,vis1,vis2,mag,Nz,Nze
  double precision,parameter::pi=datan(1)*4
  external::gravity,vis1,vis2,mag,Nz,Nze
 
  r(1:3)=x(1:3)
  v(1:3)=x(4:6)
  omg(1:3)=x(7:9)
  u(1:3)=x(10:12)
  relv(1:3)=v(1:3)-u(1:3)
 
  rkfd=0.
  if(s.le.3):
     # Differential equation of position
     # d {x,y,z}/dt = v{x,y,z}
     rkfd=v(s)
  elseif(s.le.6):
     # Differential equation of velocity
     # d v_{x,y,z}/dt = F{x,y,z}
     rkfd=gravity(s-3,m,g,r)+vis1(s-3,relv,a,eta) &
          +vis2(s-3,relv,a,eta,rho)+mag(s-3,omg,relv,a,rho)
     rkfd=rkfd/m
  elseif(s.le.9):
     # Differential equation for omega
     # d omega{x,y,z}/dt = N_{x,y,z}/I
     if(trim(omgdecay) == "yes"):
        I=0.4*m*a*a
        nv=dsqrt(relv(1)**2+relv(2)**2+relv(3)**2)
        nw=dsqrt(omg(1)**2+omg(2)**2+omg(3)**2)
        if(nw.le.1d-13):
           rkfd=0
        else
           if(ik == 0):
              rkfd=(Nz(nv,a,eta,rho,nw)/I)*x(s)/nw
           else
              rkfd=(Nze(nv,a,eta,rho,nw)/I)*x(s)/nw
           end if
        
     else
        rkfd=0
     
  else
     # Differential equation for wind
     # d u{x,y,z}/dt = Fu_{x,y,z}
     rkfd=0
  
 
  return
end def rkfd
#-----------------------------
def gravity(dir,m,g,r)
  implicit none
  integer,intent(in)::dir
  double precision::gravity
  double precision,intent(in)::m,g,r(1:3)
 
  if(dir == 3):
     gravity=-m*g
  else
     gravity=0.
  
 
  return
end def gravity
#-------------------------------
def c1(a,eta)
  implicit none
  double precision::c1
  double precision,intent(in)::eta,a
  double precision,parameter::pi=dacos(-1.)

  c1=6.*pi*eta*a
 
  return
end def c1
#------------------------------
def c2(nv,a,eta,rho)
  implicit none
  double precision,parameter::pi=dacos(-1.)
  double precision,intent(in)::nv,eta,a,rho
  double precision::Cd,Reynolds,c2
  external::Cd,Reynolds

  c2=0.5*Cd(Reynolds(nv,a,eta,rho))*rho*pi*a*a
 
  return
end def c2
#-----------------------------
def vis1(dir,v,a,eta)
  implicit none
  integer,intent(in)::dir
  double precision,intent(in)::v(1:3),a,eta
  double precision::vis1,norm,nv,c1
  external::c1
 
  nv=dsqrt(v(1)**2.+v(2)**2.+v(3)**2.)
  norm=c1(a,eta)*nv
 
  vis1=-norm*v(dir)/nv
 
  return
end def vis1
#-------------------------
def vis2(dir,v,a,eta,rho)
  implicit none
  integer,intent(in)::dir
  double precision,intent(in)::v(1:3),a,eta,rho
  double precision::vis2,norm,nv,c2
  external::c2
 
  nv=dsqrt(v(1)*v(1)+v(2)*v(2)+v(3)*v(3))
  norm=c2(nv,a,eta,rho)*nv*nv
 
  vis2=-norm*v(dir)/nv
 
  return
end def vis2
#--------------------------
def mag(dir,omg,v,a,rho)
  implicit none
  integer,intent(in)::dir
  double precision,intent(in)::omg(1:3),v(1:3),a,rho
  double precision,parameter::pi=dacos(-1.)
  double precision::mag,L(1:3),nomg,nv,nL,Cl
 
  nomg=dsqrt(omg(1)*omg(1)+omg(2)*omg(2)+omg(3)*omg(3))
  if(nomg.le.1.d-14):
     mag=0.
     return
  
 
  nv=dsqrt(v(1)*v(1)+v(2)*v(2)+v(3)*v(3))
  L(1)=v(2)*omg(3)-v(3)*omg(2)
  L(2)=v(3)*omg(1)-v(1)*omg(3)
  L(3)=v(1)*omg(2)-v(2)*omg(1)
  nL=dsqrt(L(1)*L(1)+L(2)*L(2)+L(3)*L(3))

  if(nL.le.1.d-14):
     mag=0.
     return
  
 
  # for BB bullet
  Cl=0.12
  mag=-Cl*(4./3.)*pi*(a**3)*2.*rho*nomg*nv*L(dir)/nL
 
  return
end def mag
#--------------------------
def Reynolds(nv,a,eta,rho)
  #Reynolds number
  # nv : norm of velocity of object
  #  a : radius of object
  #eta : viscosity (not Kinetic viscosity)
  #rho : density of fluid
  implicit none
  double precision,intent(in)::nv,a,eta,rho
  double precision::keta,Reynolds
 
  #keta means Kinetic viscosity
  keta=eta/rho
 
  Reynolds=nv*2.*a/keta
 
  return
end def Reynolds
#---------------------------
def Cd(Re)
  #From
  #http://www.chem.mtu.edu/~fmorriso/DataCorrelationForSphereDrag2013.pdf
  #Drag coefficient Cd,
  #Cd depend on Reynolds number,Re.
  #Fource of Drug,D is written by
  #     1
  # D= ---Cd*rho*pi*a**2*|V|**2
  #     2
  #        ^ This Cd#
  implicit none
  double precision,intent(in)::Re
  double precision::Cd,c1,c2,c3,c4
 
  c1=24./Re
  c2=Re/5.
  c2=2.6*c2/(1.+c2**(1.52))
  c3=Re/263000
  c3=0.411*c3**(-7.94)/(1.+c3**(-8))
  c4=Re**0.8/461000
 
  Cd=c1+c2+c3+c4
 
  return
end def Cd
#----------------------------
def Nz(nv,a,eta,rho,omg)
  #Moment of omg direction
  implicit none
  double precision,intent(in)::nv,a,eta,rho,omg
  double precision::Nz,Cf,Fintegral
  external::Cf,Fintegral
 
  Nz=0.5*rho*Cf(nv,a,eta,rho)*a**3
  Nz=Nz*Fintegral(nv,a,omg)
 
  return
end def Nz
def Nze(nv,a,eta,rho,omg)
  implicit none
  double precision,intent(in)::nv,a,eta,rho,omg
  double precision::Nze,Cf,Fintegral
  double precision,parameter::pi=atan(1)*4
  external::Cf,Fintegral
  double precision::pc,tc,vu,vd,t

  pc=pi/5.32065 # magic phi
  tc=pi/3.60475 # magic theta
  
  vu= nv*sin(pc)-a*omg*sin(tc)
  vd=-nv*sin(pc)-a*omg*sin(tc)
  Nze=-0.5*rho*Cf(nv,a,eta,rho)
  Nze=Nze*(4*pi*a**2)*a*0.5
  Nze=-Nze*(abs(vu)*vu+abs(vd)*vd)
  
  return
end def Nze
#-------------------
def Cf(nv,a,eta,rho)
  implicit none
  # nv : norm of velocity of object
  #  a : radius of object
  #eta : viscosity (not Kinetic viscosity)
  #rho : density of fluid
  double precision,intent(in)::nv,a,eta,rho
  double precision::Reynolds,Cf
  external::Reynolds

  # Laminar flow
  Cf=1.328/dsqrt(Reynolds(nv,a,eta,rho))
 
  # Turbulent flow Re > 10^7
  #Cf=0.455/(log10(Reynolds(nv,a,eta,rho))**(2.58))
 
  return
end def Cf
#----------------------------------
def Fintegral(u,R,omega)
  # 2016/07/23
  #   2pi     pi
  #  /       /
  #  | dphi  | |u*sin(phi)-R*omega*sin(theta)|*{u*sin(phi)-R*omega*sin(theta)}*sin^2(theta) d theta
  #  /       /
  #  0      0
  implicit none
  double precision,intent(in)::u,R,omega
  double precision,parameter::pi=datan(1.)*4
  double precision::Fintegral,s,Fphi,x(1:15),w(1:15)
  external::Fphi
  integer::i

  s=0
  call GaussKronrod15ab(0,pi,x,w)
  do i=1,15
     s=s+w(i)*Fphi(u,R,omega,x(i))
  
  call GaussKronrod15ab(pi,2*pi,x,w)
  do i=1,15
     s=s+w(i)*Fphi(u,R,omega,x(i))
      
  Fintegral=s  
 
  return
end def Fintegral
#------------------------------
def GaussKronrod15ab(a,b,x,w)
  #Gauss-Kronrod Quadrature Nodes and Weights
  #http://www.advanpix.com/2011/11/07/gauss-kronrod-quadrature-nodes-weights/
  implicit none
  double precision,intent(in)::a,b
  double precision,intent(out)::x(1:15),w(1:15)
 
  integer::i
  integer,parameter::N=15
 
  x=0; w=0
 
  x( 8) = 0
  x( 9) = 2.077849550078984676006894037732449d-1
  x(10) = 4.058451513773971669066064120769615d-1
  x(11) = 5.860872354676911302941448382587296d-1
  x(12) = 7.415311855993944398638647732807884d-1
  x(13) = 8.648644233597690727897127886409262d-1
  x(14) = 9.491079123427585245261896840478513d-1
  x(15) = 9.914553711208126392068546975263285d-1  
  do i=1,7
     x(i)=-x(N-i+1)
  
 
  w( 8) = 2.094821410847278280129991748917143d-1
  w( 9) = 2.044329400752988924141619992346491d-1
  w(10) = 1.903505780647854099132564024210137d-1
  w(11) = 1.690047266392679028265834265985503d-1
  w(12) = 1.406532597155259187451895905102379d-1
  w(13) = 1.047900103222501838398763225415180d-1
  w(14) = 6.309209262997855329070066318920429d-2
  w(15) = 2.293532201052922496373200805896959d-2
  do i=1,7
     w(i)=w(N-i+1)
  
 
  x=0.5*((b-a)*x+(a+b))
  w=0.5*(b-a)*w
 
  return
end def GaussKronrod15ab
#----------------------------------------
def Fphi(u,R,omega,phi)
  # using analysis solution.
  #
  #   pi
  #  /
  #  | |u*sin(phi)-R*omega*sin(theta)|*{u*sin(phi)-R*omega*sin(theta)}*sin^2(theta) d theta
  #  /
  #  0
  #
  implicit none
  double precision,intent(in)::u,R,omega,phi
  double precision,parameter::pi=dacos(-1)
  double precision::a,b,c,tp,Fphi,Rw

  a=(u/(R*omega))*sin(phi)
 
  if(a.le.0):
     Fphi=-(0.5*a*a*pi-8*a/3+3*pi/8)
  elseif(a.ge.1):
     Fphi=0.5*a*a*pi-8*a/3+3*pi/8
  elseif(a.gt.0.and.a.lt.1):
     tp=asin(a)
     b=-a*a*pi*0.5-8*a/3-3*pi/8+(2*a*a+3/2)*tp
     c=-(a*a+1)*dsin(2*tp)+dsin(4*tp)/8+6*a*dcos(tp)-2*a*dcos(3*tp)/3
     Fphi=b+c
  else
     print(6,*)"undefined at Fphi, program stop",a
     stop
  
 
  Rw=R*omega
  Fphi=Fphi*Rw*abs(Rw)

  return
end def Fphi

#---------------------------------
def bullet_orbit(N,x0,Nt,time,filename,ik)
  #sikinote http://slpr.sakura.ne.jp/qp/
  #Author : sikino
  #Date   : 2016/03/18 (yyyy/mm/dd)
  #       : 2016/03/20
  use GBL
  use RKmod
  implicit none
  integer,intent(in)::ik,N,Nt
  double precision,intent(in)::x0(1:N),time(0:Nt)
  character(*),intent(in)::filename
  double precision,parameter::pi=datan(1)*4
  integer::info,i,j
  double precision::x(1:N),Ene,h,t,tbound,rkfd
  external::rkfd

  x=x0
  info=0

  t=time(0)
  open(22,file=trim(filename),status='replace')
  print(22,'(A,e12.5e1,A)')"# m:  ",m,"[kg]"
  print(22,'(A,f12.5,A)')"# g:  ",g,"[m s^{-2}]"
  print(22,'(A,e12.5e1,A)')"# a:  ",a,"[m]"
  print(22,'(A,f12.5,A)')"# temperature:  ",temperature,"[degree Celsius]"
  print(22,'(A,f12.5,A)')"# moisture:  ",moisture,"[no-dimension]"
  print(22,'(A,f12.5,A)')"# pressure:  ",pressure,"[Pa]"
  print(22,'(A,e12.5e1,A)')"# eta:",eta,"[kg m^{-1} s^{-1}]"
  print(22,'(A,f12.5,A)')"# rho:",rho,"[kg m^{-3}]"
  print(22,'(A,e12.5e1)')"# tol:",tol
  print(22,'(A,e12.5e1)')"# eps:",eps
  print(22,'(A,14A13)')"#","t[s]","x[m]","y[m]","z[m]","vx[m/s]","vy[m/s]","vz[m/s]" &
   ,"wx[rot/s]","wy[rot/s]","wz[rot/s]","windx[m/s]","windy[m/s]","windz[m/s]","Energy[J]"
  Ene=0.5*m*(x(4)**2+x(5)**2+x(6)**2)
  print(22,'(14f13.6)')t,(x(i),i=1,6),(x(i)/(2*pi),i=7,9),(x(i),i=10,12),Ene
 
  call rk_preparation("rkf45")
  do j=0,Nt-1
     t=time(j); tbound=time(j+1); h=tbound-t
     call rkf451_e(rkfd,t,x,tbound,info,tol,ik)
             
     Ene=0.5*m*(x(4)**2+x(5)**2+x(6)**2)
     print(22,'(14f13.6)')t,(x(i),i=1,6),(x(i)/(2*pi),i=7,9),(x(i),i=10,12),Ene
     #Stop calculation when bullet reach ground.
     if(x(3).lt.0)exit
  
  call rk_deallocation("rkf45")
 
  close(22)
  print(6,'(A,e10.3)')"Reach ground at ",t
 
  return
end def bullet_orbit

program main
  use GBL
  implicit none
  integer::i,N,Nt,ik
  double precision,allocatable::x(:),time(:),r(:),u(:)
  double precision,parameter::pi=datan(1)*4
  double precision::rx,ry,rz,vx,vy,vz,ux,uy,uz,omgx,omgy,omgz
  double precision::energy,zeroinx,updownz,stept,exvz,exwy,theta
  character(48)::outputfile,filename,search_zeroin,search_updown
  double precision::tbound
  double precision::rho_humid,eta_air
  double precision::r0(1:3),v0(1:3),omg0(1:3),u0(1:3)
  external::rho_humid,eta_air
  
  double precision::xa,xb,za,zb
  integer::nodes,poles
  real::t1,t0
  
  namelist /input/m,a,energy,g,temperature,pressure,moisture &
       ,omgdecay,search_zeroin,zeroinx,search_updown,updownz &
       ,theta,rx,ry,rz,vy,ux,uy,uz,omgx,omgy,omgz &
       ,stept,outputfile,ik
 
  open(10,file="./input")
  read(10,nml=input)
  close(10)
  print(*,nml=input)

  N=12
  eta=eta_air(temperature)
  rho=rho_humid(temperature,pressure,moisture)
 
  vx=dsqrt((2*energy/m)-vy**2)*dcos(pi*theta/180)
  vz=dsqrt((2*energy/m)-vy**2)*dsin(pi*theta/180)
  omgx=omgx*2*pi
  omgy=omgy*2*pi
  omgz=omgz*2*pi

  allocate(x(1:N))
  x(1)=rx;  x(4)=vx;  x(7)=omgx;  x(10)=ux
  x(2)=ry;  x(5)=vy;  x(8)=omgy;  x(11)=uy
  x(3)=rz;  x(6)=vz;  x(9)=omgz;  x(12)=uz
 
  tbound=100
  Nt=nint(tbound/stept)
  allocate(time(0:Nt)); time=0
  do i=0,Nt
     time(i)=dble(i)*stept
  

  filename=trim(outputfile)//".txt"
  call bullet_orbit(N,x,Nt,time,filename,ik)
  
  call cpu_time(t0)
  if(trim(search_zeroin) == "yes"):
     print(6,'(A)')"==============================="
     print(6,'(A,f10.5)')" search vz and wy by zeroin --> ", zeroinx
     allocate(r(1:3),u(1:3)); r(1:3)=x(1:3); u(1:3)=x(10:12)
     call detwy(energy,r,x(5),x(7),x(9),u,zeroinx,exwy,exvz,ik)
     print(6,*)exwy/(2*pi),exvz
     x(4)=dsqrt((2*energy/m)-exvz**2-vy**2)
     x(5)=vy
     x(6)=exvz
     x(8)=exwy
     filename=trim(outputfile)//"_opt.txt"
     call bullet_orbit(N,x,Nt,time,filename,ik)

     deallocate(r,u)
  

  if(trim(search_updown) == "yes"):
     print(6,'(A)')"==============================="
     print(6,'(A,f10.5)')" search vz and wy by up-down height --> ", updownz
     allocate(r(1:3),u(1:3)); r(1:3)=x(1:3); u(1:3)=x(10:12)
     call detwy2(energy,r,x(5),x(7),x(9),u,updownz,exwy,exvz,ik)
     print(6,*)exwy/(2*pi),exvz
     x(4)=dsqrt((2*energy/m)-exvz**2-vy**2)
     x(5)=vy
     x(6)=exvz
     x(8)=exwy
     filename=trim(outputfile)//"_h.txt"
     call bullet_orbit(N,x,Nt,time,filename,ik)
     deallocate(r,u)
  

  call cpu_time(t1)
  print(6,'(f10.3,A)')(t1-t0),"[CPU sec]"
 
  stop
end program main
'''