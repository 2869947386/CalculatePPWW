#calculate
import numpy as np
from lhereader import LHEReader
from skhep.math import LorentzVector
import math
from scipy.integrate import tplquad,dblquad,quad
import matplotlib.pyplot as plt
#定义球谐函数pn±
def p1_plus(theta, phi):
    return np.sqrt(2) * np.sin(theta) * (5 * np.cos(theta) + 1) * np.cos(phi)

def p2_plus(theta, phi):
    return np.sqrt(2) * np.sin(theta) * (5 * np.cos(theta) + 1) * np.sin(phi)

def p3_plus(theta):
    return 0.25 * (5 + 4 * np.cos(theta) + 15 * np.cos(2 * theta))

def p4_plus(theta, phi):
    return 5 * np.sin(theta)**2 * np.cos(2 * phi)

def p5_plus(theta, phi):
    return 5 * np.sin(theta)**2 * np.sin(2 * phi)

def p6_plus(theta, phi):
    return np.sqrt(2) * np.sin(theta) * (-5 * np.cos(theta) + 1) * np.cos(phi)

def p7_plus(theta, phi):
    return np.sqrt(2) * np.sin(theta) * (-5 * np.cos(theta) + 1) * np.sin(phi)

def p8_plus(theta):
    return (1 / np.sqrt(48)) * (-5 + 12 * np.cos(theta) - 15 * np.cos(2 * theta))

def p1_minus(theta, phi):
    return np.sqrt(2) * np.sin(theta) * (5 * np.cos(theta) - 1) * np.cos(phi)

def p2_minus(theta, phi):
    return np.sqrt(2) * np.sin(theta) * (5 * np.cos(theta) - 1) * np.sin(phi)

def p3_minus(theta):
    return 0.25 * (5 - 4 * np.cos(theta) + 15 * np.cos(2 * theta))

def p4_minus(theta, phi):
    return 5 * np.sin(theta)**2 * np.cos(2 * phi)

def p5_minus(theta, phi):
    return 5 * np.sin(theta)**2 * np.sin(2 * phi)

def p6_minus(theta, phi):
    return np.sqrt(2) * np.sin(theta) * (-5 * np.cos(theta) - 1) * np.cos(phi)

def p7_minus(theta, phi):
    return np.sqrt(2)*np.sin(theta)*(-5*np.cos(theta)-1)*np.sin(phi)

def p8_minus(theta):
    return (1 / np.sqrt(48)) * (-5 - 12*np.cos(theta) - 15*np.cos(2 * theta))
#记录散射角cosΘ
cstheta=[]

I3tion=[]
#坐标系旋转和洛伦兹变换
def getxyz(pront1,k1,v,patl):
   #significant的定义用于筛选散射角
   significant1=1
   #质子速度（方向单位矢量）（质心系中质子速度均为（0,0,1））
   pront=pront1
   #求玻色子速度大小
   k2=np.power(k1,2)
   kv=math.sqrt(np.sum(k2))
   #玻色子速度方向单位矢量
   k=k1/kv
   #计算质子w玻色子散射角
   cos_Theta = np.dot(pront,k)  # 点积
   #为了方便折叠了一下
   if cos_Theta <0:
      cos_Theta=-1*cos_Theta
      pront=-1*pront
    #根据散射角筛选事件，默认为0即不筛选
   if cos_Theta<0:
      significant1 =0 
   Theta = np.arccos(cos_Theta)  #反余弦函数求夹角
   sin_Theta = np.sin(Theta)     #求正弦值
   cstheta.append(cos_Theta)     #记录散射角
   #计算新的坐标轴
   r_hat = (1/sin_Theta)*(pront-cos_Theta*k)
   n_hat = (1/sin_Theta) * np.cross(pront, k)
   #整理一下,旋转后的新坐标轴为
   kz=k
   nx=n_hat
   ry=r_hat
   #投影,得到w玻色子在新坐标轴中的速度
   a1=np.dot(v,nx)
   b1=np.dot(v,ry)
   c1=np.dot(v,kz)
   #沿z轴进行洛伦兹变换
   a2=a1*math.sqrt(1-kv*kv)/(1-c1*kv)
   b2=b1*math.sqrt(1-kv*kv)/(1-c1*kv)
   c2=(c1-kv)/(1-c1*kv)
   #得到轻子在w玻色子静止系中速度
   v_new=np.array([a2,b2,c2])
   #求二范数得到新速度的大小
   r1 = np.linalg.norm(v_new)
   #得到cosθ
   cst=np.dot(v_new,kz)/r1
   #计算横动量pta
   pbe=patl.e
   pbx=pbe*a1
   pby=pbe*b1
   pta=np.linalg.norm([pbx,pby])
   return v_new,cst,pta,significant1
#计算C2
def calculate_c2(f, g, h):
    term1 = -2/9
    g2=np.power(g,2)
    f2=np.power(f,2)
    h2=np.power(h,2)
    term2 = np.sum(g2)
    term3 =np.sum(f2)
    term4 = 4*np.sum(h2)
    return max(term1 -12*term2 + 6*term3+term4,0,term1-12*term3+6*term2+ term4)
#计算I3
def calculate_I3(h):
    return 4*(h[3,3] + h[4,4]) - (4*np.sqrt(3)/3) * (h[5,0] + h[5,5] + h[6,1] + h[6,6] + h[0,0] + h[0,5] + h[1,1] + h[1,6])
#读取lhe文件
reader = LHEReader('path to file/decay.lhe')
#W+散射出轻子的θ，φ
asind=0
afid=0

theta_p=[]
phi_p=[]
#W-散射出轻子的θ，φ
asind1=0
afid1=0

theta_m=[]
phi_m=[]

phi_all=[]

costp=[]
costm=[]
#横向动量
ptp=[]
ptm=[]
#不变质量
vmw=[]
nw=np.zeros((6,6))
i3w=np.zeros((6,6))
Htion=np.zeros((8,8))
H=np.zeros((8,8))
F=np.zeros((1,8))
G=np.zeros((1,8))

H11=np.zeros((8,8))
H12=np.zeros((8,8))
H13=np.zeros((8,8))
H14=np.zeros((8,8))
H15=np.zeros((8,8))
H16=np.zeros((8,8))
H17=np.zeros((8,8))
H18=np.zeros((8,8))
n=0
cros_sec=1.742#pb
def im(pati):#lhe文件里其实有不变质量，不过经验证这里的计算是正确的
   e2=(pati.e)**2
   px2=(pati.px)**2
   py2=(pati.py)**2
   pz2=(pati.pz)**2
   m2=e2-px2-py2-pz2
   imass=math.sqrt(m2)
   return imass
#用于计算一定不变质量和cosθ的I3
def calculaateI32(qi,cstheta,i3w,nw):
    if 0.7<cstheta[-1]<0.75:
          nw[0,qi]=nw[0,qi]+1
          for i in range(0,8):
           for j in range(0,8):
            pa1=p1(i+1,asind,afid)
            pa2=p2(j+1,asind1,afid1)
            Htion[i,j]=pa1*pa2
          i3w[0,qi]=i3w[0,qi]+calculate_I3(Htion)
    if 0.75<cstheta[-1]<0.8:
          nw[1,qi]=nw[1,qi]+1
          for i in range(0,8):
           for j in range(0,8):
            pa1=p1(i+1,asind,afid)
            pa2=p2(j+1,asind1,afid1)
            Htion[i,j]=pa1*pa2
          i3w[1,qi]=i3w[1,qi]+calculate_I3(Htion)
    if 0.8<cstheta[-1]<0.85:
          nw[2,qi]=nw[2,qi]+1
          for i in range(0,8):
           for j in range(0,8):
            pa1=p1(i+1,asind,afid)
            pa2=p2(j+1,asind1,afid1)
            Htion[i,j]=pa1*pa2
          i3w[2,qi]=i3w[2,qi]+calculate_I3(Htion)
    if 0.85<cstheta[-1]<0.9:
          nw[3,qi]=nw[3,qi]+1
          for i in range(0,8):
           for j in range(0,8):
            pa1=p1(i+1,asind,afid)
            pa2=p2(j+1,asind1,afid1)
            Htion[i,j]=pa1*pa2
          i3w[3,qi]=i3w[3,qi]+calculate_I3(Htion)
    if 0.9<cstheta[-1]<0.95:
          nw[4,qi]=nw[4,qi]+1
          for i in range(0,8):
           for j in range(0,8):
            pa1=p1(i+1,asind,afid)
            pa2=p2(j+1,asind1,afid1)
            Htion[i,j]=pa1*pa2
          i3w[4,qi]=i3w[4,qi]+calculate_I3(Htion)
    if 0.95<cstheta[-1]<1:
          nw[5,qi]=nw[5,qi]+1
          for i in range(0,8):
           for j in range(0,8):
            pa1=p1(i+1,asind,afid)
            pa2=p2(j+1,asind1,afid1)
            Htion[i,j]=pa1*pa2
          i3w[5,qi]=i3w[5,qi]+calculate_I3(Htion)
    return i3w,nw
#定义调用p+,p-函数
def p1(i,theta, phi):
    if i==1:
        return p1_plus(theta, phi)
    if i==2:
        return p2_plus(theta, phi)
    if i==3:
        return p3_plus(theta)
    if i==4:
        return p4_plus(theta, phi)
    if i==5:
        return p5_plus(theta, phi)
    if i==6:
        return p6_plus(theta, phi)
    if i==7:
        return p7_plus(theta, phi)
    if i==8:
        return p8_plus(theta)
def p2(i,theta, phi):
    if i==1:
      return p1_minus(theta, phi)
    if i==2:
      return p2_minus(theta, phi)
    if i==3:
      return p3_minus(theta)
    if i==4:
      return p4_minus(theta, phi)
    if i==5:
      return p5_minus(theta, phi)
    if i==6:
      return p6_minus(theta, phi)
    if i==7:
      return p7_minus(theta, phi)
    if i==8:
      return p8_minus(theta)
for iev, event in enumerate(reader):
    vmw1=0
    pnu=0
    inu=0
    ipp=0
    #根据pdgid寻找粒子，根据cosθ±定义选择
    aparticle = filter(lambda x: abs(x.pdgid) in(11,13,15,17), event.particles)
    #aparticle = filter(lambda x: x.pdgid in(11,12,13,14,15,16), event.particles)
    #pp为轻子，ww为w玻色子
    for pp in map(lambda x: x, aparticle):
       #排除只产生一个W玻色子的情况，与pnu有关的都是
       if pp.parent==1 or pp.parent==0:
         pnu=1
         continue
       if pnu==1:
          continue
       #定义inu，为了确保轻子与玻色子对应，这里其实并不必要，因为W+散射的一定为l+
       inu=0
       w=filter(lambda x: abs(x.pdgid)==24, event.particles)
       for ww in map(lambda x:x,w):
        #第一个w玻色子，检查轻子的parent粒子是否在第三行（pp.parent在读取时会减1）
          if pp.parent==2 and inu==0:
            #w玻色子四动量
            PA=ww.p4()
            #轻子四动量
            PB=pp.p4()
            #得到w玻色子和轻子速度 p/E
            vw=PA.boostvector
            vw1=np.array([vw.x,vw.y,vw.z])
            vl=PB.boostvector
            vl1=np.array([vl.x,vl.y,vl.z])
            #计算W玻色子不变质量
            vmw.append(im(PA))
            vmw1=vmw1+im(PA)
            #定义质子速度
            vp1=np.array([0,0,1])
            #变换到w玻色子静止系
            vl2,costheta,pt1,significantn1=getxyz(vp1,vw1,vl1,PB)
            #得到θ,φ
            atheta=np.arccos(costheta)
            aphi = np.arctan2(vl2[1],vl2[0])
            #判断由w+ 还是w-产生并记录 asind，asind1分别为记录W+W-产生轻子的cosθ
            if ww.pdgid==24:
               asind=atheta
               afid=aphi
               ptp.append(pt1)
            else :
               asind1=atheta
               afid1=aphi
               ptm.append(pt1)
        ##第二个w玻色子，检查轻子的parent粒子是否在第四行
          if pp.parent==3 and inu==1:
            #粒子四动量
            PA=ww.p4()
            PB=pp.p4()
            #vw是w玻色子速度
            vw=PA.boostvector
            vw1=np.array([vw.x,vw.y,vw.z])
            #vl轻子速度
            vl=PB.boostvector
            vl1=np.array([vl.x,vl.y,vl.z])
            #直接给出质子速度（0，0，1）
            vp1=np.array([0,0,1])
            vl2,costheta,pt1,significantn2=getxyz(vp1,vw1,vl1,PB)
            atheta = np.arccos(costheta)
            aphi = np.arctan2(vl2[1],vl2[0])
            #不变质量
            vmw.append(im(PA))
            vmw1=vmw1+im(PA)
            if ww.pdgid==24:
               asind=atheta
               afid=aphi
               ptp.append(pt1)
            else:
               asind1=atheta
               afid1=aphi
               ptm.append(pt1)
          inu=inu+1
    #计算一定不变质量和cosθ内的I3
    if 130<vmw1<140:
       i3w,nw=calculaateI32(0,cstheta,i3w,nw)
    if 140<vmw1<150:
       i3w,nw=calculaateI32(1,cstheta,i3w,nw)
    if 150<vmw1<160:
       i3w,nw=calculaateI32(2,cstheta,i3w,nw)
    if 160<vmw1<170:
       i3w,nw=calculaateI32(3,cstheta,i3w,nw)
    if 170<vmw1<180:
       i3w,nw=calculaateI32(4,cstheta,i3w,nw)
    if 180<vmw1<190:
       i3w,nw=calculaateI32(5,cstheta,i3w,nw)                     
    #避免只产生一个w玻色子的情况
    if pnu==1:
       continue
    #计算系数矩阵h,g,f
    #sig1和sig2均不为0，表明散射角均达到要求，其实此处判断一次即可，碰撞质子延z轴，对于同一事件pxpy动量守恒，散射角相等
    if significantn1*significantn2 !=0:
       #记录event数目
       n=n+1
       #计算系数矩阵
       for i in range(0,8):
        for j in range(0,8):
          pa1=p1(i+1,asind,afid)
          pa2=p2(j+1,asind1,afid1)
          H[i,j]=H[i,j]+pa1*pa2
          #H[i,j]=pa1*pa2
          Htion[i,j]=pa1*pa2
          if j==0:
            G[0,i]=G[0,i]+pa1
            #G[0,i]=pa1
          if i==0:
            F[0,j]=F[0,j]+pa2
            #F[0,j]=pa2
    #判断单个事件的C2正负
    #C21=2*calculate_c2(F,G,H)
    #if C21<0:
    #  print('c2-')
    #检验程序运行一次
    #break
    if significantn1*significantn2 !=0:
     costp.append(np.cos(asind))
     costm.append(np.cos(asind1))
     theta_p.append(asind)
     phi_p.append(afid)
     theta_m.append(asind1)
     phi_m.append(afid1)
     phi_all.append(afid)
     phi_all.append(afid1)
     I3tion.append(calculate_I3(Htion))
print(n)
F1=F/(n)
G1=G/(n)
H1=H/(n)
C2=2*calculate_c2(F1,G1,H1)
print(C2)
I3=calculate_I3(H1)
print(I3)
#画图
#散射角余弦值 1
data1=cstheta
plt.hist(data1,bins=25)
plt.title("cos(Θ)distribution")
plt.xlabel("cos(Θ)")
plt.ylabel("rate")
plt.show()
#I3分布图 2
data2=I3tion 
plt.hist(data2,bins=50)
plt.title("I3 distribution")
plt.xlabel("I3")
plt.ylabel("rate")
plt.show()
#cos + 3
data3=costp
plt.hist(data3,bins=100)
plt.title("cos(θ+) distribution")
plt.xlabel("cos(θ)")
plt.ylabel("rate")
plt.show()
#cos+ 分布2
data31=costp
bins1 = np.linspace(min(data31), max(data31),100)  # n个桶
hist,bin_edges = np.histogram(data31, bins=bins1) #分桶
cros_sec=1.742#pb
costpnu = hist
plt.step(bin_edges[:-1],costpnu,where='post')
plt.ylabel('rate')
plt.xlabel("cos(θ)")
plt.title('cos(θ+) distribution')
plt.show()
#cos - 4
data4=costm
plt.hist(data4,bins=100)
plt.title("cos(θ-) distribution")
plt.xlabel("cos(θ)")
plt.ylabel("rate")
plt.show()
#cos- 分布2
data41=costm
bins1 = np.linspace(min(data41), max(data41),100)  # 100个桶
hist,bin_edges = np.histogram(data41, bins=bins1) #分桶
costmnu = hist
plt.step(bin_edges[:-1],costmnu,where='post')
plt.ylabel('rate')
plt.xlabel("cos(θ)")
plt.title('cos(θ-) distribution')
plt.show()
#phi 正 5
data5=phi_p
plt.hist(data5,bins=50)
plt.title("φ+ distribution")
plt.xlabel("φ")
plt.ylabel("rate")          
plt.show()
#phi 负 6
data6=phi_m
plt.hist(data6,bins=50)
#plt.title("φ- distribution")
plt.title("φ- distribution")
plt.xlabel("φ")
plt.ylabel("rate")
plt.show()
#phi 总 7
data7=phi_all
plt.hist(data7,bins=50)
plt.title("φ distribution")
plt.xlabel("φ")
plt.ylabel("rate")
plt.show()
#ptp 8
data8=ptp
plt.hist(data8,bins=100)
plt.title("pt+ distribution")
plt.xlabel("pt")
plt.ylabel("rate")
plt.show()
#ptm 9
data9=ptm
plt.hist(data9,bins=100)
plt.title("pt- distribution")
plt.xlabel("pt")
plt.ylabel("rate")
plt.show()
#不变质量 10
data10=vmw
plt.hist(data10,bins=100)
plt.title("i mass  distribution")
plt.xlabel("Mw")
plt.ylabel("rate")
plt.show()
#print(H1)
#X=H1
ax = plt.matshow(H1,cmap=plt.cm.Blues)
plt.colorbar(ax.colorbar, fraction=0.025)
plt.title("matrix H")
plt.show()
#
bx = plt.matshow(G1,cmap=plt.cm.Blues)
plt.colorbar(bx.colorbar, fraction=0.025)
plt.title("matrix G")
plt.show()
#
cx = plt.matshow(F1,cmap=plt.cm.Blues)
plt.colorbar(cx.colorbar, fraction=0.025)
plt.title("matrix F")
plt.show()
#微分散射截面
data11=costp
bins1 = np.linspace(min(data11), max(data11),50)  # 50个桶
hist,bin_edges = np.histogram(data11, bins=bins1) #分桶
cros_sec=1.742#pb
diff_cross_section = hist*cros_sec/(sum(hist)*np.diff(bin_edges)*cros_sec)
plt.step(bin_edges[:-1],diff_cross_section,where='post')
plt.ylabel('dσ/dcos(θ+)*σ')
plt.xlabel("cos(θ+)")
plt.title('dσ/dcos(θ+)*σ distribution')
plt.show()
#微分散射截面2
data12=costm
bins1 = np.linspace(min(data12), max(data12),50)  # 50个桶
hist,bin_edges = np.histogram(data12, bins=bins1) #分桶
cros_sec=1.742#pb
diff_cross_section = hist*cros_sec/(sum(hist)*np.diff(bin_edges)*cros_sec)
plt.step(bin_edges[:-1],diff_cross_section,where='post')
plt.ylabel('dσ/dcos(θ-)*σ')
plt.xlabel("cos(θ-)")
plt.title('dσ/dcos(θ-)*σ distribution')
plt.show()
