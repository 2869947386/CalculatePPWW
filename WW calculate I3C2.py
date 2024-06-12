import numpy as np
from lhereader import LHEReader
from skhep.math import LorentzVector
import math
from scipy.integrate import tplquad,dblquad,quad
import matplotlib.pyplot as plt
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
#def boost
#k is for w
cstheta=[]
I3tion=[]
cros_sec=1.742#pb
#坐标系旋转和洛伦兹变换
def getxyz(pront1,k1,v,patl):
   significant1=1
   #质子方向（默认001）
   pront=pront1
   #求玻色子速度大小
   k2=np.power(k1,2)
   kv=math.sqrt(np.sum(k2))
   #玻色子速度方向单位矢量
   k=k1/kv
   #计算散射角
   cos_Theta = np.dot(pront,k)  # 点积
   if cos_Theta <0:
      cos_Theta=-1*cos_Theta
      pront=-1*pront
   if cos_Theta<0:
      significant1 =0
      
   Theta = np.arccos(cos_Theta)  # 反余弦函数求夹角
   sin_Theta = np.sin(Theta)     #sin
   cstheta.append(cos_Theta)
   #新的坐标轴
   r_hat = (1/sin_Theta)*(pront-cos_Theta*k)
   n_hat = (1/sin_Theta) * np.cross(pront, k)
   kz=k
   nx=n_hat
   ry=r_hat
   #对新坐标轴进行投影，得到w玻色子方向以z轴的坐标系速度
   a1=np.dot(v,nx)
   b1=np.dot(v,ry)
   c1=np.dot(v,kz)
   #沿z轴进行洛伦兹变换
   a2=a1*math.sqrt(1-kv*kv)/(1-c1*kv)
   b2=b1*math.sqrt(1-kv*kv)/(1-c1*kv)
   c2=(c1-kv)/(1-c1*kv)
   #得到轻子在静止系中速度
   v_new=np.array([a2,b2,c2])
   #得到新速度大小
   r1 = np.linalg.norm(v_new)
   #得到与z轴夹角的cos值
   cst=np.dot(v_new,kz)/r1
   #计算横动量，因为延z轴进行洛伦兹变换
   pbe=patl.e
   pbx=pbe*a1
   pby=pbe*b1
   pta=np.linalg.norm([pbx,pby])
   return v_new,cst,pta,significant1
def calculate_c2(f, g, h):
    term1 = -2/9
    g2=np.power(g,2)
    f2=np.power(f,2)
    h2=np.power(h,2)
    term2 = np.sum(g2)
    term3 =np.sum(f2)
    term4 = 4*np.sum(h2)
    return max(term1 -12*term2 + 6*term3+term4,0,term1-12*term3+6*term2+ term4)
def calculate_I3(h):
    return 4*(h[3,3] + h[4,4]) - (4*np.sqrt(3)/3) * (h[5,0] + h[5,5] + h[6,1] + h[6,6] + h[0,0] + h[0,5] + h[1,1] + h[1,6])
reader = LHEReader('C:/Users/Doorwaz/Downloads/decay2.lhe')
asind=0
afid=0
asind1=0
afid1=0
theta_p=[]
phi_p=[]
theta_m=[]
phi_m=[]
phi_all=[]
costp=[]
costm=[]
dcostp=[]
ptp=[]
ptm=[]
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
H11=np.zeros((8,8))
H11=np.zeros((8,8))
n=0
cros_sec=1.742#pb
def im(pati):#后来发现lhe文件里有不变质量，经验证这里的计算是正确的
   e2=(pati.e)**2
   px2=(pati.px)**2
   py2=(pati.py)**2
   pz2=(pati.pz)**2
   m2=e2-px2-py2-pz2
   imass=math.sqrt(m2)
   return imass
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
    aparticle = filter(lambda x: abs(x.pdgid) in(11,13,15,17), event.particles)
    #aparticle = filter(lambda x: x.pdgid in(11,12,13,14,15,16), event.particles)
    #pp为轻子，ww为w玻色子
    for pp in map(lambda x: x, aparticle):
       #这里是排除一个玻色子的情况，与pnu有关的都是
       if pp.parent==1 or pp.parent==0:
         pnu=1
         continue
       if pnu==1:
          continue
       #定义inu，为了确保轻子与玻色子对应
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
            #计算θ,φ
            atheta=np.arccos(costheta)
            aphi = np.arctan2(vl2[1],vl2[0])
            #选择w+ w-,记录角度
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
            vl=PB.boostvector
            vl1=np.array([vl.x,vl.y,vl.z])
            vp1=np.array([0,0,1])
            vl2,costheta,pt1,significantn2=getxyz(vp1,vw1,vl1,PB)
            point = vl2
            atheta = np.arccos(costheta)
            aphi = np.arctan2(point[1], point[0])
            vmw.append(im(PA))
            vmw1=vmw1+im(PA)
            if ww.pdgid==24:
               asind=atheta
               afid=aphi
               ptp.append(pt1)
            else :
               asind1=atheta
               afid1=aphi
               ptm.append(pt1)
          inu=inu+1
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
    #记录事件数
    #计算h,g,f
    if significantn1*significantn2 !=0:
       n=n+1
       for i in range(0,8):
        for j in range(0,8):
          pa1=p1(i+1,asind,afid)
          pa2=p2(j+1,asind1,afid1)
          H[i,j]=H[i,j]+pa1*pa2
          Htion[i,j]=pa1*pa2
          if j==0:
            G[0,i]=G[0,i]+pa1
            #G[0,i]=pa1
          if i==0:
            F[0,j]=F[0,j]+pa2
            #F[0,j]=pa2
    C21=2*calculate_c2(F,G,H)
    if C21<0:
       print('c2-')
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
plt.title("φ distribution")
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