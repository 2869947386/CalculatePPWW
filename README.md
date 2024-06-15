# Justcalculat
本程序通过pp>ww>leptons过程最终散射产物的角分布，计算可观测量I3＆C2，相应数学方法可参考：https://arxiv.org/abs/2302.00683

##模拟实验数据使用蒙特卡洛产生子Madgraph5 https://launchpad.net/mg5amcnlo 产生

![image](https://github.com/2869947386/CalculatePPWW/blob/main/image/cd.png)

generate p p > w+ w- > l+ vl l- vl~ 并output decay3文件夹

![image](https://github.com/2869947386/CalculatePPWW/blob/main/image/generate.png)

![image](https://github.com/2869947386/CalculatePPWW/blob/main/image/output.png)

检查runcard

launch decay3

![image](https://github.com/2869947386/CalculatePPWW/blob/main/image/launch.png)

得到结果与lhe文件

![image](https://github.com/2869947386/CalculatePPWW/blob/main/image/result.png)

##python程序主要由几部分组成：

1.读取lhe文件并储存W玻色子，轻子信息
2.将质心系转化为W玻色子静止系并计算相应角分布
3.计算系数矩阵与可观测量并画图

所用到的角度定义如图所示

![image](https://github.com/2869947386/CalculatePPWW/blob/main/image/图片2.png)
