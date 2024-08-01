#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# PoLoBag algorithm for d=2 polynomial features d=2 多项式特征
import pandas as pd
df = pd.read_excel('net1.xlsx') #读取文件
df.loc[:,(df==0).all()].columns #提取全为零的列

import numpy as np
from scipy import stats
from sklearn import linear_model
import warnings
from sklearn.exceptions import ConvergenceWarning


#Input expression and output sorted edgelist file names
#输入表达量 输出排序边列表的文件名
infilename = 'intersectgeneEXP.txt'
outfilename = 'CC_result_network.txt'


#Algorithm parameters 算法参数
n12 = 0.5    # controls number of linear features in each bootstrap sample 控制每个自举样本中线性特征的数量
n22 = 3.5    # controls number of nonlinear features in each bootstrap sample 控制每个自举样本中非线性特征的数量
nM = 0.6     # controls bootstrap sample size 控制自举样本大小
nB = 600     # total number of bootstrap samples in the ensemble 集合中的自举样本总数
alpha = 0.1  # the Lasso regularization parameter (Lasso 正则化参数)
#n12控制在每个bootstrap样本中线性特征的数量
#n22控制在每个bootstrap样本中非线性特征的数量
#nM控制集合中bootstrap样本的数量
#nB表示集合中bootstrap样本的总数
#alpha是Lasso正则化参数

#For repeatability
#np.random.seed(30)

#Ignore Lasso convergence warnings
#忽视LASSO收敛警告
warnings.filterwarnings("ignore", category=ConvergenceWarning)

#Algorithm

#Read input expression file
#读取输入表达文件
D = {}#循环完成后D保存的是基因及其对应的z-score表达量
genes = []
with open(infilename, 'r') as f:
    geneCount = -1
    for line in f:
        geneCount += 1 #循环完成后geneCount等于基因的总数量
        if geneCount == 0:   #header
            continue
        lineData = line.rstrip().split("\t")   #tab separated input file
        genes.append(lineData[0])#lineData是文件中的每一行,即基因名及表达量,genes保存的是所有基因的名称
        D[lineData[0]] = stats.zscore(np.array(lineData[1:]).astype(np.float))#D保存的是表达量转换成zscore后的值
genes = np.unique(genes)#保存的是174个基因的名称

#潜在调控因子即所有的基因
regs = genes #Potential regulators 潜在调控因子列表,即174个基因

edges = []#edges存储相对应的调控关系,regulator gene-target gene
w = []#w存储权重值
for t in genes:  #Regression problem for each target gene 每个靶基因的回归问题
    yt = np.transpose(D[t][np.newaxis]) #靶基因的表达量,列向量,某个基因在615个样本中的表达量
    Xt = yt[:,[]] #靶基因以外的基因的表达量,即潜在调控因子的表达量,行为样本,列为基因,615*173
    for reg in regs:#regs存储的是所有潜在调控基因的名字
        if not reg == t:  # no auto-regulation
            Xt = np.hstack((Xt,np.transpose(D[reg][np.newaxis])))
    ntR = Xt.shape[1] #ntR等于潜在调控因子的数量,173
    if ntR == 0:
        continue
    valm = np.median(np.hstack((yt,Xt))) #valm是所有基因表达量的中位数
    tau = 3
    yt = yt - tau*valm #靶基因的表达量变化,得到靶基因在所有所有样本中的表达量,615*1
    Xt = Xt - tau*valm #潜在调控因子的表达量变化,得到n-1个基因在所有样本中的表达量,615*173
    nlin = int(n12*np.sqrt(ntR)) #nlin=0.5*潜在调控因子的平方根,nlin为线性特征的数量,0.5*根号下173=6
    nnlin = int(n22*np.sqrt(ntR)) #nnlin=3.5*潜在调控因子的平方根,nnlin为非线性特征的数量,3.5*根号下173=46
    wt = np.zeros((ntR,1))  #Weight权重数组初始化为0,每个靶基因的权重数组大小为(n-1)*1,173*1
    wtM = np.zeros((ntR,1))  #Weight magnitude表示权重大小的数组初始化为0,173*1
    wtS = np.zeros((ntR,1))  #Weight sign表示权重符号的数组初始化为0,173*1
    stw = np.zeros((ntR,1)) + 1e-10  # How many times a feature was chosen in a sample 某特征在一个样本中被选择了多少次,173*1
    m = Xt.shape[0] #m表示样本总数量,615
    for b in range(nB):   #For each sample range(nB)表示从0到nB-1的整数列,nB个bootstrap
        rowindex = np.random.choice(m,int(nM*m)) #从m个整数中随机选取nM*m个整数(可以重复),从615个样本中随机选取369个样本(可以重复),369*1
        ytb = yt[rowindex,:] #靶基因在相对应随机选取的样本中的表达量,369*1
        XtbF = Xt[rowindex,:] #173个潜在调控因子在相对应随机选取的样本中的表达量,369*173
        #Separate idtb for linear and nonlinear features对于线性和非线性特征 将idtb分开
        idtblin = [] #线性特征的idtb
        idtbnlin = [] #非线性的idtb
        if nlin > 0:#如果线性特征的数量>0
            idtblin = np.random.choice(ntR,nlin,replace=False) #在潜在调控因子数量ntR中随机选取nlin个整数(不可重复),在173个潜在调控因子中随机选取6个基因作为线性特征,6*1
        Xtb = XtbF[:,idtblin]   # linear features 线性特征 随机选取的线性特征在相应随机选取的样本中的表达量,369*6
        if nnlin > 0:#如果非线性特征的数量>0
            allindex = [x for x in range(ntR) if x not in idtblin] #allindex表示随机挑出的线性基因特征以外的那些基因,173-6=167
            index = np.random.choice(allindex,2*nnlin,replace=False) #随机选取非线性特征(不可重复),从167个潜在调控基因中随机选取2*46=92个基因,92*1
            Xtb21F = XtbF[:,index[:nnlin]] #index的前半部分基因在随机样本的表达量,nnlin个非线性特征在随机样本中的表达量,369*46
            Xtb22F = XtbF[:,index[nnlin:]] #index的后半部分基因在随机样本的表达量,nnlin个非线性特征在随机样本中的表达量(与上一部分的Xtb21F不重复),369*46
            vala = Xtb21F*Xtb22F #两个非线性特征的表达矩阵点乘,即两个矩阵对应元素相乘,369*46
            vals1 = np.sign(Xtb21F) #取对应数字前的正负号,369*46
            vals12 = np.sign(vala) #取对应数字前的正负号369*46
            atbk = 0.5*np.sqrt(np.abs(vala)) #np.abs计算数组元素的绝对值,np.sqrt计算数组各元素的平方根,等式6,369*46
            ftb21 = vals1*(1+vals12)*atbk #等式6,369*46
            ftb22 = vals1*(1-vals12)*atbk #atbk ftb21 ftb22通过文中等式6计算,369*46
            Xtb = np.hstack((Xtb,ftb21[:,:int(0.5*nnlin)],ftb22[:,int(0.5*nnlin):])) # nonlinear features 线性特征和线性特征合在一起,369*52
            for l in range(int(0.5*nnlin)):
                idtbnlin.append(str(index[l]) + '#' + str(index[l+nnlin]))
            for l in range(int(0.5*nnlin),nnlin):
                idtbnlin.append('-' + str(index[l]) + '#' + str(index[l+nnlin])) # an indicator for category 分类指标
            #若非线性特征有46个,1-23对应47-69;24-46对应70-92;共46组两两对应
        clf = linear_model.Lasso(alpha=alpha,fit_intercept=False)  # Lasso,alpha为正则化参数,fit_intercept如果为FALSE,对数据进行去中心化处理
        clf.fit(Xtb, ytb)#ytb是靶基因在随机挑选的样本中的表达量,Xtb是挑选出的特征在同样的样本中的表达量,Xtb 369*52,ytb 369*1
        wtb = clf.coef_ #wtb是选出的特征的LASSO系数,Xtb是LASSO的输入,ytb是LASSO的输出,52*1
        if len(idtblin) > 0:  # for linear features 对于线性特征,如果线性特征的数量大于0
            wtM[idtblin,0] += np.abs(wtb[0:len(idtblin)]) #wtb是特征的权重,此处取wtb的线性特征并取绝对值
            wtS[idtblin,0] += wtb[0:len(idtblin)] #线性特征权重大小
            stw[idtblin,0] += 1 #特征使用的次数
        for l in range(len(idtbnlin)):  # for nonlinear features 对于非线性特征
            indi = idtbnlin[l].split("#") #将idtbnlin[l]按照"#"分割
            valai = np.sqrt(np.abs(wtb[l+len(idtblin)])) #非线性特征的权重取绝对值 然后开根号
            valsi = np.sign(wtb[l+len(idtblin)])*valai #非线性特征的权重
            if indi[0][0] == '-': #式8
                wtM[int(indi[0][1:]),0] += valai
                wtS[int(indi[0][1:]),0] += valsi
                stw[int(indi[0][1:]),0] += 1
                wtM[int(indi[1]),0] += valai
                wtS[int(indi[1]),0] += -valsi
                stw[int(indi[1]),0] += 1
            else: #式8
                wtM[int(indi[0]),0] += valai
                wtS[int(indi[0]),0] += valsi
                stw[int(indi[0]),0] += 1
                wtM[int(indi[1]),0] += valai
                wtS[int(indi[1]),0] += valsi
                stw[int(indi[1]),0] += 1
    
    wt = np.sign(wtS)*wtM/stw    # Compute weights as average 将权重计算为平均值
    #wt根据式4算出
    j = 0
    for reg in regs:
        if not reg == t:#如果调节因子不等于靶基因
            val = wt[j,0] #靶基因t的第j个调控因子的权重
            if (abs(val) > 0.0):#如果权重的绝对值大于0
                edges.append(reg+"\t"+t)#edges存储调控因子的名字+靶基因的名字
                w.append(val)#w存储权重值
            j += 1

#循环后得到w和edges
sortindex = np.flip(np.argsort(np.abs(w))[np.newaxis],1).astype(np.int)   # Sort by absolute value of edge weights 按照边权重的绝对值排序
sedgevals = [w[s] for s in sortindex[0,:]]#权重值的排序
sedgevals = sedgevals/abs(sedgevals[0])
sedges = [edges[s] for s in sortindex[0,:]]#按照权重值对基因对进行排序

ofile = open(outfilename, 'w')
for s in range(len(sedges)):
    print(sedges[s] + "\t" + str(sedgevals[s]),file=ofile)                # write to output edge list file 输出边缘列表文件
        
ofile.close()
