import numpy as np
import pandas as pd
import matplotlib.pylab as plt
from scipy.optimize import leastsq # 导入 scipy 中的最小二乘法拟合工具
from scipy.stats import linregress # 导入 scipy 中的线性回归工具
#--------------------------------------------
def df_along_Lx(df,Lz,Ly):
    if Ly == 2:
        print("---------df_res0---------")
        df_res0 = df[df["site"]%(Lz*Ly)==0]
        df_res0.sort_values(["site"],inplace=True)
        df_res0["r"] = range(1,len(df_res0)+1) 
        print(df_res0.head())
        print(df_res0.tail())
        print(len(df_res0))
        print("---------df_res1---------")
        df_res1 = df[df["site"]%(Lz*Ly)==1]
        df_res1.sort_values(["site"],inplace=True)
        df_res1["r"] = range(1,len(df_res1)+1) 
        print(df_res1.head())
        print(df_res1.tail())
        print(len(df_res1))
        print("---------df_res2---------")
        df_res2 = df[df["site"]%(Lz*Ly)==2]
        df_res2.sort_values(["site"],inplace=True)
        df_res2["r"] = range(1,len(df_res2)+1) 
        print(df_res2.head())
        print(df_res2.tail())
        print(len(df_res2))
        print("---------df_res3---------")
        df_res3 = df[df["site"]%(Lz*Ly)==3]
        df_res3.sort_values(["site"],inplace=True)
        df_res3["r"] = range(1 , len(df_res3)+1) 
        print(df_res3.head())
        print(df_res3.tail())
        print(len(df_res3))
    return df_res0, df_res1, df_res2, df_res3
#
def build_df_res1(Lx,df_res1):
    x = df_res1["r"].values
    logr = entropy_logr(Lx=Lx,x=x)
    df_res1["logr"] = logr
    df_res1.index = range(len(df_res1))
    print("ddd")
    return df_res1
#
def build_df_res2(Lx,df_res2):
    x = df_res2["r"].values
    logr = entropy_logr(Lx=Lx,x=x)
    df_res2["logr"] = logr
    df_res2.index = range(len(df_res2))
    return df_res2
#
def entropy(Lx, x, c, g):
    S = []
    for it in x:
        s = (c/6)*np.log((Lx/np.pi)*np.sin(it*np.pi/Lx)) + g
        S.append(s)    
    return S
#
def entropy_logr(Lx, x):
    logr_list = []
    for it in x:
        logr_list.append((1/6) * np.log((Lx/np.pi)*np.sin(it*np.pi/Lx)) )
    return logr_list
#
def fitfunc1(p, x): # 定义拟合函数为直线
    p0, p1 = p #拟合函数的参数
    y = p0 + p1 * x #拟合函数的表达式
    return y
#
def error1(p, x, y): # 定义观测值与拟合函数值的误差函数
    err = fitfunc1(p,x) - y # 误差
    return err
#
def fit_logr_r_leastsq(df): # 采用 leastsq 的方法做拟合
    p0 = [1, 1] # 设置拟合函数的参数初值, 截距和斜率
    x = df['logr'].values
    yObs = df['entropy'].values
    pFit, info = leastsq(error1, p0, args=(x,yObs)) # 最小二乘法拟合参数
    print("Data fitting with Scipy.optimize.leastsq")
    print("y = p[0] + p[1] * x")
    print("p[0] = {:.4f}\np[1] = {:.4f}".format(pFit[0], pFit[1]))
    intercept = pFit[0]
    slope = pFit[1]
    return intercept, slope    # 输出斜率和截距
#
def get_fit_data_by_logr_r(df, intercept, slope): # 指数型函数拟合
    x = df["logr"].values
    y = intercept + slope * x
    df["fitentropy"] = y
    df["intercept"] = intercept
    df["slope"] = slope
    return df
#
def get_fit_data(df0,df_cut):
    intercept, slope = fit_logr_r_leastsq(df_cut)  
    print("---利用 df_cut 数据点得到拟合的exponential type型函数拟合的截距和斜率---")
    print("intercept: ", intercept)
    print("slope: ", slope)
    df0 = get_fit_data_by_logr_r(df=df0, intercept=intercept, slope=slope) # 将拟合数据添加到 df      
    return df0
#
def plot_entropyfit_res1(df):
    #-----------------plot ax1-------------------------
    fig = plt.figure(figsize=(6,10)) 
    ax1 = plt.subplot(1,1,1)
    plt.sca(ax1)  ##选择对ax1进行绘图
    ax1=plt.gca() #获得坐标轴的句柄
    #
    logr = df["logr"].values
    entropy = df["entropy"].values
    fitentropy = df["fitentropy"].values
    slope = df["slope"].values[0]
    intercept = df["intercept"].values[0]
    label_data = "entropy"
    L1, = ax1.plot(logr,entropy,label=label_data,ls="-",lw=1.5,color="blue",
             marker='o',alpha=1,markersize=8,markeredgewidth=1.5, markeredgecolor="blue",
             markerfacecolor='None')
    label_fitdata = r"c = : {:.4f} g=: {:.4f}".format(slope,intercept)
    L2, = ax1.plot(logr,fitentropy,label=label_fitdata,ls="--",lw=1.5,color="k",
             marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k",
             markerfacecolor='w')
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 26, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[L1,L2], loc = 4, bbox_to_anchor=(0.48, 0.18),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    label_x = r"(1/6)log((Lx/$\pi$)sin(x$\pi$/Lx))"
    label_y = "Entropy"
    #plt.yscale("log")
    ax1.set_xlabel(label_x, size= 25)
    ax1.set_ylabel(label_y, size= 25)
    ax1.tick_params(labelsize = 25) # 设置坐标刻度对应数字的大小
    #ax1.set_xlim([0,8])
    #ax1.set_ylim([-0.1,1])
    #ax1.set_xticks([0,2,4,6,8])
    #ax1.set_yticks([-0.1,0,0.5,1])
    #ax1.text(0.3,-0.05, r'$\mathrm{(a)}$', fontsize=18)
    plt.title("dop=%.4f,Jz=%.2f,dim=%d"%(dop/(Lz*Ly*Lx),Jz,dim), fontsize=25)
    plt.show()
    return
#
def plot_entropyfit_res2(df):
    #-----------------plot ax1-------------------------
    fig = plt.figure(figsize=(6,10)) 
    ax1 = plt.subplot(1,1,1)
    plt.sca(ax1)  ##选择对ax1进行绘图
    ax1=plt.gca() #获得坐标轴的句柄
    #
    logr = df["logr"].values
    entropy = df["entropy"].values
    fitentropy = df["fitentropy"].values
    slope = df["slope"].values[0]
    intercept = df["intercept"].values[0]
    label_data = "entropy"
    L1, = ax1.plot(logr,entropy,label=label_data,ls="-",lw=1.5,color="green",
             marker='o',alpha=1,markersize=8,markeredgewidth=1.5, markeredgecolor="green",
             markerfacecolor='None')
    label_fitdata = r"c = : {:.4f} g=: {:.4f}".format(slope,intercept)
    L2, = ax1.plot(logr,fitentropy,label=label_fitdata,ls="--",lw=1.5,color="blue",
             marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k",
             markerfacecolor='w')
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 26, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[L1,L2], loc = 4, bbox_to_anchor=(0.48, 0.18),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    label_x = r"(1/6)log((Lx/$\pi$)sin(x$\pi$/Lx))"
    label_y = "Entropy"
    #plt.yscale("log")
    ax1.set_xlabel(label_x, size= 25)
    ax1.set_ylabel(label_y, size= 25)
    ax1.tick_params(labelsize = 25) # 设置坐标刻度对应数字的大小
    #ax1.set_xlim([0,8])
    #ax1.set_ylim([-0.1,1])
    #ax1.set_xticks([0,2,4,6,8])
    #ax1.set_yticks([-0.1,0,0.5,1])
    #ax1.text(0.3,-0.05, r'$\mathrm{(a)}$', fontsize=18)
    plt.title("dop=%.4f,Jz=%.2f,dim=%d"%(dop/(Lz*Ly*Lx),Jz,dim), fontsize=25)
    plt.show()
    return
#
if __name__ == "__main__":
    print()
    Lz = 2
    Ly = 2
    Lx = 48
    dop = 36
    t = 3
    J = 1
    Jz = 1.0 
    dim = 6000 # dim cutoff
    workpath = "E:\\WORK\\Work\\Project\\La3Ni2O7"
    filepath1 = "\\data_dmrgcpp\\Lz%d_Ly%d_Lx%d\\dop%g" % (Lz, Ly, Lx,dop)
    filepath2 = "\\t%d_J%d_Jz%.2f_dim%d_ent\\entanglement\\entropy" % (t, J, Jz, dim)
    filepath3 = "\\sys-ground_state"
    filename = workpath + filepath1 + filepath2 + filepath3
    print(filename)
    # load the data
    df = pd.read_csv(filename, header=None, sep='\t',encoding='utf-8')
    df.rename(columns={0: "site", 1: "entropy",}, inplace=True)
    L1 = Lz * Ly * Lx - 3
    L2 = len(df)
    df = df.iloc[L2-L1:L2]
    df.sort_values(["site"],inplace=True)
    print(df.head())
    print(df.tail())
    print(len(df))
    df_res0, df_res1, df_res2, df_res3 = df_along_Lx(df=df,Lz=Lz,Ly=Ly)
    #--------------------------------------------------------------------
    '''
    # choose df_res1
    df_res1 = build_df_res1(df_res1=df_res1, Lx=Lx)
    print('--------------------')
    print(df_res1.head())
    print(df_res1.tail())
    # cut half of Lx of df_res1
    df_res1_cut = df_res1.iloc[2:np.int32(Lx/2),:]
    df_res1_cut.index = range(len(df_res1_cut))
    print("df_res1_cut:\n", df_res1_cut)
    #         选定某些数据点
    df_cut = df_res1_cut.loc[[3,15,16,]] 
    df_fit = get_fit_data(df0=df_res1_cut, df_cut=df_cut)
    print("df_fit:\n",df_fit)
    #
    plot_entropyfit_res1(df=df_fit)
    #   
    Issave = True 
    if Issave:
        filepath4 = "\\t%d_J%d_Jz%.2f_dim%d_ent\\entanglement" % (t, J, Jz, dim)
        filepath5 = "\\measurement_entropy_fit_res1.parquet" 
        savepath = workpath + filepath1 + filepath4 + filepath5
        df_fit.to_parquet(savepath)
    '''
    #----------------------------------------------------------
    # choose df_res2
    df_res2 = build_df_res2(df_res2=df_res2,Lx=Lx)
    print('---------------------')
    print(df_res2.head())
    print(df_res2.tail())
    # cut half of Lx of df_res2
    df_res2_cut = df_res2.iloc[2:np.int32(Lx/2),:]
    df_res2_cut.index = range(len(df_res2_cut))
    print("df_res2_cut:\n", df_res2_cut)
    #             选定某些数据点
    df_cut = df_res2_cut.loc[[10,15,18,20]] 
    df_fit = get_fit_data(df0=df_res2_cut, df_cut=df_cut)
    print("df_fit:\n",df_fit)
    #
    plot_entropyfit_res2(df=df_fit)
    #   
    Issave = True 
    if Issave:
        filepath4 = "\\t%d_J%d_Jz%.2f_dim%d_ent\\entanglement" % (t, J, Jz, dim)
        filepath5 = "\\measurement_entropy_fit.parquet" 
        savepath = workpath + filepath1 + filepath4 + filepath5
        df_fit.to_parquet(savepath)