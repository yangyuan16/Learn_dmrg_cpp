import numpy as np
import pandas as pd
import matplotlib.pylab as plt
from scipy.optimize import leastsq # 导入 scipy 中的最小二乘法拟合工具
from scipy.stats import linregress # 导入 scipy 中的线性回归工具
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
    x = df['r'].values
    yObs = df['logcorre'].values
    pFit, info = leastsq(error1, p0, args=(x,yObs)) # 最小二乘法拟合参数
    print("Data fitting with Scipy.optimize.leastsq")
    print("y = p[0] + p[1] * x")
    print("p[0] = {:.4f}\np[1] = {:.4f}".format(pFit[0], pFit[1]))
    intercept = pFit[0]
    slope = pFit[1]
    return intercept, slope    # 输出斜率和截距
#
def fit_logr_r_linregress(df): # 采用 linear regress 方法做拟合
    x = df['r'].values
    yObs = df['logcorre'].values
    slope, intercept, r_value, p_value, std = linregress(x, yObs)
    print("\nLinear regress with Scipy.stats.linregress")
    print("y = p[0] + p[1] * x")
    print("p[0] = {:.4f}".format(intercept))  # 输出截距 intercept
    print("p[1] = {:.4f}".format(slope))  # 输出斜率 slope
    print("r^2_value: {:.4f}".format(r_value**2))  # 输出 r^2 值
    print("p_value: {:.4f}".format(p_value))  # 输出 p 值
    print("std: {:.4f}".format(std))  # 输出标准差 std
    return intercept, slope   # 输出斜率和截距
#
def fit_logr_logr_leastsq(df): # 采用 leastsq 的方法做拟合
    p0 = [1, 1] # 设置拟合函数的参数初值, 截距和斜率
    x = df['logr'].values
    yObs = df['logcorre'].values
    pFit, info = leastsq(error1, p0, args=(x,yObs)) # 最小二乘法拟合参数
    print("Data fitting with Scipy.optimize.leastsq")
    print("y = p[0] + p[1] * x")
    print("p[0] = {:.4f}\np[1] = {:.4f}".format(pFit[0], pFit[1]))
    intercept = pFit[0]
    slope = pFit[1]
    return intercept, slope  #输出截距和斜率
#
def fit_logr_r_linregress(df): # 采用 linear regress 方法做拟合
    x = df['logr'].values
    yObs = df['logcorre'].values
    slope, intercept, r_value, p_value, std = linregress(x, yObs)
    print("\nLinear regress with Scipy.stats.linregress")
    print("y = p[0] + p[1] * x")
    print("p[0] = {:.4f}".format(intercept))  # 输出截距 intercept
    print("p[1] = {:.4f}".format(slope))  # 输出斜率 slope
    print("r^2_value: {:.4f}".format(r_value**2))  # 输出 r^2 值
    print("p_value: {:.4f}".format(p_value))  # 输出 p 值
    print("std: {:.4f}".format(std))  # 输出标准差 std
    return intercept, slope   # 输出斜率和截距
#
def get_fit_data_by_logr_r(df, intercept, slope): # 指数型函数拟合
    x = df["r"].values
    y = intercept + slope * x
    df["fitlogcorre_exp"] = y
    df["fitcorre_exp"] = [np.power(10,i) for i in y]
    df["intercept_exp"] = intercept
    df["slope_exp"] = slope
    return df
#
def get_fit_data_by_logr_logr(df, intercept, slope): # power law 型函数拟合
    x = df["logr"].values
    y = intercept + slope * x
    df["fitlogcorre_pow"] = y
    df["fitcorre_pow"] = [np.power(10,i) for i in y]
    df["intercept_pow"] = intercept
    df["slope_pow"] = slope
    return df
#
def get_fit_data(df0, df_cut_exp, df_cut_pow):
    intercept_exp, slope_exp = fit_logr_r_leastsq(df_cut_exp) # 利用 df_cut 数据点得到拟合的指数型函数拟合的截距和斜率 
    print("---利用 df_cut 数据点得到拟合的exponential type型函数拟合的截距和斜率---")
    print("intercept_exp: ", intercept_exp)
    print("slope_exp: ", slope_exp)
    df0 = get_fit_data_by_logr_r(df=df0, intercept=intercept_exp, slope=slope_exp) # 将拟合数据添加到 df      
    #
    intercept_pow, slope_pow = fit_logr_logr_leastsq(df_cut_pow) # 利用 df_cut 数据点得到拟合的powerlaw型函数拟合的截距和斜率
    print("---利用 df_cut 数据点得到拟合的power-law type型函数拟合的截距和斜率---")
    print("intercept_pow: ", intercept_pow)
    print("slope_pow: ", slope_pow)
    df0 = get_fit_data_by_logr_logr(df=df0, intercept=intercept_pow, slope=slope_pow) # 将拟合数据添加到 df
    return df0
#
def plot_expfit(df):
    #-----------------plot ax1-------------------------
    fig = plt.figure(figsize=(6,10)) 
    ax1 = plt.subplot(1,1,1)
    plt.sca(ax1)  ##选择对ax1进行绘图
    ax1=plt.gca() #获得坐标轴的句柄
    #
    r = df["r"].values
    corre_abs = df["corre_abs"].values
    fitcorre = df["fitcorre_exp"].values
    slope = df["slope_exp"].values[0]
    xi = -1/slope 
    label_data = "< $ S_i S_j $ >"
    L1, = ax1.plot(r,corre_abs,label=label_data,ls="-",lw=1.5,color="red",
             marker='o',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="k",
             markerfacecolor='w')
    label_fitdata = r"slope: {:.4f} $\xi$: {:.4f}".format(slope,xi)
    L2, = ax1.plot(r,fitcorre,label=label_fitdata,ls="--",lw=1.5,color="blue",
             marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k",
             markerfacecolor='w')
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 26, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[L1,L2], loc = 4, bbox_to_anchor=(0.48, 0.18),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    label_x = r"|i-j|"
    label_y = "Spin Correlation"
    plt.yscale("log")
    ax1.set_xlabel(label_x, size= 25)
    ax1.set_ylabel(label_y, size= 25)
    ax1.tick_params(labelsize = 25) # 设置坐标刻度对应数字的大小
    #ax1.set_xlim([0,8])
    #ax1.set_ylim([-0.1,1])
    #ax1.set_xticks([0,2,4,6,8])
    #ax1.set_yticks([-0.1,0,0.5,1])
    #ax1.text(0.3,-0.05, r'$\mathrm{(a)}$', fontsize=18)
    plt.title("Jz=%.2f"%Jz,fontsize=25)
    plt.show()
    return
#
def plot_powfit(df):
    #-----------------plot ax1-------------------------
    fig = plt.figure(figsize=(6,10)) 
    ax1 = plt.subplot(1,1,1)
    plt.sca(ax1)  ##选择对ax1进行绘图
    ax1=plt.gca() #获得坐标轴的句柄
    #
    r = df["r"].values
    corre_abs = df["corre_abs"].values
    fitcorre = df["fitcorre_pow"].values
    slope = df["slope_pow"].values[0]
    label_data = "< $S_i S_j$ >"
    L1, = ax1.plot(r,corre_abs,label=label_data,ls="-",lw=1.5,color="red",
             marker='o',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="k",
             markerfacecolor='w')
    label_fitdata = "K:{:.4f}".format(slope)
    L2, = ax1.plot(r,fitcorre,label=label_fitdata,ls="--",lw=1.5,color="blue",
             marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k",
             markerfacecolor='w')
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 26, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[L1,L2], loc = 4, bbox_to_anchor=(0.48, 0.18),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    label_x = r"|i-j|"
    label_y = "Spin Correlation"
    plt.yscale("log")
    plt.xscale("log")
    ax1.set_xlabel(label_x, size= 25)
    ax1.set_ylabel(label_y, size= 25)
    ax1.tick_params(labelsize = 25) # 设置坐标刻度对应数字的大小
    #ax1.set_xlim([0,8])
    #ax1.set_ylim([-0.1,1])
    #ax1.set_xticks([0,2,4,6,8])
    #ax1.set_yticks([-0.1,0,0.5,1])
    #ax1.text(0.3,-0.05, r'$\mathrm{(a)}$', fontsize=18)
    plt.title("Jz=%.2f"%Jz, fontsize=25)
    plt.show()
    return
#
if __name__ == "__main__":
    print()
    Lz = 1
    Ly = 10
    Lx = 32
    dop = 40
    t = 3
    J = 1
    Jz = 0
    dim = 8000 # dim cutoff
    workpath = "E:\\WORK\\Work\\Project\\leg10-tj"
    filepath1 = "\\data_dmrgcpp\\Lz%d_Ly%d_Lx%d\\dop%g" % (Lz, Ly, Lx,dop)
    filepath2 = "\\t%d_J%d_Jz%d_dim%d" % (t, J, Jz, dim)
    filepath3 = "\\measurement_spin_correlation.dat"
    filename = workpath + filepath1 + filepath2 + filepath3
    print(filename)
    # load the data
    df = pd.read_csv(filename, header=None, sep='\t',encoding='utf-8')
    df.rename(columns={0: "site1", 1: "site2", 2: "corre"},inplace=True)
    df = df[(df["site2"]-df["site1"]) % (Lz * Ly) == 0] 
    df.sort_values(["site2"],inplace=True)
    df.index = list(range(len(df)))
    print(df.head())
    print(df.tail())
    print(len(df))
    #==============================================
    site1 = df["site1"].values
    site2 = df["site2"].values
    corre_abs = np.abs(np.array(df["corre"].values))
    r = np.array((site2 - site1) / (Lz * Ly))
    print(r)
    #
    df["r"] = r
    df["corre_abs"] = corre_abs
    df["logr"] = np.log10(r)
    df["logcorre"] = np.log10(corre_abs)
    #
    print(df.head())
    print(df.tail())
    print(len(df))
    #================选定某些数据点==========================
    df_cut_exp = df.loc[[9,13,17,21]] # exponential fit
    df_cut_pow = df.loc[[0,4,10]] # power law fit
    #========# 利用拟合得到的截距和斜率得到拟合的corre数据======
    df = get_fit_data(df0=df, df_cut_exp=df_cut_exp, df_cut_pow=df_cut_pow) 
    print(df.head())
    print(df.tail())
    print(len(df))
    #
    plot_expfit(df=df)
    plot_powfit(df=df)
    #---------------save data --------------------------
    Issave = True
    if Issave:
        filepath4 = "\\measurement_spin_correlation_fit.parquet"
        savepath = workpath + filepath1 + filepath2 + filepath4
        df.to_parquet(savepath)