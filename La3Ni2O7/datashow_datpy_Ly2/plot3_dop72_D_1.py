# 画粒子数密度关联和局域密度分布
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
#
def density_along_x_Ly2(df,Ly,Lz): # 沿着 x 方向的 electron density
    df_out = pd.DataFrame(columns = ["r0","dy0","r1","dy1","r2","dy2","r3","dy3","rmean","dymean"])
    #
    df_y0 = df[df['site'] % (Ly * Lz) ==0]
    df_out["r0"] = df_y0["site"].values
    df_out["dy0"] = df_y0["density"].values
    
    df_y1 = df[df['site'] % (Ly * Lz) ==1]
    df_out["r1"] = df_y1["site"].values
    df_out["dy1"] = df_y1["density"].values
    
    df_y2 = df[df['site'] % (Ly * Lz) ==2]
    df_out["r2"] = df_y2["site"].values
    df_out["dy2"] = df_y2["density"].values

    df_y3 = df[df['site'] % (Ly * Lz) ==3]
    df_out["r3"] = df_y3["site"].values
    df_out["dy3"] = df_y3["density"].values

    r_mean = range(len(df_y0))
    density_mean = (np.array(df_out["dy0"].values) + np.array(df_out["dy1"].values) + np.array(df_out["dy2"].values)  + 
                    np.array(df_out["dy3"].values) ) / (Lz*Ly)
    #
    df_out["rmean"] = r_mean
    df_out["dymean"] = density_mean
    #
    return df_out
#
#
def get_data_ni(Lz, Ly, Lx, dop, t, J, Jz, dim,):
    workpath = "E:\\WORK\\Work\\Project\\La3Ni2O7"
    filepath1 = "\\data_dmrgcpp\\Lz%d_Ly%d_Lx%d\\dop%g" % (Lz, Ly, Lx,dop)
    filepath2 = "\\t%d_J%d_Jz%.2f_dim%d" % (t, J, Jz, dim)
    filepath3 = "\\measurement_electron_density.dat"
    filename = workpath + filepath1 + filepath2 + filepath3
    print(filename)
    # load the data
    df = pd.read_csv(filename, header=None, sep='\t',encoding='utf-8')
    df.rename(columns={0: "site", 1: "density"},inplace=True)
    df.sort_values(['site'],inplace=True)
    df_layer1 = df[df["site"] % 2 == 0] 
    df_layer2 = df[df["site"] % 2 == 1]
    sites_layer1 = df_layer1["site"].values.reshape(-1,Ly).T
    density_layer1 = df_layer1["density"].values.reshape(-1, Ly).T
    sites_layer2 = df_layer2["site"].values.reshape(-1, Ly).T
    density_layer2 = df_layer2["density"].values.reshape(-1, Ly).T
    print(df.head())
    print(len(df))
    df_out = density_along_x_Ly2(df=df,Ly=Ly,Lz=Lz)
    print(df_out.head())
    return df_out
#
def get_data_ninj_corre(Lz, Ly, Lx, dop, t, J, Jz, dim,):
    # load data
    workpath = "E:\\WORK\\Work\\Project\\La3Ni2O7"
    filepath1 = "\\data_dmrgcpp\\Lz%d_Ly%d_Lx%d\\dop%g" % (Lz, Ly, Lx,dop)
    filepath2 = "\\t%d_J%d_Jz%.2f_dim%d" % (t, J, Jz, dim)
    filepath3 = "\\measurement_density_correlation_fit.parquet"   
    filename = workpath + filepath1 + filepath2 + filepath3
    print(filename)
    df = pd.read_parquet(filename)
    print(df.head())
    return df
#
if __name__ =="__main__":
    Lz = 2
    Ly = 2
    Lx = 48
    dop = 72
    t = 3
    J = 1
    dim = 6000
    # local ninj data
    df_ni_Jz05 = get_data_ni(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=0.5,dim=dim)
    df_ni_Jz10 = get_data_ni(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=1.0,dim=dim)
    df_ni_Jz20 = get_data_ni(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=2.0,dim=dim)
    df_ni_Jz40 = get_data_ni(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=4.0,dim=dim)
    # ninj correlation
    df_ninj_Jz02 = get_data_ninj_corre(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=0.2,dim=dim)
    df_ninj_Jz04 = get_data_ninj_corre(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=0.4,dim=dim)
    df_ninj_Jz06 = get_data_ninj_corre(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=0.6,dim=dim)
    df_ninj_Jz08 = get_data_ninj_corre(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=0.8,dim=dim)
    df_ninj_Jz10 = get_data_ninj_corre(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=1.0,dim=dim)
    df_ninj_Jz20 = get_data_ninj_corre(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=2.0,dim=dim)
    df_ninj_Jz40 = get_data_ninj_corre(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=4.0,dim=dim)
    #---------------plot logr-r fig and logr-logr fig----------------------------------
    fig = plt.figure(figsize=(10,10))
    # plt.figure(facecolor='blue',edgecolor='black') # 设置画布的颜色
    params = {
        'axes.labelsize': '30',
        'xtick.labelsize':'20',
        'ytick.labelsize':'20',
        'ytick.direction':'in',
        'xtick.direction':'in',
        'lines.linewidth':6 ,
        'legend.fontsize': '27',
        # 'figure.figsize'   : '12, 9'    # set figure size
    }
    pylab.rcParams.update(params) # set figure parameter 更新绘图的参数
    #plt.rcParams['font.family'] = 'Times New Roman'  # 设置全局字体为 Times New Roman
    # 得到子图
    ax1 = plt.axes([0.1,0.73,0.87,0.25])
    #ax2 = plt.axes([0.58,0.58,0.39,0.4])
    ax2 = plt.axes([0.1,0.1,0.39,0.55])
    ax3 = plt.axes([0.58,0.1,0.39,0.55])
    #
    # plt.subplots_adjust(left=0.1, bottom=0.1, right=0.98, top=0.96, wspace=0.32, hspace=0.26)
    # 选择子图 ax1 进行绘图
    plt.sca(ax1) ## 选择对 ax1 进行绘图
    ax1 = plt.gca()
    L, =ax1.plot([0,50],[0.625,0.625],label=" ",ls="--",lw=1.5,color="k",
             marker='o',alpha=1,markersize=0,markeredgewidth=1.5, markeredgecolor="red", markerfacecolor='w')
    label = r"$J_{\bot}$=0.5"
    L05, =ax1.plot(df_ni_Jz05["rmean"].values,df_ni_Jz05["dymean"].values,label=label,ls="--",lw=1.5,color="red",
             marker='o',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="red", markerfacecolor='w')
    label = r"$J_{\bot}$=1.0"
    L10, =ax1.plot(df_ni_Jz10["rmean"].values,df_ni_Jz10["dymean"].values,label=label,ls="--",lw=1.5,color="blue",
             marker='^',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="blue", markerfacecolor='w')
    label = r"$J_{\bot}$=2.0"
    L20, =ax1.plot(df_ni_Jz20["rmean"].values,df_ni_Jz20["dymean"].values,label=label,ls="--",lw=1.5,color="green",
             marker='v',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="green", markerfacecolor='w')
    label = r"$J_{\bot}$=4.0"
    L40, =ax1.plot(df_ni_Jz40["rmean"].values,df_ni_Jz40["dymean"].values,label=label,ls="--",lw=1.5,color="magenta",
             marker='s',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="magenta", markerfacecolor='w')
    ####图例设置
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 15, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[L05,L10,L20,L40,], loc = 4, bbox_to_anchor=(0.88, 0.01),
                       ncol = 4,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)    
    #
    label_x = r"x"
    label_y = "n(x)"
    #plt.yscale("log")
    #plt.xscale("log")
    ax1.set_xlabel(label_x, size= 16)
    ax1.set_ylabel(label_y, size= 16)
    ax1.tick_params(labelsize = 15) # 设置坐标刻度对应数字的大小
    ax1.set_xlim([0,50])
    ax1.set_xticks([0,10,20,30,40,50])
    #ax1.set_yticks([-1,-0.5,0,0.5,1]) 
    #=========================================================
    ax1.text(4,0.675,"(a)",fontsize = 16, color='black', rotation = 0)
    ax1.text(20,0.675, r'$\mathrm{\delta} = \frac{3}{8}$', fontsize = 16, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
    #=======================================================================================
    # 坐标轴设置第一层
    labels = ax1.get_xticklabels() + ax1.get_yticklabels()
    #[label.set_fontname('Times New Roman') for label in labels]###设置ticket labled的字体格式
    ax1.xaxis.set_minor_locator(MultipleLocator(5))###设置次刻度的间隔
    #ax1.yaxis.set_minor_locator(MultipleLocator(10))###设置次刻度的间隔
    ax1.xaxis.set_major_formatter(FormatStrFormatter('%1.0f'))###设置X轴标签文本格式
    #ax1.yaxis.set_major_formatter(FormatStrFormatter('%1.1f'))###设置Y轴标签文本格式
    #
    #ax1.yaxis.get_major_locator().set_params(numticks=99)
    #ax1.yaxis.get_minor_locator().set_params(numticks=99, subs=[.2,.4,.6,.8]) # 将次要刻度显示出来 
    #----将次要刻度显示出来 
    #locmin = matplotlib.ticker.LogLocator(base=10.0, subs=(0.1,0.2,0.4,0.6,0.8,1,2,4,6,8,10 )) 
    #ax1.xaxis.set_minor_locator(locmin)
    #ax1.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    #-------------------------------------------------------------------------------------
    # 选择子图 ax2 进行绘图
    plt.sca(ax2) ## 选择对 ax2 进行绘图
    ax2 = plt.gca()
    #--------- J_\bot = 0.2
    slope = df_ninj_Jz02["slope_exp"].values[0]
    xi = round(-1/slope,2)  
    label = r"$J_{\bot}$=%.1f $\xi_{CDW}$=%.2f"%(0.2,xi)
    L02, = ax2.plot(df_ninj_Jz02["r"],df_ninj_Jz02["corre_abs"],label=label,ls="-",lw=1,color="r",
             marker='o',alpha=1,markersize=4,markeredgewidth=1, markeredgecolor="k", markerfacecolor='w')
    L02_, = ax2.plot(df_ninj_Jz02["r"],df_ninj_Jz02["fitcorre_exp"],label=label,ls="--",lw=1,color="k",
             marker='o',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- J_\bot = 0.4
    slope = df_ninj_Jz04["slope_exp"].values[0]
    xi = round(-1/slope,2)  
    label = r"$J_{\bot}$=%.1f $\xi_{CDW}$=%.2f"%(0.4,xi)
    L04, = ax2.plot(df_ninj_Jz04["r"],df_ninj_Jz04["corre_abs"],label=label,ls="-",lw=1,color="blue",
             marker='s',alpha=1,markersize=4,markeredgewidth=1, markeredgecolor="k", markerfacecolor='w')
    L04_, = ax2.plot(df_ninj_Jz04["r"],df_ninj_Jz04["fitcorre_exp"],label=label,ls="--",lw=1,color="k",
             marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- J_\bot = 0.6
    slope = df_ninj_Jz06["slope_exp"].values[0]
    xi = round(-1/slope,2)  
    label = r"$J_{\bot}$=%.1f $\xi_{CDW}$=%.2f"%(0.6,xi)
    L06, = ax2.plot(df_ninj_Jz06["r"],df_ninj_Jz06["corre_abs"],label=label,ls="-",lw=1,color="green",
             marker='^',alpha=1,markersize=4,markeredgewidth=1, markeredgecolor="k", markerfacecolor='w')
    L06_, = ax2.plot(df_ninj_Jz06["r"],df_ninj_Jz06["fitcorre_exp"],label=label,ls="--",lw=1,color="k",
             marker='^',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- J_\bot = 0.8
    slope = df_ninj_Jz08["slope_exp"].values[0]
    xi = round(-1/slope,2)  
    label = r"$J_{\bot}$=%.1f $\xi_{CDW}^{yy}$=%.2f"%(0.8,xi)
    L08, = ax2.plot(df_ninj_Jz08["r"],df_ninj_Jz08["corre_abs"],label=label,ls="-",lw=1,color="magenta",
             marker='v',alpha=1,markersize=4,markeredgewidth=1, markeredgecolor="k", markerfacecolor='w')
    L08_, = ax2.plot(df_ninj_Jz08["r"],df_ninj_Jz08["fitcorre_exp"],label=label,ls="--",lw=1,color="k",
             marker='v',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 1.0
    slope = df_ninj_Jz10["slope_exp"].values[0]
    xi = round(-1/slope,2)  
    label = r"$J_{\bot}$=%.1f $\xi_{CDW}$=%.2f"%(1.0,xi)
    L10, = ax2.plot(df_ninj_Jz10["r"],df_ninj_Jz10["corre_abs"],label=label,ls="-",lw=1,color="cyan",
             marker='+',alpha=1,markersize=4,markeredgewidth=1, markeredgecolor="k", markerfacecolor='w')
    L10_, = ax2.plot(df_ninj_Jz10["r"],df_ninj_Jz10["fitcorre_exp"],label=label,ls="--",lw=1,color="k",
             marker='+',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 2.0
    slope = df_ninj_Jz20["slope_exp"].values[0]
    xi = round(-1/slope,2)  
    label = r"$J_{\bot}$=%.1f $\xi_{CDW}$=%.2f"%(2.0,xi)
    L20, = ax2.plot(df_ninj_Jz20["r"],df_ninj_Jz20["corre_abs"],label=label,ls="-",lw=1,color="pink",
             marker='*',alpha=1,markersize=4,markeredgewidth=1, markeredgecolor="k", markerfacecolor='w')
    L20_, = ax2.plot(df_ninj_Jz20["r"],df_ninj_Jz20["fitcorre_exp"],label=label,ls="--",lw=1,color="k",
             marker='*',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 4.0
    slope = df_ninj_Jz40["slope_exp"].values[0]
    xi = round(-1/slope,2)  
    label = r"$J_{\bot}$=%.1f $\xi_{CDW}$=%.2f"%(4.0,xi)
    L40, = ax2.plot(df_ninj_Jz40["r"],df_ninj_Jz40["corre_abs"],label=label,ls="-",lw=1,color="brown",
             marker='h',alpha=1,markersize=4,markeredgewidth=1, markeredgecolor="k", markerfacecolor='w')
    L40_, = ax2.plot(df_ninj_Jz40["r"],df_ninj_Jz40["fitcorre_exp"],label=label,ls="--",lw=1,color="k",
             marker='h',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #
    ####图例设置
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 14, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[L02,L04,L06,L08,L10,L20,L40], loc = 4, bbox_to_anchor=(0.98, 0.55),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    label_x = r"r"
    label_y = "|D(r)|"
    #plt.yscale("log")
    #plt.xscale("log")
    ax2.set_yscale("log",base=10,subs=[0.01,0.02,0.03])      
    ax2.set_xlabel(label_x, size= 16)
    ax2.set_ylabel(label_y, size= 16)
    ax2.tick_params(labelsize = 14) # 设置坐标刻度对应数字的大小
    ax2.set_xlim([0,40])
    #ax2.set_xticks([0,10,20,30])
    #ax2.set_yticks([-1,-0.5,0,0.5,1])  
    #=========================================================
    ax2.text(2,2e-6,"(b)",fontsize = 16, color='black', rotation = 0)
    ax2.text(17.5,2e-6, r'$\mathrm{\delta} = \frac{3}{8}$', fontsize = 16, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
    #=========================================================
    # 坐标轴设置第一层
    labels = ax2.get_xticklabels() + ax2.get_yticklabels()
    #[label.set_fontname('Times New Roman') for label in labels]###设置ticket labled的字体格式
    ax2.xaxis.set_minor_locator(MultipleLocator(5))###设置次刻度的间隔
    #ax2.yaxis.set_minor_locator(MultipleLocator(10))###设置次刻度的间隔
    ax2.xaxis.set_major_formatter(FormatStrFormatter('%1.0f'))###设置X轴标签文本格式
    #ax2.yaxis.set_major_formatter(FormatStrFormatter('%1.1f'))###设置Y轴标签文本格式
    #
    #ax2.yaxis.get_major_locator().set_params(numticks=99)
    #ax2.yaxis.get_minor_locator().set_params(numticks=99, subs=[.2,.4,.6,.8]) # 将次要刻度显示出来 
    #----将次要刻度显示出来 
    #locmin = matplotlib.ticker.LogLocator(base=10.0, subs=(0.1,0.2,0.4,0.6,0.8,1,2,4,6,8,10 )) 
    #ax2.xaxis.set_minor_locator(locmin)
    #ax2.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    #----------------------------------------------------------------------------
    #
    # 选择子图 ax3 进行绘图
    plt.sca(ax3) ## 选择对 ax1 进行绘图
    ax3 = plt.gca()
    #--------- J_\bot = 0.2
    slope = df_ninj_Jz02["slope_pow"].values[0]
    K_D = round(-slope,2)  
    label = r"$J_{\bot}$=%.1f $K_{CDW}$=%.2f"%(0.2,K_D)
    L02, = ax3.plot(df_ninj_Jz02["r"],df_ninj_Jz02["corre_abs"],label=label,ls="-",lw=1,color="r",
             marker='o',alpha=1,markersize=4,markeredgewidth=1, markeredgecolor="k", markerfacecolor='w')
    L02_, = ax3.plot(df_ninj_Jz02["r"],df_ninj_Jz02["fitcorre_pow"],label=label,ls="--",lw=1,color="k",
             marker='o',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- J_\bot = 0.4
    slope = df_ninj_Jz04["slope_pow"].values[0]
    K_F = round(-slope,2)  
    label = r"$J_{\bot}$=%.1f $K_{CDW}$=%.2f"%(0.4,K_F)
    L04, = ax3.plot(df_ninj_Jz04["r"],df_ninj_Jz04["corre_abs"],label=label,ls="-",lw=1,color="blue",
             marker='s',alpha=1,markersize=4,markeredgewidth=1, markeredgecolor="k", markerfacecolor='w')
    L04_, = ax3.plot(df_ninj_Jz04["r"],df_ninj_Jz04["fitcorre_pow"],label=label,ls="--",lw=1,color="k",
             marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- J_\bot = 0.6
    slope = df_ninj_Jz06["slope_pow"].values[0]
    K_F = round(-slope,2)  
    label = r"$J_{\bot}$=%.1f $K_{CDW}$=%.2f"%(0.6,K_F)
    L06, = ax3.plot(df_ninj_Jz06["r"],df_ninj_Jz06["corre_abs"],label=label,ls="-",lw=1,color="green",
             marker='^',alpha=1,markersize=4,markeredgewidth=1, markeredgecolor="k", markerfacecolor='w')
    L06_, = ax3.plot(df_ninj_Jz06["r"],df_ninj_Jz06["fitcorre_pow"],label=label,ls="--",lw=1,color="k",
             marker='^',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- J_\bot = 0.8
    slope = df_ninj_Jz08["slope_pow"].values[0]
    K_F = round(-slope,2)  
    label = r"$J_{\bot}$=%.1f $K_{CDW}$=%.2f"%(0.8,K_F)
    L08, = ax3.plot(df_ninj_Jz08["r"],df_ninj_Jz08["corre_abs"],label=label,ls="-",lw=1,color="magenta",
             marker='v',alpha=1,markersize=4,markeredgewidth=1, markeredgecolor="k", markerfacecolor='w')
    L08_, = ax3.plot(df_ninj_Jz08["r"],df_ninj_Jz08["fitcorre_pow"],label=label,ls="--",lw=1,color="k",
             marker='v',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 1.0
    slope = df_ninj_Jz10["slope_pow"].values[0]
    K_F = round(-slope,2)  
    label = r"$J_{\bot}$=%.1f $K_{CDW}$=%.2f"%(1.0,K_F)
    L10, = ax3.plot(df_ninj_Jz10["r"],df_ninj_Jz10["corre_abs"],label=label,ls="-",lw=1,color="cyan",
             marker='+',alpha=1,markersize=4,markeredgewidth=1, markeredgecolor="k", markerfacecolor='w')
    L10_, = ax3.plot(df_ninj_Jz10["r"],df_ninj_Jz10["fitcorre_pow"],label=label,ls="--",lw=1,color="k",
             marker='+',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 2.0
    slope = df_ninj_Jz20["slope_pow"].values[0]
    K_F = round(-slope,2)  
    label = r"$J_{\bot}$=%.1f $K_{CDW}$=%.2f"%(2.0,K_F)
    L20, = ax3.plot(df_ninj_Jz20["r"],df_ninj_Jz20["corre_abs"],label=label,ls="-",lw=1,color="pink",
             marker='*',alpha=1,markersize=4,markeredgewidth=1, markeredgecolor="k", markerfacecolor='w')
    L20_, = ax3.plot(df_ninj_Jz20["r"],df_ninj_Jz20["fitcorre_pow"],label=label,ls="--",lw=1,color="k",
             marker='*',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 4.0
    slope = df_ninj_Jz40["slope_pow"].values[0]
    K_F = round(-slope,2)  
    label = r"$J_{\bot}$=%.1f $K_{CDW}$=%.2f"%(4.0,K_F)
    L40, = ax3.plot(df_ninj_Jz40["r"],df_ninj_Jz40["corre_abs"],label=label,ls="-",lw=1,color="brown",
             marker='h',alpha=1,markersize=4,markeredgewidth=1, markeredgecolor="k", markerfacecolor='w')
    L40_, = ax3.plot(df_ninj_Jz40["r"],df_ninj_Jz40["fitcorre_pow"],label=label,ls="--",lw=1,color="k",
             marker='h',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')

    ####图例设置
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 14, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[L02,L04,L06,L08], loc = 4, bbox_to_anchor=(0.65, 0.02),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    legend2=plt.legend(handles=[L10,L20,L40], loc = 4, bbox_to_anchor=(0.98, 0.70),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    plt.gca().add_artist(legend1)
    label_x = r"r"
    label_y = "$|D(r)|$"
    #plt.yscale("log")
    #plt.xscale("log")
    ax3.set_xscale("log",base=10,subs=[0.01,0.02,0.03]) 
    ax3.set_yscale("log",base=10,subs=[0.01,0.02,0.03])      
    ax3.set_xlabel(label_x, size= 16)
    ax3.set_ylabel(label_y, size= 16)
    ax3.tick_params(labelsize = 14) # 设置坐标刻度对应数字的大小
    ax3.set_xlim([0,40])
    #ax3.set_xticks([0,10,20,30])
    #ax3.set_yticks([-1,-0.5,0,0.5,1])  
    #=========================================================
    ax3.text(28,0.05,"(c)",fontsize = 16, color='black', rotation = 0)
    ax3.text(5.1,0.05, r'$\mathrm{\delta} = \frac{3}{8}$', fontsize = 16, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
    #=========================================================
    # 坐标轴设置第一层
    labels = ax3.get_xticklabels() + ax3.get_yticklabels()
    #[label.set_fontname('Times New Roman') for label in labels]###设置ticket labled的字体格式
    ax3.xaxis.set_minor_locator(MultipleLocator(5))###设置次刻度的间隔
    #ax3.yaxis.set_minor_locator(MultipleLocator(10))###设置次刻度的间隔
    ax3.xaxis.set_major_formatter(FormatStrFormatter('%1.0f'))###设置X轴标签文本格式
    #ax3.yaxis.set_major_formatter(FormatStrFormatter('%1.1f'))###设置Y轴标签文本格式
    #
    #ax3.yaxis.get_major_locator().set_params(numticks=99)
    #ax3.yaxis.get_minor_locator().set_params(numticks=99, subs=[.2,.4,.6,.8]) # 将次要刻度显示出来 
    #----将次要刻度显示出来 
    #locmin = matplotlib.ticker.LogLocator(base=10.0, subs=(0.1,0.2,0.4,0.6,0.8,1,2,4,6,8,10 )) 
    #ax3.xaxis.set_minor_locator(locmin)
    #ax3.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    #-----------------------------------------------------------------------------------------------
    #
    plt.show()