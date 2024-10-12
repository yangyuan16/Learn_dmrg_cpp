# 对比不同 Jz 值下的density correlation function
import numpy as np
import pandas as pd
import matplotlib.pylab as plt
import matplotlib.pylab as pylab
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
def get_data_pzz(Lz, Ly, Lx, dop, t, J, Jz, dim,):
    # load data
    workpath = "E:\\WORK\\Work\\Project\\La3Ni2O7"
    filepath1 = "\\data_dmrgcpp\\Lz%d_Ly%d_Lx%d\\dop%g" % (Lz, Ly, Lx,dop)
    filepath2 = "\\t%d_J%d_Jz%.2f_dim%d" % (t, J, Jz, dim)
    filepath3 = "\\measurement_pairing_fit.parquet"
    filename = workpath + filepath1 + filepath2 + filepath3
    print(filename)
    df = pd.read_parquet(filename)
    print(df.head())
    return df
if __name__ == "__main__":
    print()
    Lz = 2
    Ly = 3
    Lx = 48
    dop = 144
    t = 3
    J = 1
    dim = 8000 # dim cutoff
    #load data
    df_pzz_Jz06 = get_data_pzz(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=0.6,dim=dim)
    df_pzz_Jz10 = get_data_pzz(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=1.0,dim=dim)
    df_pzz_Jz20 = get_data_pzz(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=2.0,dim=dim)
    df_pzz_Jz30 = get_data_pzz(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=3.0,dim=dim)
    df_pzz_Jz40 = get_data_pzz(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=4.0,dim=dim)
    #-----------------plot logr-r fig and logr-logr fig-------------------------
    fig = plt.figure(figsize=(8,10)) 
    ax1 = plt.subplot(1,1,1)
    plt.sca(ax1)  ##选择对ax2进行绘图
    ax1=plt.gca() #获得坐标轴的句柄    
    #--------- J_\bot = 0.6
    slope = df_pzz_Jz06["slope_pow"].values[0]
    K_F = round(-slope,2)  
    label = r"$J_{\bot}$=%.2f $K_{SC}^{zz}$=%.2f"%(0.6,K_F)
    L06, = ax1.plot(df_pzz_Jz06["r"],df_pzz_Jz06["corre_abs"],label=label,ls="-",lw=2,color="r",
             marker='o',alpha=1,markersize=10,markeredgewidth=1, markeredgecolor="r", markerfacecolor='None')
    L06_, = ax1.plot(df_pzz_Jz06["r"],df_pzz_Jz06["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='o',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- J_\bot = 1.0
    slope = df_pzz_Jz10["slope_pow"].values[0]
    K_F = round(-slope,2)  
    label = r"$J_{\bot}$=%.2f $K_{SC}^{zz}$=%.2f"%(1.0,K_F)
    L10, = ax1.plot(df_pzz_Jz10["r"],df_pzz_Jz10["corre_abs"],label=label,ls="-",lw=2,color="blue",
             marker='s',alpha=1,markersize=10,markeredgewidth=1, markeredgecolor="blue", markerfacecolor='None')
    L10_, = ax1.plot(df_pzz_Jz10["r"],df_pzz_Jz10["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- J_\bot = 2.0
    slope = df_pzz_Jz20["slope_pow"].values[0]
    K_F = round(-slope,2)  
    label = r"$J_{\bot}$=%.2f $K_{SC}^{zz}$=%.2f"%(2.0,K_F)
    L20, = ax1.plot(df_pzz_Jz20["r"],df_pzz_Jz20["corre_abs"],label=label,ls="-",lw=2,color="green",
             marker='^',alpha=1,markersize=10,markeredgewidth=1, markeredgecolor="green", markerfacecolor='None')
    L20_, = ax1.plot(df_pzz_Jz20["r"],df_pzz_Jz20["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='^',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- J_\bot = 3.0
    slope = df_pzz_Jz30["slope_pow"].values[0]
    K_F = round(-slope,2)  
    label = r"$J_{\bot}$=%.2f $K_{SC}^{zz}$=%.2f"%(3.0,K_F)
    L30, = ax1.plot(df_pzz_Jz30["r"],df_pzz_Jz30["corre_abs"],label=label,ls="-",lw=2,color="cyan",
             marker='v',alpha=1,markersize=10,markeredgewidth=1, markeredgecolor="cyan", markerfacecolor='None')
    L30_, = ax1.plot(df_pzz_Jz30["r"],df_pzz_Jz30["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='v',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 4.0
    slope = df_pzz_Jz40["slope_pow"].values[0]
    K_F = round(-slope,2)  
    label = r"$J_{\bot}$=%.2f $K_{SC}^{zz}$=%.2f"%(4.0,K_F)
    L40, = ax1.plot(df_pzz_Jz40["r"],df_pzz_Jz40["corre_abs"],label=label,ls="-",lw=2,color="magenta",
             marker='h',alpha=1,markersize=10,markeredgewidth=1, markeredgecolor="magenta", markerfacecolor='None')
    L40_, = ax1.plot(df_pzz_Jz40["r"],df_pzz_Jz40["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='h',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    ####图例设置
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 25, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[L06,L10,L20,L30,L40], loc = 4, bbox_to_anchor=(0.38, 0.02),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    
    label_x = r"r"
    label_y = "$\Phi^{zz}$"
    #plt.yscale("log")
    #plt.xscale("log")
    ax1.set_xscale("log",base=10,subs=[0.01,0.02,0.03]) 
    ax1.set_yscale("log",base=10,subs=[0.01,0.02,0.03])      
    ax1.set_xlabel(label_x, size= 25)
    ax1.set_ylabel(label_y, size= 25)
    ax1.tick_params(labelsize = 25) # 设置坐标刻度对应数字的大小
    ax1.set_xlim([0,40])
    #ax1.set_xticks([0,10,20,30])
    #ax1.set_yticks([-1,-0.5,0,0.5,1])  
    #=========================================================
    #ax1.text(26,0.003,"(d)",fontsize = 20, color='black', rotation = 0)
    #ax1.text(9.5,0.003, r'$\mathrm{\delta} = \frac{1}{2}$', fontsize = 20, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
    #=========================================================
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
    #=====坐标轴的第二层： 坐标轴的设置
    ax1.spines['bottom'].set_linewidth(1.5) ###设置底部坐标轴的粗细
    ax1.spines['left'].set_linewidth(1.5)   ###设置左边坐标轴的粗细
    ax1.spines['right'].set_linewidth(1.5)  ###设置右边坐标轴的粗细
    ax1.spines['top'].set_linewidth(1.5)    ###设置上部坐标轴的粗细
    #ax1.spines['right'].set_color('none')# 将右边上边的两条边颜色设置为空 其实就相当于抹掉这两条边
    #ax1.spines['top'].set_color('none')
    #====坐标轴的第三层：  主刻度线的设置
    for line in ax1.xaxis.get_ticklines():
        #line is a Line2D instance
        #line.set_color('green')
        line.set_markersize(3)####设置刻度线的长度
        line.set_markeredgewidth(1.5)####设置刻度线的宽度
    for line in ax1.yaxis.get_ticklines():
        # line is a Line2D instance
        #line.set_color('green')
        line.set_markersize(3)####设置刻度线的长度
        line.set_markeredgewidth(1.5)####设置刻度线的宽度
    plt.tick_params(axis="x", which="minor", length=2.5, width=1.5, color="k")  ### 设置次要刻度 
    plt.tick_params(axis="y", which="minor", length=2.5, width=1.5, color="k")  ### 设置次要刻度
    #------------------------------------------------------------------------------------------
    plt.show()