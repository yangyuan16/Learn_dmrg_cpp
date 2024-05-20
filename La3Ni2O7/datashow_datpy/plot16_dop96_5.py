# 画出 YY Pairing Correlation
#
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
#
def get_data_pyy(Lz, Ly, Lx, dop, t, J, Jz, dim,):
    workpath = "E:\\WORK\\Work\\Project\\La3Ni2O7"
    filepath1 = "\\data_dmrgcpp\\Lz%d_Ly%d_Lx%d\\dop%g" % (Lz, Ly, Lx,dop)
    filepath2 = "\\t%d_J%d_Jz%.2f_dim%d" % (t, J, Jz, dim)
    filepath3 = "\\measurement_pairing_yy_fit.parquet"
    filename = workpath + filepath1 + filepath2 + filepath3
    print(filename)
    df = pd.read_parquet(filename)
    print(df.head())
    return df
#
#
if __name__ =="__main__":
    Lz = 2
    Ly = 2
    Lx = 48
    dop = 96
    t = 3
    J = 1
    dim = 6000
    # load spin correlation data
    # load spin correlation data
    df_pyy_Jz025 = get_data_pyy(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=0.25,dim=dim)
    df_pyy_Jz05 = get_data_pyy(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=0.5,dim=dim)
    df_pyy_Jz075 = get_data_pyy(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=0.75,dim=dim)
    df_pyy_Jz10 = get_data_pyy(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=1.0,dim=dim)
    df_pyy_Jz125 = get_data_pyy(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=1.25,dim=dim)
    df_pyy_Jz15 = get_data_pyy(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=1.5,dim=dim)
    df_pyy_Jz175 = get_data_pyy(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=1.75,dim=dim)
    # load spin correlation data
    df_pyy_Jz20 = get_data_pyy(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=2.0,dim=dim)
    df_pyy_Jz225 = get_data_pyy(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=2.25,dim=dim)
    df_pyy_Jz25 = get_data_pyy(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=2.5,dim=dim)
    df_pyy_Jz275 = get_data_pyy(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=2.75,dim=dim)
    df_pyy_Jz30 = get_data_pyy(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=3.0,dim=dim)
    df_pyy_Jz325 = get_data_pyy(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=3.25,dim=dim)
    df_pyy_Jz35 = get_data_pyy(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=3.5,dim=dim)
    #
    #
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
    ax1 = plt.axes([0.1,0.58,0.39,0.4])
    ax2 = plt.axes([0.58,0.58,0.39,0.4])
    ax3 = plt.axes([0.1,0.1,0.39,0.4])
    ax4 = plt.axes([0.58,0.1,0.39,0.4])
    # plt.subplots_adjust(left=0.1, bottom=0.1, right=0.98, top=0.96, wspace=0.32, hspace=0.26)
    # 选择子图 ax1 进行绘图
    plt.sca(ax1) ## 选择对 ax1 进行绘图
    ax1 = plt.gca()
    #--------- J_\bot = 0.25
    slope = df_pyy_Jz025["slope_exp"].values[0]
    xi = round(-1/slope,2)  
    label = r"$J_{\bot}$=%.2f $\xi_{SC}^{yy}$=%.2f"%(0.25,xi)
    L025, = ax1.plot(df_pyy_Jz025["r"],df_pyy_Jz025["corre_abs"],label=label,ls="-",lw=2,color="r",
             marker='o',alpha=1,markersize=6,markeredgewidth=1, markeredgecolor="r", markerfacecolor='w')
    L025_, = ax1.plot(df_pyy_Jz025["r"],df_pyy_Jz025["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
             marker='o',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="r", markerfacecolor='w')
    #--------- J_\bot = 0.50
    slope = df_pyy_Jz05["slope_exp"].values[0]
    xi = round(-1/slope,2)  
    label = r"$J_{\bot}$=%.2f $\xi_{SC}^{yy}$=%.2f"%(0.5,xi)
    L05, = ax1.plot(df_pyy_Jz05["r"],df_pyy_Jz05["corre_abs"],label=label,ls="-",lw=2,color="blue",
             marker='s',alpha=1,markersize=6,markeredgewidth=1, markeredgecolor="blue", markerfacecolor='w')
    L05_, = ax1.plot(df_pyy_Jz05["r"],df_pyy_Jz05["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
             marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- J_\bot = 0.75
    slope = df_pyy_Jz075["slope_exp"].values[0]
    xi = round(-1/slope,2)  
    label = r"$J_{\bot}$=%.2f $\xi_{SC}^{yy}$=%.2f"%(0.75,xi)
    L075, = ax1.plot(df_pyy_Jz075["r"],df_pyy_Jz075["corre_abs"],label=label,ls="-",lw=2,color="green",
             marker='^',alpha=1,markersize=6,markeredgewidth=1, markeredgecolor="green", markerfacecolor='w')
    L075_, = ax1.plot(df_pyy_Jz075["r"],df_pyy_Jz075["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
             marker='^',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- J_\bot = 1.0
    slope = df_pyy_Jz10["slope_exp"].values[0]
    xi = round(-1/slope,2)  
    label = r"$J_{\bot}$=%.2f $\xi_{SC}^{yy}$=%.2f"%(1.0,xi)
    L10, = ax1.plot(df_pyy_Jz10["r"],df_pyy_Jz10["corre_abs"],label=label,ls="-",lw=2,color="magenta",
             marker='v',alpha=1,markersize=6,markeredgewidth=1, markeredgecolor="magenta", markerfacecolor='w')
    L10_, = ax1.plot(df_pyy_Jz10["r"],df_pyy_Jz10["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
             marker='v',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 1.25
    slope = df_pyy_Jz125["slope_exp"].values[0]
    xi = round(-1/slope,2)  
    label = r"$J_{\bot}$=%.2f $\xi_{SC}^{yy}$=%.2f"%(1.25,xi)
    L125, = ax1.plot(df_pyy_Jz125["r"],df_pyy_Jz125["corre_abs"],label=label,ls="-",lw=2,color="cyan",
             marker='+',alpha=1,markersize=6,markeredgewidth=1, markeredgecolor="cyan", markerfacecolor='w')
    L125_, = ax1.plot(df_pyy_Jz125["r"],df_pyy_Jz125["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
             marker='+',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 1.5
    slope = df_pyy_Jz15["slope_exp"].values[0]
    xi = round(-1/slope,2)  
    label = r"$J_{\bot}$=%.2f $\xi_{SC}^{yy}$=%.2f"%(1.5,xi)
    L15, = ax1.plot(df_pyy_Jz15["r"],df_pyy_Jz15["corre_abs"],label=label,ls="-",lw=2,color="pink",
             marker='*',alpha=1,markersize=6,markeredgewidth=1, markeredgecolor="pink", markerfacecolor='w')
    L15_, = ax1.plot(df_pyy_Jz15["r"],df_pyy_Jz15["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
             marker='*',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 1.75
    slope = df_pyy_Jz175["slope_exp"].values[0]
    xi = round(-1/slope,2)  
    label = r"$J_{\bot}$=%.2f $\xi_{SC}^{yy}$=%.2f"%(1.75,xi)
    L175, = ax1.plot(df_pyy_Jz175["r"],df_pyy_Jz175["corre_abs"],label=label,ls="-",lw=2,color="brown",
             marker='h',alpha=1,markersize=6,markeredgewidth=1, markeredgecolor="brown", markerfacecolor='w')
    L175_, = ax1.plot(df_pyy_Jz175["r"],df_pyy_Jz175["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
             marker='h',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')

    ####图例设置
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 14, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[L025,L05,L075,L10,L125,L15,L175], loc = 4, bbox_to_anchor=(0.38, 0.02),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    
    label_x = r"r"
    label_y = "$\Phi^{yy}$"
    #plt.yscale("log")
    #plt.xscale("log")
    ax1.set_yscale("log",base=10,subs=[0.01,0.02,0.03])      
    ax1.set_xlabel(label_x, size= 16)
    ax1.set_ylabel(label_y, size= 16)
    ax1.tick_params(labelsize = 14) # 设置坐标刻度对应数字的大小
    ax1.set_xlim([0,40])
    #ax1.set_xticks([0,10,20,30])
    #ax1.set_yticks([-1,-0.5,0,0.5,1])  
    #=========================================================
    ax1.text(36,0.003,"(a)",fontsize = 20, color='black', rotation = 0)
    ax1.text(17.5,0.003, r'$\mathrm{\delta} = \frac{1}{2}$', fontsize = 20, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
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
    #-----------------------------------------------------------------------------
    # 选择子图 ax2 进行绘图
    plt.sca(ax2) ## 选择对 ax1 进行绘图
    ax2 = plt.gca()
    #--------- J_\bot = 2.0
    slope = df_pyy_Jz20["slope_exp"].values[0]
    xi = round(-1/slope,2)  
    label = r"$J_{\bot}$=%.2f $\xi_{SC}^{yy}$=%.2f"%(2.0,xi)
    L20, = ax2.plot(df_pyy_Jz20["r"],df_pyy_Jz20["corre_abs"],label=label,ls="-",lw=2,color="r",
             marker='o',alpha=1,markersize=6,markeredgewidth=1, markeredgecolor="r", markerfacecolor='w')
    L20_, = ax2.plot(df_pyy_Jz20["r"],df_pyy_Jz20["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
             marker='o',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- J_\bot = 2.25
    slope = df_pyy_Jz225["slope_exp"].values[0]
    xi = round(-1/slope,2)  
    label = r"$J_{\bot}$=%.2f $\xi_{SC}^{yy}$=%.2f"%(2.25,xi)
    L225, = ax2.plot(df_pyy_Jz225["r"],df_pyy_Jz225["corre_abs"],label=label,ls="-",lw=2,color="blue",
             marker='s',alpha=1,markersize=6,markeredgewidth=1, markeredgecolor="blue", markerfacecolor='w')
    L225_, = ax2.plot(df_pyy_Jz225["r"],df_pyy_Jz225["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
             marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- J_\bot = 2.5
    slope = df_pyy_Jz25["slope_exp"].values[0]
    xi = round(-1/slope,2)  
    label = r"$J_{\bot}$=%.2f $\xi_{SC}^{yy}$=%.2f"%(2.5,xi)
    L25, = ax2.plot(df_pyy_Jz25["r"],df_pyy_Jz25["corre_abs"],label=label,ls="-",lw=2,color="green",
             marker='^',alpha=1,markersize=6,markeredgewidth=1, markeredgecolor="green", markerfacecolor='w')
    L25_, = ax2.plot(df_pyy_Jz25["r"],df_pyy_Jz25["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
             marker='^',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- J_\bot = 2.75
    slope = df_pyy_Jz275["slope_exp"].values[0]
    xi = round(-1/slope,2)  
    label = r"$J_{\bot}$=%.2f $\xi_{SC}^{yy}$=%.2f"%(2.75,xi)
    L275, = ax2.plot(df_pyy_Jz275["r"],df_pyy_Jz275["corre_abs"],label=label,ls="-",lw=2,color="magenta",
             marker='v',alpha=1,markersize=6,markeredgewidth=1, markeredgecolor="magenta", markerfacecolor='w')
    L275_, = ax2.plot(df_pyy_Jz275["r"],df_pyy_Jz275["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
             marker='v',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 3.0
    slope = df_pyy_Jz30["slope_exp"].values[0]
    xi = round(-1/slope,2)  
    label = r"$J_{\bot}$=%.2f $\xi_{SC}^{yy}$=%.2f"%(3.0,xi)
    L30, = ax2.plot(df_pyy_Jz30["r"],df_pyy_Jz30["corre_abs"],label=label,ls="-",lw=2,color="cyan",
             marker='+',alpha=1,markersize=6,markeredgewidth=1, markeredgecolor="cyan", markerfacecolor='w')
    L30_, = ax2.plot(df_pyy_Jz30["r"],df_pyy_Jz30["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
             marker='+',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 3.25
    slope = df_pyy_Jz325["slope_exp"].values[0]
    xi = round(-1/slope,2)  
    label = r"$J_{\bot}$=%.2f $\xi_{SC}^{yy}$=%.2f"%(3.25,xi)
    L325, = ax2.plot(df_pyy_Jz325["r"],df_pyy_Jz325["corre_abs"],label=label,ls="-",lw=2,color="pink",
             marker='*',alpha=1,markersize=6,markeredgewidth=1, markeredgecolor="pink", markerfacecolor='w')
    L325_, = ax2.plot(df_pyy_Jz325["r"],df_pyy_Jz325["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
             marker='*',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 3.5
    slope = df_pyy_Jz35["slope_exp"].values[0]
    xi = round(-1/slope,2)  
    label = r"$J_{\bot}$=%.2f $\xi_{SC}^{yy}$=%.2f"%(3.5,xi)
    L35, = ax2.plot(df_pyy_Jz35["r"],df_pyy_Jz35["corre_abs"],label=label,ls="-",lw=2,color="brown",
             marker='h',alpha=1,markersize=6,markeredgewidth=1, markeredgecolor="brown", markerfacecolor='w')
    L35_, = ax2.plot(df_pyy_Jz35["r"],df_pyy_Jz35["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
             marker='h',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')

    ####图例设置
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 14, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[L20,L225,L25,L275,L30,L325,L35], loc = 4, bbox_to_anchor=(0.38, 0.02),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    
    label_x = r"r"
    label_y = "$\Phi^{yy}$"
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
    ax2.text(36,0.003,"(b)",fontsize = 20, color='black', rotation = 0)
    ax2.text(17.5,0.003, r'$\mathrm{\delta} = \frac{1}{2}$', fontsize = 20, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
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
    #=====坐标轴的第二层： 坐标轴的设置
    ax2.spines['bottom'].set_linewidth(1.5) ###设置底部坐标轴的粗细
    ax2.spines['left'].set_linewidth(1.5)   ###设置左边坐标轴的粗细
    ax2.spines['right'].set_linewidth(1.5)  ###设置右边坐标轴的粗细
    ax2.spines['top'].set_linewidth(1.5)    ###设置上部坐标轴的粗细
    #ax2.spines['right'].set_color('none')# 将右边上边的两条边颜色设置为空 其实就相当于抹掉这两条边
    #ax2.spines['top'].set_color('none')
    #====坐标轴的第三层：  主刻度线的设置
    for line in ax2.xaxis.get_ticklines():
        #line is a Line2D instance
        #line.set_color('green')
        line.set_markersize(3)####设置刻度线的长度
        line.set_markeredgewidth(1.5)####设置刻度线的宽度
    for line in ax2.yaxis.get_ticklines():
        # line is a Line2D instance
        #line.set_color('green')
        line.set_markersize(3)####设置刻度线的长度
        line.set_markeredgewidth(1.5)####设置刻度线的宽度
    plt.tick_params(axis="x", which="minor", length=2.5, width=1.5, color="k")  ### 设置次要刻度 
    plt.tick_params(axis="y", which="minor", length=2.5, width=1.5, color="k")  ### 设置次要刻度  

    #----------------------------------------------------------------------------
    # 选择子图 ax3 进行绘图
    plt.sca(ax3) ## 选择对 ax1 进行绘图
    ax3 = plt.gca()
    #--------- J_\bot = 0.25
    slope = df_pyy_Jz025["slope_pow"].values[0]
    K_F = round(-slope,2)  
    label = r"$J_{\bot}$=%.2f $K_{SC}^{yy}$=%.2f"%(0.25,K_F)
    L025, = ax3.plot(df_pyy_Jz025["r"],df_pyy_Jz025["corre_abs"],label=label,ls="-",lw=2,color="r",
             marker='o',alpha=1,markersize=6,markeredgewidth=1, markeredgecolor="r", markerfacecolor='w')
    L025_, = ax3.plot(df_pyy_Jz025["r"],df_pyy_Jz025["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='o',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- J_\bot = 0.50
    slope = df_pyy_Jz05["slope_pow"].values[0]
    K_F = round(-slope,2)  
    label = r"$J_{\bot}$=%.2f $K_{SC}^{yy}$=%.2f"%(0.5,K_F)
    L05, = ax3.plot(df_pyy_Jz05["r"],df_pyy_Jz05["corre_abs"],label=label,ls="-",lw=2,color="blue",
             marker='s',alpha=1,markersize=6,markeredgewidth=1, markeredgecolor="blue", markerfacecolor='w')
    L05_, = ax3.plot(df_pyy_Jz05["r"],df_pyy_Jz05["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- J_\bot = 0.75
    slope = df_pyy_Jz075["slope_pow"].values[0]
    K_F = round(-slope,2)  
    label = r"$J_{\bot}$=%.2f $K_{SC}^{yy}$=%.2f"%(0.75,K_F)
    L075, = ax3.plot(df_pyy_Jz075["r"],df_pyy_Jz075["corre_abs"],label=label,ls="-",lw=2,color="green",
             marker='^',alpha=1,markersize=6,markeredgewidth=1, markeredgecolor="green", markerfacecolor='w')
    L075_, = ax3.plot(df_pyy_Jz075["r"],df_pyy_Jz075["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='^',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- J_\bot = 1.0
    slope = df_pyy_Jz10["slope_pow"].values[0]
    K_F = round(-slope,2)  
    label = r"$J_{\bot}$=%.2f $K_{SC}^{yy}$=%.2f"%(1.0,K_F)
    L10, = ax3.plot(df_pyy_Jz10["r"],df_pyy_Jz10["corre_abs"],label=label,ls="-",lw=2,color="magenta",
             marker='v',alpha=1,markersize=6,markeredgewidth=1, markeredgecolor="magenta", markerfacecolor='w')
    L10_, = ax3.plot(df_pyy_Jz10["r"],df_pyy_Jz10["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='v',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 1.25
    slope = df_pyy_Jz125["slope_pow"].values[0]
    K_F = round(-slope,2)  
    label = r"$J_{\bot}$=%.2f $K_{SC}^{yy}$=%.2f"%(1.25,K_F)
    L125, = ax3.plot(df_pyy_Jz125["r"],df_pyy_Jz125["corre_abs"],label=label,ls="-",lw=2,color="cyan",
             marker='+',alpha=1,markersize=6,markeredgewidth=1, markeredgecolor="cyan", markerfacecolor='w')
    L125_, = ax3.plot(df_pyy_Jz125["r"],df_pyy_Jz125["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='+',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 1.5
    slope = df_pyy_Jz15["slope_pow"].values[0]
    K_F = round(-slope,2)  
    label = r"$J_{\bot}$=%.2f $K_{SC}^{yy}$=%.2f"%(1.5,K_F)
    L15, = ax3.plot(df_pyy_Jz15["r"],df_pyy_Jz15["corre_abs"],label=label,ls="-",lw=2,color="pink",
             marker='*',alpha=1,markersize=6,markeredgewidth=1, markeredgecolor="pink", markerfacecolor='w')
    L15_, = ax3.plot(df_pyy_Jz15["r"],df_pyy_Jz15["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='*',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 1.75
    slope = df_pyy_Jz175["slope_pow"].values[0]
    K_F = round(-slope,2)  
    label = r"$J_{\bot}$=%.2f $K_{SC}^{yy}$=%.2f"%(1.75,K_F)
    L175, = ax3.plot(df_pyy_Jz175["r"],df_pyy_Jz175["corre_abs"],label=label,ls="-",lw=2,color="brown",
             marker='h',alpha=1,markersize=6,markeredgewidth=1, markeredgecolor="brown", markerfacecolor='w')
    L175_, = ax3.plot(df_pyy_Jz175["r"],df_pyy_Jz175["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='h',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')

    ####图例设置
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 14, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[L025,L05,L075,L10,L125,L15,L175], loc = 4, bbox_to_anchor=(0.38, 0.02),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    
    label_x = r"r"
    label_y = "$\Phi^{yy}$"
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
    ax3.text(26,0.003,"(c)",fontsize = 20, color='black', rotation = 0)
    ax3.text(9.5,0.003, r'$\mathrm{\delta} = \frac{1}{2}$', fontsize = 20, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
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
    #=====坐标轴的第二层： 坐标轴的设置
    ax3.spines['bottom'].set_linewidth(1.5) ###设置底部坐标轴的粗细
    ax3.spines['left'].set_linewidth(1.5)   ###设置左边坐标轴的粗细
    ax3.spines['right'].set_linewidth(1.5)  ###设置右边坐标轴的粗细
    ax3.spines['top'].set_linewidth(1.5)    ###设置上部坐标轴的粗细
    #ax3.spines['right'].set_color('none')# 将右边上边的两条边颜色设置为空 其实就相当于抹掉这两条边
    #ax3.spines['top'].set_color('none')
    #====坐标轴的第三层：  主刻度线的设置
    for line in ax3.xaxis.get_ticklines():
        #line is a Line2D instance
        #line.set_color('green')
        line.set_markersize(3)####设置刻度线的长度
        line.set_markeredgewidth(1.5)####设置刻度线的宽度
    for line in ax3.yaxis.get_ticklines():
        # line is a Line2D instance
        #line.set_color('green')
        line.set_markersize(3)####设置刻度线的长度
        line.set_markeredgewidth(1.5)####设置刻度线的宽度
    plt.tick_params(axis="x", which="minor", length=2.5, width=1.5, color="k")  ### 设置次要刻度 
    plt.tick_params(axis="y", which="minor", length=2.5, width=1.5, color="k")  ### 设置次要刻度 
    #--------------------------------------------------------------------------------------------
    # 选择子图 ax4 进行绘图
    plt.sca(ax4) ## 选择对 ax4 进行绘图
    ax4 = plt.gca()
    #--------- J_\bot = 2.0
    slope = df_pyy_Jz20["slope_pow"].values[0]
    K_F = round(-slope,2)  
    label = r"$J_{\bot}$=%.2f $K_{SC}^{yy}$=%.2f"%(2.0,K_F)
    L20, = ax4.plot(df_pyy_Jz20["r"],df_pyy_Jz20["corre_abs"],label=label,ls="-",lw=2,color="r",
             marker='o',alpha=1,markersize=6,markeredgewidth=1, markeredgecolor="r", markerfacecolor='w')
    L20_, = ax4.plot(df_pyy_Jz20["r"],df_pyy_Jz20["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='o',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- J_\bot = 2.25
    slope = df_pyy_Jz225["slope_pow"].values[0]
    K_F = round(-slope,2)  
    label = r"$J_{\bot}$=%.2f $K_{SC}^{yy}$=%.2f"%(2.25,K_F)
    L225, = ax4.plot(df_pyy_Jz225["r"],df_pyy_Jz225["corre_abs"],label=label,ls="-",lw=2,color="blue",
             marker='s',alpha=1,markersize=6,markeredgewidth=1, markeredgecolor="blue", markerfacecolor='w')
    L225_, = ax4.plot(df_pyy_Jz225["r"],df_pyy_Jz225["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- J_\bot = 2.5
    slope = df_pyy_Jz25["slope_pow"].values[0]
    K_F = round(-slope,2)  
    label = r"$J_{\bot}$=%.2f $K_{SC}^{yy}$=%.2f"%(2.5,K_F)
    L25, = ax4.plot(df_pyy_Jz25["r"],df_pyy_Jz25["corre_abs"],label=label,ls="-",lw=2,color="green",
             marker='^',alpha=1,markersize=6,markeredgewidth=1, markeredgecolor="green", markerfacecolor='w')
    L25_, = ax4.plot(df_pyy_Jz25["r"],df_pyy_Jz25["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='^',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- J_\bot = 2.75
    slope = df_pyy_Jz275["slope_pow"].values[0]
    K_F = round(-slope,2)  
    label = r"$J_{\bot}$=%.2f $K_{SC}^{yy}$=%.2f"%(2.75,K_F)
    L275, = ax4.plot(df_pyy_Jz275["r"],df_pyy_Jz275["corre_abs"],label=label,ls="-",lw=2,color="magenta",
             marker='v',alpha=1,markersize=6,markeredgewidth=1, markeredgecolor="magenta", markerfacecolor='w')
    L275_, = ax4.plot(df_pyy_Jz275["r"],df_pyy_Jz275["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='v',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 3.0
    slope = df_pyy_Jz30["slope_pow"].values[0]
    K_F = round(-slope,2)  
    label = r"$J_{\bot}$=%.2f $K_{SC}^{yy}$=%.2f"%(3.0,K_F)
    L30, = ax4.plot(df_pyy_Jz30["r"],df_pyy_Jz30["corre_abs"],label=label,ls="-",lw=2,color="cyan",
             marker='+',alpha=1,markersize=6,markeredgewidth=1, markeredgecolor="cyan", markerfacecolor='w')
    L30_, = ax4.plot(df_pyy_Jz30["r"],df_pyy_Jz30["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='+',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 3.25
    slope = df_pyy_Jz325["slope_pow"].values[0]
    K_F = round(-slope,2)  
    label = r"$J_{\bot}$=%.2f $K_{SC}^{yy}$=%.2f"%(3.25,K_F)
    L325, = ax4.plot(df_pyy_Jz325["r"],df_pyy_Jz325["corre_abs"],label=label,ls="-",lw=2,color="pink",
             marker='*',alpha=1,markersize=6,markeredgewidth=1, markeredgecolor="pink", markerfacecolor='w')
    L325_, = ax4.plot(df_pyy_Jz325["r"],df_pyy_Jz325["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='*',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 3.5
    slope = df_pyy_Jz35["slope_pow"].values[0]
    K_F = round(-slope,2)  
    label = r"$J_{\bot}$=%.2f $K_{SC}^{yy}$=%.2f"%(3.5,K_F)
    L35, = ax4.plot(df_pyy_Jz35["r"],df_pyy_Jz35["corre_abs"],label=label,ls="-",lw=2,color="brown",
             marker='h',alpha=1,markersize=6,markeredgewidth=1, markeredgecolor="brown", markerfacecolor='w')
    L35_, = ax4.plot(df_pyy_Jz35["r"],df_pyy_Jz35["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='h',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')

    ####图例设置
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 14, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[L20,L225,L25,L275,L30,L325,L35], loc = 4, bbox_to_anchor=(0.38, 0.02),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    
    label_x = r"r"
    label_y = "$\Phi^{yy}$"
    #plt.yscale("log")
    #plt.xscale("log")
    ax4.set_xscale("log",base=10,subs=[0.01,0.02,0.03]) 
    ax4.set_yscale("log",base=10,subs=[0.01,0.02,0.03])      
    ax4.set_xlabel(label_x, size= 16)
    ax4.set_ylabel(label_y, size= 16)
    ax4.tick_params(labelsize = 14) # 设置坐标刻度对应数字的大小
    ax4.set_xlim([0,40])
    #ax4.set_xticks([0,10,20,30])
    #ax4.set_yticks([-1,-0.5,0,0.5,1])  
    #=========================================================
    ax4.text(26,0.003,"(d)",fontsize = 20, color='black', rotation = 0)
    ax4.text(9.5,0.003, r'$\mathrm{\delta} = \frac{1}{2}$', fontsize = 20, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
    #=========================================================
    # 坐标轴设置第一层
    labels = ax4.get_xticklabels() + ax4.get_yticklabels()
    #[label.set_fontname('Times New Roman') for label in labels]###设置ticket labled的字体格式
    ax4.xaxis.set_minor_locator(MultipleLocator(5))###设置次刻度的间隔
    #ax4.yaxis.set_minor_locator(MultipleLocator(10))###设置次刻度的间隔
    ax4.xaxis.set_major_formatter(FormatStrFormatter('%1.0f'))###设置X轴标签文本格式
    #ax4.yaxis.set_major_formatter(FormatStrFormatter('%1.1f'))###设置Y轴标签文本格式
    #
    #ax4.yaxis.get_major_locator().set_params(numticks=99)
    #ax4.yaxis.get_minor_locator().set_params(numticks=99, subs=[.2,.4,.6,.8]) # 将次要刻度显示出来 
    #----将次要刻度显示出来 
    #locmin = matplotlib.ticker.LogLocator(base=10.0, subs=(0.1,0.2,0.4,0.6,0.8,1,2,4,6,8,10 )) 
    #ax4.xaxis.set_minor_locator(locmin)
    #ax4.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    #=====坐标轴的第二层： 坐标轴的设置
    ax4.spines['bottom'].set_linewidth(1.5) ###设置底部坐标轴的粗细
    ax4.spines['left'].set_linewidth(1.5)   ###设置左边坐标轴的粗细
    ax3.spines['right'].set_linewidth(1.5)  ###设置右边坐标轴的粗细
    ax4.spines['top'].set_linewidth(1.5)    ###设置上部坐标轴的粗细
    #ax4.spines['right'].set_color('none')# 将右边上边的两条边颜色设置为空 其实就相当于抹掉这两条边
    #ax4.spines['top'].set_color('none')
    #====坐标轴的第三层：  主刻度线的设置
    for line in ax4.xaxis.get_ticklines():
        #line is a Line2D instance
        #line.set_color('green')
        line.set_markersize(3)####设置刻度线的长度
        line.set_markeredgewidth(1.5)####设置刻度线的宽度
    for line in ax4.yaxis.get_ticklines():
        # line is a Line2D instance
        #line.set_color('green')
        line.set_markersize(3)####设置刻度线的长度
        line.set_markeredgewidth(1.5)####设置刻度线的宽度
    plt.tick_params(axis="x", which="minor", length=2.5, width=1.5, color="k")  ### 设置次要刻度 
    plt.tick_params(axis="y", which="minor", length=2.5, width=1.5, color="k")  ### 设置次要刻度
    #------------------------------------------------------------------------------------------
    plt.show()