# 将 YY pairing, ZZ pairing, Density correlation 画在一起
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

#
def get_data_cicj(Lz, Ly, Lx, dop, t, J, Jz, dim,):
    # load data
    workpath = "E:\\WORK\\Work\\Project\\La3Ni2O7"
    filepath1 = "\\data_dmrgcpp\\Lz%d_Ly%d_Lx%d\\dop%g" % (Lz, Ly, Lx,dop)
    filepath2 = "\\t%d_J%d_Jz%.2f_dim%d" % (t, J, Jz, dim)
    filepath3 = "\\measurement_green_function_fit.parquet"
    filename = workpath + filepath1 + filepath2 + filepath3
    print(filename)
    df = pd.read_parquet(filename)
    print(df.head())
    return df
#
def get_data_sisj(Lz, Ly, Lx, dop, t, J, Jz, dim,):
    workpath = "E:\\WORK\\Work\\Project\\La3Ni2O7"
    filepath1 = "\\data_dmrgcpp\\Lz%d_Ly%d_Lx%d\\dop%g" % (Lz, Ly, Lx,dop)
    filepath2 = "\\t%d_J%d_Jz%.2f_dim%d" % (t, J, Jz, dim)
    filepath3 = "\\measurement_spin_correlation_fit.parquet"
    filename = workpath + filepath1 + filepath2 + filepath3
    print(filename)
    df = pd.read_parquet(filename)
    print(df.head())
    return df
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
#
if __name__ =="__main__":
    Lz = 2
    Ly = 2
    Lx = 48
    dop = 88
    t = 3
    J = 1
    dim = 6000

    # load SiSj correlation data
    df_sisj_Jz05 = get_data_sisj(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=0.5,dim=dim)
    df_sisj_Jz15 = get_data_sisj(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=1.5,dim=dim)
    df_sisj_Jz25 = get_data_sisj(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=2.5,dim=dim)
    # load CiCj correlation data
    df_cicj_Jz05 = get_data_cicj(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=0.5,dim=dim)
    df_cicj_Jz15 = get_data_cicj(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=1.5,dim=dim)
    df_cicj_Jz25 = get_data_cicj(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=2.5,dim=dim)    
    # load YY pairing correlation data
    df_pyy_Jz05 = get_data_pyy(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=0.5,dim=dim)
    df_pyy_Jz15 = get_data_pyy(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=1.5,dim=dim)
    df_pyy_Jz25 = get_data_pyy(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=2.5,dim=dim)
    # load ZZ pairing correlation data
    df_pzz_Jz05 = get_data_pzz(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=0.5,dim=dim)
    df_pzz_Jz15 = get_data_pzz(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=1.5,dim=dim)
    df_pzz_Jz25 = get_data_pzz(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=2.5,dim=dim)
    # ninj correlation
    df_ninj_Jz05 = get_data_ninj_corre(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=0.5,dim=dim)
    df_ninj_Jz15 = get_data_ninj_corre(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=1.5,dim=dim)
    df_ninj_Jz25 = get_data_ninj_corre(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=2.5,dim=dim)
    #
    Jz_list = [0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.25,2.5,2.75,3.0,3.25,3.5]
    K_f = []
    K_g = []
    K_sc_yy = []
    K_sc_zz = []
    K_cdw = []
    for it in Jz_list:
        df_pyy = get_data_pyy(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=it,dim=dim)
        K_sc_yy.append(-df_pyy["slope_pow"].values[0])
        df_pzz = get_data_pzz(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=it,dim=dim)
        K_sc_zz.append(-df_pzz["slope_pow"].values[0])
        df_ninj = get_data_ninj_corre(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=it,dim=dim)
        K_cdw.append(-df_ninj["slope_pow"].values[0])
        df_sisj = get_data_sisj(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=it,dim=dim)
        K_f.append(-df_sisj["slope_pow"].values[0])
        df_cicj = get_data_cicj(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=it,dim=dim)
        K_g.append(-df_cicj["slope_pow"].values[0])
    #
    print("K_f:\n", K_f)
    print("K_g:\n", K_g)
    print("K_sc_yy:\n", K_sc_yy)
    print("K_sc_zz:\n",K_sc_zz)
    print("K_cdw:\n",K_cdw)
    #
    #---------------plot logr-r fig and logr-logr fig-----------------------------
    fig = plt.figure(figsize=(10,10))
    # plt.figure(facecolor='blue',edgecolor='black') # 设置画布的颜色
    params = {
        'axes.labelsize': '30',
        'xtick.labelsize':'20',
        'ytick.labelsize':'20',
        'ytick.direction':'in',
        'xtick.direction':'in',
        'lines.linewidth':10 ,
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
    #--------- J_\bot = 0.5
    slope = df_pyy_Jz05["slope_pow"].values[0]
    K_F = round(-slope,2)  
    label = r"$K_{SC}^{yy}$=%.2f"%(K_F)
    Lyy05, = ax1.plot(df_pyy_Jz05["r"],df_pyy_Jz05["corre_abs"],label=label,ls="-",lw=2,color="red",
             marker='s',alpha=1,markersize=7,markeredgewidth=1.5, markeredgecolor="red", markerfacecolor='w')
    Lyy05_, = ax1.plot(df_pyy_Jz05["r"],df_pyy_Jz05["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- J_\bot = 0.5
    slope = df_pzz_Jz05["slope_pow"].values[0]
    K_G = round(-slope,2)  
    label = r"$K_{SC}^{zz}$=%.2f"%(K_G)
    Lzz05, = ax1.plot(df_pzz_Jz05["r"],df_pzz_Jz05["corre_abs"],label=label,ls="-",lw=2,color="blue",
             marker='s',alpha=1,markersize=7,markeredgewidth=1.5, markeredgecolor="blue", markerfacecolor='w')
    Lzz05_, = ax1.plot(df_pzz_Jz05["r"],df_pzz_Jz05["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 0.5
    slope = df_ninj_Jz05["slope_pow"].values[0]
    K_F = round(-slope,2)  
    label = r"$K_{CDW}$=%.2f"%(K_F)
    Lninj05, = ax1.plot(df_ninj_Jz05["r"][1:],df_ninj_Jz05["corre_abs"][1:],label=label,ls="-",lw=2,color="green",
             marker='s',alpha=1,markersize=7,markeredgewidth=1.5, markeredgecolor="green", markerfacecolor='w')
    Lninj05_, = ax1.plot(df_ninj_Jz05["r"][1:],df_ninj_Jz05["fitcorre_pow"][1:],label=label,ls="--",lw=1.5,color="k",
             marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 0.5
    slope = df_sisj_Jz05["slope_pow"].values[0]
    K_F = round(-slope,2)  
    label = r"$K_{F}$=%.2f"%(K_F)
    Lsisj05, = ax1.plot(df_sisj_Jz05["r"][1:],df_sisj_Jz05["corre_abs"][1:],label=label,ls="-",lw=2,color="magenta",
             marker='o',alpha=1,markersize=7,markeredgewidth=1.5, markeredgecolor="magenta", markerfacecolor='w')
    Lsisj05_, = ax1.plot(df_sisj_Jz05["r"][1:],df_sisj_Jz05["fitcorre_pow"][1:],label=label,ls="--",lw=1.5,color="k",
             marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')    
    #-------- J_\bot = 0.5
    slope = df_cicj_Jz05["slope_pow"].values[0]
    K_F = round(-slope,2)  
    label = r"$K_{G}$=%.2f"%(K_F)
    Lcicj05, = ax1.plot(df_cicj_Jz05["r"][1:],df_cicj_Jz05["corre_abs"][1:],label=label,ls="-",lw=2,color="cyan",
             marker='o',alpha=1,markersize=7,markeredgewidth=1.5, markeredgecolor="cyan", markerfacecolor='w')
    Lcicj05_, = ax1.plot(df_cicj_Jz05["r"][1:],df_cicj_Jz05["fitcorre_pow"][1:],label=label,ls="--",lw=1.5,color="k",
             marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w') 
    ####图例设置
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 16, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[Lyy05,Lzz05,Lninj05,Lsisj05,Lcicj05], loc = 4, bbox_to_anchor=(0.55, 0.02),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    #legend2=plt.legend(handles=[L10,L20,L40], loc = 4, bbox_to_anchor=(0.98, 0.70),
    #                   ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    #plt.gca().add_artist(legend1)
    label_x = r"r"
    label_y = "Correlation"
    #plt.yscale("log")
    #plt.xscale("log")
    ax1.set_xscale("log",base=10,subs=[0.01,0.02,0.03]) 
    ax1.set_yscale("log",base=10,subs=[0.01,0.02,0.03])      
    ax1.set_xlabel(label_x, size= 16)
    ax1.set_ylabel(label_y, size= 16)
    ax1.tick_params(labelsize = 14) # 设置坐标刻度对应数字的大小
    ax1.set_xlim([0,40])
    #ax1.set_xticks([0,10,20,30])
    #ax1.set_yticks([-1,-0.5,0,0.5,1]) 
    #
    #=========================================================
    ax1.text(28,0.01,"(a)",fontsize = 20, color='black', rotation = 0)
    ax1.text(5.5,5e-2, r'$\mathrm{\delta} = \frac{11}{24}$', fontsize = 20, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
    ax1.text(15,5e-2, r'$J_{\bot}=0.5$', fontsize = 16, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
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
    #-----------------------------------------------------------------------------------------------
    # 选择子图 ax2 进行绘图
    plt.sca(ax2) ## 选择对 ax1 进行绘图
    ax2 = plt.gca()
    #--------- J_\bot = 1.5
    slope = df_pyy_Jz15["slope_pow"].values[0]
    K_F = round(-slope,2)  
    label = r"$K_{SC}^{yy}$=%.2f"%(K_F)
    Lyy15, = ax2.plot(df_pyy_Jz15["r"],df_pyy_Jz15["corre_abs"],label=label,ls="-",lw=2,color="red",
             marker='s',alpha=1,markersize=7,markeredgewidth=1.5, markeredgecolor="red", markerfacecolor='w')
    Lyy15_, = ax2.plot(df_pyy_Jz15["r"],df_pyy_Jz15["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- J_\bot = 1.5
    slope = df_pzz_Jz15["slope_pow"].values[0]
    K_G = round(-slope,2)  
    label = r"$K_{SC}^{zz}$=%.2f"%(K_G)
    Lzz15, = ax2.plot(df_pzz_Jz15["r"],df_pzz_Jz15["corre_abs"],label=label,ls="-",lw=2,color="blue",
             marker='s',alpha=1,markersize=7,markeredgewidth=1.5, markeredgecolor="blue", markerfacecolor='w')
    Lzz15_, = ax2.plot(df_pzz_Jz15["r"],df_pzz_Jz15["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 1.5
    slope = df_ninj_Jz15["slope_pow"].values[0]
    K_F = round(-slope,2)  
    label = r"$K_{CDW}$=%.2f"%(K_F)
    Lninj15, = ax2.plot(df_ninj_Jz15["r"][1:],df_ninj_Jz15["corre_abs"][1:],label=label,ls="-",lw=2,color="green",
             marker='s',alpha=1,markersize=7,markeredgewidth=1.5, markeredgecolor="green", markerfacecolor='w')
    Lninj15_, = ax2.plot(df_ninj_Jz15["r"][1:],df_ninj_Jz15["fitcorre_pow"][1:],label=label,ls="--",lw=1.5,color="k",
             marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- J_\bot =1.5
    slope = df_sisj_Jz15["slope_pow"].values[0]
    K_F = round(-slope,2)  
    label = r"$K_{F}$=%.2f"%(K_F)
    Lsisj15, = ax2.plot(df_sisj_Jz15["r"][1:],df_sisj_Jz15["corre_abs"][1:],label=label,ls="-",lw=2,color="magenta",
             marker='o',alpha=1,markersize=7,markeredgewidth=1.5, markeredgecolor="magenta", markerfacecolor='w')
    Lsisj15_, = ax2.plot(df_sisj_Jz15["r"][1:],df_sisj_Jz15["fitcorre_pow"][1:],label=label,ls="--",lw=1.5,color="k",
             marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')    
    #-------- J_\bot = 1.5
    slope = df_cicj_Jz15["slope_pow"].values[0]
    K_F = round(-slope,2)  
    label = r"$K_{G}$=%.2f"%(K_F)
    Lcicj15, = ax2.plot(df_cicj_Jz15["r"][1:],df_cicj_Jz15["corre_abs"][1:],label=label,ls="-",lw=2,color="cyan",
             marker='o',alpha=1,markersize=7,markeredgewidth=1.5, markeredgecolor="cyan", markerfacecolor='w')
    Lcicj15_, = ax2.plot(df_cicj_Jz15["r"][1:],df_cicj_Jz15["fitcorre_pow"][1:],label=label,ls="--",lw=1.5,color="k",
             marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w') 
    ####图例设置
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 16, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[Lyy15,Lzz15,Lninj15,Lsisj15,Lcicj15], loc = 4, bbox_to_anchor=(0.55, 0.02),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    #legend2=plt.legend(handles=[L10,L20,L40], loc = 4, bbox_to_anchor=(0.98, 0.70),
    #                   ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    #plt.gca().add_artist(legend1)
    label_x = r"r"
    label_y = "Correlation"
    #plt.yscale("log")
    #plt.xscale("log")
    ax2.set_xscale("log",base=10,subs=[0.01,0.02,0.03]) 
    ax2.set_yscale("log",base=10,subs=[0.01,0.02,0.03])      
    ax2.set_xlabel(label_x, size= 16)
    ax2.set_ylabel(label_y, size= 16)
    ax2.tick_params(labelsize = 14) # 设置坐标刻度对应数字的大小
    ax2.set_xlim([0,40])
    #ax2.set_xticks([0,10,20,30])
    #ax2.set_yticks([-1,-0.5,0,0.5,1]) 
    #
    #=========================================================
    ax2.text(28,0.015,"(b)",fontsize = 20, color='black', rotation = 0)
    ax2.text(5.5,5e-2, r'$\mathrm{\delta} = \frac{11}{24}$', fontsize = 20, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
    ax2.text(15,5e-2, r'$J_{\bot}=1.5$', fontsize = 16, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
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
    #-----------------------------------------------------------------------------------------------
    # 选择子图 ax3 进行绘图
    plt.sca(ax3) ## 选择对 ax1 进行绘图
    ax3 = plt.gca()
    #--------- J_\bot = 2.5
    slope = df_pyy_Jz25["slope_pow"].values[0]
    K_F = round(-slope,2)  
    label = r"$K_{SC}^{yy}$=%.2f"%(K_F)
    Lyy25, = ax3.plot(df_pyy_Jz25["r"],df_pyy_Jz25["corre_abs"],label=label,ls="-",lw=2,color="red",
             marker='s',alpha=1,markersize=7,markeredgewidth=1.5, markeredgecolor="red", markerfacecolor='w')
    Lyy25_, = ax3.plot(df_pyy_Jz25["r"],df_pyy_Jz25["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- J_\bot = 2.5
    slope = df_pzz_Jz25["slope_pow"].values[0]
    K_G = round(-slope,2)  
    label = r"$K_{SC}^{zz}$=%.2f"%(K_G)
    Lzz25, = ax3.plot(df_pzz_Jz25["r"],df_pzz_Jz25["corre_abs"],label=label,ls="-",lw=2,color="blue",
             marker='s',alpha=1,markersize=7,markeredgewidth=1.5, markeredgecolor="blue", markerfacecolor='w')
    Lzz25_, = ax3.plot(df_pzz_Jz25["r"],df_pzz_Jz25["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 2.5
    slope = df_ninj_Jz25["slope_pow"].values[0]
    K_F = round(-slope,2)  
    label = r"$K_{CDW}$=%.2f"%(K_F)
    Lninj25, = ax3.plot(df_ninj_Jz25["r"][1:],df_ninj_Jz25["corre_abs"][1:],label=label,ls="-",lw=2,color="green",
             marker='s',alpha=1,markersize=7,markeredgewidth=1.5, markeredgecolor="green", markerfacecolor='w')
    Lninj25_, = ax3.plot(df_ninj_Jz25["r"][1:],df_ninj_Jz25["fitcorre_pow"][1:],label=label,ls="--",lw=1.5,color="k",
             marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- J_\bot =2.5
    slope = df_sisj_Jz25["slope_pow"].values[0]
    K_F = round(-slope,2)  
    label = r"$K_{F}$=%.2f"%(K_F)
    Lsisj25, = ax3.plot(df_sisj_Jz25["r"][1:],df_sisj_Jz25["corre_abs"][1:],label=label,ls="-",lw=2,color="magenta",
             marker='o',alpha=1,markersize=7,markeredgewidth=1.5, markeredgecolor="magenta", markerfacecolor='w')
    Lsisj25_, = ax3.plot(df_sisj_Jz25["r"][1:],df_sisj_Jz25["fitcorre_pow"][1:],label=label,ls="--",lw=1.5,color="k",
             marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')    
    #-------- J_\bot = 2.5
    slope = df_cicj_Jz25["slope_pow"].values[0]
    K_F = round(-slope,2)  
    label = r"$K_{G}$=%.2f"%(K_F)
    Lcicj25, = ax3.plot(df_cicj_Jz25["r"][1:],df_cicj_Jz25["corre_abs"][1:],label=label,ls="-",lw=2,color="cyan",
             marker='o',alpha=1,markersize=7,markeredgewidth=1.5, markeredgecolor="cyan", markerfacecolor='w')
    Lcicj25_, = ax3.plot(df_cicj_Jz25["r"][1:],df_cicj_Jz25["fitcorre_pow"][1:],label=label,ls="--",lw=1.5,color="k",
             marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w') 
    ####图例设置
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 16, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[Lyy25,Lzz25,Lninj25,Lsisj25,Lcicj25], loc = 4, bbox_to_anchor=(0.55, 0.02),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    #legend2=plt.legend(handles=[L10,L20,L40], loc = 4, bbox_to_anchor=(0.98, 0.70),
    #                   ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    #plt.gca().add_artist(legend1)
    label_x = r"r"
    label_y = "Correlation"
    #plt.yscale("log")
    #plt.xscale("log")
    ax3.set_xscale("log",base=10,subs=[0.01,0.02,0.03]) 
    ax3.set_yscale("log",base=10,subs=[0.01,0.02,0.03])      
    ax3.set_xlabel(label_x, size= 16)
    ax3.set_ylabel(label_y, size= 16)
    ax3.tick_params(labelsize = 14) # 设置坐标刻度对应数字的大小
    ax3.set_xlim([0,40])
    #ax2.set_xticks([0,10,20,30])
    #ax2.set_yticks([-1,-0.5,0,0.5,1]) 
    #
    #=========================================================
    ax3.text(28,0.01,"(c)",fontsize = 20, color='black', rotation = 0)
    ax3.text(5.5,1e-1, r'$\mathrm{\delta} = \frac{11}{24}$', fontsize = 20, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
    ax3.text(15,1e-1, r'$J_{\bot}=2.5$', fontsize = 16, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
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
    #ax3.yaxis.get_minor_locator().set_params(numticks=13, subs=[.2,.4,.6,.8]) # 将次要刻度显示出来 
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
    #========================================================
    # 选择子图 ax4 进行绘图
    plt.sca(ax4) ## 选择对 ax1 进行绘图
    ax4 = plt.gca()

    LK2, = ax4.plot([0,2.5],[2,2],label=label,ls="--",lw=2,color="k",
             marker='s',alpha=1,markersize=0,markeredgewidth=1, markeredgecolor="green", markerfacecolor='w')
    label = "$K_{SC}^{yy}$"
    Lyy, = ax4.plot(Jz_list,K_sc_yy,label=label,ls="-",lw=2,color="red",
             marker='s',alpha=1,markersize=7,markeredgewidth=1.5, markeredgecolor="red", markerfacecolor='w')
    label = "$K_{SC}^{zz}$"
    Lzz, = ax4.plot(Jz_list,K_sc_zz,label=label,ls="-",lw=2,color="blue",
             marker='s',alpha=1,markersize=7,markeredgewidth=1.5, markeredgecolor="blue", markerfacecolor='w')
    label = "$K_{CDW}$"
    Lninj, = ax4.plot(Jz_list,K_cdw,label=label,ls="-",lw=2,color="green",
             marker='s',alpha=1,markersize=7,markeredgewidth=1.5, markeredgecolor="green", markerfacecolor='w')
    label = "$K_{F}$"
    Lsisj, = ax4.plot(Jz_list,K_f,label=label,ls="-",lw=2,color="magenta",
             marker='o',alpha=1,markersize=7,markeredgewidth=1.5, markeredgecolor="magenta", markerfacecolor='w')    
    label = "$K_{G}$"
    Lcicj, = ax4.plot(Jz_list,K_g,label=label,ls="-",lw=2,color="cyan",
             marker='o',alpha=1,markersize=7,markeredgewidth=1.5, markeredgecolor="cyan", markerfacecolor='w')       
    #
    ####图例设置
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 16, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[Lyy,Lzz,Lninj,Lsisj,Lcicj], loc = 4, bbox_to_anchor=(0.52, 0.38),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    label_x = r"$J_{\bot}$"
    label_y = "K"
    #plt.yscale("log")
    #plt.xscale("log")
    #ax4.set_xscale("log",base=10,subs=[0.01,0.02,0.03]) 
    #ax4.set_yscale("log",base=10,subs=[0.01,0.02,0.03])      
    ax4.set_xlabel(label_x, size= 18)
    ax4.set_ylabel(label_y, size= 18)
    ax4.tick_params(labelsize = 16) # 设置坐标刻度对应数字的大小
    ax4.set_xlim([0,2.5])
    ax4.set_ylim([0,7])
    #ax4.set_xticks([0,10,20,30])
    #ax4.set_yticks([-1,-0.5,0,0.5,1]) 
    #
    #=========================================================
    ax4.text(2.1,6,"(d)",fontsize = 20, color='black', rotation = 0)
    ax4.text(0.8,6, r'$\mathrm{\delta} = \frac{11}{24}$', fontsize = 20, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
    #=========================================================
    #======== 坐标轴设置第一层
    labels = ax4.get_xticklabels() + ax4.get_yticklabels()
    #[label.set_fontname('Times New Roman') for label in labels]###设置ticket labled的字体格式
    ax4.xaxis.set_minor_locator(MultipleLocator(0.25))###设置次刻度的间隔
    ax4.yaxis.set_minor_locator(MultipleLocator(0.5))###设置次刻度的间隔
    ax4.xaxis.set_major_formatter(FormatStrFormatter('%1.1f'))###设置X轴标签文本格式
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
    ax4.spines['right'].set_linewidth(1.5)  ###设置右边坐标轴的粗细
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
    #========================================================
    plt.show()