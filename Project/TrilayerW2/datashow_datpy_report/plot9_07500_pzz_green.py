import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib.colors import ListedColormap
import seaborn as sns
#
def get_data_cicj_layer1(Lz, Ly, Lx, dop, t, J, Jz, dim,):
    workpath = "E:\\WORK\\Work\\Project\\TrilayerW2-tj-Jz-openz"
    filepath1 = "\\data_dmrgcpp\\Lz%d_Ly%d_Lx%d\\dop%g" % (Lz, Ly, Lx,dop)
    filepath2 = "\\t%.2f_J%.2f_Jz%.2f_dim%d" % (t, J, Jz, dim)
    filepath3 = "\\measurement_green_function_s2_fit.parquet"
    filename = workpath + filepath1 + filepath2 + filepath3
    print(filename)
    df = pd.read_parquet(filename)
    print(df.head())
    return df
#
#
def get_data_cicj_layer2(Lz, Ly, Lx, dop, t, J, Jz, dim,):
    workpath = "E:\\WORK\\Work\\Project\\TrilayerW2-tj-Jz-openz"
    filepath1 = "\\data_dmrgcpp\\Lz%d_Ly%d_Lx%d\\dop%g" % (Lz, Ly, Lx,dop)
    filepath2 = "\\t%.2f_J%.2f_Jz%.2f_dim%d" % (t, J, Jz, dim)
    filepath3 = "\\measurement_green_function_fit.parquet"
    filename = workpath + filepath1 + filepath2 + filepath3
    print(filename)
    df = pd.read_parquet(filename)
    print(df.head())
    return df
#
def get_data_cicj_layer3(Lz, Ly, Lx, dop, t, J, Jz, dim,):
    workpath = "E:\\WORK\\Work\\Project\\TrilayerW2-tj-Jz-openz"
    filepath1 = "\\data_dmrgcpp\\Lz%d_Ly%d_Lx%d\\dop%g" % (Lz, Ly, Lx,dop)
    filepath2 = "\\t%.2f_J%.2f_Jz%.2f_dim%d" % (t, J, Jz, dim)
    filepath3 = "\\measurement_green_function_s1_fit.parquet"
    filename = workpath + filepath1 + filepath2 + filepath3
    print(filename)
    df = pd.read_parquet(filename)
    print(df.head())
    return df
#
def get_data_sisj_layer1(Lz, Ly, Lx, dop, t, J, Jz, dim,):
    workpath = "E:\\WORK\\Work\\Project\\TrilayerW2-tj-Jz-openz"
    filepath1 = "\\data_dmrgcpp\\Lz%d_Ly%d_Lx%d\\dop%g" % (Lz, Ly, Lx,dop)
    filepath2 = "\\t%.2f_J%.2f_Jz%.2f_dim%d" % (t, J, Jz, dim)
    filepath3 = "\\measurement_spin_correlation_s2_fit.parquet"
    filename = workpath + filepath1 + filepath2 + filepath3
    print(filename)
    df = pd.read_parquet(filename)
    print(df.head())
    return df
#
def get_data_sisj_layer2(Lz, Ly, Lx, dop, t, J, Jz, dim,):
    workpath = "E:\\WORK\\Work\\Project\\TrilayerW2-tj-Jz-openz"
    filepath1 = "\\data_dmrgcpp\\Lz%d_Ly%d_Lx%d\\dop%g" % (Lz, Ly, Lx,dop)
    filepath2 = "\\t%.2f_J%.2f_Jz%.2f_dim%d" % (t, J, Jz, dim)
    filepath3 = "\\measurement_spin_correlation_fit.parquet"
    filename = workpath + filepath1 + filepath2 + filepath3
    print(filename)
    df = pd.read_parquet(filename)
    print(df.head())
    return df
#
def get_data_sisj_layer3(Lz, Ly, Lx, dop, t, J, Jz, dim,):
    workpath = "E:\\WORK\\Work\\Project\\TrilayerW2-tj-Jz-openz"
    filepath1 = "\\data_dmrgcpp\\Lz%d_Ly%d_Lx%d\\dop%g" % (Lz, Ly, Lx,dop)
    filepath2 = "\\t%.2f_J%.2f_Jz%.2f_dim%d" % (t, J, Jz, dim)
    filepath3 = "\\measurement_spin_correlation_s1_fit.parquet"
    filename = workpath + filepath1 + filepath2 + filepath3
    print(filename)
    df = pd.read_parquet(filename)
    print(df.head())
    return df
#
def get_data_ninj_corre_layer1(Lz, Ly, Lx, dop, t, J, Jz, dim,):
    workpath = "E:\\WORK\\Work\\Project\\TrilayerW2-tj-Jz-openz"
    filepath1 = "\\data_dmrgcpp\\Lz%d_Ly%d_Lx%d\\dop%g" % (Lz, Ly, Lx,dop)
    filepath2 = "\\t%.2f_J%.2f_Jz%.2f_dim%d" % (t, J, Jz, dim)
    filepath3 = "\\measurement_density_correlation_s2_fit.parquet"
    filename = workpath + filepath1 + filepath2 + filepath3
    print(filename)
    df = pd.read_parquet(filename)
    print(df.head())
    return df
#
#
def get_data_ninj_corre_layer2(Lz, Ly, Lx, dop, t, J, Jz, dim,):
    workpath = "E:\\WORK\\Work\\Project\\TrilayerW2-tj-Jz-openz"
    filepath1 = "\\data_dmrgcpp\\Lz%d_Ly%d_Lx%d\\dop%g" % (Lz, Ly, Lx,dop)
    filepath2 = "\\t%.2f_J%.2f_Jz%.2f_dim%d" % (t, J, Jz, dim)
    filepath3 = "\\measurement_density_correlation_fit.parquet"
    filename = workpath + filepath1 + filepath2 + filepath3
    print(filename)
    df = pd.read_parquet(filename)
    print(df.head())
    return df
#
def get_data_ninj_corre_layer3(Lz, Ly, Lx, dop, t, J, Jz, dim,):
    workpath = "E:\\WORK\\Work\\Project\\TrilayerW2-tj-Jz-openz"
    filepath1 = "\\data_dmrgcpp\\Lz%d_Ly%d_Lx%d\\dop%g" % (Lz, Ly, Lx,dop)
    filepath2 = "\\t%.2f_J%.2f_Jz%.2f_dim%d" % (t, J, Jz, dim)
    filepath3 = "\\measurement_density_correlation_s1_fit.parquet"
    filename = workpath + filepath1 + filepath2 + filepath3
    print(filename)
    df = pd.read_parquet(filename)
    print(df.head())
    return df
#
def get_data_pzz_layer12(Lz, Ly, Lx, dop, t, J, Jz, dim,):
    workpath = "E:\\WORK\\Work\\Project\\TrilayerW2-tj-Jz-openz"
    filepath1 = "\\data_dmrgcpp\\Lz%d_Ly%d_Lx%d\\dop%g" % (Lz, Ly, Lx,dop)
    filepath2 = "\\t%.2f_J%.2f_Jz%.2f_dim%d" % (t, J, Jz, dim)
    filepath3 = "\\measurement_pairing_zz_la23_fit.parquet"
    filename = workpath + filepath1 + filepath2 + filepath3
    print(filename)
    df = pd.read_parquet(filename)
    print(df.head())
    return df
#
def get_data_pzz_layer23(Lz, Ly, Lx, dop, t, J, Jz, dim,):
    workpath = "E:\\WORK\\Work\\Project\\TrilayerW2-tj-Jz-openz"
    filepath1 = "\\data_dmrgcpp\\Lz%d_Ly%d_Lx%d\\dop%g" % (Lz, Ly, Lx,dop)
    filepath2 = "\\t%.2f_J%.2f_Jz%.2f_dim%d" % (t, J, Jz, dim)
    filepath3 = "\\measurement_pairing_zz_fit.parquet"
    filename = workpath + filepath1 + filepath2 + filepath3
    print(filename)
    df = pd.read_parquet(filename)
    print(df.head())
    return df
#
def get_data_pzz_layer13(Lz, Ly, Lx, dop, t, J, Jz, dim,):
    workpath = "E:\\WORK\\Work\\Project\\TrilayerW2-tj-Jz-openz"
    filepath1 = "\\data_dmrgcpp\\Lz%d_Ly%d_Lx%d\\dop%g" % (Lz, Ly, Lx,dop)
    filepath2 = "\\t%.2f_J%.2f_Jz%.2f_dim%d" % (t, J, Jz, dim)
    filepath3 = "\\measurement_pairing_zz_la13_fit.parquet"
    filename = workpath + filepath1 + filepath2 + filepath3
    print(filename)
    df = pd.read_parquet(filename)
    print(df.head())
    return df
#
if __name__ =="__main__":
    Lz = 3
    Ly = 2
    Lx = 48
    dop = 216
    t = 3
    J = 1
    #---------------plot logr-r fig and logr-logr fig-----------------------------
    fig = plt.figure(figsize=(20,10))
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
    ax1 = plt.axes([0.1,0.55,0.24,0.4])
    ax2 = plt.axes([0.4,0.55,0.24,0.4])
    ax3 = plt.axes([0.7,0.55,0.24,0.4])
    ax4 = plt.axes([0.1,0.1,0.24,0.4])
    ax5 = plt.axes([0.4,0.1,0.24,0.4])
    ax6 = plt.axes([0.7,0.1,0.24,0.4])
    # plt.subplots_adjust(left=0.1, bottom=0.1, right=0.98, top=0.96, wspace=0.32, hspace=0.26)
    # 选择子图 ax2 进行绘图
    plt.sca(ax2) ## 选择对 ax2 进行绘图
    ax2 = plt.gca()
    #
    df_pzz_Jz120_mu12_10000 = get_data_pzz_layer12(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=12.0,dim=10000)
    df_pzz_Jz120_mu23_10000 = get_data_pzz_layer23(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=12.0,dim=10000)
    #
    df_cicj_Jz120_mu1_10000 = get_data_cicj_layer1(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=12.0,dim=10000)
    df_cicj_Jz120_mu2_10000 = get_data_cicj_layer2(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=12.0,dim=10000)
    df_cicj_Jz120_mu3_10000 = get_data_cicj_layer3(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=12.0,dim=10000)
    #----------------- plot pzz^{12}-----------------------------------------------------------------------
    slope = df_pzz_Jz120_mu12_10000["slope_pow"].values[0]
    Ks = round(-slope,2) 
    label = r"$P_{zz}^{1,2}$, $K_{P_{zz}}^{1,2}$=%.2f"%(Ks)
    L120_pzzmu12_10000, = ax2.plot(df_pzz_Jz120_mu12_10000["r"],df_pzz_Jz120_mu12_10000["corre_abs"],label=label,ls="-",lw=1.5,color="red",
             marker='o',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="red", markerfacecolor='None')
    L120_pzzmu12_10000_, = ax2.plot(df_pzz_Jz120_mu12_10000["r"],df_pzz_Jz120_mu12_10000["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='o',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #----------------- plot pzz^{23}-----------------------------------------------------------------------
    slope = df_pzz_Jz120_mu23_10000["slope_pow"].values[0]
    Ks = round(-slope,2) 
    label = r"$P_{zz}^{2,3}$, $K_{P_{zz}}^{2,3}$=%.2f"%(Ks)
    L120_pzzmu23_10000, = ax2.plot(df_pzz_Jz120_mu23_10000["r"],df_pzz_Jz120_mu23_10000["corre_abs"],label=label,ls="-",lw=1.5,color="magenta",
             marker='o',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="magenta", markerfacecolor='None')
    L120_pzzmu23_10000_, = ax2.plot(df_pzz_Jz120_mu23_10000["r"],df_pzz_Jz120_mu23_10000["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='o',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #----------------- plot G(r)^{\mu=1} * G(r)^{\mu=2} -------------------------------------------------------------------
    #slope = df_cicj_Jz120_mu1_10000['slope_pow'].values[0]
    #Ks = round(-slope,2)
    label = r"$G(r)^{\mu=1}G(r)^{\mu=2}$"
    L120_Gmu1mu2_10000, = ax2.plot(df_cicj_Jz120_mu1_10000["r"],df_cicj_Jz120_mu1_10000["corre_abs"]*df_cicj_Jz120_mu2_10000["corre_abs"],
                                   label=label,ls="-",lw=1.5,color="blue", marker='^',alpha=1,markersize=8,markeredgewidth=1, 
                                   markeredgecolor="blue", markerfacecolor='None')
    #----------------- plot G(r)^{\mu=2} * G(r)^{\mu=3} -------------------------------------------------------------------
    #slope = df_cicj_Jz120_mu2_10000['slope_pow'].values[0]
    #Ks = round(-slope,2)
    label = r"$G(r)^{\mu=2}G(r)^{\mu=3}$"
    L120_Gmu2mu3_10000, = ax2.plot(df_cicj_Jz120_mu2_10000["r"],df_cicj_Jz120_mu2_10000["corre_abs"]*df_cicj_Jz120_mu3_10000["corre_abs"],
                                   label=label,ls="-",lw=1.5,color="green", marker='v',alpha=1,markersize=8,markeredgewidth=1, 
                                   markeredgecolor="green", markerfacecolor='None')
    #---------------------------------------------------------------------------------------------------------------------
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 14, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[L120_pzzmu12_10000,L120_pzzmu23_10000,L120_Gmu1mu2_10000,L120_Gmu2mu3_10000,], 
                       loc = 4, bbox_to_anchor=(0.55, -0.05),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    #legend2=plt.legend(handles=[], loc = 4, bbox_to_anchor=(0.65, -0.04),
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
    ax2.set_xlim([0,35])
    ax2.set_xticks([5,10,15,20,30,])
    #ax2.set_yticks([-1,-0.5,0,0.5,1]) 
    #=========================================================
    ax2.text(25,0.03,"(b)",fontsize = 20, color='black', rotation = 0)
    ax2.text(5,0.03, r'$\mathrm{\delta} = \frac{3}{4}$', fontsize = 16, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
    ax2.text(8,0.03, r'$J_{\bot}=12.0$', fontsize = 16, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
    ax2.text(8,0.01, r'D=10000', fontsize = 16, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
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
    #ax2.yaxis.get_minor_locator().set_params(numticks=13, subs=[.2,.4,.6,.8]) # 将次要刻度显示出来 
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
    #------------------------------------------------------------------------------------------

    #---------------------------------------------------------------------------------------------
    # 选择子图 ax6 进行绘图
    plt.sca(ax6) ## 选择对 ax6 进行绘图
    ax6 = plt.gca()
    #
    df_pzz_Jz240_mu12_6000 = get_data_pzz_layer12(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=24.0,dim=6000)
    df_pzz_Jz240_mu23_6000 = get_data_pzz_layer23(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=24.0,dim=6000)
    #
    df_cicj_Jz240_mu1_6000 = get_data_cicj_layer1(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=24.0,dim=6000)
    df_cicj_Jz240_mu2_6000 = get_data_cicj_layer2(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=24.0,dim=6000)
    df_cicj_Jz240_mu3_6000 = get_data_cicj_layer3(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=24.0,dim=6000)
    #----------------- plot pzz^{12}-----------------------------------------------------------------------
    slope = df_pzz_Jz240_mu12_6000["slope_pow"].values[0]
    Ks = round(-slope,2) 
    label = r"$P_{zz}^{1,2}$, $K_{P_{zz}}^{1,2}$=%.2f"%(Ks)
    L240_pzzmu12_6000, = ax6.plot(df_pzz_Jz240_mu12_6000["r"],df_pzz_Jz240_mu12_6000["corre_abs"],label=label,ls="-",lw=1.5,color="red",
             marker='o',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="red", markerfacecolor='None')
    L240_pzzmu12_6000_, = ax6.plot(df_pzz_Jz240_mu12_6000["r"],df_pzz_Jz240_mu12_6000["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='o',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #----------------- plot pzz^{23}-----------------------------------------------------------------------
    slope = df_pzz_Jz240_mu23_6000["slope_pow"].values[0]
    Ks = round(-slope,2) 
    label = r"$P_{zz}^{2,3}$, $K_{P_{zz}}^{2,3}$=%.2f"%(Ks)
    L240_pzzmu23_6000, = ax6.plot(df_pzz_Jz240_mu23_6000["r"],df_pzz_Jz240_mu23_6000["corre_abs"],label=label,ls="-",lw=1.5,color="magenta",
             marker='o',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="magenta", markerfacecolor='None')
    L240_pzzmu23_6000_, = ax6.plot(df_pzz_Jz240_mu23_6000["r"],df_pzz_Jz240_mu23_6000["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='o',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #----------------- plot G(r)^{\mu=1} * G(r)^{\mu=2} -------------------------------------------------------------------
    #slope = df_cicj_Jz120_mu1_10000['slope_pow'].values[0]
    #Ks = round(-slope,2)
    label = r"$G(r)^{\mu=1}G(r)^{\mu=2}$"
    L240_Gmu1mu2_6000, = ax6.plot(df_cicj_Jz240_mu1_6000["r"],df_cicj_Jz240_mu1_6000["corre_abs"]*df_cicj_Jz240_mu2_6000["corre_abs"],
                                   label=label,ls="-",lw=1.5,color="blue", marker='^',alpha=1,markersize=8,markeredgewidth=1, 
                                   markeredgecolor="blue", markerfacecolor='None')
    #----------------- plot G(r)^{\mu=2} * G(r)^{\mu=3} -------------------------------------------------------------------
    #slope = df_cicj_Jz120_mu2_10000['slope_pow'].values[0]
    #Ks = round(-slope,2)
    label = r"$G(r)^{\mu=2}G(r)^{\mu=3}$"
    L240_Gmu2mu3_6000, = ax6.plot(df_cicj_Jz240_mu2_6000["r"],df_cicj_Jz240_mu2_6000["corre_abs"]*df_cicj_Jz240_mu3_6000["corre_abs"],
                                   label=label,ls="-",lw=1.5,color="green", marker='v',alpha=1,markersize=8,markeredgewidth=1, 
                                   markeredgecolor="green", markerfacecolor='None')
    #---------------------------------------------------------------------------------------------------------------------
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 14, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[L240_pzzmu12_6000,L240_pzzmu23_6000,L240_Gmu1mu2_6000,L240_Gmu2mu3_6000,], 
                       loc = 4, bbox_to_anchor=(0.55, -0.05),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    #legend2=plt.legend(handles=[], loc = 4, bbox_to_anchor=(0.65, -0.04),
    #                   ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    #plt.gca().add_artist(legend1)
    label_x = r"r"
    label_y = "Correlation"
    #plt.yscale("log")
    #plt.xscale("log")
    ax6.set_xscale("log",base=10,subs=[0.01,0.02,0.03]) 
    ax6.set_yscale("log",base=10,subs=[0.01,0.02,0.03])      
    ax6.set_xlabel(label_x, size= 16)
    ax6.set_ylabel(label_y, size= 16)
    ax6.tick_params(labelsize = 14) # 设置坐标刻度对应数字的大小
    ax6.set_xlim([0,35])
    ax6.set_xticks([5,10,15,20,30,])
    #ax6.set_yticks([-1,-0.5,0,0.5,1]) 
    #=========================================================
    ax6.text(25,0.03,"(f)",fontsize = 20, color='black', rotation = 0)
    ax6.text(5,0.03, r'$\mathrm{\delta} = \frac{3}{4}$', fontsize = 16, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
    ax6.text(8,0.03, r'$J_{\bot}=24.0$', fontsize = 16, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
    ax6.text(8,0.01, r'D=6000', fontsize = 16, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
    #=========================================================
    # 坐标轴设置第一层
    labels = ax6.get_xticklabels() + ax6.get_yticklabels()
    #[label.set_fontname('Times New Roman') for label in labels]###设置ticket labled的字体格式
    ax6.xaxis.set_minor_locator(MultipleLocator(5))###设置次刻度的间隔
    #ax6.yaxis.set_minor_locator(MultipleLocator(10))###设置次刻度的间隔
    ax6.xaxis.set_major_formatter(FormatStrFormatter('%1.0f'))###设置X轴标签文本格式
    #ax6.yaxis.set_major_formatter(FormatStrFormatter('%1.1f'))###设置Y轴标签文本格式
    #
    #ax6.yaxis.get_major_locator().set_params(numticks=99)
    #ax6.yaxis.get_minor_locator().set_params(numticks=13, subs=[.2,.4,.6,.8]) # 将次要刻度显示出来 
    #----将次要刻度显示出来 
    #locmin = matplotlib.ticker.LogLocator(base=10.0, subs=(0.1,0.2,0.4,0.6,0.8,1,2,4,6,8,10 )) 
    #ax6.xaxis.set_minor_locator(locmin)
    #ax6.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    #=====坐标轴的第二层： 坐标轴的设置
    ax6.spines['bottom'].set_linewidth(1.5) ###设置底部坐标轴的粗细
    ax6.spines['left'].set_linewidth(1.5)   ###设置左边坐标轴的粗细
    ax6.spines['right'].set_linewidth(1.5)  ###设置右边坐标轴的粗细
    ax6.spines['top'].set_linewidth(1.5)    ###设置上部坐标轴的粗细
    #ax6.spines['right'].set_color('none')# 将右边上边的两条边颜色设置为空 其实就相当于抹掉这两条边
    #ax6.spines['top'].set_color('none')
    #====坐标轴的第三层：  主刻度线的设置
    for line in ax6.xaxis.get_ticklines():
        #line is a Line2D instance
        #line.set_color('green')
        line.set_markersize(3)####设置刻度线的长度
        line.set_markeredgewidth(1.5)####设置刻度线的宽度
    for line in ax6.yaxis.get_ticklines():
        # line is a Line2D instance
        #line.set_color('green')
        line.set_markersize(3)####设置刻度线的长度
        line.set_markeredgewidth(1.5)####设置刻度线的宽度
    plt.tick_params(axis="x", which="minor", length=2.5, width=1.5, color="k")  ### 设置次要刻度 
    plt.tick_params(axis="y", which="minor", length=2.5, width=1.5, color="k")  ### 设置次要刻度 
    #-----------------------------------------------------------------------------------------------
    plt.show() 
