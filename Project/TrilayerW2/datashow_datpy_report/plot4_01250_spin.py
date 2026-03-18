import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib.colors import ListedColormap
import seaborn as sns
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
def get_sign_value_matrix(df_layer1,df_layer2,df_layer3):
    layer1 = list(df_layer1["corre"].values)
    layer2 = list(df_layer2["corre"].values)
    layer3 = list(df_layer3["corre"].values)
    M = np.sign(np.array([layer3,layer2,layer1]))
    return M

if __name__ =="__main__":
    Lz = 3
    Ly = 2
    Lx = 48
    t = 3
    J = 1
    #
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
    ax1_sub = plt.axes([0.11,0.56,0.18,0.06])
    ax2 = plt.axes([0.4,0.55,0.24,0.4])
    ax2_sub = plt.axes([0.41,0.56,0.18,0.06])
    ax3 = plt.axes([0.7,0.55,0.24,0.4])
    ax3_sub = plt.axes([0.71,0.56,0.18,0.06])
    ax4 = plt.axes([0.1,0.1,0.24,0.4])
    ax4_sub = plt.axes([0.11,0.11,0.18,0.06])
    ax5 = plt.axes([0.4,0.1,0.24,0.4])
    ax5_sub = plt.axes([0.47,0.41,0.18,0.06])
    # plt.subplots_adjust(left=0.1, bottom=0.1, right=0.98, top=0.96, wspace=0.32, hspace=0.26)
    # 选择子图 ax1 进行绘图
    plt.sca(ax1) ## 选择对 ax1 进行绘图
    ax1 = plt.gca()
    #
    df_sisj_Jz10_mu1 = get_data_sisj_layer1(Lz=Lz,Ly=Ly,Lx=Lx,dop=36,t=t,J=J,Jz=1.0,dim=10000)
    df_sisj_Jz10_mu2 = get_data_sisj_layer2(Lz=Lz,Ly=Ly,Lx=Lx,dop=36,t=t,J=J,Jz=1.0,dim=10000)
    df_sisj_Jz10_mu3 = get_data_sisj_layer3(Lz=Lz,Ly=Ly,Lx=Lx,dop=36,t=t,J=J,Jz=1.0,dim=10000)
    #df_sisj_Jz20 = get_data_sisj(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=2.0,dim=dim)
    #--------- layer 1
    slope = df_sisj_Jz10_mu1["slope_exp"].values[0]
    Ks = round(-1/slope,2)  
    label = r"$\mu=1$, $\xi_F$=%.2f"%(Ks)
    L10_mu1, = ax1.plot(df_sisj_Jz10_mu1["r"],df_sisj_Jz10_mu1["corre_abs"],label=label,ls="-",lw=1.5,color="red",
             marker='o',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="red", markerfacecolor='None')
    L10_mu1_, = ax1.plot(df_sisj_Jz10_mu1["r"],df_sisj_Jz10_mu1["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
             marker='o',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- layer 2
    slope = df_sisj_Jz10_mu2['slope_exp'].values[0]
    Ks = round(-1/slope,2)
    label = r"$\mu=2$"
    L10_mu2, = ax1.plot(df_sisj_Jz10_mu2["r"],df_sisj_Jz10_mu2["corre_abs"],label=label,ls="-",lw=1.5,color="blue",
             marker='o',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="blue", markerfacecolor='None')
    #L10_mu2_, = ax1.plot(df_sisj_Jz10_mu2["r"],df_sisj_Jz10_mu2["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
    #         marker='o',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- layer 3
    slope = df_sisj_Jz10_mu3['slope_exp'].values[0]
    Ks = round(-1/slope,2)
    label = r"$\mu=3$"
    L10_mu3, = ax1.plot(df_sisj_Jz10_mu3["r"],df_sisj_Jz10_mu3["corre_abs"],label=label,ls="-",lw=1.5,color="magenta",
             marker='o',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="magenta", markerfacecolor='None')
    #L10_mu3_, = ax1.plot(df_sisj_Jz10_mu3["r"],df_sisj_Jz10_mu3["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
    #         marker='o',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    # 
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 14, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[L10_mu1,L10_mu2,L10_mu3], loc = 4, bbox_to_anchor=(0.5, 0.2),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    #legend2=plt.legend(handles=[], loc = 4, bbox_to_anchor=(0.65, -0.04),
    #                   ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    #plt.gca().add_artist(legend1)
    label_x = r"r"
    label_y = "$F(r)$"
    #plt.yscale("log")
    #plt.xscale("log")
    #ax1.set_xscale("log",base=10,subs=[0.01,0.02,0.03]) 
    ax1.set_yscale("log",base=10,subs=[0.01,0.02,0.03])      
    ax1.set_xlabel(label_x, size= 16)
    ax1.set_ylabel(label_y, size= 16)
    ax1.tick_params(labelsize = 14) # 设置坐标刻度对应数字的大小
    ax1.set_xlim([0,35])
    ax1.set_xticks([5,10,15,20,30,])
    #ax1.set_yticks([-1,-0.5,0,0.5,1]) 
    #
    #=========================================================
    ax1.text(32,0.15,"(a)",fontsize = 20, color='black', rotation = 0)
    ax1.text(4,1e-1, r'$\mathrm{\delta} = 0.125$', fontsize = 16, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
    ax1.text(15,1e-1, r'$J_{\bot}=1.0$', fontsize = 16, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
    ax1.text(15,0.4e-1, r'$D=10000$', fontsize = 16, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
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
    #ax1.yaxis.get_minor_locator().set_params(numticks=13, subs=[.2,.4,.6,.8]) # 将次要刻度显示出来 
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
    #----------------------------------------------------------------------------------------------------
    # plot ax1_sub
    sign_matrix_10 = get_sign_value_matrix(df_layer1=df_sisj_Jz10_mu1,df_layer2=df_sisj_Jz10_mu2,df_layer3=df_sisj_Jz10_mu3)
    print("sign_matrix_10:\n",sign_matrix_10)
    cmap = ListedColormap(["lightcoral","lightblue"]) # -1 -> light coral, +1 -> light blue
    #plt.sca(ax1_sub) ## 选择对 ax1 进行绘图
    #ax1_sub = plt.gca()
    hm = sns.heatmap(sign_matrix_10, ax=ax1_sub, cmap=cmap, cbar=True,
                cbar_kws={"ticks":[-1,1]},   # no "0"
                linewidths=0.5, linecolor='black',  # grid lines
                xticklabels=False, yticklabels=False, # hide tick labels
                square=False, vmin=-1, vmax=1)  # square -> false, 每一个色块儿的大小会自动调整为适应 ax1_sub size 比例的大小
    # Remove tick lines on colorbar
    cbar = hm.collections[0].colorbar
    cbar.ax.tick_params(size=0, labelsize=14)  # no tick lines, set label font size=14
    ax1_sub.set_title("Sign of $F(r)$")
    #-----------------------------------------------------------------------------------------------------
    plt.sca(ax2) ## 选择对 ax2 进行绘图
    ax2 = plt.gca()
    #
    df_sisj_Jz40_mu1 = get_data_sisj_layer1(Lz=Lz,Ly=Ly,Lx=Lx,dop=36,t=t,J=J,Jz=4.0,dim=8000)
    df_sisj_Jz40_mu2 = get_data_sisj_layer2(Lz=Lz,Ly=Ly,Lx=Lx,dop=36,t=t,J=J,Jz=4.0,dim=8000)
    df_sisj_Jz40_mu3 = get_data_sisj_layer3(Lz=Lz,Ly=Ly,Lx=Lx,dop=36,t=t,J=J,Jz=4.0,dim=8000)
    #--------- layer 1
    slope = df_sisj_Jz40_mu1["slope_pow"].values[0]
    Ks = round(-slope,2)  
    label = r"$\mu=1$, $K_F$=%.2f"%(Ks)
    L40_mu1, = ax2.plot(df_sisj_Jz40_mu1["r"],df_sisj_Jz40_mu1["corre_abs"],label=label,ls="-",lw=1.5,color="red",
             marker='o',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="red", markerfacecolor='None')
    L40_mu1_, = ax2.plot(df_sisj_Jz40_mu1["r"],df_sisj_Jz40_mu1["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='o',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- layer 2
    slope = df_sisj_Jz40_mu2['slope_pow'].values[0]
    Ks = round(-slope,2)
    label = r"$\mu=2$"
    L40_mu2, = ax2.plot(df_sisj_Jz40_mu2["r"],df_sisj_Jz40_mu2["corre_abs"],label=label,ls="-",lw=1.5,color="blue",
             marker='o',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="blue", markerfacecolor='None')
    #L40_mu2_, = ax2.plot(df_sisj_Jz40_mu2["r"],df_sisj_Jz40_mu2["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
    #         marker='o',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')    
    #--------- layer 3
    slope = df_sisj_Jz40_mu3['slope_pow'].values[0]
    Ks = round(-slope,2)
    label = r"$\mu=3$"
    L40_mu3, = ax2.plot(df_sisj_Jz40_mu3["r"],df_sisj_Jz40_mu3["corre_abs"],label=label,ls="-",lw=1.5,color="magenta",
             marker='o',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="magenta", markerfacecolor='None')
    #L40_mu3_, = ax2.plot(df_sisj_Jz40_mu3["r"],df_sisj_Jz40_mu3["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
    #         marker='o',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    # 
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 14, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[L40_mu1,L40_mu2,L40_mu3], loc = 4, bbox_to_anchor=(0.5, 0.2),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    #legend2=plt.legend(handles=[], loc = 4, bbox_to_anchor=(0.65, -0.04),
    #                   ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    #plt.gca().add_artist(legend1)
    label_x = r"r"
    label_y = "$F(r)$"
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
    #
    #=========================================================
    ax2.text(25,0.1,"(b)",fontsize = 20, color='black', rotation = 0)
    ax2.text(3,0.8e-1, r'$\mathrm{\delta} = 0.125$', fontsize = 16, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
    ax2.text(8,0.8e-1, r'$J_{\bot}=4.0$', fontsize = 16, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
    ax2.text(8,0.4e-1, r'$D=8000$', fontsize = 16, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
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
    #-------------------------------------------------------------------------------------------------
    # plot ax2_sub
    sign_matrix_40 = get_sign_value_matrix(df_layer1=df_sisj_Jz40_mu1,df_layer2=df_sisj_Jz40_mu2,df_layer3=df_sisj_Jz40_mu3)
    print("sign_matrix_40:\n",sign_matrix_40)
    cmap = ListedColormap(["lightcoral","lightblue"]) # -1 -> light coral, +1 -> light blue
    hm_ax2 = sns.heatmap(sign_matrix_40, ax=ax2_sub, cmap=cmap, cbar=True,
                cbar_kws={"ticks":[-1,1]},   # no "0"
                linewidths=0.5, linecolor='black',  # grid lines
                xticklabels=False, yticklabels=False, # hide tick labels
                square=False, vmin=-1, vmax=1)  # square -> false, 每一个色块儿的大小会自动调整为适应 ax1_sub size 比例的大小
    # Remove tick lines on colorbar
    cbar = hm_ax2.collections[0].colorbar
    cbar.ax.tick_params(size=0, labelsize=14)  # no tick lines, set label font size=14
    ax2_sub.set_title("Sign of $F(r)$")
    #-------------------------------------------------------------------------------------------------
    plt.sca(ax3) ## 选择对 ax3 进行绘图
    ax3 = plt.gca()
    #
    df_sisj_Jz80_mu1 = get_data_sisj_layer1(Lz=Lz,Ly=Ly,Lx=Lx,dop=36,t=t,J=J,Jz=8.0,dim=6000)
    df_sisj_Jz80_mu2 = get_data_sisj_layer2(Lz=Lz,Ly=Ly,Lx=Lx,dop=36,t=t,J=J,Jz=8.0,dim=6000)
    df_sisj_Jz80_mu3 = get_data_sisj_layer3(Lz=Lz,Ly=Ly,Lx=Lx,dop=36,t=t,J=J,Jz=8.0,dim=6000)
    #--------- layer 1
    slope = df_sisj_Jz80_mu1["slope_pow"].values[0]
    Ks = round(-slope,2)  
    label = r"$\mu=1$, $K_F$=%.2f"%(Ks)
    L80_mu1, = ax3.plot(df_sisj_Jz80_mu1["r"],df_sisj_Jz80_mu1["corre_abs"],label=label,ls="-",lw=1.5,color="red",
             marker='o',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="red", markerfacecolor='None')
    L80_mu1_, = ax3.plot(df_sisj_Jz80_mu1["r"],df_sisj_Jz80_mu1["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='o',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- layer 2
    slope = df_sisj_Jz80_mu2['slope_pow'].values[0]
    Ks = round(-slope,2)
    label = r"$\mu=2$"
    L80_mu2, = ax3.plot(df_sisj_Jz80_mu2["r"],df_sisj_Jz80_mu2["corre_abs"],label=label,ls="-",lw=1.5,color="blue",
             marker='o',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="blue", markerfacecolor='None')
    #L80_mu2_, = ax3.plot(df_sisj_Jz80_mu2["r"],df_sisj_Jz80_mu2["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
    #         marker='o',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- layer 3
    slope = df_sisj_Jz80_mu3['slope_pow'].values[0]
    Ks = round(-slope,2)
    label = r"$\mu=3$"
    L80_mu3, = ax3.plot(df_sisj_Jz80_mu3["r"],df_sisj_Jz80_mu3["corre_abs"],label=label,ls="-",lw=1.5,color="magenta",
             marker='o',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="magenta", markerfacecolor='None')
    #L80_mu3_, = ax3.plot(df_sisj_Jz80_mu3["r"],df_sisj_Jz80_mu3["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
    #         marker='o',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    # 
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 14, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[L80_mu1,L80_mu2,L80_mu3], loc = 4, bbox_to_anchor=(0.40, 0.2),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    #legend2=plt.legend(handles=[], loc = 4, bbox_to_anchor=(0.65, -0.04),
    #                   ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    #plt.gca().add_artist(legend1)
    label_x = r"r"
    label_y = "$F(r)$"
    #plt.yscale("log")
    #plt.xscale("log")
    ax3.set_xscale("log",base=10,subs=[0.01,0.02,0.03]) 
    ax3.set_yscale("log",base=10,subs=[0.01,0.02,0.03])      
    ax3.set_xlabel(label_x, size= 16)
    ax3.set_ylabel(label_y, size= 16)
    ax3.tick_params(labelsize = 14) # 设置坐标刻度对应数字的大小
    ax3.set_xlim([0,35])
    ax3.set_xticks([5,10,15,20,30,])
    #ax3.set_yticks([-1,-0.5,0,0.5,1]) 
    #
    #=========================================================
    ax3.text(25,0.07,"(c)",fontsize = 20, color='black', rotation = 0)
    ax3.text(3,0.7e-1, r'$\mathrm{\delta} = 0.125$', fontsize = 16, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
    ax3.text(8,0.7e-1, r'$J_{\bot}=8.0$', fontsize = 16, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
    ax3.text(8,0.3e-1, r'$D=6000$', fontsize = 16, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
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
    #-------------------------------------------------------------------------------------------------------
    # plot ax3_sub
    sign_matrix_80 = get_sign_value_matrix(df_layer1=df_sisj_Jz80_mu1,df_layer2=df_sisj_Jz80_mu2,df_layer3=df_sisj_Jz80_mu3)
    print("sign_matrix_80:\n",sign_matrix_80)
    cmap = ListedColormap(["lightcoral","lightblue"]) # -1 -> light coral, +1 -> light blue
    hm_ax3 = sns.heatmap(sign_matrix_80, ax=ax3_sub, cmap=cmap, cbar=True,
                cbar_kws={"ticks":[-1,1]},   # no "0"
                linewidths=0.5, linecolor='black',  # grid lines
                xticklabels=False, yticklabels=False, # hide tick labels
                square=False, vmin=-1, vmax=1)  # square -> false, 每一个色块儿的大小会自动调整为适应 ax_sub size 比例的大小
    # Remove tick lines on colorbar
    cbar = hm_ax3.collections[0].colorbar
    cbar.ax.tick_params(size=0, labelsize=14)  # no tick lines, set label font size=14
    ax3_sub.set_title("Sign of $F(r)$")
    #-------------------------------------------------------------------------------------------------------
    plt.sca(ax4) ## 选择对 ax4 进行绘图
    ax4 = plt.gca()
    #
    df_sisj_Jz120_mu1 = get_data_sisj_layer1(Lz=Lz,Ly=Ly,Lx=Lx,dop=36,t=t,J=J,Jz=12.0,dim=8000)
    df_sisj_Jz120_mu2 = get_data_sisj_layer2(Lz=Lz,Ly=Ly,Lx=Lx,dop=36,t=t,J=J,Jz=12.0,dim=8000)
    df_sisj_Jz120_mu3 = get_data_sisj_layer3(Lz=Lz,Ly=Ly,Lx=Lx,dop=36,t=t,J=J,Jz=12.0,dim=8000)
    #--------- layer 1
    slope = df_sisj_Jz120_mu1["slope_pow"].values[0]
    Ks = round(-slope,2)  
    label = r"$\mu=1$, $K_F$=%.2f"%(Ks)
    L120_mu1, = ax4.plot(df_sisj_Jz120_mu1["r"],df_sisj_Jz120_mu1["corre_abs"],label=label,ls="-",lw=1.5,color="red",
             marker='o',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="red", markerfacecolor='None')
    L120_mu1_, = ax4.plot(df_sisj_Jz120_mu1["r"],df_sisj_Jz120_mu1["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='o',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- layer 2
    slope = df_sisj_Jz120_mu2['slope_pow'].values[0]
    Ks = round(-slope,2)
    label = r"$\mu=2$"
    L120_mu2, = ax4.plot(df_sisj_Jz120_mu2["r"],df_sisj_Jz120_mu2["corre_abs"],label=label,ls="-",lw=1.5,color="blue",
             marker='o',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="blue", markerfacecolor='None')
    #L120_mu2_, = ax4.plot(df_sisj_Jz120_mu2["r"],df_sisj_Jz120_mu2["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
    #         marker='o',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- layer 3
    slope = df_sisj_Jz120_mu3['slope_pow'].values[0]
    Ks = round(-slope,2)
    label = r"$\mu=3$"
    L120_mu3, = ax4.plot(df_sisj_Jz120_mu3["r"],df_sisj_Jz120_mu3["corre_abs"],label=label,ls="-",lw=1.5,color="magenta",
             marker='o',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="magenta", markerfacecolor='None')
    #L120_mu3_, = ax4.plot(df_sisj_Jz120_mu3["r"],df_sisj_Jz120_mu3["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
    #         marker='o',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    # 
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 14, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[L120_mu1,L120_mu2,L120_mu3], loc = 4, bbox_to_anchor=(0.5, 0.2),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    #legend2=plt.legend(handles=[], loc = 4, bbox_to_anchor=(0.65, -0.04),
    #                   ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    #plt.gca().add_artist(legend1)
    label_x = r"r"
    label_y = "$F(r)$"
    #plt.yscale("log")
    #plt.xscale("log")
    ax4.set_xscale("log",base=10,subs=[0.01,0.02,0.03]) 
    ax4.set_yscale("log",base=10,subs=[0.01,0.02,0.03])      
    ax4.set_xlabel(label_x, size= 16)
    ax4.set_ylabel(label_y, size= 16)
    ax4.tick_params(labelsize = 14) # 设置坐标刻度对应数字的大小
    ax4.set_xlim([0,35])
    ax4.set_xticks([5,10,15,20,30,])
    #ax4.set_yticks([-1,-0.5,0,0.5,1]) 
    #
    #=========================================================
    ax4.text(25,0.05,"(d)",fontsize = 20, color='black', rotation = 0)
    ax4.text(3,0.5e-1, r'$\mathrm{\delta} = 0.125$', fontsize = 16, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
    ax4.text(8,0.5e-1, r'$J_{\bot}=12.0$', fontsize = 16, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
    ax4.text(8,0.2e-1, r'$D=8000$', fontsize = 16, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
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
    #ax4.yaxis.get_minor_locator().set_params(numticks=13, subs=[.2,.4,.6,.8]) # 将次要刻度显示出来 
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
    #---------------------------------------------------------------------------------------------
    # plot ax4_sub
    sign_matrix_120 = get_sign_value_matrix(df_layer1=df_sisj_Jz120_mu1,df_layer2=df_sisj_Jz120_mu2,df_layer3=df_sisj_Jz120_mu3)
    print("sign_matrix_120:\n",sign_matrix_120)
    cmap = ListedColormap(["lightcoral","lightblue"]) # -1 -> light coral, +1 -> light blue
    hm_ax4 = sns.heatmap(sign_matrix_120, ax=ax4_sub, cmap=cmap, cbar=True,
                cbar_kws={"ticks":[-1,1]},   # no "0"
                linewidths=0.5, linecolor='black',  # grid lines
                xticklabels=False, yticklabels=False, # hide tick labels
                square=False, vmin=-1, vmax=1)  # square -> false, 每一个色块儿的大小会自动调整为适应 ax_sub size 比例的大小
    # Remove tick lines on colorbar
    cbar = hm_ax4.collections[0].colorbar
    cbar.ax.tick_params(size=0, labelsize=14)  # no tick lines, set label font size=14
    ax4_sub.set_title("Sign of $F(r)$")
    #----------------------------------------------------------------------------------------------
    plt.sca(ax5) ## 选择对 ax5 进行绘图
    ax5 = plt.gca()
    #
    df_sisj_Jz240_mu1 = get_data_sisj_layer1(Lz=Lz,Ly=Ly,Lx=Lx,dop=36,t=t,J=J,Jz=24.0,dim=6000)
    df_sisj_Jz240_mu2 = get_data_sisj_layer2(Lz=Lz,Ly=Ly,Lx=Lx,dop=36,t=t,J=J,Jz=24.0,dim=6000)
    df_sisj_Jz240_mu3 = get_data_sisj_layer3(Lz=Lz,Ly=Ly,Lx=Lx,dop=36,t=t,J=J,Jz=24.0,dim=6000)
    #--------- layer 1
    slope = df_sisj_Jz240_mu1["slope_pow"].values[0]
    Ks = round(-slope,2)  
    label = r"$\mu=1$, $K_F$=%.2f"%(Ks)
    L240_mu1, = ax5.plot(df_sisj_Jz240_mu1["r"],df_sisj_Jz240_mu1["corre_abs"],label=label,ls="-",lw=1.5,color="red",
             marker='o',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="red", markerfacecolor='None')
    L240_mu1_, = ax5.plot(df_sisj_Jz240_mu1["r"],df_sisj_Jz240_mu1["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='o',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- layer 2
    slope = df_sisj_Jz240_mu2['slope_pow'].values[0]
    Ks = round(-slope,2)
    label = r"$\mu=2$"
    L240_mu2, = ax5.plot(df_sisj_Jz240_mu2["r"],df_sisj_Jz240_mu2["corre_abs"],label=label,ls="-",lw=1.5,color="blue",
             marker='o',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="blue", markerfacecolor='None')
    #L240_mu2_, = ax5.plot(df_sisj_Jz240_mu2["r"],df_sisj_Jz240_mu2["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
    #         marker='o',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- layer 3
    slope = df_sisj_Jz240_mu3['slope_pow'].values[0]
    Ks = round(-slope,2)
    label = r"$\mu=3$"
    L240_mu3, = ax5.plot(df_sisj_Jz240_mu3["r"],df_sisj_Jz240_mu3["corre_abs"],label=label,ls="-",lw=1.5,color="magenta",
             marker='o',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="magenta", markerfacecolor='None')
    #L240_mu3_, = ax5.plot(df_sisj_Jz240_mu3["r"],df_sisj_Jz240_mu3["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
    #         marker='o',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    # 
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 14, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[L240_mu1,L240_mu2,L240_mu3], loc = 4, bbox_to_anchor=(0.75, -0.04),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    #legend2=plt.legend(handles=[], loc = 4, bbox_to_anchor=(0.65, -0.04),
    #                   ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    #plt.gca().add_artist(legend1)
    label_x = r"r"
    label_y = "$F(r)$"
    #plt.yscale("log")
    #plt.xscale("log")
    ax5.set_xscale("log",base=10,subs=[0.01,0.02,0.03]) 
    ax5.set_yscale("log",base=10,subs=[0.01,0.02,0.03])      
    ax5.set_xlabel(label_x, size= 16)
    ax5.set_ylabel(label_y, size= 16)
    ax5.tick_params(labelsize = 14) # 设置坐标刻度对应数字的大小
    ax5.set_xlim([0,35])
    ax5.set_xticks([5,10,15,20,30,])
    #ax5.set_yticks([-1,-0.5,0,0.5,1]) 
    #
    #=========================================================
    ax5.text(31,0.5e-5,"(e)",fontsize = 20, color='black', rotation = 0)
    ax5.text(1,0.1e-2, r'$\mathrm{\delta} = 0.125$', fontsize = 16, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
    ax5.text(1,0.1e-3, r'$J_{\bot}=24.0$', fontsize = 16, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
    ax5.text(1,0.1e-4, r'$D=6000$', fontsize = 16, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
    #=========================================================
    # 坐标轴设置第一层
    labels = ax5.get_xticklabels() + ax5.get_yticklabels()
    #[label.set_fontname('Times New Roman') for label in labels]###设置ticket labled的字体格式
    ax5.xaxis.set_minor_locator(MultipleLocator(5))###设置次刻度的间隔
    #ax5.yaxis.set_minor_locator(MultipleLocator(10))###设置次刻度的间隔
    ax5.xaxis.set_major_formatter(FormatStrFormatter('%1.0f'))###设置X轴标签文本格式
    #ax5.yaxis.set_major_formatter(FormatStrFormatter('%1.1f'))###设置Y轴标签文本格式
    #
    #ax5.yaxis.get_major_locator().set_params(numticks=99)
    #ax5.yaxis.get_minor_locator().set_params(numticks=13, subs=[.2,.4,.6,.8]) # 将次要刻度显示出来 
    #----将次要刻度显示出来 
    #locmin = matplotlib.ticker.LogLocator(base=10.0, subs=(0.1,0.2,0.4,0.6,0.8,1,2,4,6,8,10 )) 
    #ax5.xaxis.set_minor_locator(locmin)
    #ax5.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    #=====坐标轴的第二层： 坐标轴的设置
    ax5.spines['bottom'].set_linewidth(1.5) ###设置底部坐标轴的粗细
    ax5.spines['left'].set_linewidth(1.5)   ###设置左边坐标轴的粗细
    ax5.spines['right'].set_linewidth(1.5)  ###设置右边坐标轴的粗细
    ax5.spines['top'].set_linewidth(1.5)    ###设置上部坐标轴的粗细
    #ax5.spines['right'].set_color('none')# 将右边上边的两条边颜色设置为空 其实就相当于抹掉这两条边
    #ax5.spines['top'].set_color('none')
    #====坐标轴的第三层：  主刻度线的设置
    for line in ax5.xaxis.get_ticklines():
        #line is a Line2D instance
        #line.set_color('green')
        line.set_markersize(3)####设置刻度线的长度
        line.set_markeredgewidth(1.5)####设置刻度线的宽度
    for line in ax5.yaxis.get_ticklines():
        # line is a Line2D instance
        #line.set_color('green')
        line.set_markersize(3)####设置刻度线的长度
        line.set_markeredgewidth(1.5)####设置刻度线的宽度
    plt.tick_params(axis="x", which="minor", length=2.5, width=1.5, color="k")  ### 设置次要刻度 
    plt.tick_params(axis="y", which="minor", length=2.5, width=1.5, color="k")  ### 设置次要刻度  
    #------------------------------------------------------------------------------------------
    # plot ax5_sub
    sign_matrix_240 = get_sign_value_matrix(df_layer1=df_sisj_Jz240_mu1,df_layer2=df_sisj_Jz240_mu2,df_layer3=df_sisj_Jz240_mu3)
    print("sign_matrix_120:\n",sign_matrix_240)
    cmap = ListedColormap(["lightcoral","lightblue"]) # -1 -> light coral, +1 -> light blue
    hm_ax5 = sns.heatmap(sign_matrix_240, ax=ax5_sub, cmap=cmap, cbar=True,
                cbar_kws={"ticks":[-1,1]},   # no "0"
                linewidths=0.5, linecolor='black',  # grid lines
                xticklabels=False, yticklabels=False, # hide tick labels
                square=False, vmin=-1, vmax=1)  # square -> false, 每一个色块儿的大小会自动调整为适应 ax_sub size 比例的大小
    # Remove tick lines on colorbar
    cbar = hm_ax5.collections[0].colorbar
    cbar.ax.tick_params(size=0, labelsize=14)  # no tick lines, set label font size=14
    ax5_sub.set_title("Sign of $F(r)$")
    #------------------------------------------------------------------------------------------
    #fig.tight_layout() # 自动调整 subplot 间的间隙参数
    #plt.savefig("E:\\WORK\\Work\\Project\\TrilayerW2-tj-Jz-openz\\datashow_datpy_report\\figs\\fig_01250_spin.eps",
    #            dpi=300, format='eps',bbox_inches='tight') # 白边紧凑型
    #
    plt.show()