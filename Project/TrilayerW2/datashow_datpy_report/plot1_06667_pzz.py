import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib.colors import ListedColormap
import seaborn as sns
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
def get_sign_value_matrix(df_layer12,df_layer23,df_layer13):
    layer12 = list(df_layer12["corre"].values)
    layer23 = list(df_layer23["corre"].values)
    layer13 = list(df_layer13["corre"].values)

    layer_list = [layer13, layer23, layer12]
    # find the maximum length
    max_len = max(len(lst) for lst in layer_list)
    # Pad with zeros at the top
    padded_lists = [[0]*(max_len - len(lst)) + lst for lst in layer_list]
    M = np.sign(np.array(padded_lists))
    return M
#
if __name__ =="__main__":
    Lz = 3
    Ly = 2
    Lx = 48
    dop = 192
    t = 3
    J = 1
    dim = 6000
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
    ax1_sub = plt.axes([0.12,0.56,0.18,0.06])
    ax2 = plt.axes([0.4,0.55,0.24,0.4])
    ax2_sub = plt.axes([0.42,0.56,0.18,0.06])
    ax3 = plt.axes([0.7,0.55,0.24,0.4])
    ax3_sub = plt.axes([0.72,0.56,0.18,0.06])
    ax4 = plt.axes([0.1,0.1,0.24,0.4])
    ax4_sub = plt.axes([0.12,0.11,0.18,0.06])
    ax5 = plt.axes([0.4,0.1,0.24,0.4])
    ax5_sub = plt.axes([0.42,0.11,0.18,0.06])
    ax6 = plt.axes([0.7,0.1,0.24,0.4])
    ax6_sub = plt.axes([0.72,0.11,0.18,0.06])
    # plt.subplots_adjust(left=0.1, bottom=0.1, right=0.98, top=0.96, wspace=0.32, hspace=0.26)
    # 选择子图 ax1 进行绘图
    plt.sca(ax1) ## 选择对 ax1 进行绘图
    ax1 = plt.gca()
    #
    df_pzz_Jz10_mu12 = get_data_pzz_layer12(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=1.0,dim=dim)
    df_pzz_Jz10_mu23 = get_data_pzz_layer23(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=1.0,dim=dim)
    df_pzz_Jz10_mu13 = get_data_pzz_layer13(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=1.0,dim=dim)
    #--------- layer 12
    slope = df_pzz_Jz10_mu12["slope_pow"].values[0]
    Ks = round(-slope,2)  
    label = r"$K_{P_{zz}}^{1,2}$=%.2f"%(Ks)
    L10_mu12, = ax1.plot(df_pzz_Jz10_mu12["r"],df_pzz_Jz10_mu12["corre_abs"],label=label,ls="-",lw=1.5,color="red",
             marker='o',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="red", markerfacecolor='None')
    L10_mu12_, = ax1.plot(df_pzz_Jz10_mu12["r"],df_pzz_Jz10_mu12["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='o',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- layer 23
    slope = df_pzz_Jz10_mu23['slope_pow'].values[0]
    Ks = round(-slope,2)
    label = r"$K_{P_{zz}}^{2,3}$=%.2f"%(Ks)
    L10_mu23, = ax1.plot(df_pzz_Jz10_mu23["r"],df_pzz_Jz10_mu23["corre_abs"],label=label,ls="-",lw=1.5,color="magenta",
             marker='o',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="magenta", markerfacecolor='None')
    #L10_mu23_, = ax1.plot(df_pyy_Jz10_mu23["r"],df_pyy_Jz10_mu23["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
    #         marker='o',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- layer 13
    slope = df_pzz_Jz10_mu13['slope_pow'].values[0]
    Ks = round(-slope,2)
    label = r"$K_{P_{zz}}^{1,3}$=%.2f"%(Ks)
    L10_mu13, = ax1.plot(df_pzz_Jz10_mu13["r"],df_pzz_Jz10_mu13["corre_abs"],label=label,ls="-",lw=1.5,color="blue",
             marker='o',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="blue", markerfacecolor='None')
    #L10_mu13_, = ax1.plot(df_pzz_Jz10_mu13["r"],df_pzz_Jz10_mu13["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
    #         marker='o',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    # 
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 14, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[L10_mu12,L10_mu23,L10_mu13], loc = 4, bbox_to_anchor=(0.35, 0.18),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    #legend2=plt.legend(handles=[], loc = 4, bbox_to_anchor=(0.65, -0.04),
    #                   ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    #plt.gca().add_artist(legend1)
    label_x = r"r"
    label_y = "$P_{zz}(r)$"
    #plt.yscale("log")
    #plt.xscale("log")
    ax1.set_xscale("log",base=10,subs=[0.01,0.02,0.03]) 
    ax1.set_yscale("log",base=10,subs=[0.01,0.02,0.03])      
    ax1.set_xlabel(label_x, size= 16)
    ax1.set_ylabel(label_y, size= 16)
    ax1.tick_params(labelsize = 14) # 设置坐标刻度对应数字的大小
    ax1.set_xlim([0,35])
    ax1.set_xticks([5,10,15,20,30,])
    #ax1.set_yticks([-1,-0.5,0,0.5,1]) 
    #
    #=========================================================
    ax1.text(25,0.01,"(a)",fontsize = 20, color='black', rotation = 0)
    ax1.text(5,0.003, r'$\mathrm{\delta} = \frac{2}{3}$', fontsize = 16, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
    ax1.text(8,0.003, r'$J_{\bot}=1.0$', fontsize = 16, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
    ax1.text(8,0.001, r'$D=6000$', fontsize = 16, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
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
    sign_matrix_10 = get_sign_value_matrix(df_layer12=df_pzz_Jz10_mu12,
                                                                df_layer23=df_pzz_Jz10_mu23,
                                                                df_layer13=df_pzz_Jz10_mu13)
    print("sign_matrix_10_n:\n",sign_matrix_10)
    cmap = ListedColormap(["lightcoral", "white", "lightblue"]) # -1 -> light coral, +1 -> light blue
    row_labels = ["$\mu^{1,3}$", "$\mu^{2,3}$", "$\mu^{1,2}$"]
    hm = sns.heatmap(sign_matrix_10, ax=ax1_sub, cmap=cmap, cbar=True,
                cbar_kws={"ticks":[-1,1]},   # no "0"
                linewidths=0.5, linecolor='black',  # grid lines
                xticklabels=False,  # hide tick labels
                yticklabels=row_labels,
                square=False, vmin=-1, vmax=1)  # square -> false, 每一个色块儿的大小会自动调整为适应 ax1_sub size 比例的大小
    # Hide tick marks but keep labels
    ax1_sub.tick_params(axis="y", length=0) 
    # Control row label size via ax
    ax1_sub.set_yticklabels(ax1_sub.get_yticklabels(), fontsize=10)
    # Remove tick lines on colorbar
    cbar = hm.collections[0].colorbar
    cbar.ax.tick_params(size=0, labelsize=14)  # no tick lines, set label font size=14
    ax1_sub.set_title("Sign of $P_{zz}(r)$")
    #-------------------------------------------------------------------------------------------------------
    plt.sca(ax2) ## 选择对 ax2 进行绘图
    ax2 = plt.gca()
    #
    df_pzz_Jz40_mu12 = get_data_pzz_layer12(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=4.0,dim=dim)
    df_pzz_Jz40_mu23 = get_data_pzz_layer23(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=4.0,dim=dim)
    df_pzz_Jz40_mu13 = get_data_pzz_layer13(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=4.0,dim=dim)
    #--------- layer 12
    slope = df_pzz_Jz40_mu12["slope_pow"].values[0]
    Ks = round(-slope,2)  
    label = r"$K_{P_{zz}}^{1,2}$=%.2f"%(Ks)
    L40_mu12, = ax2.plot(df_pzz_Jz40_mu12["r"],df_pzz_Jz40_mu12["corre_abs"],label=label,ls="-",lw=1.5,color="red",
             marker='o',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="red", markerfacecolor='None')
    L40_mu12_, = ax2.plot(df_pzz_Jz40_mu12["r"],df_pzz_Jz40_mu12["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='o',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- layer 23
    slope = df_pzz_Jz40_mu23['slope_pow'].values[0]
    Ks = round(-slope,2)
    label = r"$K_{P_{zz}}^{2,3}$=%.2f"%(Ks)
    L40_mu23, = ax2.plot(df_pzz_Jz40_mu23["r"],df_pzz_Jz40_mu23["corre_abs"],label=label,ls="-",lw=1.5,color="magenta",
             marker='o',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="magenta", markerfacecolor='None')
    #L40_mu23_, = ax2.plot(df_pzz_Jz40_mu23["r"],df_pzz_Jz40_mu23["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
    #         marker='o',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- layer 13
    slope = df_pzz_Jz40_mu13['slope_pow'].values[0]
    Ks = round(-slope,2)
    label = r"$K_{P_{zz}}^{1,3}$=%.2f"%(Ks)
    L40_mu13, = ax2.plot(df_pzz_Jz40_mu13["r"],df_pzz_Jz40_mu13["corre_abs"],label=label,ls="-",lw=1.5,color="blue",
             marker='o',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="blue", markerfacecolor='None')
    #L40_mu13_, = ax2.plot(df_pzz_Jz40_mu13["r"],df_pzz_Jz40_mu13["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
    #         marker='o',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    # 
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 14, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[L40_mu12,L40_mu23,L40_mu13], loc = 4, bbox_to_anchor=(0.35, 0.18),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    #legend2=plt.legend(handles=[], loc = 4, bbox_to_anchor=(0.65, -0.04),
    #                   ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    #plt.gca().add_artist(legend1)
    label_x = r"r"
    label_y = "$P_{zz}(r)$"
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
    ax2.text(25,0.02,"(b)",fontsize = 20, color='black', rotation = 0)
    ax2.text(5,0.006, r'$\mathrm{\delta} = \frac{2}{3}$', fontsize = 16, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
    ax2.text(8,0.006, r'$J_{\bot}=4.0$', fontsize = 16, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
    ax2.text(8,0.002, r'$D=6000$', fontsize = 16, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
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
    #----------------------------------------------------------------------------------------------------
    # plot ax2_sub
    sign_matrix_40 = get_sign_value_matrix(df_layer12=df_pzz_Jz40_mu12,
                                                                df_layer23=df_pzz_Jz40_mu23,
                                                                df_layer13=df_pzz_Jz40_mu13)
    print("sign_matrix_40_n:\n",sign_matrix_40)
    cmap = ListedColormap(["lightcoral", "white", "lightblue"]) # -1 -> light coral, +1 -> light blue
    row_labels = ["$\mu^{1,3}$", "$\mu^{2,3}$", "$\mu^{1,2}$"]
    hm = sns.heatmap(sign_matrix_40, ax=ax2_sub, cmap=cmap, cbar=True,
                cbar_kws={"ticks":[-1,1]},   # no "0"
                linewidths=0.5, linecolor='black',  # grid lines
                xticklabels=False,  # hide tick labels
                yticklabels=row_labels,
                square=False, vmin=-1, vmax=1)  # square -> false, 每一个色块儿的大小会自动调整为适应 ax1_sub size 比例的大小
    # Hide tick marks but keep labels
    ax2_sub.tick_params(axis="y", length=0) 
    # Control row label size via ax
    ax2_sub.set_yticklabels(ax2_sub.get_yticklabels(), fontsize=10)
    # Remove tick lines on colorbar
    cbar = hm.collections[0].colorbar
    cbar.ax.tick_params(size=0, labelsize=14)  # no tick lines, set label font size=14
    ax2_sub.set_title("Sign of $P_{zz}(r)$")
    #---------------------------------------------------------------------------------------------
    plt.sca(ax3) ## 选择对 ax2 进行绘图
    ax3 = plt.gca()
    #
    df_pzz_Jz80_mu12 = get_data_pzz_layer12(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=8.0,dim=dim)
    df_pzz_Jz80_mu23 = get_data_pzz_layer23(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=8.0,dim=dim)
    df_pzz_Jz80_mu13 = get_data_pzz_layer13(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=8.0,dim=dim)
    #--------- layer 12
    slope = df_pzz_Jz80_mu12["slope_pow"].values[0]
    Ks = round(-slope,2)  
    label = r"$K_{P_{zz}}^{1,2}$=%.2f"%(Ks)
    L80_mu12, = ax3.plot(df_pzz_Jz80_mu12["r"],df_pzz_Jz80_mu12["corre_abs"],label=label,ls="-",lw=1.5,color="red",
             marker='o',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="red", markerfacecolor='None')
    L80_mu12_, = ax3.plot(df_pzz_Jz80_mu12["r"],df_pzz_Jz80_mu12["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='o',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- layer 23
    slope = df_pzz_Jz80_mu23['slope_pow'].values[0]
    Ks = round(-slope,2)
    label = r"$K_{P_{zz}}^{2,3}$=%.2f"%(Ks)
    L80_mu23, = ax3.plot(df_pzz_Jz80_mu23["r"],df_pzz_Jz80_mu23["corre_abs"],label=label,ls="-",lw=1.5,color="magenta",
             marker='o',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="magenta", markerfacecolor='None')
    #L80_mu23_, = ax3.plot(df_pzz_Jz80_mu23["r"],df_pzz_Jz80_mu23["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
    #         marker='o',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- layer 13
    slope = df_pzz_Jz80_mu13['slope_pow'].values[0]
    Ks = round(-slope,2)
    label = r"$K_{P_{zz}}^{1,3}$=%.2f"%(Ks)
    L80_mu13, = ax3.plot(df_pzz_Jz80_mu13["r"],df_pzz_Jz80_mu13["corre_abs"],label=label,ls="-",lw=1.5,color="blue",
             marker='o',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="blue", markerfacecolor='None')
    #L80_mu13_, = ax3.plot(df_pzz_Jz80_mu13["r"],df_pzz_Jz80_mu13["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
    #         marker='o',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    # 
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 14, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[L80_mu12,L80_mu23,L80_mu13], loc = 4, bbox_to_anchor=(0.35, 0.18),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    #legend2=plt.legend(handles=[], loc = 4, bbox_to_anchor=(0.65, -0.04),
    #                   ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    #plt.gca().add_artist(legend1)
    label_x = r"r"
    label_y = "$P_{zz}(r)$"
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
    ax3.text(25,0.035,"(c)",fontsize = 20, color='black', rotation = 0)
    ax3.text(5,0.035, r'$\mathrm{\delta} = \frac{2}{3}$', fontsize = 16, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
    ax3.text(8,0.035, r'$J_{\bot}=8.0$', fontsize = 16, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
    ax3.text(8,0.008, r'$D=6000$', fontsize = 16, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
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
    #----------------------------------------------------------------------------------------------------
    # plot ax3_sub
    sign_matrix_80 = get_sign_value_matrix(df_layer12=df_pzz_Jz80_mu12,
                                                                df_layer23=df_pzz_Jz80_mu23,
                                                                df_layer13=df_pzz_Jz80_mu13)
    print("sign_matrix_80_n:\n",sign_matrix_80)
    cmap = ListedColormap(["lightcoral", "white", "lightblue"]) # -1 -> light coral, +1 -> light blue
    row_labels = ["$\mu^{1,3}$", "$\mu^{2,3}$", "$\mu^{1,2}$"]
    hm = sns.heatmap(sign_matrix_80, ax=ax3_sub, cmap=cmap, cbar=True,
                cbar_kws={"ticks":[-1,1]},   # no "0"
                linewidths=0.5, linecolor='black',  # grid lines
                xticklabels=False,  # hide tick labels
                yticklabels=row_labels,
                square=False, vmin=-1, vmax=1)  # square -> false, 每一个色块儿的大小会自动调整为适应 ax1_sub size 比例的大小
    # Hide tick marks but keep labels
    ax3_sub.tick_params(axis="y", length=0) 
    # Control row label size via ax
    ax3_sub.set_yticklabels(ax3_sub.get_yticklabels(), fontsize=10)
    # Remove tick lines on colorbar
    cbar = hm.collections[0].colorbar
    cbar.ax.tick_params(size=0, labelsize=14)  # no tick lines, set label font size=14
    ax3_sub.set_title("Sign of $P_{zz}(r)$")
    #--------------------------------------------------------------------------------------------
    plt.sca(ax4) ## 选择对 ax4 进行绘图
    ax4 = plt.gca()
    #
    df_pzz_Jz120_mu12 = get_data_pzz_layer12(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=12.0,dim=dim)
    df_pzz_Jz120_mu23 = get_data_pzz_layer23(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=12.0,dim=dim)
    df_pzz_Jz120_mu13 = get_data_pzz_layer13(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=12.0,dim=dim)
    #--------- layer 12
    slope = df_pzz_Jz120_mu12["slope_pow"].values[0]
    Ks = round(-slope,2)  
    label = r"$K_{P_{zz}}^{1,2}$=%.2f"%(Ks)
    L120_mu12, = ax4.plot(df_pzz_Jz120_mu12["r"],df_pzz_Jz120_mu12["corre_abs"],label=label,ls="-",lw=1.5,color="red",
             marker='o',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="red", markerfacecolor='None')
    L120_mu12_, = ax4.plot(df_pzz_Jz120_mu12["r"],df_pzz_Jz120_mu12["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='o',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- layer 23
    slope = df_pzz_Jz120_mu23['slope_pow'].values[0]
    Ks = round(-slope,2)
    label = r"$K_{P_{zz}}^{2,3}$=%.2f"%(Ks)
    L120_mu23, = ax4.plot(df_pzz_Jz120_mu23["r"],df_pzz_Jz120_mu23["corre_abs"],label=label,ls="-",lw=1.5,color="magenta",
             marker='o',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="magenta", markerfacecolor='None')
    #L120_mu23_, = ax4.plot(df_pzz_Jz120_mu23["r"],df_pzz_Jz80_mu23["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
    #         marker='o',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- layer 13
    slope = df_pzz_Jz120_mu13['slope_pow'].values[0]
    Ks = round(-slope,2)
    label = r"$K_{P_{zz}}^{1,3}$=%.2f"%(Ks)
    L120_mu13, = ax4.plot(df_pzz_Jz120_mu13["r"],df_pzz_Jz120_mu13["corre_abs"],label=label,ls="-",lw=1.5,color="blue",
             marker='o',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="blue", markerfacecolor='None')
    #L120_mu13_, = ax4.plot(df_pzz_Jz120_mu13["r"],df_pzz_Jz120_mu13["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
    #         marker='o',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    # 
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 14, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[L120_mu12,L120_mu23,L120_mu13], loc = 4, bbox_to_anchor=(0.35, 0.18),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    #legend2=plt.legend(handles=[], loc = 4, bbox_to_anchor=(0.65, -0.04),
    #                   ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    #plt.gca().add_artist(legend1)
    label_x = r"r"
    label_y = "$P_{zz}(r)$"
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
    ax4.text(5,0.05, r'$\mathrm{\delta} = \frac{2}{3}$', fontsize = 16, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
    ax4.text(8,0.05, r'$J_{\bot}=12.0$', fontsize = 16, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
    ax4.text(8,0.01, r'$D=6000$', fontsize = 16, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
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
    #----------------------------------------------------------------------------------------------------
    # plot ax4_sub
    sign_matrix_120 = get_sign_value_matrix(df_layer12=df_pzz_Jz120_mu12,
                                                                df_layer23=df_pzz_Jz120_mu23,
                                                                df_layer13=df_pzz_Jz120_mu13)
    print("sign_matrix_80_n:\n",sign_matrix_120)
    cmap = ListedColormap(["lightcoral", "white", "lightblue"]) # -1 -> light coral, +1 -> light blue
    row_labels = ["$\mu^{1,3}$", "$\mu^{2,3}$", "$\mu^{1,2}$"]
    hm = sns.heatmap(sign_matrix_120, ax=ax4_sub, cmap=cmap, cbar=True,
                cbar_kws={"ticks":[-1,1]},   # no "0"
                linewidths=0.5, linecolor='black',  # grid lines
                xticklabels=False,  # hide tick labels
                yticklabels=row_labels,
                square=False, vmin=-1, vmax=1)  # square -> false, 每一个色块儿的大小会自动调整为适应 ax1_sub size 比例的大小
    # Hide tick marks but keep labels
    ax4_sub.tick_params(axis="y", length=0) 
    # Control row label size via ax
    ax4_sub.set_yticklabels(ax4_sub.get_yticklabels(), fontsize=10)
    # Remove tick lines on colorbar
    cbar = hm.collections[0].colorbar
    cbar.ax.tick_params(size=0, labelsize=14)  # no tick lines, set label font size=14
    ax4_sub.set_title("Sign of $P_{zz}(r)$")
    #--------------------------------------------------------------------------------------------------
    plt.sca(ax5) ## 选择对 ax4 进行绘图
    ax5 = plt.gca()
    #
    df_pzz_Jz240_mu12 = get_data_pzz_layer12(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=24.0,dim=dim)
    df_pzz_Jz240_mu23 = get_data_pzz_layer23(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=24.0,dim=dim)
    df_pzz_Jz240_mu13 = get_data_pzz_layer13(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=24.0,dim=dim)
    #--------- layer 12
    slope = df_pzz_Jz240_mu12["slope_pow"].values[0]
    Ks = round(-slope,2)  
    label = r"$K_{P_{zz}}^{1,2}$=%.2f"%(Ks)
    L240_mu12, = ax5.plot(df_pzz_Jz240_mu12["r"],df_pzz_Jz240_mu12["corre_abs"],label=label,ls="-",lw=1.5,color="red",
             marker='o',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="red", markerfacecolor='None')
    L240_mu12_, = ax5.plot(df_pzz_Jz240_mu12["r"],df_pzz_Jz240_mu12["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='o',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- layer 23
    slope = df_pzz_Jz240_mu23['slope_pow'].values[0]
    Ks = round(-slope,2)
    label = r"$K_{P_{zz}}^{2,3}$=%.2f"%(Ks)
    L240_mu23, = ax5.plot(df_pzz_Jz240_mu23["r"],df_pzz_Jz240_mu23["corre_abs"],label=label,ls="-",lw=1.5,color="magenta",
             marker='o',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="magenta", markerfacecolor='None')
    #L240_mu23_, = ax5.plot(df_pzz_Jz240_mu23["r"],df_pzz_Jz240_mu23["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
    #         marker='o',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- layer 13
    slope = df_pzz_Jz240_mu13['slope_pow'].values[0]
    Ks = round(-slope,2)
    label = r"$K_{P_{zz}}^{1,3}$=%.2f"%(Ks)
    L240_mu13, = ax5.plot(df_pzz_Jz240_mu13["r"],df_pzz_Jz240_mu13["corre_abs"],label=label,ls="-",lw=1.5,color="blue",
             marker='o',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="blue", markerfacecolor='None')
    #L240_mu13_, = ax5.plot(df_pzz_Jz240_mu13["r"],df_pzz_Jz240_mu13["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
    #         marker='o',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    # 
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 14, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[L240_mu12,L240_mu23,L240_mu13], loc = 4, bbox_to_anchor=(0.35, 0.18),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    #legend2=plt.legend(handles=[], loc = 4, bbox_to_anchor=(0.65, -0.04),
    #                   ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    #plt.gca().add_artist(legend1)
    label_x = r"r"
    label_y = "$P_{zz}(r)$"
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
    ax5.text(25,0.05,"(e)",fontsize = 20, color='black', rotation = 0)
    ax5.text(5,0.0001, r'$\mathrm{\delta} = \frac{2}{3}$', fontsize = 16, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
    ax5.text(8,0.0001, r'$J_{\bot}=24.0$', fontsize = 16, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
    ax5.text(8,0.00001, r'$D=6000$', fontsize = 16, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
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
    #----------------------------------------------------------------------------------------------------
    # plot ax5_sub
    sign_matrix_240 = get_sign_value_matrix(df_layer12=df_pzz_Jz240_mu12,
                                                                df_layer23=df_pzz_Jz240_mu23,
                                                                df_layer13=df_pzz_Jz240_mu13)
    print("sign_matrix_240_n:\n",sign_matrix_240)
    cmap = ListedColormap(["lightcoral", "white", "lightblue"]) # -1 -> light coral, +1 -> light blue
    row_labels = ["$\mu^{1,3}$", "$\mu^{2,3}$", "$\mu^{1,2}$"]
    hm = sns.heatmap(sign_matrix_240, ax=ax5_sub, cmap=cmap, cbar=True,
                cbar_kws={"ticks":[-1,1]},   # no "0"
                linewidths=0.5, linecolor='black',  # grid lines
                xticklabels=False,  # hide tick labels
                yticklabels=row_labels,
                square=False, vmin=-1, vmax=1)  # square -> false, 每一个色块儿的大小会自动调整为适应 ax1_sub size 比例的大小
    # Hide tick marks but keep labels
    ax5_sub.tick_params(axis="y", length=0) 
    # Control row label size via ax
    ax5_sub.set_yticklabels(ax5_sub.get_yticklabels(), fontsize=10)
    # Remove tick lines on colorbar
    cbar = hm.collections[0].colorbar
    cbar.ax.tick_params(size=0, labelsize=14)  # no tick lines, set label font size=14
    ax5_sub.set_title("Sign of $P_{zz}(r)$")
    #-----------------------------------------------------------------------------------------------------------
    plt.sca(ax6) ## 选择对 ax4 进行绘图
    ax6 = plt.gca()
    #
    df_pzz_Jz480_mu12 = get_data_pzz_layer12(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=48.0,dim=dim)
    df_pzz_Jz480_mu23 = get_data_pzz_layer23(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=48.0,dim=dim)
    df_pzz_Jz480_mu13 = get_data_pzz_layer13(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=48.0,dim=dim)
    #--------- layer 12
    slope = df_pzz_Jz480_mu12["slope_pow"].values[0]
    Ks = round(-slope,2)  
    label = r"$K_{P_{zz}}^{1,2}$=%.2f"%(Ks)
    L480_mu12, = ax6.plot(df_pzz_Jz480_mu12["r"],df_pzz_Jz480_mu12["corre_abs"],label=label,ls="-",lw=1.5,color="red",
             marker='o',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="red", markerfacecolor='None')
    L480_mu12_, = ax6.plot(df_pzz_Jz480_mu12["r"],df_pzz_Jz480_mu12["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='o',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- layer 23
    slope = df_pzz_Jz480_mu23['slope_pow'].values[0]
    Ks = round(-slope,2)
    label = r"$K_{P_{zz}}^{2,3}$=%.2f"%(Ks)
    L480_mu23, = ax6.plot(df_pzz_Jz480_mu23["r"],df_pzz_Jz480_mu23["corre_abs"],label=label,ls="-",lw=1.5,color="magenta",
             marker='o',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="magenta", markerfacecolor='None')
    #L480_mu23_, = ax6.plot(df_pzz_Jz480_mu23["r"],df_pzz_Jz480_mu23["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
    #         marker='o',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- layer 13
    slope = df_pzz_Jz480_mu13['slope_pow'].values[0]
    Ks = round(-slope,2)
    label = r"$K_{P_{zz}}^{1,3}$=%.2f"%(Ks)
    L480_mu13, = ax6.plot(df_pzz_Jz480_mu13["r"],df_pzz_Jz480_mu13["corre_abs"],label=label,ls="-",lw=1.5,color="blue",
             marker='o',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="blue", markerfacecolor='None')
    #L480_mu13_, = ax6.plot(df_pzz_Jz480_mu13["r"],df_pzz_Jz480_mu13["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
    #         marker='o',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    # 
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 14, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[L480_mu12,L480_mu23,L480_mu13], loc = 4, bbox_to_anchor=(0.35, 0.18),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    #legend2=plt.legend(handles=[], loc = 4, bbox_to_anchor=(0.65, -0.04),
    #                   ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    #plt.gca().add_artist(legend1)
    label_x = r"r"
    label_y = "$P_{zz}(r)$"
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
    #
    #=========================================================
    ax6.text(25,0.03,"(f)",fontsize = 20, color='black', rotation = 0)
    ax6.text(5,0.0001, r'$\mathrm{\delta} = \frac{2}{3}$', fontsize = 16, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
    ax6.text(8,0.0001, r'$J_{\bot}=48.0$', fontsize = 16, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
    ax6.text(8,0.00001, r'$D=6000$', fontsize = 16, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
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
    #----------------------------------------------------------------------------------------------------
    # plot ax6_sub
    sign_matrix_480 = get_sign_value_matrix(df_layer12=df_pzz_Jz480_mu12,
                                                                df_layer23=df_pzz_Jz480_mu23,
                                                                df_layer13=df_pzz_Jz480_mu13)
    print("sign_matrix_480_n:\n",sign_matrix_480)
    cmap = ListedColormap(["lightcoral", "white", "lightblue"]) # -1 -> light coral, +1 -> light blue
    row_labels = ["$\mu^{1,3}$", "$\mu^{2,3}$", "$\mu^{1,2}$"]
    hm = sns.heatmap(sign_matrix_480, ax=ax6_sub, cmap=cmap, cbar=True,
                cbar_kws={"ticks":[-1,1]},   # no "0"
                linewidths=0.5, linecolor='black',  # grid lines
                xticklabels=False,  # hide tick labels
                yticklabels=row_labels,
                square=False, vmin=-1, vmax=1)  # square -> false, 每一个色块儿的大小会自动调整为适应 ax1_sub size 比例的大小
    # Hide tick marks but keep labels
    ax6_sub.tick_params(axis="y", length=0) 
    # Control row label size via ax
    ax6_sub.set_yticklabels(ax6_sub.get_yticklabels(), fontsize=10)
    # Remove tick lines on colorbar
    cbar = hm.collections[0].colorbar
    cbar.ax.tick_params(size=0, labelsize=14)  # no tick lines, set label font size=14
    ax6_sub.set_title("Sign of $P_{zz}(r)$")
    #---------------------------------------------------
    #fig.tight_layout() # 自动调整 subplot 间的间隙参数
    #plt.savefig("E:\\WORK\\Work\\Project\\TrilayerW2-tj-Jz-openz\\datashow_datpy_report\\figs\\fig_06667_pzz.eps",
    #            dpi=300, format='eps',bbox_inches='tight') # 白边紧凑型
    #
    plt.show()