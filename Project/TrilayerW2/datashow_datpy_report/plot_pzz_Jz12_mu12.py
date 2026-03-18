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
    #---------------plot logr-logr fig-----------------------------
    fig = plt.figure(figsize=(4,6))
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
    #ax1 = plt.axes([0.12,0.57,0.8,0.4])
    ax2 = plt.axes([0.1,0.1,0.8,0.8])
    #ax3 = plt.axes([0.1,0.1,0.37,0.375])
    #ax4 = plt.axes([0.58,0.1,0.37,0.375])
    # plt.subplots_adjust(left=0.1, bottom=0.1, right=0.98, top=0.96, wspace=0.32, hspace=0.26)
    #-------------------------------------------------------------------------------------------
    # 选择子图 ax2 进行绘图
    plt.sca(ax2) ## 选择对 ax2 进行绘图
    ax2 = plt.gca()
    #------- J_\bot = 12.0, delta = 0.3750
    df_pzz_0375_Jz120_mu12 = get_data_pzz_layer12(Lz=3,Ly=2,Lx=48,dop=108,t=3,J=1,Jz=12.0,dim=14000)
    df_pzz_0375_Jz120_mu23 = get_data_pzz_layer23(Lz=3,Ly=2,Lx=48,dop=108,t=3,J=1,Jz=12.0,dim=14000)
    df_pzz_0375_Jz120_mu13 = get_data_pzz_layer13(Lz=3,Ly=2,Lx=48,dop=108,t=3,J=1,Jz=12.0,dim=14000)
    #------ J_\bot = 12.0, delta = 0.50
    df_pzz_05_Jz120_mu12 = get_data_pzz_layer12(Lz=3,Ly=2,Lx=48,dop=144,t=3,J=1,Jz=12.0,dim=14000)
    df_pzz_05_Jz120_mu23 = get_data_pzz_layer23(Lz=3,Ly=2,Lx=48,dop=144,t=3,J=1,Jz=12.0,dim=14000)
    df_pzz_05_Jz120_mu13 = get_data_pzz_layer13(Lz=3,Ly=2,Lx=48,dop=144,t=3,J=1,Jz=12.0,dim=14000)
    #----- J_\bot = 12.0, delta = 0.6667
    df_pzz_0667_Jz120_mu12 = get_data_pzz_layer12(Lz=3,Ly=2,Lx=48,dop=192,t=3,J=1,Jz=12.0,dim=6000)
    df_pzz_0667_Jz120_mu23 = get_data_pzz_layer23(Lz=3,Ly=2,Lx=48,dop=192,t=3,J=1,Jz=12.0,dim=6000)
    df_pzz_0667_Jz120_mu13 = get_data_pzz_layer13(Lz=3,Ly=2,Lx=48,dop=192,t=3,J=1,Jz=12.0,dim=6000)
    #
    #--------- layer 12  J_\bot = 12.0, delta = 0.3750
    slope = df_pzz_0375_Jz120_mu12["slope_pow"].values[0]
    Ks = round(-slope,2)  
    label = r"$\delta=0.375,K_{P_{zz}}^{1,2}$=%.2f,D=%d"%(Ks,14000)
    L0375_120_mu12, = ax2.plot(df_pzz_0375_Jz120_mu12["r"],df_pzz_0375_Jz120_mu12["corre_abs"],label=label,ls="-",lw=1.5,color="red",
             marker='o',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="red", markerfacecolor='None')
    L0375_120_mu12_, = ax2.plot(df_pzz_0375_Jz120_mu12["r"],df_pzz_0375_Jz120_mu12["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='o',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- layer 23  J_\bot = 12.0, delta = 0.3750
    #slope = df_pzz_0375_Jz120_mu23['slope_pow'].values[0]
    #Ks = round(-slope,2)
    #label = "$K_{P_{zz}}^{2,3}$=%.2f"%(Ks)
    #L0375_120_mu23, = ax2.plot(df_pzz_0375_Jz120_mu23["r"],df_pzz_0375_Jz120_mu23["corre_abs"],label=label,ls="-",lw=1.5,color="magenta",
    #         marker='o',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="magenta", markerfacecolor='None')
    #--------------------------------------------------------------------------------------
    #--------- layer 12  J_\bot = 12.0, delta = 0.50
    slope = df_pzz_05_Jz120_mu12["slope_pow"].values[0]
    Ks = round(-slope,2)  
    label = r"$\delta=0.5,K_{P_{zz}}^{1,2}$=%.2f,D=%d"%(Ks,14000)
    L05_120_mu12, = ax2.plot(df_pzz_05_Jz120_mu12["r"],df_pzz_05_Jz120_mu12["corre_abs"],label=label,ls="-",lw=1.5,color="blue",
             marker='^',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="blue", markerfacecolor='None')
    L05_120_mu12_, = ax2.plot(df_pzz_05_Jz120_mu12["r"],df_pzz_05_Jz120_mu12["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='^',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- layer 23  J_\bot = 12.0, delta = 0.50
    #slope = df_pzz_05_Jz120_mu23['slope_pow'].values[0]
    #Ks = round(-slope,2)
    #label = "$K_{P_{zz}}^{2,3}$=%.2f"%(Ks)
    #L05_120_mu23, = ax2.plot(df_pzz_05_Jz120_mu23["r"],df_pzz_05_Jz120_mu23["corre_abs"],label=label,ls="-",lw=1.5,color="cyan",
    #         marker='^',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="cyan", markerfacecolor='None')
    #-------------------------------------------------------------------------------------------
    #--------- layer 12  J_\bot = 12.0, delta = 0.6667
    slope = df_pzz_0667_Jz120_mu12["slope_pow"].values[0]
    Ks = round(-slope,2)  
    label = r"$\delta=0.6667,K_{P_{zz}}^{1,2}$=%.2f,D=%d"%(Ks,6000)
    L0667_120_mu12, = ax2.plot(df_pzz_0667_Jz120_mu12["r"],df_pzz_0667_Jz120_mu12["corre_abs"],label=label,ls="-",lw=1.5,color="green",
             marker='s',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="green", markerfacecolor='None')
    L0667_120_mu12_, = ax2.plot(df_pzz_0667_Jz120_mu12["r"],df_pzz_0667_Jz120_mu12["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- layer 23  J_\bot = 12.0, delta = 0.6667
    #slope = df_pzz_0667_Jz120_mu23['slope_pow'].values[0]
    #Ks = round(-slope,2)
    #label = "$K_{P_{zz}}^{2,3}$=%.2f"%(Ks)
    #L0667_120_mu23, = ax2.plot(df_pzz_0667_Jz120_mu23["r"],df_pzz_0667_Jz120_mu23["corre_abs"],label=label,ls="-",lw=1.5,color="olive",
    #         marker='s',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="olive", markerfacecolor='None')
    #-------------------------------------------------------------------------------------------
    ####图例设置
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 15, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[L0375_120_mu12,L05_120_mu12,L0667_120_mu12], 
                       loc = 4, bbox_to_anchor=(0.67, 0.05),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    #legend2=plt.legend(handles=[Lyy10_36,Lyy10_48,Lyy10_72], loc = 4, bbox_to_anchor=(0.64, -0.01),
    #                   ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    #plt.gca().add_artist(legend1)#####把图例legend1重新加载回来
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
    ax2.text(25,0.05,"(a)",fontsize = 20, color='black', rotation = 0)
    #ax2.text(3,0.05, r'$\mathrm{\delta} = 0.3750$', fontsize = 16, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
    ax2.text(8,0.05, r'$J_{\bot}=12.0$', fontsize = 16, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
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
    #
    #fig.tight_layout() # 自动调整 subplot 间的间隙参数
    #plt.savefig("E:\\WORK\\Work\\Project\\TrilayerW2-tj-Jz-openz\\datashow_datpy_report\\figs\\fig_pzz_J12.eps",
    #            dpi=300, format='eps',bbox_inches='tight') # 白边紧凑型
    plt.show()
    