# 将 YY pairing, ZZ pairing, Density correlation 画在一起
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
#
def entropy(Lx, x, c, g):
    S = []
    for it in x:
        s = (c/6)*np.log((Lx/np.pi)*np.sin(it*np.pi/Lx)) + g
        S.append(s)    
    return S
#
def df_along_Lx(df,Lz,Ly):
    if Ly == 2:
        print("---------df_res0---------")
        df_res0 = df[df["site"]%(Lz*Ly)==0]
        df_res0["r"] = range(1,len(df_res0)+1) 
        print(df_res0.head())
        print(df_res0.tail())
        print(len(df_res0))
        print("---------df_res1---------")
        df_res1 = df[df["site"]%(Lz*Ly)==1]
        df_res1["r"] = range(1,len(df_res1)+1) 
        print(df_res1.head())
        print(df_res1.tail())
        print(len(df_res1))
        print("---------df_res2---------")
        df_res2 = df[df["site"]%(Lz*Ly)==2]
        df_res2["r"] = range(1,len(df_res2)+1) 
        print(df_res2.head())
        print(df_res2.tail())
        print(len(df_res2))
        print("---------df_res3---------")
        df_res3 = df[df["site"]%(Lz*Ly)==3]
        df_res3["r"] = range(1 , len(df_res3)+1) 
        print(df_res3.head())
        print(df_res3.tail())
        print(len(df_res3))
    return df_res0, df_res1, df_res2, df_res3
#
def get_data_entropy(Lz,Ly,Lx,dop,t,J,Jz,dim):
    #Lz = 2
    #Ly = 2
    #Lx = 48
    #dop = 72
    #t = 3
    #J = 1
    #Jz = 0.4
    #dim = 6000 # dim cutoff
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
    return df_res0
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
    dop = 96
    t = 3
    J = 1
    dim = 6000
    # load spin correlation data
    df_sisj_Jz025 = get_data_sisj(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=0.25,dim=dim)
    df_sisj_Jz075 = get_data_sisj(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=0.75,dim=dim)
    df_sisj_Jz175 = get_data_sisj(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=1.75,dim=dim)
    # load spin correlation data
    df_sisj_Jz25 = get_data_sisj(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=2.5,dim=dim)
    df_sisj_Jz30 = get_data_sisj(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=3.0,dim=dim)
    df_sisj_Jz35 = get_data_sisj(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=3.5,dim=dim)
    #
    #
    fig = plt.figure(figsize=(20,10))
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
    ax1 = plt.axes([0.1,0.58,0.15,0.4])
    ax2 = plt.axes([0.30,0.58,0.15,0.4])
    ax5 = plt.axes([0.1,0.1,0.15,0.4])
    ax6 = plt.axes([0.30,0.1,0.15,0.4])

    ax3 = plt.axes([0.50,0.58,0.15,0.4])   
    ax4 = plt.axes([0.70,0.58,0.15,0.4])

    ax7 = plt.axes([0.50,0.1,0.15,0.4])
    ax8 = plt.axes([0.70,0.1,0.15,0.4])
    # plt.subplots_adjust(left=0.1, bottom=0.1, right=0.98, top=0.96, wspace=0.32, hspace=0.26)
    #--------------------------------------------------------------------------------------------
    # 选择子图 ax1 进行绘图
    plt.sca(ax1) ## 选择对 ax1 进行绘图
    ax1 = plt.gca()
    
    df_ni_Jz05 = get_data_ni(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=0.5,dim=dim)
    df_ni_Jz15 = get_data_ni(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=1.5,dim=dim)
    df_ni_Jz25 = get_data_ni(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=2.5,dim=dim)
    
    label = r"$J_{\bot}$=0.5"
    L05, =ax1.plot(df_ni_Jz05["rmean"].values,df_ni_Jz05["dymean"].values,label=label,ls="--",lw=1.5,color="red",
             marker='o',alpha=1,markersize=9,markeredgewidth=1.5, markeredgecolor="red", markerfacecolor='None')
    label = r"$J_{\bot}$=1.5"
    L15, =ax1.plot(df_ni_Jz15["rmean"].values,df_ni_Jz15["dymean"].values,label=label,ls="--",lw=1.5,color="blue",
             marker='^',alpha=1,markersize=9,markeredgewidth=1.5, markeredgecolor="blue", markerfacecolor='None')
    label = r"$J_{\bot}$=2.5"
    L25, =ax1.plot(df_ni_Jz25["rmean"].values,df_ni_Jz25["dymean"].values,label=label,ls="--",lw=1.5,color="green",
             marker='v',alpha=1,markersize=9,markeredgewidth=1.5, markeredgecolor="green", markerfacecolor='None')
    ####图例设置
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 15, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[L05,L15,L25,], loc = 4, bbox_to_anchor=(0.75, 0.01),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)    
    #
    label_x = r"x"
    label_y = "n(x)"
    #plt.yscale("log")
    #plt.xscale("log")
    ax1.set_xlabel(label_x, size= 16)
    ax1.set_ylabel(label_y, size= 16)
    ax1.tick_params(labelsize = 15) # 设置坐标刻度对应数字的大小
    ax1.set_xlim([0,48])
    ax1.set_xticks([0,10,20,30,40,48])
    #ax1.set_yticks([-1,-0.5,0,0.5,1]) 
    #=========================================================
    ax1.text(2,0.5027,"(a)",fontsize = 20, color='black', rotation = 0)
    ax1.text(20,0.5025, r'$\mathrm{\delta} = 0.5$', fontsize = 18, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
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
    #
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
    #-------------------------------------------------------------------------------
    # 选择子图 ax2 进行绘图
    plt.sca(ax2) ## 选择对 ax1 进行绘图
    ax2 = plt.gca()
    #--------- J_\bot = 0.25
    slope = df_sisj_Jz025["slope_pow"].values[0]
    K_F = round(-slope,2)  
    label = r"$J_{\bot}$=%.2f $K_F$=%.2f"%(0.25,K_F)
    L025, = ax2.plot(df_sisj_Jz025["r"],df_sisj_Jz025["corre_abs"],label=label,ls="-",lw=1.5,color="r",
             marker='o',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="r", markerfacecolor='None')
    L025_, = ax2.plot(df_sisj_Jz025["r"],df_sisj_Jz025["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='o',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- J_\bot = 0.75
    slope = df_sisj_Jz075["slope_pow"].values[0]
    K_F = round(-slope,2)  
    label = r"$J_{\bot}$=%.2f $K_F$=%.2f"%(0.75,K_F)
    L075, = ax2.plot(df_sisj_Jz075["r"],df_sisj_Jz075["corre_abs"],label=label,ls="-",lw=1.5,color="blue",
             marker='s',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="blue", markerfacecolor='None')
    L075_, = ax2.plot(df_sisj_Jz075["r"],df_sisj_Jz075["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 1.75
    slope = df_sisj_Jz175["slope_pow"].values[0]
    K_F = round(-slope,2)  
    label = r"$J_{\bot}$=%.2f $K_F$=%.2f"%(1.75,K_F)
    L175, = ax2.plot(df_sisj_Jz175["r"],df_sisj_Jz175["corre_abs"],label=label,ls="-",lw=1.5,color="green",
             marker='^',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="green", markerfacecolor='None')
    L175_, = ax2.plot(df_sisj_Jz175["r"],df_sisj_Jz175["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='^',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- J_\bot = 2.5
    slope = df_sisj_Jz25["slope_pow"].values[0]
    K_F = round(-slope,2)  
    label = r"$J_{\bot}$=%.1f"%(2.5,)
    L25, = ax2.plot(df_sisj_Jz25["r"],df_sisj_Jz25["corre_abs"],label=label,ls="-",lw=1.5,color="magenta",
             marker='v',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="magenta", markerfacecolor='None')
    #L25_, = ax2.plot(df_sisj_Jz25["r"],df_sisj_Jz25["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
    #         marker='v',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 3.0
    slope = df_sisj_Jz30["slope_pow"].values[0]
    K_F = round(-slope,2)  
    label = r"$J_{\bot}$=%.1f"%(3.0,)
    L30, = ax2.plot(df_sisj_Jz30["r"],df_sisj_Jz30["corre_abs"],label=label,ls="-",lw=1.5,color="cyan",
             marker='>',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="cyan", markerfacecolor='None')
    #L30_, = ax2.plot(df_sisj_Jz30["r"],df_sisj_Jz30["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
    #         marker='>',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 3.5
    #slope = df_sisj_Jz35["slope_pow"].values[0]
    #K_F = round(-slope,2)  
    #label = r"$J_{\bot}$=%.2f"%(3.5)
    #L35, = ax2.plot(df_sisj_Jz35["r"],df_sisj_Jz35["corre_abs"],label=label,ls="-",lw=1.5,color="magenta",
    #         marker='<',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="magenta", markerfacecolor='None')
    #L35_, = ax2.plot(df_sisj_Jz35["r"],df_sisj_Jz35["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
    #         marker='<',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #
    ####图例设置
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 15, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[L025,L075,L175,L25,L30,], loc = 4, bbox_to_anchor=(0.80, -0.01),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    #legend2=plt.legend(handles=[L25,L30,L35], loc = 4, bbox_to_anchor=(0.38, 0.02),
    #                   ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    #plt.gca().add_artist(legend1)
    #
    label_x = r"r"
    label_y = "|F(r)|"
    #plt.yscale("log")
    #plt.xscale("log")
    ax2.set_xscale("log",base=10,subs=[0.01,0.02,0.03]) 
    ax2.set_yscale("log",base=10,subs=[1])      
    ax2.set_xlabel(label_x, size= 16)
    ax2.set_ylabel(label_y, size= 16)
    ax2.tick_params(labelsize = 14) # 设置坐标刻度对应数字的大小
    ax2.set_xlim([0,40])
    #ax2.set_ylim([1e-6,1e-1])
    ax2.set_xticks([5,10,20,30])
    #ax2.set_yticks([-1,-0.5,0,0.5,1])  
    #=========================================================
    ax2.text(23,0.027,"(b)",fontsize = 20, color='black', rotation = 0)
    ax2.text(5,0.019, r'$\mathrm{\delta} = 0.5$', fontsize = 18, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
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
    plt.sca(ax3) ## 选择对 ax3 进行绘图
    ax3 = plt.gca()
    #
    df_cicj_Jz025 = get_data_cicj(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=0.25,dim=dim)
    df_cicj_Jz075 = get_data_cicj(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=0.75,dim=dim)
    df_cicj_Jz175 = get_data_cicj(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=1.75,dim=dim)
    df_cicj_Jz25 = get_data_cicj(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=2.5,dim=dim)
    df_cicj_Jz30 = get_data_cicj(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=3.0,dim=dim)
    df_cicj_Jz35 = get_data_cicj(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=3.5,dim=dim) 
    #--------- J_\bot = 0.25
    #slope = df_cicj_Jz025["slope_exp"].values[0]
    #K_G = round(-1/slope,2)  
    #label = r"$J_{\bot}$=%.2f $\xi_G$=%.2f"%(0.25,K_G)
    #L025, = ax3.plot(df_cicj_Jz025["r"],df_cicj_Jz025["corre_abs"],label=label,ls="-",lw=1.5,color="r",
    #         marker='o',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="r", markerfacecolor='None')
    #L025_, = ax3.plot(df_cicj_Jz025["r"],df_cicj_Jz025["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
    #         marker='o',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- J_\bot = 0.75
    slope = df_cicj_Jz075["slope_exp"].values[0]
    K_G = round(-1/slope,2)  
    label = r"$J_{\bot}$=%.2f $\xi_G$=%.2f"%(0.75,K_G)
    L075, = ax3.plot(df_cicj_Jz075["r"],df_cicj_Jz075["corre_abs"],label=label,ls="-",lw=1.5,color="blue",
             marker='s',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="blue", markerfacecolor='None')
    L075_, = ax3.plot(df_cicj_Jz075["r"],df_cicj_Jz075["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
             marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 1.75
    slope = df_cicj_Jz175["slope_exp"].values[0]
    K_G = round(-1/slope,2)  
    label = r"$J_{\bot}$=%.2f $\xi_G$=%.2f"%(1.75,K_G)
    L175, = ax3.plot(df_cicj_Jz175["r"],df_cicj_Jz175["corre_abs"],label=label,ls="-",lw=1.5,color="green",
             marker='^',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="green", markerfacecolor='None')
    L175_, = ax3.plot(df_cicj_Jz175["r"],df_cicj_Jz175["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
             marker='^',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- J_\bot = 2.5
    slope = df_cicj_Jz25["slope_exp"].values[0]
    K_G = round(-1/slope,2)  
    label = r"$J_{\bot}$=%.1f $\xi_G$=%.2f"%(2.5,K_G)
    L25, = ax3.plot(df_cicj_Jz25["r"],df_cicj_Jz25["corre_abs"],label=label,ls="-",lw=1.5,color="magenta",
             marker='v',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="magenta", markerfacecolor='None')
    L25_, = ax3.plot(df_cicj_Jz25["r"],df_cicj_Jz25["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
             marker='v',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 3.0
    slope = df_cicj_Jz30["slope_exp"].values[0]
    K_G = round(-1/slope,2)  
    label = r"$J_{\bot}$=%.1f $\xi_G$=%.2f"%(3.0,K_G)
    L30, = ax3.plot(df_cicj_Jz30["r"],df_cicj_Jz30["corre_abs"],label=label,ls="-",lw=1.5,color="cyan",
             marker='>',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="cyan", markerfacecolor='None')
    L30_, = ax3.plot(df_cicj_Jz30["r"],df_cicj_Jz30["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
             marker='>',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 3.5
    #slope = df_cicj_Jz35["slope_exp"].values[0]
    #K_G = round(-1/slope,2)  
    #label = r"$J_{\bot}$=%.2f $\xi_G$=%.2f"%(3.5,K_G)
    #L35, = ax3.plot(df_cicj_Jz35["r"],df_cicj_Jz35["corre_abs"],label=label,ls="-",lw=1.5,color="magenta",
    #         marker='h',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="magenta", markerfacecolor='None')
    #L35_, = ax3.plot(df_cicj_Jz35["r"],df_cicj_Jz35["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
    #         marker='h',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    ####图例设置
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 14, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[L075,L175,L25,L30,], loc = 4, bbox_to_anchor=(1.05, 0.63),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    #legend2=plt.legend(handles=[L075,L175,L25,L30,L35], loc = 4, bbox_to_anchor=(0.48, -0.02),
    #                   ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    #plt.gca().add_artist(legend1)
    label_x = r"r"
    label_y = "|G(r)|"
    #plt.yscale("log")
    #plt.xscale("log")
    #ax3.set_xscale("log",base=10,subs=[0.01,0.02,0.03]) 
    ax3.set_yscale("log",base=10,subs=[1])      
    ax3.set_xlabel(label_x, size= 16)
    ax3.set_ylabel(label_y, size= 16)
    ax3.tick_params(labelsize = 14) # 设置坐标刻度对应数字的大小
    ax3.set_xlim([0,37])
    #ax3.set_ylim([0.5e-6,5e-1])
    ax3.set_xticks([0,5,10,15,20,25,30,35])
    #ax3.set_yticks([-1,-0.5,0,0.5,1])  
    #=========================================================
    ax3.text(1,0.5e-5,"(c)",fontsize = 20, color='black', rotation = 0)
    ax3.text(1,0.3e-4, r'$\mathrm{\delta} = 0.5$', fontsize = 18, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
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
    #
    #----------------------------------------------------------------------------
    # 选择子图 ax4 进行绘图
    plt.sca(ax4) ## 选择对 ax4 进行绘图
    ax4 = plt.gca()
    #
    df_ninj_Jz025 = get_data_ninj_corre(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=0.25,dim=dim)
    df_ninj_Jz075 = get_data_ninj_corre(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=0.75,dim=dim)
    df_ninj_Jz175 = get_data_ninj_corre(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=1.75,dim=dim)
    # load spin correlation data
    df_ninj_Jz25 = get_data_ninj_corre(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=2.5,dim=dim)
    df_ninj_Jz30 = get_data_ninj_corre(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=3.0,dim=dim)
    df_ninj_Jz35 = get_data_ninj_corre(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=3.5,dim=dim)
    #--------- J_\bot = 0.25
    #slope = df_ninj_Jz025["slope_exp"].values[0]
    #xi = round(-1/slope,2)  
    #label = r"$J_{\bot}$=%.2f $\xi_{D}$=%.2f"%(0.25,xi)
    #L025, = ax4.plot(df_ninj_Jz025["r"],df_ninj_Jz025["corre_abs"],label=label,ls="-",lw=1.5,color="r",
    #         marker='o',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="r", markerfacecolor='None')
    #L025_, = ax4.plot(df_ninj_Jz025["r"],df_ninj_Jz025["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
    #         marker='o',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="r", markerfacecolor='None')
    #--------- J_\bot = 0.75
    slope = df_ninj_Jz075["slope_exp"].values[0]
    xi = round(-1/slope,2)  
    label = r"$J_{\bot}$=%.2f $\xi_{D}$=%.2f"%(0.75,xi)
    L075, = ax4.plot(df_ninj_Jz075["r"],df_ninj_Jz075["corre_abs"],label=label,ls="-",lw=1.5,color="blue",
             marker='s',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="blue", markerfacecolor='None')
    L075_, = ax4.plot(df_ninj_Jz075["r"],df_ninj_Jz075["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
             marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 1.75
    slope = df_ninj_Jz175["slope_exp"].values[0]
    xi = round(-1/slope,2)  
    label = r"$J_{\bot}$=%.2f $\xi_{D}$=%.2f"%(1.75,xi)
    L175, = ax4.plot(df_ninj_Jz175["r"],df_ninj_Jz175["corre_abs"],label=label,ls="-",lw=1.5,color="green",
             marker='^',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="green", markerfacecolor='None')
    L175_, = ax4.plot(df_ninj_Jz175["r"],df_ninj_Jz175["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
             marker='^',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- J_\bot = 2.5
    slope = df_ninj_Jz25["slope_exp"].values[0]
    xi = round(-1/slope,2)  
    label = r"$J_{\bot}$=%.1f $\xi_{D}$=%.2f"%(2.5,xi)
    L25, = ax4.plot(df_ninj_Jz25["r"],df_ninj_Jz25["corre_abs"],label=label,ls="-",lw=1.5,color="magenta",
             marker='v',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="magenta", markerfacecolor='None')
    L25_, = ax4.plot(df_ninj_Jz25["r"],df_ninj_Jz25["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
             marker='v',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 3.0
    slope = df_ninj_Jz30["slope_exp"].values[0]
    xi = round(-1/slope,2)  
    label = r"$J_{\bot}$=%.1f $\xi_{D}$=%.2f"%(3.0,xi)
    L30, = ax4.plot(df_ninj_Jz30["r"],df_ninj_Jz30["corre_abs"],label=label,ls="-",lw=1.5,color="cyan",
             marker='>',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="cyan", markerfacecolor='None')
    L30_, = ax4.plot(df_ninj_Jz30["r"],df_ninj_Jz30["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
             marker='>',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 3.5
    #slope = df_ninj_Jz35["slope_exp"].values[0]
    #xi = round(-1/slope,2)  
    #label = r"$J_{\bot}$=%.2f $\xi_{D}$=%.2f"%(3.5,xi)
    #L35, = ax4.plot(df_ninj_Jz35["r"],df_ninj_Jz35["corre_abs"],label=label,ls="-",lw=1.5,color="magenta",
    #         marker='h',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="magenta", markerfacecolor='None')
    #L35_, = ax4.plot(df_ninj_Jz35["r"],df_ninj_Jz35["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
    #         marker='h',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    ####图例设置
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 15, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[L075,L175,L25,L30,], loc = 4, bbox_to_anchor=(1.05, 0.62),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    #legend2=plt.legend(handles=[L35], loc = 4, bbox_to_anchor=(0.80, -0.02),
    #                   ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    #plt.gca().add_artist(legend1)
    label_x = r"r"
    label_y = "|D(r)|"
    #plt.yscale("log")
    #plt.xscale("log")
    ax4.set_yscale("log",base=10,subs=[1])      
    ax4.set_xlabel(label_x, size= 16)
    ax4.set_ylabel(label_y, size= 16)
    ax4.tick_params(labelsize = 14) # 设置坐标刻度对应数字的大小
    ax4.set_xlim([0,37])
    ax4.set_xticks([0,5,10,15,20,25,30,35])
    #ax4.set_yticks([-1,-0.5,0,0.5,1])  
    #=========================================================
    ax4.text(1,0.22e-8,"(d)",fontsize = 20, color='black', rotation = 0)
    ax4.text(1,0.4e-7, r'$\mathrm{\delta} = 0.5$', fontsize = 18, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
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
    #
    #----------------------------------------------------------------------------
    # 选择子图 ax5 进行绘图
    plt.sca(ax5) ## 选择对 ax5 进行绘图
    ax5 = plt.gca()
    
    df_pyy_Jz025 = get_data_pyy(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=0.25,dim=dim)
    df_pyy_Jz075 = get_data_pyy(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=0.75,dim=dim)
    df_pyy_Jz175 = get_data_pyy(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=1.75,dim=dim)
    df_pyy_Jz25 = get_data_pyy(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=2.5,dim=dim)
    df_pyy_Jz30 = get_data_pyy(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=3.0,dim=dim)
    df_pyy_Jz35 = get_data_pyy(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=3.5,dim=dim)
    #--------- J_\bot = 0.25
    #slope = df_pyy_Jz025["slope_exp"].values[0]
    #xi = round(-1/slope,2)  
    #label = r"$J_{\bot}$=%.2f $\xi_{SC}^{yy}$=%.2f"%(0.25,xi)
    #L025, = ax5.plot(df_pyy_Jz025["r"],df_pyy_Jz025["corre_abs"],label=label,ls="-",lw=1.5,color="r",
    #         marker='o',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="r", markerfacecolor='w')
    #L025_, = ax5.plot(df_pyy_Jz025["r"],df_pyy_Jz025["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
    #         marker='o',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="r", markerfacecolor='w')
    #--------- J_\bot = 0.75
    slope = df_pyy_Jz075["slope_exp"].values[0]
    xi = round(-1/slope,2)  
    label = r"$J_{\bot}$=%.2f $\xi_{SC}^{yy}$=%.2f"%(0.75,xi)
    L075, = ax5.plot(df_pyy_Jz075["r"],df_pyy_Jz075["corre_abs"],label=label,ls="-",lw=1.5,color="blue",
             marker='s',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="blue", markerfacecolor='None')
    L075_, = ax5.plot(df_pyy_Jz075["r"],df_pyy_Jz075["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
             marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 1.75
    slope = df_pyy_Jz175["slope_exp"].values[0]
    xi = round(-1/slope,2)  
    label = r"$J_{\bot}$=%.2f $\xi_{SC}^{yy}$=%.2f"%(1.75,xi)
    L175, = ax5.plot(df_pyy_Jz175["r"],df_pyy_Jz175["corre_abs"],label=label,ls="-",lw=1.5,color="green",
             marker='^',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="green", markerfacecolor='None')
    L175_, = ax5.plot(df_pyy_Jz175["r"],df_pyy_Jz175["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
             marker='^',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- J_\bot = 2.5
    slope = df_pyy_Jz25["slope_exp"].values[0]
    xi = round(-1/slope,2)  
    label = r"$J_{\bot}$=%.1f $\xi_{SC}^{yy}$=%.2f"%(2.5,xi)
    L25, = ax5.plot(df_pyy_Jz25["r"],df_pyy_Jz25["corre_abs"],label=label,ls="-",lw=1.5,color="magenta",
             marker='v',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="magenta", markerfacecolor='None')
    L25_, = ax5.plot(df_pyy_Jz25["r"],df_pyy_Jz25["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
             marker='v',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 3.0
    slope = df_pyy_Jz30["slope_exp"].values[0]
    xi = round(-1/slope,2)  
    label = r"$J_{\bot}$=%.1f $\xi_{SC}^{yy}$=%.2f"%(3.0,xi)
    L30, = ax5.plot(df_pyy_Jz30["r"],df_pyy_Jz30["corre_abs"],label=label,ls="-",lw=1.5,color="cyan",
             marker='>',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="cyan", markerfacecolor='None')
    L30_, = ax5.plot(df_pyy_Jz30["r"],df_pyy_Jz30["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
             marker='>',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 3.5
    #slope = df_pyy_Jz35["slope_exp"].values[0]
    #xi = round(-1/slope,2)  
    #label = r"$J_{\bot}$=%.2f $\xi_{SC}^{yy}$=%.2f"%(3.5,xi)
    #L35, = ax5.plot(df_pyy_Jz35["r"],df_pyy_Jz35["corre_abs"],label=label,ls="-",lw=1.5,color="magenta",
    #         marker='h',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="magenta", markerfacecolor='w')
    #L35_, = ax5.plot(df_pyy_Jz35["r"],df_pyy_Jz35["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
    #         marker='h',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    ####图例设置
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 15, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[L075,L175,], loc = 4, bbox_to_anchor=(1.05, 0.75),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    legend2=plt.legend(handles=[L25,L30], loc = 4, bbox_to_anchor=(0.82, -0.02),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    plt.gca().add_artist(legend1)
    label_x = r"r"
    label_y = "$P^{yy}(r)$"
    #plt.yscale("log")
    #plt.xscale("log")
    ax5.set_yscale("log",base=10,subs=[1])      
    ax5.set_xlabel(label_x, size= 16)
    ax5.set_ylabel(label_y, size= 16)
    ax5.tick_params(labelsize = 14) # 设置坐标刻度对应数字的大小
    ax5.set_xlim([1,37])
    ax5.set_xticks([5,10,15,20,25,30,35])
    #ax5.set_yticks([-1,-0.5,0,0.5,1])  
    #=========================================================
    ax5.text(2.8,0.004,"(e)",fontsize = 20, color='black', rotation = 0)
    ax5.text(2,1e-6, r'$\mathrm{\delta} = 0.5$', fontsize = 18, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
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
    #ax5.yaxis.get_minor_locator().set_params(numticks=99, subs=[.2,.4,.6,.8]) # 将次要刻度显示出来 
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
    #
    #----------------------------------------------------------------------------
    # 选择子图 ax6 进行绘图
    plt.sca(ax6) ## 选择对 ax6 进行绘图
    ax6 = plt.gca()
    
    df_pzz_Jz025 = get_data_pzz(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=0.25,dim=dim)
    df_pzz_Jz05 = get_data_pzz(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=0.5,dim=dim)
    df_pzz_Jz075 = get_data_pzz(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=0.75,dim=dim)
    df_pzz_Jz10 = get_data_pzz(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=1.0,dim=dim)
    df_pzz_Jz125 = get_data_pzz(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=1.25,dim=dim)
    df_pzz_Jz15 = get_data_pzz(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=1.5,dim=dim)
    df_pzz_Jz175 = get_data_pzz(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=1.75,dim=dim)
    # load spin correlation data
    df_pzz_Jz20 = get_data_pzz(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=2.0,dim=dim)
    df_pzz_Jz225 = get_data_pzz(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=2.25,dim=dim)
    df_pzz_Jz25 = get_data_pzz(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=2.5,dim=dim)
    df_pzz_Jz275 = get_data_pzz(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=2.75,dim=dim)
    df_pzz_Jz30 = get_data_pzz(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=3.0,dim=dim)
    df_pzz_Jz325 = get_data_pzz(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=3.25,dim=dim)
    df_pzz_Jz35 = get_data_pzz(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=3.5,dim=dim)
    
    #--------- J_\bot = 0.25
    #slope = df_pzz_Jz025["slope_exp"].values[0]
    #xi = round(-1/slope,2)  
    #label = r"$J_{\bot}$=%.2f $\xi_{SC}^{zz}$=%.2f"%(0.25,xi)
    #L025, = ax6.plot(df_pzz_Jz025["r"],df_pzz_Jz025["corre_abs"],label=label,ls="-",lw=1.5,color="r",
    #         marker='o',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="r", markerfacecolor='None')
    #L025_, = ax6.plot(df_pzz_Jz025["r"],df_pzz_Jz025["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
    #         marker='o',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="r", markerfacecolor='w')
    #--------- J_\bot = 0.75
    slope = df_pzz_Jz075["slope_exp"].values[0]
    xi = round(-1/slope,2)  
    label = r"$J_{\bot}$=%.2f $\xi_{SC}^{zz}$=%.2f"%(0.75,xi)
    L075, = ax6.plot(df_pzz_Jz075["r"],df_pzz_Jz075["corre_abs"],label=label,ls="-",lw=1.5,color="blue",
             marker='s',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="blue", markerfacecolor='None')
    L075_, = ax6.plot(df_pzz_Jz075["r"],df_pzz_Jz075["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
             marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 1.75
    slope = df_pzz_Jz175["slope_exp"].values[0]
    xi = round(-1/slope,2)  
    label = r"$J_{\bot}$=%.2f $\xi_{SC}^{zz}$=%.2f"%(1.75,xi)
    L175, = ax6.plot(df_pzz_Jz175["r"],df_pzz_Jz175["corre_abs"],label=label,ls="-",lw=1.5,color="green",
             marker='^',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="green", markerfacecolor='None')
    L175_, = ax6.plot(df_pzz_Jz175["r"],df_pzz_Jz175["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
             marker='^',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- J_\bot = 2.5
    slope = df_pzz_Jz25["slope_exp"].values[0]
    xi = round(-1/slope,2)  
    label = r"$J_{\bot}$=%.1f $\xi_{SC}^{zz}$=%.2f"%(2.5,xi)
    L25, = ax6.plot(df_pzz_Jz25["r"],df_pzz_Jz25["corre_abs"],label=label,ls="-",lw=1.5,color="magenta",
             marker='v',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="magenta", markerfacecolor='None')
    L25_, = ax6.plot(df_pzz_Jz25["r"],df_pzz_Jz25["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
             marker='v',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 3.0
    slope = df_pzz_Jz30["slope_exp"].values[0]
    xi = round(-1/slope,2)  
    label = r"$J_{\bot}$=%.1f $\xi_{SC}^{zz}$=%.2f"%(3.0,xi)
    L30, = ax6.plot(df_pzz_Jz30["r"],df_pzz_Jz30["corre_abs"],label=label,ls="-",lw=1.5,color="cyan",
             marker='>',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="cyan", markerfacecolor='None')
    L30_, = ax6.plot(df_pzz_Jz30["r"],df_pzz_Jz30["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
             marker='>',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 3.5
    #slope = df_pzz_Jz35["slope_exp"].values[0]
    #xi = round(-1/slope,2)  
    #label = r"$J_{\bot}$=%.2f $\xi_{SC}^{zz}$=%.2f"%(3.5,xi)
    #L35, = ax6.plot(df_pzz_Jz35["r"],df_pzz_Jz35["corre_abs"],label=label,ls="-",lw=1.5,color="magenta",
    #         marker='h',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="magenta", markerfacecolor='None')
    #L35_, = ax6.plot(df_pzz_Jz35["r"],df_pzz_Jz35["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
    #         marker='h',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    ####图例设置
    ####图例设置
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 14, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[L075,L175,], loc = 4, bbox_to_anchor=(1.05, 0.74),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    legend2=plt.legend(handles=[L25,L30,], loc = 4, bbox_to_anchor=(0.78, -0.02),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    plt.gca().add_artist(legend1)
    
    label_x = r"r"
    label_y = "$P^{zz}(r)$"
    #plt.yscale("log")
    #plt.xscale("log")
    ax6.set_yscale("log",base=10,subs=[1])      
    ax6.set_xlabel(label_x, size= 16)
    ax6.set_ylabel(label_y, size= 16)
    ax6.tick_params(labelsize = 14) # 设置坐标刻度对应数字的大小
    ax6.set_xlim([2,37])
    ax6.set_xticks([5,10,15,20,25,30,35])
    #ax6.set_yticks([-1,-0.5,0,0.5,1])  
    #=========================================================
    ax6.text(32,0.008,"(f)",fontsize = 20, color='black', rotation = 0)
    ax6.text(3.5,3e-6, r'$\mathrm{\delta} = 0.5$', fontsize = 18, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
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
    #ax6.yaxis.get_minor_locator().set_params(numticks=99, subs=[.2,.4,.6,.8]) # 将次要刻度显示出来 
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
    #----------------------------------------------------------------------------
    # 选择子图 ax7 进行绘图
    plt.sca(ax7) ## 选择对 ax7 进行绘图
    ax7 = plt.gca()
    #--------- J_\bot = 0.25
    #slope = df_pzz_Jz025["slope_pow"].values[0]
    #K_F = round(-slope,2)  
    #label = r"$J_{\bot}$=%.2f $K_{SC}^{zz}$=%.2f"%(0.25,K_F)
    #L025, = ax7.plot(df_pzz_Jz025["r"],df_pzz_Jz025["corre_abs"],label=label,ls="-",lw=1.5,color="r",
    #         marker='o',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="r", markerfacecolor='None')
    #L025_, = ax7.plot(df_pzz_Jz025["r"],df_pzz_Jz025["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
    #         marker='o',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- J_\bot = 0.75
    slope = df_pzz_Jz075["slope_pow"].values[0]
    K_F = round(-slope,2)  
    label = r"$J_{\bot}$=%.2f"%(0.75)
    L075, = ax7.plot(df_pzz_Jz075["r"],df_pzz_Jz075["corre_abs"],label=label,ls="-",lw=1.5,color="blue",
             marker='s',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="blue", markerfacecolor='None')
    #L075_, = ax7.plot(df_pzz_Jz075["r"],df_pzz_Jz075["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
    #         marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 1.75
    slope = df_pzz_Jz175["slope_pow"].values[0]
    K_F = round(-slope,2)  
    label = r"$J_{\bot}$=%.2f"%(1.75,)
    L175, = ax7.plot(df_pzz_Jz175["r"],df_pzz_Jz175["corre_abs"],label=label,ls="-",lw=1.5,color="green",
             marker='^',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="green", markerfacecolor='None')
    #L175_, = ax7.plot(df_pzz_Jz175["r"],df_pzz_Jz175["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
    #      marker='^',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- J_\bot = 2.5
    slope = df_pzz_Jz25["slope_pow"].values[0]
    K_F = round(-slope,2)  
    label = r"$J_{\bot}$=%.1f"%(2.5,)
    L25, = ax7.plot(df_pzz_Jz25["r"],df_pzz_Jz25["corre_abs"],label=label,ls="-",lw=1.5,color="magenta",
             marker='v',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="magenta", markerfacecolor='None')
    #L25_, = ax7.plot(df_pzz_Jz25["r"],df_pzz_Jz25["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
    #         marker='v',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 3.0
    slope = df_pzz_Jz30["slope_pow"].values[0]
    K_F = round(-slope,2)  
    label = r"$J_{\bot}$=%.1f"%(3.0,)
    L30, = ax7.plot(df_pzz_Jz30["r"],df_pzz_Jz30["corre_abs"],label=label,ls="-",lw=1.5,color="cyan",
             marker='>',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="cyan", markerfacecolor='None')
    #L30_, = ax7.plot(df_pzz_Jz30["r"],df_pzz_Jz30["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
    #         marker='>',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 3.5
    #slope = df_pzz_Jz35["slope_pow"].values[0]
    #K_F = round(-slope,2)  
    #label = r"$J_{\bot}$=%.2f $K_{SC}^{zz}$=%.2f"%(3.5,K_F)
    #L35, = ax7.plot(df_pzz_Jz35["r"],df_pzz_Jz35["corre_abs"],label=label,ls="-",lw=1.5,color="magenta",
    #         marker='h',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="magenta", markerfacecolor='None')
    #L35_, = ax7.plot(df_pzz_Jz35["r"],df_pzz_Jz35["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
    #         marker='h',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    ####图例设置
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 14, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    #legend1=plt.legend(handles=[], loc = 4, bbox_to_anchor=(1.05, 0.88),
    #                   ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    legend2=plt.legend(handles=[L075,L175,L25,L30,], loc = 4, bbox_to_anchor=(0.55, -0.03),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    #plt.gca().add_artist(legend1)
    label_x = r"r"
    label_y = "$P^{zz}(r)$"
    #plt.yscale("log")
    #plt.xscale("log")
    ax7.set_xscale("log",base=10,subs=[0.01,0.02,0.03]) 
    ax7.set_yscale("log",base=10,subs=[1])      
    ax7.set_xlabel(label_x, size= 16)
    ax7.set_ylabel(label_y, size= 16)
    ax7.tick_params(labelsize = 14) # 设置坐标刻度对应数字的大小
    ax7.set_xlim([0,40])
    ax7.set_xticks([5,10,20,30])
    #ax7.set_yticks([-1,-0.5,0,0.5,1])  
    #=========================================================
    ax7.text(25,0.008,"(g)",fontsize = 20, color='black', rotation = 0)
    ax7.text(2,0.000008, r'$\mathrm{\delta} = 0.5$', fontsize = 18, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
    #=========================================================
    # 坐标轴设置第一层
    labels = ax7.get_xticklabels() + ax7.get_yticklabels()
    #[label.set_fontname('Times New Roman') for label in labels]###设置ticket labled的字体格式
    ax7.xaxis.set_minor_locator(MultipleLocator(5))###设置次刻度的间隔
    #ax7.yaxis.set_minor_locator(MultipleLocator(10))###设置次刻度的间隔
    ax7.xaxis.set_major_formatter(FormatStrFormatter('%1.0f'))###设置X轴标签文本格式
    #ax7.yaxis.set_major_formatter(FormatStrFormatter('%1.1f'))###设置Y轴标签文本格式
    #
    #ax7.yaxis.get_major_locator().set_params(numticks=99)
    #ax7.yaxis.get_minor_locator().set_params(numticks=99, subs=[.2,.4,.6,.8]) # 将次要刻度显示出来 
    #----将次要刻度显示出来 
    #locmin = matplotlib.ticker.LogLocator(base=10.0, subs=(0.1,0.2,0.4,0.6,0.8,1,2,4,6,8,10 )) 
    #ax7.xaxis.set_minor_locator(locmin)
    #ax7.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    #=====坐标轴的第二层： 坐标轴的设置
    ax7.spines['bottom'].set_linewidth(1.5) ###设置底部坐标轴的粗细
    ax7.spines['left'].set_linewidth(1.5)   ###设置左边坐标轴的粗细
    ax7.spines['right'].set_linewidth(1.5)  ###设置右边坐标轴的粗细
    ax7.spines['top'].set_linewidth(1.5)    ###设置上部坐标轴的粗细
    #ax7.spines['right'].set_color('none')# 将右边上边的两条边颜色设置为空 其实就相当于抹掉这两条边
    #ax7.spines['top'].set_color('none')
    #====坐标轴的第三层：  主刻度线的设置
    for line in ax7.xaxis.get_ticklines():
        #line is a Line2D instance
        #line.set_color('green')
        line.set_markersize(3)####设置刻度线的长度
        line.set_markeredgewidth(1.5)####设置刻度线的宽度
    for line in ax7.yaxis.get_ticklines():
        # line is a Line2D instance
        #line.set_color('green')
        line.set_markersize(3)####设置刻度线的长度
        line.set_markeredgewidth(1.5)####设置刻度线的宽度
    plt.tick_params(axis="x", which="minor", length=2.5, width=1.5, color="k")  ### 设置次要刻度 
    plt.tick_params(axis="y", which="minor", length=2.5, width=1.5, color="k")  ### 设置次要刻度 
    #--------------------------------------------------------------------------------------------
    # 选择子图 ax8 进行绘图
    plt.sca(ax8) ## 选择对 ax6 进行绘图
    ax8 = plt.gca()
    c_01 = 2.00
    g_01 = 2.55
    df_ent_01 = get_data_entropy(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=0.1,dim=6000)
    S_res01 = entropy(Lx=Lx,x=df_ent_01["r"].values,c=c_01,g=g_01)
    #
    #c_40 = 0.0
    #g_40 = 2.52
    #df_ent_40 = get_data_entropy(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=4.00,dim=6000)
    #S_res40 = entropy(Lx=Lx,x=df_ent_40["r"].values,c=c_40,g=g_40)
    #
    c_30 = 0.0
    g_30 = 2.70
    df_ent_30 = get_data_entropy(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=3.00,dim=6000)
    S_res30 = entropy(Lx=Lx,x=df_ent_30["r"].values,c=c_30,g=g_30)
    #
    L01, = ax8.plot(df_ent_01['r'].values,df_ent_01["entropy"],label=r"$J_{\bot}$=%.1f, c=%g"%(0.1,c_01),ls="-",lw=1.5,color='r',
             marker='o',alpha=1,markersize=8,markeredgewidth=1.5, markeredgecolor='r',
             markerfacecolor='None')
    L01_fit, = ax8.plot(df_ent_01['r'].values,S_res01,label="c=%.2f"%(c_01),ls="--",lw=1.5,color="k",
             marker='o',alpha=1,markersize=0,markeredgewidth=1.5, markeredgecolor="k",
             markerfacecolor='None')
    #
    #L40, = ax8.plot(df_ent_40['r'].values,df_ent_40["entropy"],label=r"$J_{\bot}$=%.2f"%(4.0),ls="-",lw=1.5,color='blue',
    #         marker='o',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor='blue',
    #         markerfacecolor='None')
    #L40_fit, = ax8.plot(df_ent_40['r'].values,S_res40,label="c=%.2f,g=%.2f"%(c_40,g_40),ls="--",lw=2.5,color="green",
    #         marker='o',alpha=1,markersize=0,markeredgewidth=1.5, markeredgecolor="k",
    #         markerfacecolor='None')
    #
    L30, = ax8.plot(df_ent_30['r'].values,df_ent_30["entropy"],label=r"$J_{\bot}$=%.1f, c=%g"%(3.0,c_30),ls="-",lw=1.5,color='magenta',
             marker='s',alpha=1,markersize=8,markeredgewidth=1.5, markeredgecolor='magenta',
             markerfacecolor='None')
    L30_fit, = ax8.plot(df_ent_30['r'].values,S_res30,label="c=%.2f"%(c_30),ls="--",lw=2.5,color="green",
             marker='s',alpha=1,markersize=0,markeredgewidth=1.5, markeredgecolor="k",
             markerfacecolor='None')
    #
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 16, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[L01,L30], loc = 4, bbox_to_anchor=(0.86, 0.01),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    #legend2=plt.legend(handles=[L30,L30_fit], loc = 4, bbox_to_anchor=(0.86, 0.10),
    #                   ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    #plt.gca().add_artist(legend1)#####把图例legend1重新加载回来

    label_x = r"x"
    label_y = "S(x)"
    #plt.yscale("log")
    ax8.set_xlabel(label_x, size= 16)
    ax8.set_ylabel(label_y, size= 16)
    ax8.tick_params(labelsize = 16) # 设置坐标刻度对应数字的大小
    ax8.set_xlim([0,47])
    #ax8.set_ylim([-0.1,1])
    ax8.set_xticks([0,10,20,30,40])
    #ax8.set_yticks([-0.1,0,0.5,1])
    #ax8.text(0.3,-0.05, r'$\mathrm{(a)}$', fontsize=18)
    ax8.text(2,3.435,"(h)",fontsize = 20, color='black', rotation = 0)
    ax8.text(10,2.4,"$L_x$=48",fontsize = 18, color='black', rotation = 0)
    ax8.text(25,2.4,"$\delta$=0.5",fontsize = 18, color='black', rotation = 0)
    #plt.title("Jz=%.2f"%Jz,fontsize=25)
    #
    #=======================================================================================
    # 坐标轴设置第一层
    labels = ax8.get_xticklabels() + ax8.get_yticklabels()
    #[label.set_fontname('Times New Roman') for label in labels]###设置ticket labled的字体格式
    ax8.xaxis.set_minor_locator(MultipleLocator(5))###设置次刻度的间隔
    #ax8.yaxis.set_minor_locator(MultipleLocator(10))###设置次刻度的间隔
    ax8.xaxis.set_major_formatter(FormatStrFormatter('%1.0f'))###设置X轴标签文本格式
    #ax8.yaxis.set_major_formatter(FormatStrFormatter('%1.1f'))###设置Y轴标签文本格式
    #
    #ax8.yaxis.get_major_locator().set_params(numticks=99)
    #ax8.yaxis.get_minor_locator().set_params(numticks=99, subs=[.2,.4,.6,.8]) # 将次要刻度显示出来 
    #----将次要刻度显示出来 
    #locmin = matplotlib.ticker.LogLocator(base=10.0, subs=(0.1,0.2,0.4,0.6,0.8,1,2,4,6,8,10 )) 
    #ax8.xaxis.set_minor_locator(locmin)
    #ax8.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    #
    #=====坐标轴的第二层： 坐标轴的设置
    ax8.spines['bottom'].set_linewidth(1.5) ###设置底部坐标轴的粗细
    ax8.spines['left'].set_linewidth(1.5)   ###设置左边坐标轴的粗细
    ax8.spines['right'].set_linewidth(1.5)  ###设置右边坐标轴的粗细
    ax8.spines['top'].set_linewidth(1.5)    ###设置上部坐标轴的粗细
    #ax8.spines['right'].set_color('none')# 将右边上边的两条边颜色设置为空 其实就相当于抹掉这两条边
    #ax8.spines['top'].set_color('none')
    #====坐标轴的第三层：  主刻度线的设置
    for line in ax8.xaxis.get_ticklines():
        #line is a Line2D instance
        #line.set_color('green')
        line.set_markersize(3)####设置刻度线的长度
        line.set_markeredgewidth(1.5)####设置刻度线的宽度
    for line in ax8.yaxis.get_ticklines():
        # line is a Line2D instance
        #line.set_color('green')
        line.set_markersize(3)####设置刻度线的长度
        line.set_markeredgewidth(1.5)####设置刻度线的宽度
    plt.tick_params(axis="x", which="minor", length=2.5, width=1.5, color="k")  ### 设置次要刻度 
    plt.tick_params(axis="y", which="minor", length=2.5, width=1.5, color="k")  ### 设置次要刻度 
    #
    fig.tight_layout() # 自动调整 subplot 间的间隙参数
    plt.savefig("E:\\WORK\\Work\\Project\\La3Ni2O7\\datashow_datpy_paper\\figs\\fig_05_1_v2.eps",
                dpi=300, format='eps',bbox_inches='tight') # 白边紧凑型
    #
    plt.show()