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
def get_entropy_Lz2(Lz, Ly, Lx, dop, t, J, Jz, dim):
    workpath = "E:\\WORK\\Work\\Project\\La3Ni2O7"
    filepath1 = "\\data_dmrgcpp\\Lz%d_Ly%d_Lx%d\\dop%g" % (Lz, Ly, Lx,dop)
    filepath2 = "\\t%d_J%d_Jz%.2f_dim%d_ent\\entanglement" % (t, J, Jz, dim)
    filepath3 = "\\measurement_entropy_fit.parquet"
    filename = workpath + filepath1 + filepath2 + filepath3
    print(filename)
    df = pd.read_parquet(filename)
    print(df.head())
    return df
#
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
def get_data_pyy_leg2(Lz, Ly, Lx, dop, t, J, Jz, dim,):  # single-layer 2-leg ladder
    workpath = "E:\\WORK\\Work\\Project\\leg2-tj"
    filepath1 = "\\data_dmrgcpp\\Lz%d_Ly%d_Lx%d\\dop%g" % (Lz, Ly, Lx,dop)
    filepath2 = "\\t%d_J%d_Jz%d_dim%d" % (t, J, Jz, dim)
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
    #
    #---------------plot logr-r fig and logr-logr fig-----------------------------
    fig = plt.figure(figsize=(16,10))
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
    #-----------------------------------------------------------------------------------------------
    plt.sca(ax1) ## 选择对 ax1 进行绘图 CDW profile
    ax1 = plt.gca()
    df_ni_Jz05 = get_data_ni(Lz=Lz,Ly=Ly,Lx=48,dop=88,t=3,J=1,Jz=0.5,dim=6000)
    df_ni_Jz10 = get_data_ni(Lz=Lz,Ly=Ly,Lx=48,dop=88,t=3,J=1,Jz=1.0,dim=6000)
    df_ni_Jz20 = get_data_ni(Lz=Lz,Ly=Ly,Lx=48,dop=88,t=3,J=1,Jz=2.0,dim=6000)
    df_ni_Jz30 = get_data_ni(Lz=Lz,Ly=Ly,Lx=48,dop=88,t=3,J=1,Jz=3.0,dim=6000)
    #
    L, =ax1.plot([0,50],[0.5417,0.5417],label=" ",ls="--",lw=1.5,color="k",
             marker='o',alpha=1,markersize=0,markeredgewidth=1.5, markeredgecolor="red", markerfacecolor='None')
    label = r"$J_{\bot}$=0.5"
    L05, =ax1.plot(df_ni_Jz05["rmean"].values,df_ni_Jz05["dymean"].values,label=label,ls="--",lw=1.5,color="red",
             marker='o',alpha=1,markersize=9,markeredgewidth=1.5, markeredgecolor="red", markerfacecolor='None')
    label = r"$J_{\bot}$=1.0"
    L10, =ax1.plot(df_ni_Jz10["rmean"].values,df_ni_Jz10["dymean"].values,label=label,ls="--",lw=1.5,color="blue",
             marker='s',alpha=1,markersize=9,markeredgewidth=1.5, markeredgecolor="blue", markerfacecolor='None')
    label = r"$J_{\bot}$=2.0"
    L20, =ax1.plot(df_ni_Jz20["rmean"].values,df_ni_Jz20["dymean"].values,label=label,ls="--",lw=1.5,color="green",
             marker='^',alpha=1,markersize=9,markeredgewidth=1.5, markeredgecolor="green", markerfacecolor='None')
    label = r"$J_{\bot}$=3.0"
    L30, =ax1.plot(df_ni_Jz30["rmean"].values,df_ni_Jz30["dymean"].values,label=label,ls="--",lw=1.5,color="magenta",
             marker='v',alpha=1,markersize=9,markeredgewidth=1.5, markeredgecolor="magenta", markerfacecolor='None')
    ####图例设置
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 15, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[L05,L10,L20,L30,], loc = 4, bbox_to_anchor=(0.92, 0.01),
                       ncol = 2,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)    
    #legend2=plt.legend(handles=[L04g,L10g,L20g], loc = 4, bbox_to_anchor=(0.64, 0.02),
    #                   ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    #plt.gca().add_artist(legend1)#####把图例legend1重新加载回来
    #
    label_x = r"x"
    label_y = "n(x)"
    #plt.yscale("log")
    #plt.xscale("log")
    ax1.set_xlabel(label_x, size= 16)
    ax1.set_ylabel(label_y, size= 16)
    ax1.tick_params(labelsize = 15) # 设置坐标刻度对应数字的大小
    ax1.set_xlim([0,48])
    #ax1.set_ylim([])
    ax1.set_xticks([0,10,20,30,40,48])
    #ax1.set_yticks([-1,-0.5,0,0.5,1]) 
    #=========================================================
    ax1.text(7,0.567,"(a)",fontsize = 20, color='black', rotation = 0)
    ax1.text(20,0.564, r'$\mathrm{\delta} = 0.4583$', fontsize = 16, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
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
    #-----------------------------------------------------------------------------------------------
    # 选择子图 ax2 进行绘图
    plt.sca(ax2) ## 选择对 ax4 进行绘图
    ax2 = plt.gca()
    #------------------------
    df_cicj_Jz01 = get_data_cicj(Lz=Lz,Ly=Ly,Lx=48,dop=88,t=3,J=1,Jz=0.1,dim=6000)
    df_cicj_Jz06 = get_data_cicj(Lz=Lz,Ly=Ly,Lx=48,dop=88,t=3,J=1,Jz=0.6,dim=6000)
    df_cicj_Jz10 = get_data_cicj(Lz=Lz,Ly=Ly,Lx=48,dop=88,t=3,J=1,Jz=1.0,dim=6000)
    df_cicj_Jz15 = get_data_cicj(Lz=Lz,Ly=Ly,Lx=48,dop=88,t=3,J=1,Jz=1.5,dim=6000)
    df_cicj_Jz175 = get_data_cicj(Lz=Lz,Ly=Ly,Lx=48,dop=88,t=3,J=1,Jz=1.75,dim=6000)
    df_cicj_Jz20 = get_data_cicj(Lz=Lz,Ly=Ly,Lx=48,dop=88,t=3,J=1,Jz=2.0,dim=6000)
    df_cicj_Jz225 = get_data_cicj(Lz=Lz,Ly=Ly,Lx=48,dop=88,t=3,J=1,Jz=2.25,dim=6000)
    df_cicj_Jz30 = get_data_cicj(Lz=Lz,Ly=Ly,Lx=48,dop=88,t=3,J=1,Jz=3.0,dim=6000)
    #--------- J_\bot = 0.1
    slope = df_cicj_Jz01["slope_pow"].values[0]
    Kg = round(-slope,2)  
    label = r"$J_{\bot}$=%.1f $K_G$=%.2f"%(0.1,Kg)
    L01, = ax2.plot(df_cicj_Jz01["r"],df_cicj_Jz01["corre_abs"],label=label,ls="-",lw=1.5,color="red",
             marker='o',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="red", markerfacecolor='None')
    L01_, = ax2.plot(df_cicj_Jz01["r"],df_cicj_Jz01["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='o',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- J_\bot = 0.6
    #slope = df_cicj_Jz06["slope_pow"].values[0]
    #Kg = round(-slope,2)  
    #label = r"$J_{\bot}$=%.1f $K_G$=%.2f"%(0.6,Kg)
    #L06, = ax2.plot(df_cicj_Jz06["r"],df_cicj_Jz06["corre_abs"],label=label,ls="-",lw=1.5,color="blue",
    #         marker='s',alpha=1,markersize=6,markeredgewidth=1, markeredgecolor="blue", markerfacecolor='None')
    #L06_, = ax2.plot(df_cicj_Jz06["r"],df_cicj_Jz06["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
    #         marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- J_\bot = 1.0
    slope = df_cicj_Jz10["slope_pow"].values[0]
    Kg = round(-slope,2)  
    label = r"$J_{\bot}$=%.1f $K_G$=%.2f"%(1.0,Kg)
    L10, = ax2.plot(df_cicj_Jz10["r"],df_cicj_Jz10["corre_abs"],label=label,ls="-",lw=1.5,color="blue",
             marker='s',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="blue", markerfacecolor='None')
    L10_, = ax2.plot(df_cicj_Jz10["r"],df_cicj_Jz10["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- J_\bot = 1.5
    slope = df_cicj_Jz15["slope_pow"].values[0]
    Kg = round(-slope,2)  
    label = r"$J_{\bot}$=%.1f $K_G$=%.2f"%(1.5,Kg)
    L15, = ax2.plot(df_cicj_Jz15["r"],df_cicj_Jz15["corre_abs"],label=label,ls="-",lw=1.5,color="green",
             marker='v',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="green", markerfacecolor='None')
    L15_, = ax2.plot(df_cicj_Jz15["r"],df_cicj_Jz15["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='v',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 1.75
    #slope = df_cicj_Jz175["slope_pow"].values[0]
    #Kg = round(-slope,2)  
    #label = r"$J_{\bot}$=%.2f $K_G$=%.2f"%(1.75,Kg)
    #L175, = ax2.plot(df_cicj_Jz175["r"],df_cicj_Jz175["corre_abs"],label=label,ls="-",lw=1.5,color="cyan",
    #         marker='<',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="cyan", markerfacecolor='None')
    #L175_, = ax2.plot(df_cicj_Jz175["r"],df_cicj_Jz175["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
    #         marker='<',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 2.0
    slope = df_cicj_Jz20["slope_pow"].values[0]
    Kg = round(-slope,2)  
    #label = r"$J_{\bot}$=%.1f $K_G$=%.2f"%(2.0,Kg)
    label = r"$J_{\bot}$=%.1f"%(2.0)
    L20, = ax2.plot(df_cicj_Jz20["r"],df_cicj_Jz20["corre_abs"],label=label,ls="-",lw=1.5,color="magenta",
             marker='>',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="magenta", markerfacecolor='None')
    #L20_, = ax2.plot(df_cicj_Jz20["r"],df_cicj_Jz20["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
    #         marker='>',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 2.25
    slope = df_cicj_Jz225["slope_pow"].values[0]
    Kg = round(-slope,2)  
    #label = r"$J_{\bot}$=%.2f $K_G$=%.2f"%(2.25,Kg)
    label = r"$J_{\bot}$=%.2f"%(2.25,)
    L225, = ax2.plot(df_cicj_Jz225["r"],df_cicj_Jz225["corre_abs"],label=label,ls="-",lw=1.5,color="cyan",
             marker='^',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="cyan", markerfacecolor='None')
    #L225_, = ax2.plot(df_cicj_Jz225["r"],df_cicj_Jz225["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
    #         marker='^',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 3.0
    slope = df_cicj_Jz30["slope_pow"].values[0]
    Kg = round(-slope,2)  
    #label = r"$J_{\bot}$=%.2f $K_G$=%.2f"%(3.0,Kg)
    label = r"$J_{\bot}$=%.1f"%(3.0)
    L30, = ax2.plot(df_cicj_Jz30["r"],df_cicj_Jz30["corre_abs"],label=label,ls="-",lw=1.5,color="brown",
             marker='<',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="brown", markerfacecolor='None')
    #L30_, = ax2.plot(df_cicj_Jz30["r"],df_cicj_Jz30["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
    #         marker='<',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------
    ####图例设置
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 15, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    #legend1=plt.legend(handles=[], loc = 4, bbox_to_anchor=(0.99, 0.72),
    #                   ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    legend2=plt.legend(handles=[L01,L10,L15,L20,L225,L30], loc = 4, bbox_to_anchor=(0.64, -0.01),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    #plt.gca().add_artist(legend1)#####把图例legend1重新加载回来
    label_x = r"r"
    label_y = "|G(r)|"
    #plt.yscale("log")
    #plt.xscale("log")
    ax2.set_yscale("log",base=10,subs=[1])   
    ax2.set_xscale("log",base=10,subs=[0.01,0.02,0.03])    
    ax2.set_xlabel(label_x, size= 14)
    ax2.set_ylabel(label_y, size= 14)
    ax2.tick_params(labelsize = 12) # 设置坐标刻度对应数字的大小
    ax2.set_xlim([1.1,40])
    ax2.set_xticks([1,5,10,20,30])
    #ax2.set_yticks([-1,-0.5,0,0.5,1])  
    #=========================================================
    ax2.text(27,0.22,"(b)",fontsize = 20, color='black', rotation = 0)
    ax2.text(5,0.11, r'$\mathrm{\delta} = 0.4583$', fontsize = 16, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
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
    #----------------------------------------------------------------------------------------------
    # 选择子图 ax3 进行绘图
    plt.sca(ax3) ## 选择对 ax3 进行绘图
    ax3 = plt.gca()
    df_sisj_Jz01 = get_data_sisj(Lz=Lz,Ly=Ly,Lx=48,dop=88,t=3,J=1,Jz=0.1,dim=6000)
    df_sisj_Jz06 = get_data_sisj(Lz=Lz,Ly=Ly,Lx=48,dop=88,t=3,J=1,Jz=0.6,dim=6000)
    df_sisj_Jz10 = get_data_sisj(Lz=Lz,Ly=Ly,Lx=48,dop=88,t=3,J=1,Jz=1.0,dim=6000)
    df_sisj_Jz15 = get_data_sisj(Lz=Lz,Ly=Ly,Lx=48,dop=88,t=3,J=1,Jz=1.5,dim=6000)
    df_sisj_Jz175 = get_data_sisj(Lz=Lz,Ly=Ly,Lx=48,dop=88,t=3,J=1,Jz=1.75,dim=6000)
    df_sisj_Jz20 = get_data_sisj(Lz=Lz,Ly=Ly,Lx=48,dop=88,t=3,J=1,Jz=2.0,dim=6000)
    df_sisj_Jz225 = get_data_sisj(Lz=Lz,Ly=Ly,Lx=48,dop=88,t=3,J=1,Jz=2.25,dim=6000)
    df_sisj_Jz30 = get_data_sisj(Lz=Lz,Ly=Ly,Lx=48,dop=88,t=3,J=1,Jz=3.0,dim=6000)
    #--------- J_\bot = 0.1
    slope = df_sisj_Jz01["slope_pow"].values[0]
    Kf = round(-slope,2)  
    label = r"$J_{\bot}$=%.1f $K_F$=%.2f"%(0.1,Kf)
    L01, = ax3.plot(df_sisj_Jz01["r"],df_sisj_Jz01["corre_abs"],label=label,ls="-",lw=1.5,color="red",
             marker='o',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="red", markerfacecolor='None')
    L01_, = ax3.plot(df_sisj_Jz01["r"],df_sisj_Jz01["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='o',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- J_\bot = 0.6
    #slope = df_sisj_Jz06["slope_pow"].values[0]
    #Kf = round(-slope,2)  
    #label = r"$J_{\bot}$=%.1f $K_F$=%.2f"%(0.6,Kf)
    #L06, = ax3.plot(df_sisj_Jz06["r"],df_sisj_Jz06["corre_abs"],label=label,ls="-",lw=1.5,color="blue",
    #         marker='s',alpha=1,markersize=6,markeredgewidth=1, markeredgecolor="blue", markerfacecolor='None')
    #L06_, = ax3.plot(df_sisj_Jz06["r"],df_sisj_Jz06["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
    #         marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- J_\bot = 1.0
    slope = df_sisj_Jz10["slope_pow"].values[0]
    Kf = round(-slope,2)  
    label = r"$J_{\bot}$=%.1f $K_F$=%.2f"%(1.0,Kf)
    L10, = ax3.plot(df_sisj_Jz10["r"],df_sisj_Jz10["corre_abs"],label=label,ls="-",lw=1.5,color="blue",
             marker='s',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="blue", markerfacecolor='None')
    L10_, = ax3.plot(df_sisj_Jz10["r"],df_sisj_Jz10["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- J_\bot = 1.5
    slope = df_sisj_Jz15["slope_pow"].values[0]
    Kf = round(-slope,2)  
    label = r"$J_{\bot}$=%.1f $K_F$=%.2f"%(1.5,Kf)
    L15, = ax3.plot(df_sisj_Jz15["r"],df_sisj_Jz15["corre_abs"],label=label,ls="-",lw=1.5,color="green",
             marker='v',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="green", markerfacecolor='None')
    L15_, = ax3.plot(df_sisj_Jz15["r"],df_sisj_Jz15["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='v',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 1.75
    #slope = df_sisj_Jz175["slope_pow"].values[0]
    #Kf = round(-slope,2)  
    #label = r"$J_{\bot}$=%.2f $K_F$=%.2f"%(1.75,Kf)
    #L175, = ax3.plot(df_sisj_Jz175["r"],df_sisj_Jz175["corre_abs"],label=label,ls="-",lw=1.5,color="cyan",
    #         marker='<',alpha=1,markersize=6,markeredgewidth=1, markeredgecolor="cyan", markerfacecolor='None')
    #L175_, = ax3.plot(df_sisj_Jz175["r"],df_sisj_Jz175["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
    #         marker='<',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 2.0
    slope = df_sisj_Jz20["slope_pow"].values[0]
    Kf = round(-slope,2)  
    #label = r"$J_{\bot}$=%.1f $K_F$=%.2f"%(2.0,Kf)
    label = r"$J_{\bot}$=%.1f"%(2.0,)
    L20, = ax3.plot(df_sisj_Jz20["r"],df_sisj_Jz20["corre_abs"],label=label,ls="-",lw=1.5,color="magenta",
             marker='>',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="magenta", markerfacecolor='None')
    #L20_, = ax3.plot(df_sisj_Jz20["r"],df_sisj_Jz20["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
    #         marker='>',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 2.25
    slope = df_sisj_Jz225["slope_pow"].values[0]
    Kf = round(-slope,2)  
    #label = r"$J_{\bot}$=%.2f $K_F$=%.2f"%(2.25,Kf)
    label = r"$J_{\bot}$=%.2f"%(2.25)
    L225, = ax3.plot(df_sisj_Jz225["r"],df_sisj_Jz225["corre_abs"],label=label,ls="-",lw=1.5,color="cyan",
             marker='^',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="cyan", markerfacecolor='None')
    #L225_, = ax3.plot(df_sisj_Jz225["r"],df_sisj_Jz225["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
    #         marker='^',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 3.0
    slope = df_sisj_Jz30["slope_pow"].values[0]
    Kf = round(-slope,2)  
    #label = r"$J_{\bot}$=%.2f $K_F$=%.2f"%(3.0,Kf)
    label = r"$J_{\bot}$=%.1f"%(3.0,)
    L30, = ax3.plot(df_sisj_Jz30["r"],df_sisj_Jz30["corre_abs"],label=label,ls="-",lw=1.5,color="brown",
             marker='<',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="brown", markerfacecolor='None')
    #L30_, = ax3.plot(df_sisj_Jz30["r"],df_sisj_Jz30["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
    #         marker='<',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------
    ####图例设置
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 15, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    #legend1=plt.legend(handles=[], loc = 4, bbox_to_anchor=(0.99, 0.678),
    #                   ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    legend2=plt.legend(handles=[L01,L10,L15,L20,L225,L30], loc = 4, bbox_to_anchor=(0.58, -0.01),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    #plt.gca().add_artist(legend1)#####把图例legend1重新加载回来
    label_x = r"r"
    label_y = "|F(r)|"
    #plt.yscale("log")
    #plt.xscale("log")
    ax3.set_yscale("log",base=10,subs=[1])   
    ax3.set_xscale("log",base=10,subs=[0.01,0.02,0.03])    
    ax3.set_xlabel(label_x, size= 14)
    ax3.set_ylabel(label_y, size= 14)
    ax3.tick_params(labelsize = 12) # 设置坐标刻度对应数字的大小
    ax3.set_xlim([1.1,40])
    ax3.set_xticks([1,5,10,20,30])
    #ax3.set_yticks([-1,-0.5,0,0.5,1])  
    #=========================================================
    ax3.text(27, 0.105,"(c)",fontsize = 20, color='black', rotation = 0)
    ax3.text(5, 4.0e-2, r'$\mathrm{\delta} = 0.4583$', fontsize = 16, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
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
    #-------------------------------------------------------------------------------------------
    # 选择子图 ax4 进行绘图
    plt.sca(ax4) ## 选择对 ax4 进行绘图
    ax4 = plt.gca()
    df_ninj_Jz01 = get_data_ninj_corre(Lz=Lz,Ly=Ly,Lx=Lx,dop=88,t=3,J=1,Jz=0.1,dim=6000)
    df_ninj_Jz10 = get_data_ninj_corre(Lz=Lz,Ly=Ly,Lx=Lx,dop=88,t=3,J=1,Jz=1.0,dim=6000)
    df_ninj_Jz15 = get_data_ninj_corre(Lz=Lz,Ly=Ly,Lx=Lx,dop=88,t=3,J=1,Jz=1.5,dim=6000)
    df_ninj_Jz225 = get_data_ninj_corre(Lz=Lz,Ly=Ly,Lx=Lx,dop=88,t=3,J=1,Jz=2.25,dim=6000)
    df_ninj_Jz30 = get_data_ninj_corre(Lz=Lz,Ly=Ly,Lx=Lx,dop=88,t=3,J=1,Jz=3.0,dim=6000)
    #
    #--------- J_\bot = 0.1
    slope = df_ninj_Jz01["slope_pow"].values[0]
    xi = round(-slope,2)  
    label = r"$J_{\bot}$=%.1f $K_D$=%.2f"%(0.1,xi)
    L01, = ax4.plot(df_ninj_Jz01["r"],df_ninj_Jz01["corre_abs"],label=label,ls="-",lw=1.5,color="red",
             marker='o',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="red", markerfacecolor='None')
    L01_, = ax4.plot(df_ninj_Jz01["r"],df_ninj_Jz01["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='o',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- J_\bot = 1.0
    #slope = df_ninj_Jz10["slope_pow"].values[0]
    #xi = round(-slope,2)  
    #label = r"$J_{\bot}$=%.1f $K_D$=%.2f"%(0.6,xi)
    #L10, = ax4.plot(df_ninj_Jz10["r"],df_ninj_Jz10["corre_abs"],label=label,ls="-",lw=1.5,color="blue",
    #         marker='s',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="blue", markerfacecolor='None')
    #L10_, = ax4.plot(df_ninj_Jz10["r"],df_ninj_Jz10["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
    #         marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- J_\bot = 0.8
    #slope = df_ninj_Jz08["slope_pow"].values[0]
    #xi = round(-slope,2)  
    #label = r"$J_{\bot}$=%.1f $K_D$=%.2f"%(0.8,xi)
    #L08, = ax4.plot(df_ninj_Jz08["r"],df_ninj_Jz08["corre_abs"],label=label,ls="-",lw=1.5,color="green",
    #         marker='^',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="green", markerfacecolor='None')
    #L08_, = ax4.plot(df_ninj_Jz08["r"],df_ninj_Jz08["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
    #         marker='^',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- J_\bot = 1.5
    slope = df_ninj_Jz15["slope_pow"].values[0]
    xi = round(-slope,2)  
    label = r"$J_{\bot}$=%.1f $K_D$=%.2f"%(1.5,xi)
    L15, = ax4.plot(df_ninj_Jz15["r"],df_ninj_Jz15["corre_abs"],label=label,ls="-",lw=1.5,color="blue",
             marker='s',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="blue", markerfacecolor='None')
    L15_, = ax4.plot(df_ninj_Jz15["r"],df_ninj_Jz15["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- J_\bot = 2.25
    #slope = df_ninj_Jz225["slope_pow"].values[0]
    #xi = round(-slope,2)  
    #label = r"$J_{\bot}$=%.2f $K_D$=%.2f"%(2.25,xi)
    #L225, = ax4.plot(df_ninj_Jz225["r"],df_ninj_Jz225["corre_abs"],label=label,ls="-",lw=1.5,color="magenta",
    #         marker='<',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="magenta", markerfacecolor='None')
    #L225_, = ax4.plot(df_ninj_Jz225["r"],df_ninj_Jz225["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
    #         marker='<',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- J_\bot = 3.0
    slope = df_ninj_Jz30["slope_pow"].values[0]
    xi = round(-slope,2)  
    label = r"$J_{\bot}$=%.1f $K_D$=%.2f"%(3.0,xi)
    L30, = ax4.plot(df_ninj_Jz30["r"],df_ninj_Jz30["corre_abs"],label=label,ls="-",lw=1.5,color="green",
             marker='^',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="green", markerfacecolor='None')
    L30_, = ax4.plot(df_ninj_Jz30["r"],df_ninj_Jz30["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='^',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #
    ####图例设置
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 15, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    #legend1=plt.legend(handles=[L04,L06,L08,], loc = 4, bbox_to_anchor=(0.99, 0.778),
    #                   ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    legend2=plt.legend(handles=[L01,L15,L30], loc = 4, bbox_to_anchor=(0.62, -0.03),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    #plt.gca().add_artist(legend1)#####把图例legend1重新加载回来
    label_x = r"r"
    label_y = "D(r)"
    #plt.yscale("log")
    #plt.xscale("log")
    #ax4.set_yscale("log",base=10,subs=[0.01,0.02,0.03])
    ax4.set_yscale("log",base=10,subs=[1])
    ax4.set_xscale("log",base=10,subs=[0.01,0.02,0.03])       
    ax4.set_xlabel(label_x, size= 14)
    ax4.set_ylabel(label_y, size= 14)
    ax4.tick_params(labelsize = 12) # 设置坐标刻度对应数字的大小
    ax4.set_xlim([1.1,40])
    ax4.set_xticks([1,5,10,20,30])
    #ax4.set_yticks([-1,-0.5,0,0.5,1])  
    #=========================================================
    ax4.text(25,0.06,"(d)",fontsize = 20, color='black', rotation = 0)
    ax4.text(5,0.03, r'$\mathrm{\delta} = 0.4583$', fontsize = 16, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
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
    #----------------------------------------------------------------------------------------------
    plt.sca(ax5) ## 选择对 ax5 进行绘图
    ax5 = plt.gca()
    #------------------------
    df_pyy_Jz01 = get_data_pyy(Lz=Lz,Ly=Ly,Lx=48,dop=88,t=3,J=1,Jz=0.1,dim=6000)
    df_pyy_Jz06 = get_data_pyy(Lz=Lz,Ly=Ly,Lx=48,dop=88,t=3,J=1,Jz=0.6,dim=6000)
    df_pyy_Jz10 = get_data_pyy(Lz=Lz,Ly=Ly,Lx=48,dop=88,t=3,J=1,Jz=1.0,dim=6000)
    df_pyy_Jz15 = get_data_pyy(Lz=Lz,Ly=Ly,Lx=48,dop=88,t=3,J=1,Jz=1.5,dim=6000)
    df_pyy_Jz175 = get_data_pyy(Lz=Lz,Ly=Ly,Lx=48,dop=88,t=3,J=1,Jz=1.75,dim=6000)
    df_pyy_Jz20 = get_data_pyy(Lz=Lz,Ly=Ly,Lx=48,dop=88,t=3,J=1,Jz=2.0,dim=6000)
    df_pyy_Jz225 = get_data_pyy(Lz=Lz,Ly=Ly,Lx=48,dop=88,t=3,J=1,Jz=2.25,dim=6000)
    df_pyy_Jz30 = get_data_pyy(Lz=Lz,Ly=Ly,Lx=48,dop=88,t=3,J=1,Jz=3.0,dim=6000)
    #--------- J_\bot = 0.1
    slope = df_pyy_Jz01["slope_pow"].values[0]
    Kp = round(-slope,2)  
    label = r"$J_{\bot}$=%.1f $K_{SC}^{yy}$=%.2f"%(0.1,Kp)
    L01, = ax5.plot(df_pyy_Jz01["r"],df_pyy_Jz01["corre_abs"],label=label,ls="-",lw=1.5,color="red",
             marker='o',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="red", markerfacecolor='None')
    L01_, = ax5.plot(df_pyy_Jz01["r"],df_pyy_Jz01["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='o',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- J_\bot = 0.6
    #slope = df_pyy_Jz06["slope_pow"].values[0]
    #Kp = round(-slope,2)  
    #label = r"$J_{\bot}$=%.1f $K_{SC}^{yy}$=%.2f"%(0.6,Kp)
    #L06, = ax5.plot(df_pyy_Jz06["r"],df_pyy_Jz06["corre_abs"],label=label,ls="-",lw=1.5,color="blue",
    #         marker='s',alpha=1,markersize=6,markeredgewidth=1, markeredgecolor="blue", markerfacecolor='None')
    #L06_, = ax5.plot(df_pyy_Jz06["r"],df_pyy_Jz06["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
    #         marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- J_\bot = 1.0
    slope = df_pyy_Jz10["slope_pow"].values[0]
    Kp = round(-slope,2)  
    label = r"$J_{\bot}$=%.1f $K_{SC}^{yy}$=%.2f"%(1.0,Kp)
    L10, = ax5.plot(df_pyy_Jz10["r"],df_pyy_Jz10["corre_abs"],label=label,ls="-",lw=1.5,color="blue",
             marker='s',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="blue", markerfacecolor='None')
    L10_, = ax5.plot(df_pyy_Jz10["r"],df_pyy_Jz10["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- J_\bot = 1.5
    slope = df_pyy_Jz15["slope_pow"].values[0]
    Kp = round(-slope,2)  
    label = r"$J_{\bot}$=%.1f $K_{SC}^{yy}$=%.2f"%(1.5,Kp)
    L15, = ax5.plot(df_pyy_Jz15["r"],df_pyy_Jz15["corre_abs"],label=label,ls="-",lw=1.5,color="green",
             marker='v',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="green", markerfacecolor='None')
    L15_, = ax5.plot(df_pyy_Jz15["r"],df_pyy_Jz15["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='v',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 1.75
    #slope = df_pyy_Jz175["slope_pow"].values[0]
    #Kp = round(-slope,2)  
    #label = r"$J_{\bot}$=%.2f $K_{SC}^{yy}$=%.2f"%(1.75,Kp)
    #L175, = ax5.plot(df_pyy_Jz175["r"],df_pyy_Jz175["corre_abs"],label=label,ls="-",lw=1.5,color="cyan",
    #         marker='<',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="cyan", markerfacecolor='None')
    #L175_, = ax5.plot(df_pyy_Jz175["r"],df_pyy_Jz175["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
    #         marker='<',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 2.0
    slope = df_pyy_Jz20["slope_pow"].values[0]
    Kp = round(-slope,2)  
    #label = r"$J_{\bot}$=%.1f $K_{SC}^{yy}$=%.2f"%(2.0,Kp)
    label = r"$J_{\bot}$=%.1f"%(2.0,)
    L20, = ax5.plot(df_pyy_Jz20["r"],df_pyy_Jz20["corre_abs"],label=label,ls="-",lw=1.5,color="cyan",
             marker='>',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="cyan", markerfacecolor='None')
    #L20_, = ax5.plot(df_pyy_Jz20["r"],df_pyy_Jz20["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
    #         marker='>',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 2.25
    #slope = df_pyy_Jz225["slope_pow"].values[0]
    #Kf = round(-slope,2)  
    #label = r"$J_{\bot}$=%.2f $K_{SC}^{yy}$=%.2f"%(2.25,Kf)
    #L225, = ax5.plot(df_pyy_Jz225["r"],df_pyy_Jz225["corre_abs"],label=label,ls="-",lw=1.5,color="olive",
    #         marker='h',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="olive", markerfacecolor='None')
    #L225_, = ax5.plot(df_pyy_Jz225["r"],df_pyy_Jz225["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
    #         marker='h',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 3.0
    slope = df_pyy_Jz30["slope_pow"].values[0]
    Kf = round(-slope,2)  
    #label = r"$J_{\bot}$=%.2f $K_{SC}^{yy}$=%.2f"%(3.0,Kf)
    label = r"$J_{\bot}$=%.1f"%(3.0,)
    L30, = ax5.plot(df_pyy_Jz30["r"],df_pyy_Jz30["corre_abs"],label=label,ls="-",lw=1.5,color="magenta",
             marker='<',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="magenta", markerfacecolor='None')
    #L30_, = ax5.plot(df_pyy_Jz30["r"],df_pyy_Jz30["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
    #         marker='<',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    ####图例设置
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 15, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    #legend1=plt.legend(handles=[L01,L06,], loc = 4, bbox_to_anchor=(0.99, 0.75),
    #                   ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    legend2=plt.legend(handles=[L01,L10,L15,L20,L30], loc = 4, bbox_to_anchor=(0.62, -0.01),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    #plt.gca().add_artist(legend1)#####把图例legend1重新加载回来
    label_x = r"r"
    label_y = "$P^{yy}(r)$"
    #plt.yscale("log")
    #plt.xscale("log")
    ax5.set_yscale("log",base=10,subs=[1])   
    ax5.set_xscale("log",base=10,subs=[0.01,0.02,0.03])    
    ax5.set_xlabel(label_x, size= 14)
    ax5.set_ylabel(label_y, size= 14)
    ax5.tick_params(labelsize = 12) # 设置坐标刻度对应数字的大小
    ax5.set_xlim([1.7,40])
    ax5.set_xticks([5,10,20,30])
    #ax5.set_yticks([-1,-0.5,0,0.5,1])  
    #=========================================================
    ax5.text(28, 0.012,"(e)",fontsize = 20, color='black', rotation = 0)
    ax5.text(7.5, 3.0e-3, r'$\mathrm{\delta} = 0.4583$', fontsize = 16, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
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
    #-----------------------------------------------------------------------------------------------
    # 选择子图 ax6 进行绘图
    plt.sca(ax6) ## 选择对 ax6 进行绘图
    ax6 = plt.gca()
    #------------------------
    df_pzz_Jz01 = get_data_pzz(Lz=Lz,Ly=Ly,Lx=48,dop=88,t=3,J=1,Jz=0.1,dim=8000)
    df_pzz_Jz06 = get_data_pzz(Lz=Lz,Ly=Ly,Lx=48,dop=88,t=3,J=1,Jz=0.6,dim=8000)
    df_pzz_Jz10 = get_data_pzz(Lz=Lz,Ly=Ly,Lx=48,dop=88,t=3,J=1,Jz=1.0,dim=8000)
    df_pzz_Jz15 = get_data_pzz(Lz=Lz,Ly=Ly,Lx=48,dop=88,t=3,J=1,Jz=1.5,dim=6000)
    df_pzz_Jz175 = get_data_pzz(Lz=Lz,Ly=Ly,Lx=48,dop=88,t=3,J=1,Jz=1.75,dim=6000)
    df_pzz_Jz20 = get_data_pzz(Lz=Lz,Ly=Ly,Lx=48,dop=88,t=3,J=1,Jz=2.0,dim=6000)
    df_pzz_Jz225 = get_data_pzz(Lz=Lz,Ly=Ly,Lx=48,dop=88,t=3,J=1,Jz=2.25,dim=6000)
    df_pzz_Jz30 = get_data_pzz(Lz=Lz,Ly=Ly,Lx=48,dop=88,t=3,J=1,Jz=3.0,dim=6000)
    #--------- J_\bot = 0.1
    slope = df_pzz_Jz01["slope_pow"].values[0]
    Kp = round(-slope,2)  
    label = r"$J_{\bot}$=%.1f $K_{SC}^{zz}$=%.2f"%(0.1,Kp)
    L01, = ax6.plot(df_pzz_Jz01["r"],df_pzz_Jz01["corre_abs"],label=label,ls="-",lw=1.5,color="red",
             marker='o',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="red", markerfacecolor='None')
    L01_, = ax6.plot(df_pzz_Jz01["r"],df_pzz_Jz01["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='o',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- J_\bot = 0.6
    slope = df_pzz_Jz06["slope_pow"].values[0]
    Kp = round(-slope,2)  
    label = r"$J_{\bot}$=%.1f $K_{SC}^{zz}$=%.2f"%(0.6,Kp)
    L06, = ax6.plot(df_pzz_Jz06["r"],df_pzz_Jz06["corre_abs"],label=label,ls="-",lw=1.5,color="blue",
             marker='s',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="blue", markerfacecolor='None')
    L06_, = ax6.plot(df_pzz_Jz06["r"],df_pzz_Jz06["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- J_\bot = 1.0
    slope = df_pzz_Jz10["slope_pow"].values[0]
    Kg = round(-slope,2)  
    label = r"$J_{\bot}$=%.1f $K_{SC}^{zz}$=%.2f"%(1.0,Kg)
    L10, = ax6.plot(df_pzz_Jz10["r"],df_pzz_Jz10["corre_abs"],label=label,ls="-",lw=1.5,color="green",
             marker='^',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="green", markerfacecolor='None')
    L10_, = ax6.plot(df_pzz_Jz10["r"],df_pzz_Jz10["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='^',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- J_\bot = 1.5
    slope = df_pzz_Jz15["slope_pow"].values[0]
    Kg = round(-slope,2)  
    label = r"$J_{\bot}$=%.1f $K_{SC}^{zz}$=%.2f"%(1.5,Kg)
    L15, = ax6.plot(df_pzz_Jz15["r"],df_pzz_Jz15["corre_abs"],label=label,ls="-",lw=1.5,color="magenta",
             marker='v',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="magenta", markerfacecolor='None')
    L15_, = ax6.plot(df_pzz_Jz15["r"],df_pzz_Jz15["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='v',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 1.75
    #slope = df_pzz_Jz175["slope_pow"].values[0]
    #Kg = round(-slope,2)  
    #label = r"$J_{\bot}$=%.2f $K_{SC}^{zz}$=%.2f"%(1.75,Kg)
    #L175, = ax6.plot(df_pzz_Jz175["r"],df_pzz_Jz175["corre_abs"],label=label,ls="-",lw=1.5,color="cyan",
    #         marker='<',alpha=1,markersize=6,markeredgewidth=1, markeredgecolor="cyan", markerfacecolor='None')
    #L175_, = ax6.plot(df_pzz_Jz175["r"],df_pzz_Jz175["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
    #         marker='<',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 2.0
    slope = df_pzz_Jz20["slope_pow"].values[0]
    Kg = round(-slope,2)  
    label = r"$J_{\bot}$=%.1f $K_{SC}^{zz}$=%.2f"%(2.0,Kg)
    L20, = ax6.plot(df_pzz_Jz20["r"],df_pzz_Jz20["corre_abs"],label=label,ls="-",lw=1.5,color="brown",
             marker='<',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="brown", markerfacecolor='None')
    L20_, = ax6.plot(df_pzz_Jz20["r"],df_pzz_Jz20["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='<',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 2.25
    #slope = df_pzz_Jz225["slope_pow"].values[0]
    #Kg = round(-slope,2)  
    #label = r"$J_{\bot}$=%.2f $K_{SC}^{zz}$=%.2f"%(2.25,Kg)
    #L225, = ax6.plot(df_pzz_Jz225["r"],df_pzz_Jz225["corre_abs"],label=label,ls="-",lw=1.5,color="olive",
    #         marker='h',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="olive", markerfacecolor='None')
    #L225_, = ax6.plot(df_pzz_Jz225["r"],df_pzz_Jz225["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
    #         marker='h',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 3.0
    slope = df_pzz_Jz30["slope_pow"].values[0]
    Kg = round(-slope,2)  
    label = r"$J_{\bot}$=%.1f $K_{SC}^{zz}$=%.2f"%(3.0,Kg)
    L30, = ax6.plot(df_pzz_Jz30["r"],df_pzz_Jz30["corre_abs"],label=label,ls="-",lw=1.5,color="olive",
             marker='>',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="olive", markerfacecolor='None')
    L30_, = ax6.plot(df_pzz_Jz30["r"],df_pzz_Jz30["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='>',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------
    ####图例设置
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 13.6, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    #legend1=plt.legend(handles=[L01,], loc = 4, bbox_to_anchor=(0.99, 0.89),
    #                   ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    legend2=plt.legend(handles=[L01,L06,L10,L15,L20,L30], loc = 4, bbox_to_anchor=(0.56, -0.045),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    #plt.gca().add_artist(legend1)#####把图例legend1重新加载回来
    label_x = r"r"
    label_y = "$P^{zz}(r)$"
    #plt.yscale("log")
    #plt.xscale("log")
    ax6.set_yscale("log",base=10,subs=[1])   
    ax6.set_xscale("log",base=10,subs=[0.01,0.02,0.03])    
    ax6.set_xlabel(label_x, size= 14)
    ax6.set_ylabel(label_y, size= 14)
    ax6.tick_params(labelsize = 12) # 设置坐标刻度对应数字的大小
    ax6.set_xlim([0,40])
    ax6.set_ylim([1.0e-7,2.0e-2])
    ax6.set_xticks([5,10,20,30])
    #ax6.set_yticks([-1,-0.5,0,0.5,1])  
    #=========================================================
    ax6.text(30,1.0e-2,"(f)",fontsize = 20, color='black', rotation = 0)
    ax6.text(8,1.0e-2, r'$\mathrm{\delta} = 0.4583$', fontsize = 16, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
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
    #
    #===========================================================================================
    fig.tight_layout() # 自动调整 subplot 间的间隙参数
    plt.savefig("E:\\WORK\\Work\\Project\\La3Ni2O7\\datashow_datpy_paper\\figs\\fig_04583_1_v2.eps",
                dpi=300, format='eps',bbox_inches='tight') # 白边紧凑型
    plt.show()