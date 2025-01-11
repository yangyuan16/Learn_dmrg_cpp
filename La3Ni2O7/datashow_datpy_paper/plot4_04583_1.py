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
def get_data_Lz1(Lz, Ly, Lx, dop, t, J, Jz, dim):
    workpath = "E:\\WORK\\Work\\Project\\leg2-tj"
    filepath1 = "\\data_dmrgcpp\\Lz%d_Ly%d_Lx%d\\dop%g" % (Lz, Ly, Lx,dop)
    filepath2 = "\\t%d_J%d_Jz%d_dim%d_ent\\entanglement" % (t, J, Jz, dim)
    filepath3 = "\\measurement_entropy_fit.parquet"
    filename = workpath + filepath1 + filepath2 + filepath3
    print(filename)
    df = pd.read_parquet(filename)
    print(df.head())
    return df
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
def get_data_pxx(Lz, Ly, Lx, dop, t, J, Jz, dim,):
    workpath = "E:\\WORK\\Work\\Project\\La3Ni2O7"
    filepath1 = "\\data_dmrgcpp\\Lz%d_Ly%d_Lx%d\\dop%g" % (Lz, Ly, Lx,dop)
    filepath2 = "\\t%d_J%d_Jz%.2f_dim%d" % (t, J, Jz, dim)
    filepath3 = "\\measurement_pairing_xx_fit.parquet"
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
    t = 3
    J = 1
    dim = 6000
    #
    #---------------plot logr-logr fig-----------------------------
    fig = plt.figure(figsize=(11,10))
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
    ax1 = plt.axes([0.1,0.58,0.37,0.375])
    ax2 = plt.axes([0.58,0.58,0.37,0.375])
    ax3 = plt.axes([0.1,0.1,0.37,0.375])
    ax4 = plt.axes([0.58,0.1,0.37,0.375])
    # plt.subplots_adjust(left=0.1, bottom=0.1, right=0.98, top=0.96, wspace=0.32, hspace=0.26)
    #--------------------------------------------------------------------------
    # 选择子图 ax1 进行绘图
    plt.sca(ax1) ## 选择对 ax1 进行绘图
    ax1 = plt.gca()
    #------------------------
    df_sisj_Jz01 = get_data_sisj(Lz=Lz,Ly=Ly,Lx=48,dop=88,t=t,J=J,Jz=0.1,dim=dim)
    df_sisj_Jz06 = get_data_sisj(Lz=Lz,Ly=Ly,Lx=48,dop=88,t=t,J=J,Jz=0.6,dim=dim)
    df_sisj_Jz10 = get_data_sisj(Lz=Lz,Ly=Ly,Lx=48,dop=88,t=t,J=J,Jz=1.0,dim=dim)
    df_sisj_Jz15 = get_data_sisj(Lz=Lz,Ly=Ly,Lx=48,dop=88,t=t,J=J,Jz=1.5,dim=dim)
    df_sisj_Jz175 = get_data_sisj(Lz=Lz,Ly=Ly,Lx=48,dop=88,t=t,J=J,Jz=1.75,dim=dim)
    df_sisj_Jz20 = get_data_sisj(Lz=Lz,Ly=Ly,Lx=48,dop=88,t=t,J=J,Jz=2.0,dim=dim)
    df_sisj_Jz225 = get_data_sisj(Lz=Lz,Ly=Ly,Lx=48,dop=88,t=t,J=J,Jz=2.25,dim=dim)
    #--------- J_\bot = 0.1
    slope = df_sisj_Jz01["slope_exp"].values[0]
    xi = round(-1/slope,2)  
    label = r"$J_{\bot}$=%.1f $\xi_F$=%.2f"%(0.1,xi)
    L01, = ax1.plot(df_sisj_Jz01["r"],df_sisj_Jz01["corre_abs"],label=label,ls="-",lw=1.5,color="red",
             marker='o',alpha=1,markersize=6,markeredgewidth=1, markeredgecolor="red", markerfacecolor='None')
    L01_, = ax1.plot(df_sisj_Jz01["r"],df_sisj_Jz01["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
             marker='o',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- J_\bot = 0.6
    slope = df_sisj_Jz06["slope_exp"].values[0]
    xi = round(-1/slope,2)  
    label = r"$J_{\bot}$=%.1f $\xi_F$=%.2f"%(0.6,xi)
    L06, = ax1.plot(df_sisj_Jz06["r"],df_sisj_Jz06["corre_abs"],label=label,ls="-",lw=1.5,color="blue",
             marker='s',alpha=1,markersize=6,markeredgewidth=1, markeredgecolor="blue", markerfacecolor='None')
    L06_, = ax1.plot(df_sisj_Jz06["r"],df_sisj_Jz06["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
             marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- J_\bot = 1.0
    slope = df_sisj_Jz10["slope_exp"].values[0]
    xi = round(-1/slope,2)  
    label = r"$J_{\bot}$=%.1f $\xi_F$=%.2f"%(1.0,xi)
    L10, = ax1.plot(df_sisj_Jz10["r"],df_sisj_Jz10["corre_abs"],label=label,ls="-",lw=1.5,color="green",
             marker='^',alpha=1,markersize=6,markeredgewidth=1, markeredgecolor="green", markerfacecolor='None')
    L10_, = ax1.plot(df_sisj_Jz10["r"],df_sisj_Jz10["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
             marker='^',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- J_\bot = 1.5
    slope = df_sisj_Jz15["slope_exp"].values[0]
    xi = round(-1/slope,2)  
    label = r"$J_{\bot}$=%.1f $\xi_F$=%.2f"%(1.5,xi)
    L15, = ax1.plot(df_sisj_Jz15["r"],df_sisj_Jz15["corre_abs"],label=label,ls="-",lw=1.5,color="magenta",
             marker='v',alpha=1,markersize=6,markeredgewidth=1, markeredgecolor="magenta", markerfacecolor='None')
    L15_, = ax1.plot(df_sisj_Jz15["r"],df_sisj_Jz15["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
             marker='v',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 1.75
    slope = df_sisj_Jz175["slope_exp"].values[0]
    xi = round(-1/slope,2)  
    label = r"$J_{\bot}$=%.2f $\xi_F$=%.2f"%(1.75,xi)
    L175, = ax1.plot(df_sisj_Jz175["r"],df_sisj_Jz175["corre_abs"],label=label,ls="-",lw=1.5,color="cyan",
             marker='<',alpha=1,markersize=6,markeredgewidth=1, markeredgecolor="cyan", markerfacecolor='None')
    L175_, = ax1.plot(df_sisj_Jz175["r"],df_sisj_Jz175["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
             marker='<',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 2.0
    slope = df_sisj_Jz20["slope_exp"].values[0]
    xi = round(-1/slope,2)  
    label = r"$J_{\bot}$=%.1f $\xi_F$=%.2f"%(2.0,xi)
    L20, = ax1.plot(df_sisj_Jz20["r"],df_sisj_Jz20["corre_abs"],label=label,ls="-",lw=1.5,color="brown",
             marker='>',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="brown", markerfacecolor='None')
    L20_, = ax1.plot(df_sisj_Jz20["r"],df_sisj_Jz20["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
             marker='>',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 2.25
    slope = df_sisj_Jz225["slope_exp"].values[0]
    xi = round(-1/slope,2)  
    label = r"$J_{\bot}$=%.2f $\xi_F$=%.2f"%(2.25,xi)
    L225, = ax1.plot(df_sisj_Jz225["r"],df_sisj_Jz225["corre_abs"],label=label,ls="-",lw=1.5,color="olive",
             marker='h',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="olive", markerfacecolor='None')
    L225_, = ax1.plot(df_sisj_Jz225["r"],df_sisj_Jz225["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
             marker='h',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------
    ####图例设置
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 12, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[L01,L06,L10,L15], loc = 4, bbox_to_anchor=(0.99, 0.678),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    legend2=plt.legend(handles=[L175,L20,L225], loc = 4, bbox_to_anchor=(0.54, -0.01),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    plt.gca().add_artist(legend1)#####把图例legend1重新加载回来
    label_x = r"r"
    label_y = "|F(r)|"
    #plt.yscale("log")
    #plt.xscale("log")
    ax1.set_yscale("log",base=10,subs=[0.01,0.02,0.03])      
    ax1.set_xlabel(label_x, size= 14)
    ax1.set_ylabel(label_y, size= 14)
    ax1.tick_params(labelsize = 12) # 设置坐标刻度对应数字的大小
    ax1.set_xlim([0,40])
    #ax1.set_xticks([0,10,20,30])
    #ax1.set_yticks([-1,-0.5,0,0.5,1])  
    #=========================================================
    ax1.text(4,0.02,"(a)",fontsize = 20, color='black', rotation = 0)
    ax1.text(6,5.0e-3, r'$\mathrm{\delta} = 0.4583$', fontsize = 16, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
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
    plt.sca(ax2) ## 选择对 ax2 进行绘图
    ax2 = plt.gca()
    #--------- J_\bot = 0.1
    slope = df_sisj_Jz01["slope_pow"].values[0]
    Kf = round(-slope,2)  
    label = r"$J_{\bot}$=%.1f $K_F$=%.2f"%(0.1,Kf)
    L01, = ax2.plot(df_sisj_Jz01["r"],df_sisj_Jz01["corre_abs"],label=label,ls="-",lw=1.5,color="red",
             marker='o',alpha=1,markersize=6,markeredgewidth=1, markeredgecolor="red", markerfacecolor='None')
    L01_, = ax2.plot(df_sisj_Jz01["r"],df_sisj_Jz01["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='o',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- J_\bot = 0.6
    slope = df_sisj_Jz06["slope_pow"].values[0]
    Kf = round(-slope,2)  
    label = r"$J_{\bot}$=%.1f $K_F$=%.2f"%(0.6,Kf)
    L06, = ax2.plot(df_sisj_Jz06["r"],df_sisj_Jz06["corre_abs"],label=label,ls="-",lw=1.5,color="blue",
             marker='s',alpha=1,markersize=6,markeredgewidth=1, markeredgecolor="blue", markerfacecolor='None')
    L06_, = ax2.plot(df_sisj_Jz06["r"],df_sisj_Jz06["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- J_\bot = 1.0
    slope = df_sisj_Jz10["slope_pow"].values[0]
    Kf = round(-slope,2)  
    label = r"$J_{\bot}$=%.1f $K_F$=%.2f"%(1.0,Kf)
    L10, = ax2.plot(df_sisj_Jz10["r"],df_sisj_Jz10["corre_abs"],label=label,ls="-",lw=1.5,color="green",
             marker='^',alpha=1,markersize=6,markeredgewidth=1, markeredgecolor="green", markerfacecolor='None')
    L10_, = ax2.plot(df_sisj_Jz10["r"],df_sisj_Jz10["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='^',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- J_\bot = 1.5
    slope = df_sisj_Jz15["slope_pow"].values[0]
    Kf = round(-slope,2)  
    label = r"$J_{\bot}$=%.1f $K_F$=%.2f"%(1.5,Kf)
    L15, = ax2.plot(df_sisj_Jz15["r"],df_sisj_Jz15["corre_abs"],label=label,ls="-",lw=1.5,color="magenta",
             marker='v',alpha=1,markersize=6,markeredgewidth=1, markeredgecolor="magenta", markerfacecolor='None')
    L15_, = ax2.plot(df_sisj_Jz15["r"],df_sisj_Jz15["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='v',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 1.75
    slope = df_sisj_Jz175["slope_pow"].values[0]
    Kf = round(-slope,2)  
    label = r"$J_{\bot}$=%.2f $K_F$=%.2f"%(1.75,Kf)
    L175, = ax2.plot(df_sisj_Jz175["r"],df_sisj_Jz175["corre_abs"],label=label,ls="-",lw=1.5,color="cyan",
             marker='<',alpha=1,markersize=6,markeredgewidth=1, markeredgecolor="cyan", markerfacecolor='None')
    L175_, = ax2.plot(df_sisj_Jz175["r"],df_sisj_Jz175["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='<',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 2.0
    slope = df_sisj_Jz20["slope_pow"].values[0]
    Kf = round(-slope,2)  
    label = r"$J_{\bot}$=%.1f $K_F$=%.2f"%(2.0,Kf)
    L20, = ax2.plot(df_sisj_Jz20["r"],df_sisj_Jz20["corre_abs"],label=label,ls="-",lw=1.5,color="brown",
             marker='>',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="brown", markerfacecolor='None')
    L20_, = ax2.plot(df_sisj_Jz20["r"],df_sisj_Jz20["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='>',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 2.25
    slope = df_sisj_Jz225["slope_pow"].values[0]
    Kf = round(-slope,2)  
    label = r"$J_{\bot}$=%.2f $K_F$=%.2f"%(2.25,Kf)
    L225, = ax2.plot(df_sisj_Jz225["r"],df_sisj_Jz225["corre_abs"],label=label,ls="-",lw=1.5,color="olive",
             marker='h',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="olive", markerfacecolor='None')
    L225_, = ax2.plot(df_sisj_Jz225["r"],df_sisj_Jz225["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='h',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------
    ####图例设置
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 12, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[L01,L06,L10,], loc = 4, bbox_to_anchor=(0.99, 0.678),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    legend2=plt.legend(handles=[L15,L175,L20,L225], loc = 4, bbox_to_anchor=(0.54, -0.01),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    plt.gca().add_artist(legend1)#####把图例legend1重新加载回来
    label_x = r"r"
    label_y = "|F(r)|"
    #plt.yscale("log")
    #plt.xscale("log")
    ax2.set_yscale("log",base=10,subs=[0.01,0.02,0.03])   
    ax2.set_xscale("log",base=10,subs=[0.01,0.02,0.03])    
    ax2.set_xlabel(label_x, size= 14)
    ax2.set_ylabel(label_y, size= 14)
    ax2.tick_params(labelsize = 12) # 设置坐标刻度对应数字的大小
    ax2.set_xlim([0,40])
    ax2.set_xticks([1,5,10,20,30])
    #ax2.set_yticks([-1,-0.5,0,0.5,1])  
    #=========================================================
    ax2.text(1.1, 0.125,"(b)",fontsize = 20, color='black', rotation = 0)
    ax2.text(1, 1.0e-3, r'$\mathrm{\delta} = 0.4583$', fontsize = 16, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
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
    #--------------------------------------------------------------------------------------------
    # 选择子图 ax3 进行绘图
    plt.sca(ax3) ## 选择对 ax2 进行绘图
    ax3 = plt.gca()
    #------------------------
    df_cicj_Jz01 = get_data_cicj(Lz=Lz,Ly=Ly,Lx=48,dop=88,t=t,J=J,Jz=0.1,dim=dim)
    df_cicj_Jz06 = get_data_cicj(Lz=Lz,Ly=Ly,Lx=48,dop=88,t=t,J=J,Jz=0.6,dim=dim)
    df_cicj_Jz10 = get_data_cicj(Lz=Lz,Ly=Ly,Lx=48,dop=88,t=t,J=J,Jz=1.0,dim=dim)
    df_cicj_Jz15 = get_data_cicj(Lz=Lz,Ly=Ly,Lx=48,dop=88,t=t,J=J,Jz=1.5,dim=dim)
    df_cicj_Jz175 = get_data_cicj(Lz=Lz,Ly=Ly,Lx=48,dop=88,t=t,J=J,Jz=1.75,dim=dim)
    df_cicj_Jz20 = get_data_cicj(Lz=Lz,Ly=Ly,Lx=48,dop=88,t=t,J=J,Jz=2.0,dim=dim)
    df_cicj_Jz225 = get_data_cicj(Lz=Lz,Ly=Ly,Lx=48,dop=88,t=t,J=J,Jz=2.25,dim=dim)
    #--------- J_\bot = 0.1
    slope = df_cicj_Jz01["slope_exp"].values[0]
    xi = round(-1/slope,2)  
    label = r"$J_{\bot}$=%.1f $\xi_G$=%.2f"%(0.1,xi)
    L01, = ax3.plot(df_cicj_Jz01["r"],df_cicj_Jz01["corre_abs"],label=label,ls="-",lw=1.5,color="red",
             marker='o',alpha=1,markersize=6,markeredgewidth=1, markeredgecolor="red", markerfacecolor='None')
    L01_, = ax3.plot(df_cicj_Jz01["r"],df_cicj_Jz01["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
             marker='o',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- J_\bot = 0.6
    slope = df_cicj_Jz06["slope_exp"].values[0]
    xi = round(-1/slope,2)  
    label = r"$J_{\bot}$=%.1f $\xi_G$=%.2f"%(0.6,xi)
    L06, = ax3.plot(df_cicj_Jz06["r"],df_cicj_Jz06["corre_abs"],label=label,ls="-",lw=1.5,color="blue",
             marker='s',alpha=1,markersize=6,markeredgewidth=1, markeredgecolor="blue", markerfacecolor='None')
    L06_, = ax3.plot(df_cicj_Jz06["r"],df_cicj_Jz06["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
             marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- J_\bot = 1.0
    slope = df_cicj_Jz10["slope_exp"].values[0]
    xi = round(-1/slope,2)  
    label = r"$J_{\bot}$=%.1f $\xi_G$=%.2f"%(1.0,xi)
    L10, = ax3.plot(df_cicj_Jz10["r"],df_cicj_Jz10["corre_abs"],label=label,ls="-",lw=1.5,color="green",
             marker='^',alpha=1,markersize=6,markeredgewidth=1, markeredgecolor="green", markerfacecolor='None')
    L10_, = ax3.plot(df_cicj_Jz10["r"],df_cicj_Jz10["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
             marker='^',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- J_\bot = 1.5
    slope = df_cicj_Jz15["slope_exp"].values[0]
    xi = round(-1/slope,2)  
    label = r"$J_{\bot}$=%.1f $\xi_G$=%.2f"%(1.5,xi)
    L15, = ax3.plot(df_cicj_Jz15["r"],df_cicj_Jz15["corre_abs"],label=label,ls="-",lw=1.5,color="magenta",
             marker='v',alpha=1,markersize=6,markeredgewidth=1, markeredgecolor="magenta", markerfacecolor='None')
    L15_, = ax3.plot(df_cicj_Jz15["r"],df_cicj_Jz15["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
             marker='v',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 1.75
    slope = df_cicj_Jz175["slope_exp"].values[0]
    xi = round(-1/slope,2)  
    label = r"$J_{\bot}$=%.2f $\xi_G$=%.2f"%(1.75,xi)
    L175, = ax3.plot(df_cicj_Jz175["r"],df_cicj_Jz175["corre_abs"],label=label,ls="-",lw=1.5,color="cyan",
             marker='<',alpha=1,markersize=6,markeredgewidth=1, markeredgecolor="cyan", markerfacecolor='None')
    L175_, = ax3.plot(df_cicj_Jz175["r"],df_cicj_Jz175["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
             marker='<',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 2.0
    slope = df_cicj_Jz20["slope_exp"].values[0]
    xi = round(-1/slope,2)  
    label = r"$J_{\bot}$=%.1f $\xi_G$=%.2f"%(2.0,xi)
    L20, = ax3.plot(df_cicj_Jz20["r"],df_cicj_Jz20["corre_abs"],label=label,ls="-",lw=1.5,color="brown",
             marker='>',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="brown", markerfacecolor='None')
    L20_, = ax3.plot(df_cicj_Jz20["r"],df_cicj_Jz20["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
             marker='>',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 2.25
    slope = df_cicj_Jz225["slope_exp"].values[0]
    xi = round(-1/slope,2)  
    label = r"$J_{\bot}$=%.2f $\xi_G$=%.2f"%(2.25,xi)
    L225, = ax3.plot(df_cicj_Jz225["r"],df_cicj_Jz225["corre_abs"],label=label,ls="-",lw=1.5,color="olive",
             marker='h',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="olive", markerfacecolor='None')
    L225_, = ax3.plot(df_cicj_Jz225["r"],df_cicj_Jz225["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
             marker='h',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------
    ####图例设置
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 12, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[L01,L06,L10,L15], loc = 4, bbox_to_anchor=(0.99, 0.678),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    legend2=plt.legend(handles=[L175,L20,L225], loc = 4, bbox_to_anchor=(0.54, -0.01),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    plt.gca().add_artist(legend1)#####把图例legend1重新加载回来
    label_x = r"r"
    label_y = "|G(r)|"
    #plt.yscale("log")
    #plt.xscale("log")
    ax3.set_yscale("log",base=10,subs=[0.01,0.02,0.03])      
    ax3.set_xlabel(label_x, size= 14)
    ax3.set_ylabel(label_y, size= 14)
    ax3.tick_params(labelsize = 12) # 设置坐标刻度对应数字的大小
    ax3.set_xlim([0,40])
    #ax3.set_xticks([0,10,20,30])
    #ax3.set_yticks([-1,-0.5,0,0.5,1])  
    #=========================================================
    ax3.text(4,0.2,"(c)",fontsize = 20, color='black', rotation = 0)
    ax3.text(6,5.0e-2, r'$\mathrm{\delta} = 0.4583$', fontsize = 16, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
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
    #-----------------------------------------------------------------------------------------------
    # 选择子图 ax4 进行绘图
    plt.sca(ax4) ## 选择对 ax4 进行绘图
    ax4 = plt.gca()
    #--------- J_\bot = 0.1
    slope = df_cicj_Jz01["slope_pow"].values[0]
    Kg = round(-slope,2)  
    label = r"$J_{\bot}$=%.1f $K_G$=%.2f"%(0.1,Kg)
    L01, = ax4.plot(df_cicj_Jz01["r"],df_cicj_Jz01["corre_abs"],label=label,ls="-",lw=1.5,color="red",
             marker='o',alpha=1,markersize=6,markeredgewidth=1, markeredgecolor="red", markerfacecolor='None')
    L01_, = ax4.plot(df_cicj_Jz01["r"],df_cicj_Jz01["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='o',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- J_\bot = 0.6
    slope = df_cicj_Jz06["slope_pow"].values[0]
    Kg = round(-slope,2)  
    label = r"$J_{\bot}$=%.1f $K_G$=%.2f"%(0.6,Kg)
    L06, = ax4.plot(df_cicj_Jz06["r"],df_cicj_Jz06["corre_abs"],label=label,ls="-",lw=1.5,color="blue",
             marker='s',alpha=1,markersize=6,markeredgewidth=1, markeredgecolor="blue", markerfacecolor='None')
    L06_, = ax4.plot(df_cicj_Jz06["r"],df_cicj_Jz06["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- J_\bot = 1.0
    slope = df_cicj_Jz10["slope_pow"].values[0]
    Kg = round(-slope,2)  
    label = r"$J_{\bot}$=%.1f $K_G$=%.2f"%(1.0,Kg)
    L10, = ax4.plot(df_cicj_Jz10["r"],df_cicj_Jz10["corre_abs"],label=label,ls="-",lw=1.5,color="green",
             marker='^',alpha=1,markersize=6,markeredgewidth=1, markeredgecolor="green", markerfacecolor='None')
    L10_, = ax4.plot(df_cicj_Jz10["r"],df_cicj_Jz10["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='^',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- J_\bot = 1.5
    slope = df_cicj_Jz15["slope_pow"].values[0]
    Kg = round(-slope,2)  
    label = r"$J_{\bot}$=%.1f $K_G$=%.2f"%(1.5,Kg)
    L15, = ax4.plot(df_cicj_Jz15["r"],df_cicj_Jz15["corre_abs"],label=label,ls="-",lw=1.5,color="magenta",
             marker='v',alpha=1,markersize=6,markeredgewidth=1, markeredgecolor="magenta", markerfacecolor='None')
    L15_, = ax4.plot(df_cicj_Jz15["r"],df_cicj_Jz15["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='v',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 1.75
    slope = df_cicj_Jz175["slope_pow"].values[0]
    Kg = round(-slope,2)  
    label = r"$J_{\bot}$=%.2f $K_G$=%.2f"%(1.75,Kg)
    L175, = ax4.plot(df_cicj_Jz175["r"],df_cicj_Jz175["corre_abs"],label=label,ls="-",lw=1.5,color="cyan",
             marker='<',alpha=1,markersize=6,markeredgewidth=1, markeredgecolor="cyan", markerfacecolor='None')
    L175_, = ax4.plot(df_cicj_Jz175["r"],df_cicj_Jz175["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='<',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 2.0
    slope = df_cicj_Jz20["slope_pow"].values[0]
    Kg = round(-slope,2)  
    label = r"$J_{\bot}$=%.1f $K_G$=%.2f"%(2.0,Kg)
    L20, = ax4.plot(df_cicj_Jz20["r"],df_cicj_Jz20["corre_abs"],label=label,ls="-",lw=1.5,color="brown",
             marker='>',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="brown", markerfacecolor='None')
    L20_, = ax4.plot(df_cicj_Jz20["r"],df_cicj_Jz20["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='>',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 2.25
    slope = df_cicj_Jz225["slope_pow"].values[0]
    Kg = round(-slope,2)  
    label = r"$J_{\bot}$=%.2f $K_G$=%.2f"%(2.25,Kg)
    L225, = ax4.plot(df_cicj_Jz225["r"],df_cicj_Jz225["corre_abs"],label=label,ls="-",lw=1.5,color="olive",
             marker='h',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="olive", markerfacecolor='None')
    L225_, = ax4.plot(df_cicj_Jz225["r"],df_cicj_Jz225["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='h',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------
    ####图例设置
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 12, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[L01,L06,L10,], loc = 4, bbox_to_anchor=(0.99, 0.72),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    legend2=plt.legend(handles=[L15,L175,L20,L225], loc = 4, bbox_to_anchor=(0.54, -0.01),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    plt.gca().add_artist(legend1)#####把图例legend1重新加载回来
    label_x = r"r"
    label_y = "|G(r)|"
    #plt.yscale("log")
    #plt.xscale("log")
    ax4.set_yscale("log",base=10,subs=[0.01,0.02,0.03])   
    ax4.set_xscale("log",base=10,subs=[0.01,0.02,0.03])    
    ax4.set_xlabel(label_x, size= 14)
    ax4.set_ylabel(label_y, size= 14)
    ax4.tick_params(labelsize = 12) # 设置坐标刻度对应数字的大小
    ax4.set_xlim([0,40])
    ax4.set_xticks([1,5,10,20,30])
    #ax4.set_yticks([-1,-0.5,0,0.5,1])  
    #=========================================================
    ax4.text(1.1,0.2,"(d)",fontsize = 20, color='black', rotation = 0)
    ax4.text(1,1.0e-3, r'$\mathrm{\delta} = 0.4583$', fontsize = 16, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
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
    fig.tight_layout() # 自动调整 subplot 间的间隙参数
    plt.savefig("E:\\WORK\\Work\\Project\\La3Ni2O7\\datashow_datpy_paper\\figs\\fig_04583_1.eps",
                dpi=300, format='eps',bbox_inches='tight') # 白边紧凑型
    #
    plt.show()