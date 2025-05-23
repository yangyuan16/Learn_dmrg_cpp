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
    dop = 36
    t = 3
    J = 1
    dim = 6000
    #
    # ninj correlation
    '''
    df_ninj_Jz04 = get_data_ninj_corre(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=0.4,dim=dim)
    df_ninj_Jz10 = get_data_ninj_corre(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=1.0,dim=dim)
    df_ninj_Jz20 = get_data_ninj_corre(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=2.0,dim=dim)
    #
    Jz_list = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,]
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
    #
    print("K_sc_yy:\n", K_sc_yy)
    print("K_sc_zz:\n",K_sc_zz)
    print("K_cdw:\n",K_cdw)
    '''
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
    ax1 = plt.axes([0.1,0.58,0.15,0.4])
    ax2 = plt.axes([0.30,0.58,0.15,0.4])
    ax3 = plt.axes([0.1,0.1,0.15,0.4])
    ax4 = plt.axes([0.30,0.1,0.15,0.4])
    ax5 = plt.axes([0.50,0.58,0.37,0.4])   
    ax6 = plt.axes([0.50,0.1,0.15,0.4])
    ax7 = plt.axes([0.685,0.1,0.185,0.4])
    # plt.subplots_adjust(left=0.1, bottom=0.1, right=0.98, top=0.96, wspace=0.32, hspace=0.26)
    # 选择子图 ax1 进行绘图
    plt.sca(ax1) ## 选择对 ax1 进行绘图
    ax1 = plt.gca()
    #
    df_sisj_Jz01 = get_data_sisj(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=0.1,dim=dim)
    df_sisj_Jz02 = get_data_sisj(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=0.2,dim=dim)
    df_sisj_Jz04 = get_data_sisj(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=0.4,dim=dim)
    df_sisj_Jz06 = get_data_sisj(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=0.6,dim=dim)
    df_sisj_Jz08 = get_data_sisj(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=0.8,dim=dim)
    df_sisj_Jz10 = get_data_sisj(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=1.0,dim=dim)
    df_sisj_Jz12 = get_data_sisj(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=1.2,dim=dim)
    #df_sisj_Jz20 = get_data_sisj(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=2.0,dim=dim)
    #--------- J_\bot = 0.1
    slope = df_sisj_Jz01["slope_exp"].values[0]
    xi = round(-1/slope,2)  
    label = r"$J_{\bot}$=%.1f $\xi_F$=%.2f"%(0.1,xi)
    L01, = ax1.plot(df_sisj_Jz01["r"],df_sisj_Jz01["corre_abs"],label=label,ls="-",lw=1.5,color="red",
             marker='o',alpha=1,markersize=6,markeredgewidth=1, markeredgecolor="red", markerfacecolor='None')
    L01_, = ax1.plot(df_sisj_Jz01["r"],df_sisj_Jz01["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
             marker='o',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- J_\bot = 0.2
    slope = df_sisj_Jz02["slope_exp"].values[0]
    xi = round(-1/slope,2)  
    label = r"$J_{\bot}$=%.1f $\xi_F$=%.2f"%(0.2,xi)
    L02, = ax1.plot(df_sisj_Jz02["r"],df_sisj_Jz02["corre_abs"],label=label,ls="-",lw=1.5,color="blue",
             marker='s',alpha=1,markersize=6,markeredgewidth=1, markeredgecolor="blue", markerfacecolor='None')
    L02_, = ax1.plot(df_sisj_Jz02["r"],df_sisj_Jz02["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
             marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- J_\bot = 0.4
    slope = df_sisj_Jz04["slope_exp"].values[0]
    xi = round(-1/slope,2)  
    label = r"$J_{\bot}$=%.1f $\xi_F$=%.2f"%(0.4,xi)
    L04, = ax1.plot(df_sisj_Jz04["r"],df_sisj_Jz04["corre_abs"],label=label,ls="-",lw=1.5,color="green",
             marker='^',alpha=1,markersize=6,markeredgewidth=1, markeredgecolor="green", markerfacecolor='None')
    L04_, = ax1.plot(df_sisj_Jz04["r"],df_sisj_Jz04["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
             marker='^',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- J_\bot = 0.6
    slope = df_sisj_Jz06["slope_exp"].values[0]
    xi = round(-1/slope,2)  
    label = r"$J_{\bot}$=%.1f $\xi_F$=%.2f"%(0.6,xi)
    L06, = ax1.plot(df_sisj_Jz06["r"],df_sisj_Jz06["corre_abs"],label=label,ls="-",lw=1.5,color="magenta",
             marker='v',alpha=1,markersize=6,markeredgewidth=1, markeredgecolor="magenta", markerfacecolor='None')
    L06_, = ax1.plot(df_sisj_Jz06["r"],df_sisj_Jz06["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
             marker='v',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 0.8
    slope = df_sisj_Jz08["slope_exp"].values[0]
    xi = round(-1/slope,2)  
    label = r"$J_{\bot}$=%.1f $\xi_F$=%.2f"%(0.8,xi)
    L08, = ax1.plot(df_sisj_Jz08["r"],df_sisj_Jz08["corre_abs"],label=label,ls="-",lw=1.5,color="cyan",
             marker='<',alpha=1,markersize=6,markeredgewidth=1, markeredgecolor="cyan", markerfacecolor='None')
    L08_, = ax1.plot(df_sisj_Jz08["r"],df_sisj_Jz08["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
             marker='<',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 1.0
    slope = df_sisj_Jz10["slope_exp"].values[0]
    xi = round(-1/slope,2)  
    label = r"$J_{\bot}$=%.1f $\xi_F$=%.2f"%(1.0,xi)
    L10, = ax1.plot(df_sisj_Jz10["r"],df_sisj_Jz10["corre_abs"],label=label,ls="-",lw=1.5,color="brown",
             marker='>',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="brown", markerfacecolor='None')
    L10_, = ax1.plot(df_sisj_Jz10["r"],df_sisj_Jz10["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
             marker='>',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 1.2
    slope = df_sisj_Jz12["slope_exp"].values[0]
    xi = round(-1/slope,2)  
    label = r"$J_{\bot}$=%.1f $\xi_F$=%.2f"%(1.2,xi)
    L12, = ax1.plot(df_sisj_Jz12["r"],df_sisj_Jz12["corre_abs"],label=label,ls="-",lw=1.5,color="olive",
             marker='h',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="olive", markerfacecolor='None')
    L12_, = ax1.plot(df_sisj_Jz12["r"],df_sisj_Jz12["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
             marker='h',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------
    ####图例设置
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 12, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[L01,L02,], loc = 4, bbox_to_anchor=(0.99, 0.778),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    legend2=plt.legend(handles=[L04,L06,L08,L10,L12], loc = 4, bbox_to_anchor=(0.64, -0.01),
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
    ax1.text(4,0.06,"(a)",fontsize = 20, color='black', rotation = 0)
    ax1.text(12,2.0e-1, r'$\mathrm{\delta} = 0.1875$', fontsize = 16, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
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
    #
    df_cicj_Jz01 = get_data_cicj(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=0.1,dim=dim)
    df_cicj_Jz02 = get_data_cicj(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=0.2,dim=dim)
    df_cicj_Jz04 = get_data_cicj(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=0.4,dim=dim)
    df_cicj_Jz06 = get_data_cicj(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=0.6,dim=dim)
    df_cicj_Jz08 = get_data_cicj(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=0.8,dim=dim)
    df_cicj_Jz10 = get_data_cicj(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=1.0,dim=dim)
    df_cicj_Jz12 = get_data_cicj(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=1.2,dim=dim)
    #df_cicj_Jz20 = get_data_cicj(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=2.0,dim=dim)
    #--------- J_\bot = 0.1
    slope = df_cicj_Jz01["slope_exp"].values[0]
    xi = round(-1/slope,2)  
    label = r"$J_{\bot}$=%.1f $\xi_G$=%.2f"%(0.1,xi)
    L01, = ax2.plot(df_cicj_Jz01["r"],df_cicj_Jz01["corre_abs"],label=label,ls="-",lw=1.5,color="red",
             marker='o',alpha=1,markersize=6,markeredgewidth=1, markeredgecolor="red", markerfacecolor='None')
    L01_, = ax2.plot(df_cicj_Jz01["r"],df_cicj_Jz01["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
             marker='o',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- J_\bot = 0.2
    slope = df_cicj_Jz02["slope_exp"].values[0]
    xi = round(-1/slope,2)  
    label = r"$J_{\bot}$=%.1f $\xi_G$=%.2f"%(0.2,xi)
    L02, = ax2.plot(df_cicj_Jz02["r"],df_cicj_Jz02["corre_abs"],label=label,ls="-",lw=1.5,color="blue",
             marker='s',alpha=1,markersize=6,markeredgewidth=1, markeredgecolor="blue", markerfacecolor='None')
    L02_, = ax2.plot(df_cicj_Jz02["r"],df_cicj_Jz02["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
             marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- J_\bot = 0.4
    slope = df_cicj_Jz04["slope_exp"].values[0]
    xi = round(-1/slope,2)  
    label = r"$J_{\bot}$=%.1f $\xi_G$=%.2f"%(0.4,xi)
    L04, = ax2.plot(df_cicj_Jz04["r"],df_cicj_Jz04["corre_abs"],label=label,ls="-",lw=1.5,color="green",
             marker='^',alpha=1,markersize=6,markeredgewidth=1, markeredgecolor="green", markerfacecolor='None')
    L04_, = ax2.plot(df_cicj_Jz04["r"],df_cicj_Jz04["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
             marker='^',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------- J_\bot = 0.6
    slope = df_cicj_Jz06["slope_exp"].values[0]
    xi = round(-1/slope,2)  
    label = r"$J_{\bot}$=%.1f $\xi_G$=%.2f"%(0.6,xi)
    L06, = ax2.plot(df_cicj_Jz06["r"],df_cicj_Jz06["corre_abs"],label=label,ls="-",lw=1.5,color="magenta",
             marker='v',alpha=1,markersize=6,markeredgewidth=1, markeredgecolor="magenta", markerfacecolor='None')
    L06_, = ax2.plot(df_cicj_Jz06["r"],df_cicj_Jz06["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
             marker='v',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 0.8
    slope = df_cicj_Jz08["slope_exp"].values[0]
    xi = round(-1/slope,2)  
    label = r"$J_{\bot}$=%.1f $\xi_G$=%.2f"%(0.8,xi)
    L08, = ax2.plot(df_cicj_Jz08["r"],df_cicj_Jz08["corre_abs"],label=label,ls="-",lw=1.5,color="cyan",
             marker='<',alpha=1,markersize=6,markeredgewidth=1, markeredgecolor="cyan", markerfacecolor='None')
    L08_, = ax2.plot(df_cicj_Jz08["r"],df_cicj_Jz08["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
             marker='<',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 1.0
    slope = df_cicj_Jz10["slope_exp"].values[0]
    xi = round(-1/slope,2)  
    label = r"$J_{\bot}$=%.1f $\xi_G$=%.2f"%(1.0,xi)
    L10, = ax2.plot(df_cicj_Jz10["r"],df_cicj_Jz10["corre_abs"],label=label,ls="-",lw=1.5,color="brown",
             marker='>',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="brown", markerfacecolor='None')
    L10_, = ax2.plot(df_cicj_Jz10["r"],df_cicj_Jz10["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
             marker='>',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 1.2
    slope = df_cicj_Jz12["slope_exp"].values[0]
    xi = round(-1/slope,2)  
    label = r"$J_{\bot}$=%.1f $\xi_G$=%.2f"%(1.2,xi)
    L12, = ax2.plot(df_cicj_Jz12["r"],df_cicj_Jz12["corre_abs"],label=label,ls="-",lw=1.5,color="olive",
             marker='h',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="olive", markerfacecolor='None')
    L12_, = ax2.plot(df_cicj_Jz12["r"],df_cicj_Jz12["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
             marker='h',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #--------
    ####图例设置
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 12, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[L01,L02,L04,], loc = 4, bbox_to_anchor=(0.99, 0.68),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    
    legend2=plt.legend(handles=[L06,L08,L10,L12], loc = 4, bbox_to_anchor=(0.64, -0.01),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    
    plt.gca().add_artist(legend1)#####把图例legend1重新加载回来

    label_x = r"r"
    label_y = "|G(r)|"
    #plt.yscale("log")
    #plt.xscale("log")
    ax2.set_yscale("log",base=10,subs=[0.01,0.02,0.03])      
    ax2.set_xlabel(label_x, size= 14)
    ax2.set_ylabel(label_y, size= 14)
    ax2.tick_params(labelsize = 12) # 设置坐标刻度对应数字的大小
    ax2.set_xlim([0,40])
    #ax2.set_xticks([0,10,20,30])
    #ax2.set_yticks([-1,-0.5,0,0.5,1])  
    #=========================================================
    ax2.text(4,0.06,"(b)",fontsize = 20, color='black', rotation = 0)
    ax2.text(12,1.0e-1, r'$\mathrm{\delta} = 0.1875$', fontsize = 16, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
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
    # load YY pairing correlation data
    df_pyy_Jz01 = get_data_pyy(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=0.1,dim=dim)
    df_pyy_Jz02 = get_data_pyy(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=0.2,dim=dim)
    df_pyy_Jz04 = get_data_pyy(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=0.4,dim=dim)
    df_pyy_Jz06 = get_data_pyy(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=0.6,dim=dim)
    df_pyy_Jz08 = get_data_pyy(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=0.8,dim=dim)
    df_pyy_Jz10 = get_data_pyy(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=1.0,dim=dim)
    df_pyy_Jz12 = get_data_pyy(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=1.2,dim=dim)
    #df_pyy_Jz20 = get_data_pyy(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=2.0,dim=dim)
    #
    #--------- J_\bot = 0.1
    slope = df_pyy_Jz01["slope_exp"].values[0]
    xi = round(-1/slope,2)
    label = r"$J_{\bot}$=%.1f $\xi_{SC}^{yy}$=%.2f"%(0.1,xi)
    Lyy01, = ax3.plot(df_pyy_Jz01["r"],df_pyy_Jz01["corre_abs"],label=label,ls="-",lw=1.5,color="red",
             marker='o',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="red", markerfacecolor='None')
    Lyy01_, = ax3.plot(df_pyy_Jz01["r"],df_pyy_Jz01["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
             marker='o',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='None')
    #--------- J_\bot = 0.2
    slope = df_pyy_Jz02["slope_exp"].values[0]
    xi = round(-1/slope,2)
    label = r"$J_{\bot}$=%.1f $\xi_{SC}^{yy}$=%.2f"%(0.2,xi)
    Lyy02, = ax3.plot(df_pyy_Jz02["r"],df_pyy_Jz02["corre_abs"],label=label,ls="-",lw=1.5,color="blue",
             marker='s',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="blue", markerfacecolor='None')
    Lyy02_, = ax3.plot(df_pyy_Jz02["r"],df_pyy_Jz02["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
             marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='None')
    #--------- J_\bot = 0.4
    slope = df_pyy_Jz04["slope_exp"].values[0]
    xi = round(-1/slope,2)
    label = r"$J_{\bot}$=%.1f $\xi_{SC}^{yy}$=%.2f"%(0.4,xi)
    Lyy04, = ax3.plot(df_pyy_Jz04["r"],df_pyy_Jz04["corre_abs"],label=label,ls="-",lw=1.5,color="green",
             marker='^',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="green", markerfacecolor='None')
    Lyy04_, = ax3.plot(df_pyy_Jz04["r"],df_pyy_Jz04["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
             marker='^',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='None')
    #--------- J_\bot = 0.6
    slope = df_pyy_Jz06["slope_exp"].values[0]
    xi = round(-1/slope,2)
    label = r"$J_{\bot}$=%.1f $\xi_{SC}^{yy}$=%.2f"%(0.6,xi)
    Lyy06, = ax3.plot(df_pyy_Jz06["r"],df_pyy_Jz06["corre_abs"],label=label,ls="-",lw=1.5,color="magenta",
             marker='v',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="magenta", markerfacecolor='None')
    Lyy06_, = ax3.plot(df_pyy_Jz06["r"],df_pyy_Jz06["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
             marker='v',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='None')
    #--------- J_\bot = 0.8
    slope = df_pyy_Jz08["slope_exp"].values[0]
    xi = round(-1/slope,2)
    label = r"$J_{\bot}$=%.1f $\xi_{SC}^{yy}$=%.2f"%(0.8,xi)
    Lyy08, = ax3.plot(df_pyy_Jz08["r"],df_pyy_Jz08["corre_abs"],label=label,ls="-",lw=1.5,color="cyan",
             marker='<',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="cyan", markerfacecolor='None')
    Lyy08_, = ax3.plot(df_pyy_Jz08["r"],df_pyy_Jz08["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
             marker='<',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='None')
    #--------- J_\bot = 1.0
    slope = df_pyy_Jz10["slope_exp"].values[0]
    xi = round(-1/slope,2)
    label = r"$J_{\bot}$=%.1f $\xi_{SC}^{yy}$=%.2f"%(1.0,xi)
    Lyy10, = ax3.plot(df_pyy_Jz10["r"],df_pyy_Jz10["corre_abs"],label=label,ls="-",lw=1.5,color="brown",
             marker='>',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="brown", markerfacecolor='None')
    Lyy10_, = ax3.plot(df_pyy_Jz10["r"],df_pyy_Jz10["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
             marker='>',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='None')
    #--------- J_\bot = 1.2
    slope = df_pyy_Jz12["slope_exp"].values[0]
    xi = round(-1/slope,2)
    label = r"$J_{\bot}$=%.1f $\xi_{SC}^{yy}$=%.2f"%(1.2,xi)
    Lyy12, = ax3.plot(df_pyy_Jz12["r"],df_pyy_Jz12["corre_abs"],label=label,ls="-",lw=1.5,color="olive",
             marker='h',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="olive", markerfacecolor='None')
    Lyy12_, = ax3.plot(df_pyy_Jz12["r"],df_pyy_Jz12["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
             marker='h',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='None')
    #
    ####图例设置
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 12, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[Lyy01,Lyy02,], loc = 4, bbox_to_anchor=(1.05, 0.82),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    legend2=plt.legend(handles=[Lyy04,Lyy06,Lyy08,Lyy10,Lyy12], loc = 4, bbox_to_anchor=(0.65, -0.02),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    plt.gca().add_artist(legend1)
    label_x = r"r"
    label_y = "$P^{yy}(r)$"
    #plt.yscale("log")
    #plt.xscale("log")
    #ax3.set_xscale("log",base=10,subs=[0.01,0.02,0.03]) 
    ax3.set_yscale("log",base=10,subs=[0.01,0.02,0.03])      
    ax3.set_xlabel(label_x, size= 16)
    ax3.set_ylabel(label_y, size= 16)
    ax3.tick_params(labelsize = 14) # 设置坐标刻度对应数字的大小
    ax3.set_xlim([0,40])
    #ax3.set_xticks([0,10,20,30])
    #ax3.set_yticks([-1,-0.5,0,0.5,1]) 
    #
    #=========================================================
    ax3.text(8,0.004,"(d)",fontsize = 20, color='black', rotation = 0)
    ax3.text(0,0.15e-10, r'$\mathrm{\delta} = 0.1875$', fontsize = 16, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
    #ax3.text(2,0.3e-6, r'$J_{\bot}=2.0$', fontsize = 20, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
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
    # load ZZ pairing correlation data
    df_pzz_Jz01 = get_data_pzz(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=0.1,dim=dim)
    df_pzz_Jz02 = get_data_pzz(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=0.2,dim=dim)
    df_pzz_Jz04 = get_data_pzz(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=0.4,dim=dim)
    df_pzz_Jz06 = get_data_pzz(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=0.6,dim=dim)
    df_pzz_Jz08 = get_data_pzz(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=0.8,dim=dim)
    df_pzz_Jz10 = get_data_pzz(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=1.0,dim=dim)
    df_pzz_Jz12 = get_data_pzz(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=1.2,dim=dim)
    #df_pzz_Jz20 = get_data_pzz(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=2.0,dim=dim)
    #--------- J_\bot = 0.1
    slope = df_pzz_Jz01["slope_exp"].values[0]
    xi = round(-1/slope,2)
    label = r"$J_{\bot}$=%.1f $\xi_{SC}^{zz}$=%.2f"%(0.1,xi)
    Lzz01, = ax4.plot(df_pzz_Jz01["r"],df_pzz_Jz01["corre_abs"],label=label,ls="-",lw=1.5,color="red",
             marker='o',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="red", markerfacecolor='None')
    Lzz01_, = ax4.plot(df_pzz_Jz01["r"],df_pzz_Jz01["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
             marker='o',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='None')
    #--------- J_\bot = 0.2
    slope = df_pzz_Jz02["slope_exp"].values[0]
    xi = round(-1/slope,2)
    label = r"$J_{\bot}$=%.1f $\xi_{SC}^{zz}$=%.2f"%(0.2,xi)
    Lzz02, = ax4.plot(df_pzz_Jz02["r"],df_pzz_Jz02["corre_abs"],label=label,ls="-",lw=1.5,color="blue",
             marker='s',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="blue", markerfacecolor='None')
    Lzz02_, = ax4.plot(df_pzz_Jz02["r"],df_pzz_Jz02["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
             marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='None')
    #--------- J_\bot = 0.4
    slope = df_pzz_Jz04["slope_exp"].values[0]
    xi = round(-1/slope,2)
    label = r"$J_{\bot}$=%.1f $\xi_{SC}^{zz}$=%.2f"%(0.4,xi)
    Lzz04, = ax4.plot(df_pzz_Jz04["r"],df_pzz_Jz04["corre_abs"],label=label,ls="-",lw=1.5,color="green",
             marker='^',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="green", markerfacecolor='None')
    Lzz04_, = ax4.plot(df_pzz_Jz04["r"],df_pzz_Jz04["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
             marker='^',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='None')
    #--------- J_\bot = 0.6
    slope = df_pzz_Jz06["slope_exp"].values[0]
    xi = round(-1/slope,2)
    label = r"$J_{\bot}$=%.1f $\xi_{SC}^{zz}$=%.2f"%(0.6,xi)
    Lzz06, = ax4.plot(df_pzz_Jz06["r"],df_pzz_Jz06["corre_abs"],label=label,ls="-",lw=1.5,color="magenta",
             marker='v',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="magenta", markerfacecolor='None')
    Lzz06_, = ax4.plot(df_pzz_Jz06["r"],df_pzz_Jz06["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
             marker='v',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='None')
    #--------- J_\bot = 0.8
    slope = df_pzz_Jz08["slope_exp"].values[0]
    xi = round(-1/slope,2)
    #label = r"$J_{\bot}$=%.1f $\xi_{SC}^{zz}$=%.2f"%(0.8,xi)
    label = r"$J_{\bot}$=%.1f"%(0.8)
    Lzz08, = ax4.plot(df_pzz_Jz08["r"],df_pzz_Jz08["corre_abs"],label=label,ls="-",lw=1.5,color="cyan",
             marker='<',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="cyan", markerfacecolor='None')
    #Lzz08_, = ax4.plot(df_pzz_Jz08["r"],df_pzz_Jz08["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
    #         marker='<',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='None')
    #--------- J_\bot = 1.0
    slope = df_pzz_Jz10["slope_exp"].values[0]
    xi = round(-1/slope,2)
    #label = r"$J_{\bot}$=%.1f $\xi_{SC}^{zz}$=%.2f"%(1.0,xi)
    label = r"$J_{\bot}$=%.1f"%(1.0)
    Lzz10, = ax4.plot(df_pzz_Jz10["r"],df_pzz_Jz10["corre_abs"],label=label,ls="-",lw=1.5,color="brown",
             marker='>',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="brown", markerfacecolor='None')
    #Lzz10_, = ax4.plot(df_pzz_Jz10["r"],df_pzz_Jz10["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
    #         marker='>',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='None')
    #--------- J_\bot = 1.2
    slope = df_pzz_Jz12["slope_exp"].values[0]
    xi = round(-1/slope,2)
    #label = r"$J_{\bot}$=%.1f $\xi_{SC}^{zz}$=%.2f"%(1.2,xi)
    label = r"$J_{\bot}$=%.1f"%(1.2)
    Lzz12, = ax4.plot(df_pzz_Jz12["r"],df_pzz_Jz12["corre_abs"],label=label,ls="-",lw=1.5,color="olive",
             marker='h',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="olive", markerfacecolor='None')
    #Lzz12_, = ax4.plot(df_pzz_Jz12["r"],df_pzz_Jz12["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
    #         marker='h',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='None')
    #
    ####图例设置
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 12, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[Lzz01,Lzz02,Lzz04,Lzz06], loc = 4, bbox_to_anchor=(0.65, -0.02),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    legend2=plt.legend(handles=[Lzz08,Lzz10,Lzz12], loc = 4, bbox_to_anchor=(0.4, 0.28),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    plt.gca().add_artist(legend1)
    label_x = r"r"
    label_y = "$P^{zz}(r)$"
    #plt.yscale("log")
    #plt.xscale("log")
    #ax4.set_xscale("log",base=10,subs=[0.01,0.02,0.03]) 
    ax4.set_yscale("log",base=10,subs=[0.01,0.02,0.03])      
    ax4.set_xlabel(label_x, size= 16)
    ax4.set_ylabel(label_y, size= 16)
    ax4.tick_params(labelsize = 14) # 设置坐标刻度对应数字的大小
    ax4.set_xlim([0,40])
    #ax4.set_xticks([0,10,20,30])
    #ax4.set_yticks([-1,-0.5,0,0.5,1]) 
    #
    #=========================================================
    ax4.text(8,0.01,"(e)",fontsize = 20, color='black', rotation = 0)
    ax4.text(20,1e-2, r'$\mathrm{\delta} = 0.1875$', fontsize = 16, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
    #ax4.text(2,0.3e-6, r'$J_{\bot}=2.0$', fontsize = 20, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
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
    #=============================================================================================
    # 选择子图 ax6 进行绘图
    plt.sca(ax6) ## 选择对 ax1 进行绘图
    ax6 = plt.gca()
    #--------- J_\bot = 0.1
    slope = df_pzz_Jz01["slope_pow"].values[0]
    Ksc = round(-slope,2) 
    label = r"$J_{\bot}$=%.1f $K_{SC}^{zz}$=%.2f"%(0.1,Ksc)
    Lpzz01, = ax6.plot(df_pzz_Jz01["r"],df_pzz_Jz01["corre_abs"],label=label,ls="-",lw=1.5,color="red",
             marker='o',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="red", markerfacecolor='None')
    Lpzz01_, = ax6.plot(df_pzz_Jz01["r"],df_pzz_Jz01["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='o',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='None')
    #--------- J_\bot = 0.2
    slope = df_pzz_Jz02["slope_pow"].values[0]
    Ksc = round(-slope,2)
    label = r"$J_{\bot}$=%.1f $K_{SC}^{zz}$=%.2f"%(0.2,Ksc)
    Lpzz02, = ax6.plot(df_pzz_Jz02["r"],df_pzz_Jz02["corre_abs"],label=label,ls="-",lw=1.5,color="blue",
             marker='s',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="blue", markerfacecolor='None')
    Lpzz02_, = ax6.plot(df_pzz_Jz02["r"],df_pzz_Jz02["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='None')
    #--------- J_\bot = 0.4
    slope = df_pzz_Jz04["slope_pow"].values[0]
    Ksc = round(-slope,2)
    label = r"$J_{\bot}$=%.1f $K_{SC}^{zz}$=%.2f"%(0.4,Ksc)
    Lpzz04, = ax6.plot(df_pzz_Jz04["r"],df_pzz_Jz04["corre_abs"],label=label,ls="-",lw=1.5,color="green",
             marker='^',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="green", markerfacecolor='None')
    Lpzz04_, = ax6.plot(df_pzz_Jz04["r"],df_pzz_Jz04["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='^',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='None')
    #--------- J_\bot = 0.6
    slope = df_pzz_Jz06["slope_pow"].values[0]
    Ksc = round(-slope,2)
    label = r"$J_{\bot}$=%.1f $K_{SC}^{zz}$=%.2f"%(0.6,Ksc)
    Lpzz06, = ax6.plot(df_pzz_Jz06["r"],df_pzz_Jz06["corre_abs"],label=label,ls="-",lw=1.5,color="magenta",
             marker='v',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="magenta", markerfacecolor='None')
    Lpzz06_, = ax6.plot(df_pzz_Jz06["r"],df_pzz_Jz06["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='v',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='None')
    #--------- J_\bot = 0.8
    slope = df_pzz_Jz08["slope_pow"].values[0]
    Ksc = round(-slope,2)
    label = r"$J_{\bot}$=%.1f $K_{SC}^{zz}$=%.2f"%(0.8,Ksc)
    Lpzz08, = ax6.plot(df_pzz_Jz08["r"],df_pzz_Jz08["corre_abs"],label=label,ls="-",lw=1.5,color="cyan",
             marker='<',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="cyan", markerfacecolor='None')
    Lpzz08_, = ax6.plot(df_pzz_Jz08["r"],df_pzz_Jz08["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='<',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='None')
    #--------- J_\bot = 1.0
    slope = df_pzz_Jz10["slope_pow"].values[0]
    Ksc = round(-slope,2)
    label = r"$J_{\bot}$=%.1f $K_{SC}^{zz}$=%.2f"%(1.0,Ksc)
    Lpzz10, = ax6.plot(df_pzz_Jz10["r"],df_pzz_Jz10["corre_abs"],label=label,ls="-",lw=1.5,color="brown",
             marker='>',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="brown", markerfacecolor='None')
    Lpzz10_, = ax6.plot(df_pzz_Jz10["r"],df_pzz_Jz10["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='>',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='None')
    #--------- J_\bot = 1.2
    slope = df_pzz_Jz12["slope_pow"].values[0]
    Ksc = round(-slope,2)
    label = r"$J_{\bot}$=%.1f $K_{SC}^{zz}$=%.2f"%(1.2,Ksc)
    Lpzz12, = ax6.plot(df_pzz_Jz12["r"],df_pzz_Jz12["corre_abs"],label=label,ls="-",lw=1.5,color="olive",
             marker='h',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="olive", markerfacecolor='None')
    Lpzz12_, = ax6.plot(df_pzz_Jz12["r"],df_pzz_Jz12["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='h',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='None')
    #
    ####图例设置
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 12, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[Lpzz01,Lpzz02,Lpzz04,Lpzz06], loc = 4, bbox_to_anchor=(0.65, 0.18),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    legend2=plt.legend(handles=[Lpzz08,Lpzz10,Lpzz12], loc = 4, bbox_to_anchor=(0.65, -0.04),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    plt.gca().add_artist(legend1)
    label_x = r"r"
    label_y = "$P^{zz}(r)$"
    #plt.yscale("log")
    #plt.xscale("log")
    ax6.set_xscale("log",base=10,subs=[0.01,0.02,0.03]) 
    ax6.set_yscale("log",base=10,subs=[0.01,0.02,0.03])      
    ax6.set_xlabel(label_x, size= 16)
    ax6.set_ylabel(label_y, size= 16)
    ax6.tick_params(labelsize = 14) # 设置坐标刻度对应数字的大小
    ax6.set_xlim([0,40])
    ax6.set_xticks([5,10,15,20,30,])
    #ax6.set_yticks([-1,-0.5,0,0.5,1]) 
    #
    #=========================================================
    ax6.text(25,0.01,"(f)",fontsize = 20, color='black', rotation = 0)
    ax6.text(4,1e-2, r'$\mathrm{\delta} = 0.1875$', fontsize = 16, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
    #ax6.text(2,0.3e-6, r'$J_{\bot}=2.0$', fontsize = 20, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
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
    #
    # 选择子图 ax7 进行绘图
    plt.sca(ax7) ## 选择对 ax7 进行绘图
    ax7 = plt.gca()
    #--------------------------------------------------------------------------
    # load the entropy data
    df_ent_01 = get_entropy_Lz2(Lz=2,Ly=2,Lx=48,dop=36,t=3,J=1,Jz= 0.1,dim=6000)
    df_ent_02 = get_entropy_Lz2(Lz=2,Ly=2,Lx=48,dop=36,t=3,J=1,Jz= 0.2,dim=6000)
    df_ent_04 = get_entropy_Lz2(Lz=2,Ly=2,Lx=48,dop=36,t=3,J=1,Jz= 0.4,dim=6000)
    df_ent_06 = get_entropy_Lz2(Lz=2,Ly=2,Lx=48,dop=36,t=3,J=1,Jz= 0.6,dim=6000)
    df_ent_08 = get_entropy_Lz2(Lz=2,Ly=2,Lx=48,dop=36,t=3,J=1,Jz= 0.8,dim=6000)
    df_ent_10 = get_entropy_Lz2(Lz=2,Ly=2,Lx=48,dop=36,t=3,J=1,Jz= 1.0,dim=6000)
    df_ent_12 = get_entropy_Lz2(Lz=2,Ly=2,Lx=48,dop=36,t=3,J=1,Jz= 1.2,dim=6000)
    #df_14 = get_entropy_Lz2(Lz=2,Ly=2,Lx=48,dop=36,t=3,J=1,Jz= 1.4,dim=6000)
    # J_\bot = 0.1
    slope = df_ent_01["slope"].values[0]
    intercept = df_ent_01["intercept"].values[0]
    label_fitdata = r"$J_{\bot}$=%.1f c = %.2f g= %.2f"%(0.1,slope,intercept)
    Lent01, = ax7.plot(df_ent_01["logr"].values,df_ent_01["entropy"].values,label=label_fitdata,ls="-",lw=1.5,color="red",
             marker='o',alpha=1,markersize=8,markeredgewidth=1.0, markeredgecolor="magenta", markerfacecolor='None')
    Lent01_fit, = ax7.plot(df_ent_01["logr"].values[6:-1],df_ent_01["fitentropy"].values[6:-1], label=label_fitdata,ls="--",lw=1.5,color="k",
             marker='o',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k",markerfacecolor='None')
    # J_\bot = 0.2
    slope = df_ent_02["slope"].values[0]
    intercept = df_ent_02["intercept"].values[0]
    label_fitdata = r"$J_{\bot}$=%.1f c = %.2f g= %.2f"%(0.2,slope,intercept)
    Lent02, = ax7.plot(df_ent_02["logr"].values,df_ent_02["entropy"].values,label=label_fitdata,ls="-",lw=1.5,color="blue",
             marker='s',alpha=1,markersize=8,markeredgewidth=1.0, markeredgecolor="blue", markerfacecolor='None')
    Lent02_fit, = ax7.plot(df_ent_02["logr"].values[6:-1],df_ent_02["fitentropy"].values[6:-1], label=label_fitdata,ls="--",lw=1.5,color="k",
             marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k",markerfacecolor='None')
    # J_\bot = 0.4
    slope = df_ent_04["slope"].values[0]
    intercept = df_ent_04["intercept"].values[0]
    label_fitdata = r"$J_{\bot}$=%.1f c = %.2f g= %.2f"%(0.4,slope,intercept)
    Lent04, = ax7.plot(df_ent_04["logr"].values,df_ent_04["entropy"].values,label=label_fitdata,ls="-",lw=1.5,color="green",
             marker='^',alpha=1,markersize=8,markeredgewidth=1.0, markeredgecolor="green", markerfacecolor='None')
    Lent04_fit, = ax7.plot(df_ent_04["logr"].values[6:-1],df_ent_04["fitentropy"].values[6:-1], label=label_fitdata,ls="--",lw=1.5,color="k",
             marker='^',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k",markerfacecolor='None')
    # J_\bot = 0.6
    slope = df_ent_06["slope"].values[0]
    intercept = df_ent_06["intercept"].values[0]
    label_fitdata = r"$J_{\bot}$=%.1f c = %.2f g= %.2f"%(0.6,slope,intercept)
    Lent06, = ax7.plot(df_ent_06["logr"].values,df_ent_06["entropy"].values,label=label_fitdata,ls="-",lw=1.5,color="magenta",
             marker='v',alpha=1,markersize=8,markeredgewidth=1.0, markeredgecolor="magenta", markerfacecolor='None')
    Lent06_fit, = ax7.plot(df_ent_06["logr"].values[6:-1],df_ent_06["fitentropy"].values[6:-1], label=label_fitdata,ls="--",lw=1.5,color="k",
             marker='v',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k",markerfacecolor='None')
    # J_\bot = 0.8
    slope = df_ent_08["slope"].values[0]
    intercept = df_ent_08["intercept"].values[0]
    label_fitdata = r"$J_{\bot}$=%.1f c = %.2f g= %.2f"%(0.8,slope,intercept)
    Lent08, = ax7.plot(df_ent_08["logr"].values,df_ent_08["entropy"].values,label=label_fitdata,ls="-",lw=1.5,color="cyan",
             marker='<',alpha=1,markersize=8,markeredgewidth=1.0, markeredgecolor="cyan", markerfacecolor='None')
    Lent08_fit, = ax7.plot(df_ent_08["logr"].values[6:-1],df_ent_08["fitentropy"].values[6:-1], label=label_fitdata,ls="--",lw=1.5,color="k",
             marker='<',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k",markerfacecolor='None')
    # J_\bot = 1.0
    slope = df_ent_10["slope"].values[0]
    intercept = df_ent_10["intercept"].values[0]
    label_fitdata = r"$J_{\bot}$=%.1f c = %.2f g= %.2f"%(1.0,slope,intercept)
    Lent10, = ax7.plot(df_ent_10["logr"].values,df_ent_10["entropy"].values,label=label_fitdata,ls="-",lw=1.5,color="brown",
             marker='>',alpha=1,markersize=8,markeredgewidth=1.0, markeredgecolor="brown", markerfacecolor='None')
    Lent10_fit, = ax7.plot(df_ent_10["logr"].values,df_ent_10["fitentropy"].values, label=label_fitdata,ls="--",lw=1.5,color="k",
             marker='>',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k",markerfacecolor='None')
    # J_\bot = 1.2
    slope = df_ent_12["slope"].values[0]
    intercept = df_ent_12["intercept"].values[0]
    label_fitdata = r"$J_{\bot}$=%.1f c = %.2f g= %.2f"%(1.2,slope,intercept)
    Lent12, = ax7.plot(df_ent_12["logr"].values,df_ent_12["entropy"].values,label=label_fitdata,ls="-",lw=1.5,color="olive",
             marker='h',alpha=1,markersize=8,markeredgewidth=1.0, markeredgecolor="olive", markerfacecolor='None')
    Lent12_fit, = ax7.plot(df_ent_12["logr"].values,df_ent_12["fitentropy"].values, label=label_fitdata,ls="--",lw=1.5,color="k",
             marker='h',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k",markerfacecolor='None')
    #
    ####图例设置
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 12, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[Lent01,Lent02,Lent04,Lent06], loc = 4, bbox_to_anchor=(0.65, 0.18),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    legend2=plt.legend(handles=[Lent08,Lent10,Lent12], loc = 4, bbox_to_anchor=(0.65, -0.04),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    plt.gca().add_artist(legend1)
    label_x = r"(1/6)log((Lx/$\pi$)sin(x$\pi$/Lx))"
    label_y = "S(x)"
    #plt.yscale("log")
    #plt.xscale("log")
    #ax7.set_xscale("log",base=10,subs=[0.01,0.02,0.03]) 
    #ax7.set_yscale("log",base=10,subs=[0.01,0.02,0.03])      
    ax7.set_xlabel(label_x, size= 16)
    ax7.set_ylabel(label_y, size= 16)
    ax7.tick_params(labelsize = 14) # 设置坐标刻度对应数字的大小
    #ax7.set_xlim([0,40])
    #ax7.set_xticks([5,10,15,20,30,])
    #ax6.set_yticks([-1,-0.5,0,0.5,1]) 
    #
    #=========================================================
    ax7.text(0.2,5.15,"(g)",fontsize = 20, color='black', rotation = 0)
    ax7.text(0.2,4.9, r'$\mathrm{\delta} = 0.1875$', fontsize = 16, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
    #ax7.text(2,0.3e-6, r'$J_{\bot}=2.0$', fontsize = 20, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
    #=========================================================
    # 坐标轴设置第一层
    labels = ax7.get_xticklabels() + ax7.get_yticklabels()
    #[label.set_fontname('Times New Roman') for label in labels]###设置ticket labled的字体格式
    ax7.xaxis.set_minor_locator(MultipleLocator(5))###设置次刻度的间隔
    #ax7.yaxis.set_minor_locator(MultipleLocator(10))###设置次刻度的间隔
    ax7.xaxis.set_major_formatter(FormatStrFormatter('%1.2f'))###设置X轴标签文本格式
    #ax7.yaxis.set_major_formatter(FormatStrFormatter('%1.1f'))###设置Y轴标签文本格式
    #
    #ax7.yaxis.get_major_locator().set_params(numticks=99)
    #ax7.yaxis.get_minor_locator().set_params(numticks=13, subs=[.2,.4,.6,.8]) # 将次要刻度显示出来 
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
    #-----------------------------------------------------------------------------------------------
    plt.sca(ax5) ## 选择对 ax5 进行绘图
    ax5 = plt.gca()
    df_ni_Jz01 = get_data_ni(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=0.1,dim=dim)
    df_ni_Jz02 = get_data_ni(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=0.2,dim=dim)
    df_ni_Jz04 = get_data_ni(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=0.4,dim=dim)
    df_ni_Jz06 = get_data_ni(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=0.6,dim=dim)
    df_ni_Jz08 = get_data_ni(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=0.8,dim=dim)
    df_ni_Jz10 = get_data_ni(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=1.0,dim=dim)
    df_ni_Jz12 = get_data_ni(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=1.2,dim=dim)
    #
    L, =ax5.plot([0,50],[0.8125,0.8125],label=" ",ls="--",lw=1.5,color="k",
             marker='o',alpha=1,markersize=0,markeredgewidth=1.5, markeredgecolor="red", markerfacecolor='None')
    label = r"$J_{\bot}$=0.1"
    L01, =ax5.plot(df_ni_Jz01["rmean"].values,df_ni_Jz01["dymean"].values,label=label,ls="--",lw=1.5,color="red",
             marker='o',alpha=1,markersize=9,markeredgewidth=1.5, markeredgecolor="red", markerfacecolor='None')
    label = r"$J_{\bot}$=0.2"
    L02, =ax5.plot(df_ni_Jz02["rmean"].values,df_ni_Jz02["dymean"].values,label=label,ls="--",lw=1.5,color="blue",
             marker='s',alpha=1,markersize=9,markeredgewidth=1.5, markeredgecolor="blue", markerfacecolor='None')
    label = r"$J_{\bot}$=0.4"
    L04, =ax5.plot(df_ni_Jz04["rmean"].values,df_ni_Jz04["dymean"].values,label=label,ls="--",lw=1.5,color="green",
             marker='^',alpha=1,markersize=9,markeredgewidth=1.5, markeredgecolor="green", markerfacecolor='None')
    label = r"$J_{\bot}$=0.6"
    L06, =ax5.plot(df_ni_Jz06["rmean"].values,df_ni_Jz06["dymean"].values,label=label,ls="--",lw=1.5,color="magenta",
             marker='v',alpha=1,markersize=9,markeredgewidth=1.5, markeredgecolor="magenta", markerfacecolor='None')
    label = r"$J_{\bot}$=0.8"
    L08, =ax5.plot(df_ni_Jz08["rmean"].values,df_ni_Jz08["dymean"].values,label=label,ls="--",lw=1.5,color="cyan",
             marker='<',alpha=1,markersize=9,markeredgewidth=1.5, markeredgecolor="cyan", markerfacecolor='None')
    label = r"$J_{\bot}$=1.0"
    L10, =ax5.plot(df_ni_Jz10["rmean"].values,df_ni_Jz10["dymean"].values,label=label,ls="--",lw=1.5,color="brown",
             marker='>',alpha=1,markersize=9,markeredgewidth=1.5, markeredgecolor="brown", markerfacecolor='None')
    label = r"$J_{\bot}$=1.2"
    L12, =ax5.plot(df_ni_Jz12["rmean"].values,df_ni_Jz12["dymean"].values,label=label,ls="--",lw=1.5,color="olive",
             marker='h',alpha=1,markersize=9,markeredgewidth=1.5, markeredgecolor="olive", markerfacecolor='None')
    #label = r"$J_{\bot}$=2.0"
    #L20, =ax5.plot(df_ni_Jz20["rmean"].values,df_ni_Jz20["dymean"].values,label=label,ls="--",lw=1.5,color="green",
    #         marker='v',alpha=1,markersize=9,markeredgewidth=1.5, markeredgecolor="green", markerfacecolor='None')
    ####图例设置
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 15, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[L01,L02,L04,L06,L08,L10,L12], loc = 4, bbox_to_anchor=(0.96, 0.51),
                       ncol = 4,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)    
    #legend2=plt.legend(handles=[L04g,L10g,L20g], loc = 4, bbox_to_anchor=(0.64, 0.02),
    #                   ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    #plt.gca().add_artist(legend1)#####把图例legend1重新加载回来
    #
    label_x = r"x"
    label_y = "n(x)"
    #plt.yscale("log")
    #plt.xscale("log")
    ax5.set_xlabel(label_x, size= 16)
    ax5.set_ylabel(label_y, size= 16)
    ax5.tick_params(labelsize = 15) # 设置坐标刻度对应数字的大小
    ax5.set_xlim([0,48])
    #ax5.set_ylim([])
    ax5.set_xticks([0,10,20,30,40,48])
    #ax5.set_yticks([-1,-0.5,0,0.5,1]) 
    #=========================================================
    ax5.text(4,0.9375,"(c)",fontsize = 20, color='black', rotation = 0)
    ax5.text(20,0.9375, r'$\mathrm{\delta} = 0.1875$', fontsize = 16, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
    #=======================================================================================
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
    #
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
    #=========================================================================
    fig.tight_layout() # 自动调整 subplot 间的间隙参数
    plt.savefig("E:\\WORK\\Work\\Project\\La3Ni2O7\\datashow_datpy_paper\\figs\\fig_01875_1.eps",
                dpi=300, format='eps',bbox_inches='tight') # 白边紧凑型
    plt.show()