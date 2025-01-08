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
    #---------------plot logr-r fig and logr-logr fig-----------------------------
    fig = plt.figure(figsize=(16,6))
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
    ax1 = plt.axes([0.1,0.1,0.25,0.8])
    ax2 = plt.axes([0.4,0.1,0.25,0.8])
    ax3 = plt.axes([0.7,0.1,0.25,0.8])
    # plt.subplots_adjust(left=0.1, bottom=0.1, right=0.98, top=0.96, wspace=0.32, hspace=0.26)
    #--------------------------------------------------------------------------
    # 选择子图 ax1 进行绘图
    plt.sca(ax1) ## 选择对 ax1 进行绘图
    ax1 = plt.gca()
    # load the single layer entropy data
    df_SL = get_data_Lz1(Lz=1,Ly=2,Lx=48,dop=36,t=3,J=1,Jz= 0,dim=6000)
    # load the entropy data
    df_ent_001 = get_entropy_Lz2(Lz=2,Ly=2,Lx=48,dop=72,t=3,J=1,Jz= 0.01,dim=6000)
    df_ent_01 = get_entropy_Lz2(Lz=2,Ly=2,Lx=48,dop=72,t=3,J=1,Jz= 0.1,dim=6000)
    df_ent_04 = get_entropy_Lz2(Lz=2,Ly=2,Lx=48,dop=72,t=3,J=1,Jz= 0.4,dim=6000)
    df_ent_08 = get_entropy_Lz2(Lz=2,Ly=2,Lx=48,dop=72,t=3,J=1,Jz= 0.8,dim=6000)
    df_ent_10 = get_entropy_Lz2(Lz=2,Ly=2,Lx=48,dop=72,t=3,J=1,Jz= 1.0,dim=6000)
    df_ent_26 = get_entropy_Lz2(Lz=2,Ly=2,Lx=48,dop=72,t=3,J=1,Jz= 2.6,dim=6000)
    #df_14 = get_entropy_Lz2(Lz=2,Ly=2,Lx=48,dop=36,t=3,J=1,Jz= 1.4,dim=6000)
    # plot single layer 
    slope = df_SL["slope"].values[0]
    intercept = df_SL["intercept"].values[0]
    label_fitdata = "Single Layer t-J, c = %.2f g=%.2f"%(slope,intercept)
    LSL, = ax1.plot(df_SL["logr"].values,df_SL["entropy"].values,label=label_fitdata,ls="-",lw=1.5,color="red",
             marker='o',alpha=1,markersize=8,markeredgewidth=1.0, markeredgecolor="red",markerfacecolor='k')
    LSL_fit, = ax1.plot(df_SL["logr"].values,df_SL["fitentropy"].values,label=label_fitdata,ls="--",lw=1.5,color="red",
             marker='o',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="red",
             markerfacecolor='w')
    # plot Bilayer layer
    # J_\bot = 0.1
    slope = df_ent_001["slope"].values[0]
    intercept = df_ent_001["intercept"].values[0]
    label_fitdata = r"$J_{\bot}$=%.1f c = %.2f g= %.2f"%(0.01,slope,intercept)
    Lent001, = ax1.plot(df_ent_001["logr"].values,df_ent_001["entropy"].values,label=label_fitdata,ls="-",lw=1.5,color="magenta",
             marker='D',alpha=1,markersize=8,markeredgewidth=1.0, markeredgecolor="magenta", markerfacecolor='None')
    Lent001_fit, = ax1.plot(df_ent_001["logr"].values[2:-1],df_ent_001["fitentropy"].values[2:-1], label=label_fitdata,ls="--",lw=1.5,color="k",
             marker='D',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k",markerfacecolor='None')
    # J_\bot = 0.1
    slope = df_ent_01["slope"].values[0]
    intercept = df_ent_01["intercept"].values[0]
    label_fitdata = r"$J_{\bot}$=%.1f c = %.2f g= %.2f"%(0.1,slope,intercept)
    Lent01, = ax1.plot(df_ent_01["logr"].values,df_ent_01["entropy"].values,label=label_fitdata,ls="-",lw=1.5,color="blue",
             marker='s',alpha=1,markersize=8,markeredgewidth=1.0, markeredgecolor="blue", markerfacecolor='None')
    Lent01_fit, = ax1.plot(df_ent_01["logr"].values[2:-1],df_ent_01["fitentropy"].values[2:-1], label=label_fitdata,ls="--",lw=1.5,color="k",
             marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k",markerfacecolor='None')
    # J_\bot = 0.4
    slope = df_ent_04["slope"].values[0]
    intercept = df_ent_04["intercept"].values[0]
    label_fitdata = r"$J_{\bot}$=%.1f c = %.2f g= %.2f"%(0.4,slope,intercept)
    Lent04, = ax1.plot(df_ent_04["logr"].values,df_ent_04["entropy"].values,label=label_fitdata,ls="-",lw=1.5,color="green",
             marker='^',alpha=1,markersize=8,markeredgewidth=1.0, markeredgecolor="green", markerfacecolor='None')
    Lent04_fit, = ax1.plot(df_ent_04["logr"].values[2:-1],df_ent_04["fitentropy"].values[2:-1], label=label_fitdata,ls="--",lw=1.5,color="k",
             marker='^',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k",markerfacecolor='None')
    # J_\bot = 0.8
    slope = df_ent_08["slope"].values[0]
    intercept = df_ent_08["intercept"].values[0]
    label_fitdata = r"$J_{\bot}$=%.1f c = %.2f g= %.2f"%(0.8,slope,intercept)
    Lent08, = ax1.plot(df_ent_08["logr"].values,df_ent_08["entropy"].values,label=label_fitdata,ls="-",lw=1.5,color="cyan",
             marker='<',alpha=1,markersize=8,markeredgewidth=1.0, markeredgecolor="cyan", markerfacecolor='None')
    Lent08_fit, = ax1.plot(df_ent_08["logr"].values[2:-1],df_ent_08["fitentropy"].values[2:-1], label=label_fitdata,ls="--",lw=1.5,color="k",
             marker='<',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k",markerfacecolor='None')
    # J_\bot = 1.0
    slope = df_ent_10["slope"].values[0]
    intercept = df_ent_10["intercept"].values[0]
    label_fitdata = r"$J_{\bot}$=%.1f c = %.2f g= %.2f"%(1.0,slope,intercept)
    Lent10, = ax1.plot(df_ent_10["logr"].values,df_ent_10["entropy"].values,label=label_fitdata,ls="-",lw=1.5,color="brown",
             marker='>',alpha=1,markersize=8,markeredgewidth=1.0, markeredgecolor="brown", markerfacecolor='None')
    Lent10_fit, = ax1.plot(df_ent_10["logr"].values,df_ent_10["fitentropy"].values, label=label_fitdata,ls="--",lw=1.5,color="k",
             marker='>',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k",markerfacecolor='None')
    # J_\bot = 2.6
    slope = df_ent_26["slope"].values[0]
    intercept = df_ent_26["intercept"].values[0]
    label_fitdata = r"$J_{\bot}$=%.1f c = %.2f g= %.2f"%(1.2,slope,intercept)
    Lent26, = ax1.plot(df_ent_26["logr"].values,df_ent_26["entropy"].values,label=label_fitdata,ls="-",lw=1.5,color="olive",
             marker='h',alpha=1,markersize=8,markeredgewidth=1.0, markeredgecolor="olive", markerfacecolor='None')
    Lent26_fit, = ax1.plot(df_ent_26["logr"].values,df_ent_26["fitentropy"].values, label=label_fitdata,ls="--",lw=1.5,color="k",
             marker='h',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k",markerfacecolor='None')
    #
    ####图例设置
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 14, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[Lent001,Lent01,Lent04,Lent08,Lent10,Lent26,], loc = 4, bbox_to_anchor=(0.75, 0.1),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    legend2=plt.legend(handles=[LSL], loc = 4, bbox_to_anchor=(0.95, 0.48),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    plt.gca().add_artist(legend1)
    label_x = r"(1/6)log((Lx/$\pi$)sin(x$\pi$/Lx))"
    label_y = "S(x)"
    #plt.yscale("log")
    #plt.xscale("log")
    #ax1.set_xscale("log",base=10,subs=[0.01,0.02,0.03]) 
    #ax1.set_yscale("log",base=10,subs=[0.01,0.02,0.03])      
    ax1.set_xlabel(label_x, size= 16)
    ax1.set_ylabel(label_y, size= 16)
    ax1.tick_params(labelsize = 14) # 设置坐标刻度对应数字的大小
    #ax1.set_xlim([0,40])
    #ax1.set_xticks([5,10,15,20,30,])
    #ax1.set_yticks([-1,-0.5,0,0.5,1]) 
    #
    #=========================================================
    ax1.text(0.2,4.5,"(a)",fontsize = 20, color='black', rotation = 0)
    ax1.text(0.25,3.5, r'$\mathrm{\delta} = 0.375$', fontsize = 16, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
    #ax1.text(2,0.3e-6, r'$J_{\bot}=2.0$', fontsize = 20, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
    #=========================================================
    # 坐标轴设置第一层
    labels = ax1.get_xticklabels() + ax1.get_yticklabels()
    #[label.set_fontname('Times New Roman') for label in labels]###设置ticket labled的字体格式
    ax1.xaxis.set_minor_locator(MultipleLocator(5))###设置次刻度的间隔
    #ax1.yaxis.set_minor_locator(MultipleLocator(10))###设置次刻度的间隔
    ax1.xaxis.set_major_formatter(FormatStrFormatter('%1.2f'))###设置X轴标签文本格式
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
    #-----------------------------------------------------------------------------------------------
    # 选择子图 ax2 进行绘图
    plt.sca(ax2) ## 选择对 ax2 进行绘图
    ax2 = plt.gca()
    df_ent_01_6000 = get_entropy_Lz2(Lz=2,Ly=2,Lx=48,dop=72,t=3,J=1,Jz= 0.1,dim=6000)
    df_ent_01_8000 = get_entropy_Lz2(Lz=2,Ly=2,Lx=48,dop=72,t=3,J=1,Jz= 0.1,dim=8000)
    df_ent_01_10000 = get_entropy_Lz2(Lz=2,Ly=2,Lx=48,dop=72,t=3,J=1,Jz= 0.1,dim=10000)
    df_ent_01_12000 = get_entropy_Lz2(Lz=2,Ly=2,Lx=48,dop=72,t=3,J=1,Jz= 0.1,dim=12000)
    #
    slope = df_SL["slope"].values[0]
    intercept = df_SL["intercept"].values[0]
    label_fitdata = "$2\cdot S_{SL}$, c = %.2f g=%.2f"%(slope * 2, 2*intercept)
    LSL_2, = ax2.plot(df_SL["logr"].values,2*df_SL["entropy"].values,label=label_fitdata,ls="-",lw=1.5,color="red",
             marker='o',alpha=1,markersize=8,markeredgewidth=1.0, markeredgecolor="red",markerfacecolor='k')
    LSL_fit_2, = ax2.plot(df_SL["logr"].values,2*df_SL["fitentropy"].values,label=label_fitdata,ls="--",lw=1.5,color="red",
             marker='o',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="red",
             markerfacecolor='w')
    #
    # J_\bot = 0.1  dim = 6000
    slope = df_ent_01_6000["slope"].values[0]
    intercept = df_ent_01_6000["intercept"].values[0]
    label_fitdata = r"Dim=%d c = %.2f g= %.2f"%(6000,slope,intercept)
    Lent01_6000, = ax2.plot(df_ent_01_6000["logr"].values,df_ent_01_6000["entropy"].values,label=label_fitdata,ls="-",lw=1.5,color="blue",
             marker='s',alpha=1,markersize=8,markeredgewidth=1.0, markeredgecolor="blue", markerfacecolor='None')
    Lent01_fit_6000, = ax2.plot(df_ent_01_6000["logr"].values[2:-1],df_ent_01_6000["fitentropy"].values[2:-1], label=label_fitdata,ls="--",lw=1.5,color="k",
             marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k",markerfacecolor='None')
    # J_\bot = 0.1  dim = 8000
    slope = df_ent_01_8000["slope"].values[0]
    intercept = df_ent_01_8000["intercept"].values[0]
    label_fitdata = r"Dim=%d c = %.2f g= %.2f"%(8000,slope,intercept)
    Lent01_8000, = ax2.plot(df_ent_01_8000["logr"].values,df_ent_01_8000["entropy"].values,label=label_fitdata,ls="-",lw=1.5,color="green",
             marker='^',alpha=1,markersize=8,markeredgewidth=1.0, markeredgecolor="green", markerfacecolor='None')
    Lent01_fit_8000, = ax2.plot(df_ent_01_8000["logr"].values[2:-1],df_ent_01_8000["fitentropy"].values[2:-1], label=label_fitdata,ls="--",lw=1.5,color="k",
             marker='^',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k",markerfacecolor='None')
    # J_\bot = 0.1  dim = 10000
    slope = df_ent_01_10000["slope"].values[0]
    intercept = df_ent_01_10000["intercept"].values[0]
    label_fitdata = r"Dim=%d c = %.2f g= %.2f"%(10000,slope,intercept)
    Lent01_10000, = ax2.plot(df_ent_01_10000["logr"].values,df_ent_01_10000["entropy"].values,label=label_fitdata,ls="-",lw=1.5,color="cyan",
             marker='s',alpha=1,markersize=8,markeredgewidth=1.0, markeredgecolor="cyan", markerfacecolor='None')
    Lent01_fit_10000, = ax2.plot(df_ent_01_10000["logr"].values[2:-1],df_ent_01_10000["fitentropy"].values[2:-1], label=label_fitdata,ls="--",lw=1.5,color="k",
             marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k",markerfacecolor='None')
    # J_\bot = 0.1  dim = 12000
    slope = df_ent_01_12000["slope"].values[0]
    intercept = df_ent_01_12000["intercept"].values[0]
    label_fitdata = r"Dim=%d c = %.2f g= %.2f"%(12000,slope,intercept)
    Lent01_12000, = ax2.plot(df_ent_01_12000["logr"].values,df_ent_01_12000["entropy"].values,label=label_fitdata,ls="-",lw=1.5,color="brown",
             marker='s',alpha=1,markersize=8,markeredgewidth=1.0, markeredgecolor="brown", markerfacecolor='None')
    Lent01_fit_12000, = ax2.plot(df_ent_01_12000["logr"].values[2:-1],df_ent_01_12000["fitentropy"].values[2:-1], label=label_fitdata,ls="--",lw=1.5,color="k",
             marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k",markerfacecolor='None')
    #
    ####图例设置
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 14, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[Lent01_6000,Lent01_8000,Lent01_10000,Lent01_12000,LSL_2], loc = 4, bbox_to_anchor=(0.95, 0.05),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    #legend2=plt.legend(handles=[LSL_2], loc = 4, bbox_to_anchor=(0.9, 0.09),
    #                   ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    #plt.gca().add_artist(legend1)
    label_x = r"(1/6)log((Lx/$\pi$)sin(x$\pi$/Lx))"
    label_y = "S(x)"
    #plt.yscale("log")
    #plt.xscale("log")
    #ax2.set_xscale("log",base=10,subs=[0.01,0.02,0.03]) 
    #ax2.set_yscale("log",base=10,subs=[0.01,0.02,0.03])      
    ax2.set_xlabel(label_x, size= 16)
    ax2.set_ylabel(label_y, size= 16)
    ax2.tick_params(labelsize = 14) # 设置坐标刻度对应数字的大小
    #ax2.set_xlim([0,40])
    #ax2.set_xticks([5,10,15,20,30,])
    #ax2.set_yticks([-1,-0.5,0,0.5,1]) 
    #
    #=========================================================
    ax2.text(0.2,4.75,"(b)",fontsize = 20, color='black', rotation = 0)
    ax2.text(0.2,4.6, r'$\mathrm{\delta} = 0.375$', fontsize = 16, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
    ax2.text(0.2,4.5, r'$J_{\bot} = 0.1$', fontsize = 16, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
    #ax2.text(2,0.3e-6, r'$J_{\bot}=2.0$', fontsize = 20, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
    #=========================================================
    # 坐标轴设置第一层
    labels = ax2.get_xticklabels() + ax2.get_yticklabels()
    #[label.set_fontname('Times New Roman') for label in labels]###设置ticket labled的字体格式
    ax2.xaxis.set_minor_locator(MultipleLocator(5))###设置次刻度的间隔
    #ax2.yaxis.set_minor_locator(MultipleLocator(10))###设置次刻度的间隔
    ax2.xaxis.set_major_formatter(FormatStrFormatter('%1.2f'))###设置X轴标签文本格式
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
    #----------------------------------------------------------------------------------------------------------
    # 选择子图 ax3 进行绘图
    plt.sca(ax3) ## 选择对 ax3 进行绘图
    ax3 = plt.gca()
    # load the single layer entropy data
    df18_SL = get_data_Lz1(Lz=1,Ly=2,Lx=48,dop=18,t=3,J=1,Jz= 0,dim=6000)
    #
    df36_ent_01_6000 = get_entropy_Lz2(Lz=2,Ly=2,Lx=48,dop=36,t=3,J=1,Jz= 0.1,dim=6000)
    df36_ent_01_8000 = get_entropy_Lz2(Lz=2,Ly=2,Lx=48,dop=36,t=3,J=1,Jz= 0.1,dim=8000)
    df36_ent_01_10000 = get_entropy_Lz2(Lz=2,Ly=2,Lx=48,dop=36,t=3,J=1,Jz= 0.1,dim=10000)
    #df36_ent_01_12000 = get_entropy_Lz2(Lz=2,Ly=2,Lx=48,dop=36,t=3,J=1,Jz= 0.1,dim=12000)
    #
    slope = df18_SL["slope"].values[0]
    intercept = df18_SL["intercept"].values[0]
    label_fitdata = "$2\cdot S_{SL}$, c = %.2f g=%.2f"%(slope * 2, 2*intercept)
    L18SL_2, = ax3.plot(df18_SL["logr"].values,2*df18_SL["entropy"].values,label=label_fitdata,ls="-",lw=1.5,color="red",
             marker='o',alpha=1,markersize=8,markeredgewidth=1.0, markeredgecolor="red",markerfacecolor='k')
    L18SL_fit_2, = ax3.plot(df18_SL["logr"].values,2*df18_SL["fitentropy"].values,label=label_fitdata,ls="--",lw=1.5,color="red",
             marker='o',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="red",
             markerfacecolor='w')
    #
    # J_\bot = 0.1  dim = 6000
    slope = df36_ent_01_6000["slope"].values[0]
    intercept = df36_ent_01_6000["intercept"].values[0]
    label_fitdata = r"Dim=%d c = %.2f g= %.2f"%(6000,slope,intercept)
    L36ent01_6000, = ax3.plot(df36_ent_01_6000["logr"].values,df36_ent_01_6000["entropy"].values,label=label_fitdata,ls="-",lw=1.5,color="blue",
             marker='s',alpha=1,markersize=8,markeredgewidth=1.0, markeredgecolor="blue", markerfacecolor='None')
    L36ent01_fit_6000, = ax3.plot(df36_ent_01_6000["logr"].values[2:-1],df36_ent_01_6000["fitentropy"].values[2:-1], label=label_fitdata,ls="--",lw=1.5,color="k",
             marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k",markerfacecolor='None')
    # J_\bot = 0.1  dim = 8000
    slope = df36_ent_01_8000["slope"].values[0]
    intercept = df36_ent_01_8000["intercept"].values[0]
    label_fitdata = r"Dim=%d c = %.2f g= %.2f"%(8000,slope,intercept)
    L36ent01_8000, = ax3.plot(df36_ent_01_8000["logr"].values,df36_ent_01_8000["entropy"].values,label=label_fitdata,ls="-",lw=1.5,color="green",
             marker='^',alpha=1,markersize=8,markeredgewidth=1.0, markeredgecolor="green", markerfacecolor='None')
    L36ent01_fit_8000, = ax3.plot(df36_ent_01_8000["logr"].values[2:-1],df36_ent_01_8000["fitentropy"].values[2:-1], label=label_fitdata,ls="--",lw=1.5,color="k",
             marker='^',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k",markerfacecolor='None')
    # J_\bot = 0.1  dim = 10000
    slope = df36_ent_01_10000["slope"].values[0]
    intercept = df36_ent_01_10000["intercept"].values[0]
    label_fitdata = r"Dim=%d c = %.2f g= %.2f"%(10000,slope,intercept)
    L36ent01_10000, = ax3.plot(df36_ent_01_10000["logr"].values,df36_ent_01_10000["entropy"].values,label=label_fitdata,ls="-",lw=1.5,color="cyan",
             marker='s',alpha=1,markersize=8,markeredgewidth=1.0, markeredgecolor="cyan", markerfacecolor='None')
    L36ent01_fit_10000, = ax3.plot(df36_ent_01_10000["logr"].values[2:-1],df36_ent_01_10000["fitentropy"].values[2:-1], label=label_fitdata,ls="--",lw=1.5,color="k",
             marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k",markerfacecolor='None')
    #
    ####图例设置
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 14, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[L36ent01_6000,L36ent01_8000,L36ent01_10000,L18SL_2], loc = 4, bbox_to_anchor=(0.8, 0.74),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    #legend2=plt.legend(handles=[L18SL_2], loc = 4, bbox_to_anchor=(0.9, 0.09),
    #                   ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    #plt.gca().add_artist(legend1)
    label_x = r"(1/6)log((Lx/$\pi$)sin(x$\pi$/Lx))"
    label_y = "S(x)"
    #plt.yscale("log")
    #plt.xscale("log")
    #ax3.set_xscale("log",base=10,subs=[0.01,0.02,0.03]) 
    #ax3.set_yscale("log",base=10,subs=[0.01,0.02,0.03])      
    ax3.set_xlabel(label_x, size= 16)
    ax3.set_ylabel(label_y, size= 16)
    ax3.tick_params(labelsize = 14) # 设置坐标刻度对应数字的大小
    #ax3.set_xlim([0,40])
    #ax3.set_xticks([5,10,15,20,30,])
    #ax3.set_yticks([-1,-0.5,0,0.5,1]) 
    #
    #=========================================================
    ax3.text(0.42,4.3,"(c)",fontsize = 20, color='black', rotation = 0)
    ax3.text(0.2,5.0, r'$\mathrm{\delta} = 0.1875$', fontsize = 16, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
    ax3.text(0.2,4.9, r'$J_{\bot} = 0.1$', fontsize = 16, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
    #ax3.text(2,0.3e-6, r'$J_{\bot}=2.0$', fontsize = 20, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
    #=========================================================
    # 坐标轴设置第一层
    labels = ax3.get_xticklabels() + ax3.get_yticklabels()
    #[label.set_fontname('Times New Roman') for label in labels]###设置ticket labled的字体格式
    ax3.xaxis.set_minor_locator(MultipleLocator(5))###设置次刻度的间隔
    #ax3.yaxis.set_minor_locator(MultipleLocator(10))###设置次刻度的间隔
    ax3.xaxis.set_major_formatter(FormatStrFormatter('%1.2f'))###设置X轴标签文本格式
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
    #
    fig.tight_layout() # 自动调整 subplot 间的间隙参数
    plt.savefig("E:\\WORK\\Work\\Project\\La3Ni2O7\\datashow_datpy_paper\\figs\\fig_0375_3.eps",
                dpi=300, format='eps',bbox_inches='tight') # 白边紧凑型
    #
    plt.show()