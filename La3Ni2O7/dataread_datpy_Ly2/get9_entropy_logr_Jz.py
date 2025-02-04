# 对比不同 Jz 值下的自旋关联
import numpy as np
import pandas as pd
import matplotlib.pylab as plt
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

def get_data_Lz2(Lz, Ly, Lx, dop, t, J, Jz, dim):
    workpath = "E:\\WORK\\Work\\Project\\La3Ni2O7"
    filepath1 = "\\data_dmrgcpp\\Lz%d_Ly%d_Lx%d\\dop%g" % (Lz, Ly, Lx,dop)
    filepath2 = "\\t%d_J%d_Jz%.2f_dim%d_ent\\entanglement" % (t, J, Jz, dim)
    filepath3 = "\\measurement_entropy_fit.parquet"
    filename = workpath + filepath1 + filepath2 + filepath3
    print(filename)
    df = pd.read_parquet(filename)
    print(df.head())
    return df

if __name__ == "__main__":
    print()
    
    #-------- Lz = 1, Jz = 0, dim = 6000
    #df_SL = get_data_Lz1(Lz=1,Ly=2,Lx=48,dop=36,t=3,J=1,Jz= 0,dim=6000)
    #--------- Lz = 2, Jz = [0.01,2.6], dim = 6000
    #df_001 = get_data_Lz2(Lz=2,Ly=2,Lx=48,dop=72,t=3,J=1,Jz= 0.01,dim=6000)
    #df_002 = get_data_Lz2(Lz=2,Ly=2,Lx=48,dop=72,t=3,J=1,Jz= 0.02,dim=6000)
    #df_01 = get_data_Lz2(Lz=2,Ly=2,Lx=48,dop=72,t=3,J=1,Jz= 0.1,dim=6000)
    #df_02 = get_data_Lz2(Lz=2,Ly=2,Lx=48,dop=72,t=3,J=1,Jz= 0.2,dim=6000)
    #df_04 = get_data_Lz2(Lz=2,Ly=2,Lx=48,dop=72,t=3,J=1,Jz= 0.4,dim=6000)
    #df_06 = get_data_Lz2(Lz=2,Ly=2,Lx=48,dop=72,t=3,J=1,Jz= 0.6,dim=6000)
    #df_08 = get_data_Lz2(Lz=2,Ly=2,Lx=48,dop=72,t=3,J=1,Jz= 0.8,dim=6000)
    #df_10 = get_data_Lz2(Lz=2,Ly=2,Lx=48,dop=72,t=3,J=1,Jz= 1.0,dim=6000)
    #df_26 = get_data_Lz2(Lz=2,Ly=2,Lx=48,dop=72,t=3,J=1,Jz= 2.6,dim=6000)
    #--------------------------
    
    df_SL = get_data_Lz1(Lz=1,Ly=2,Lx=48,dop=18,t=3,J=1,Jz= 0,dim=6000)
    #--------- Lz = 2, Jz = [0.01,2.6], dim = 6000
    df_01 = get_data_Lz2(Lz=2,Ly=2,Lx=48,dop=36,t=3,J=1,Jz= 0.1,dim=6000)
    df_02 = get_data_Lz2(Lz=2,Ly=2,Lx=48,dop=36,t=3,J=1,Jz= 0.2,dim=6000)
    df_04 = get_data_Lz2(Lz=2,Ly=2,Lx=48,dop=36,t=3,J=1,Jz= 0.4,dim=6000)
    df_06 = get_data_Lz2(Lz=2,Ly=2,Lx=48,dop=36,t=3,J=1,Jz= 0.6,dim=6000)
    df_08 = get_data_Lz2(Lz=2,Ly=2,Lx=48,dop=36,t=3,J=1,Jz= 0.8,dim=6000)
    df_10 = get_data_Lz2(Lz=2,Ly=2,Lx=48,dop=36,t=3,J=1,Jz= 1.0,dim=6000)
    df_12 = get_data_Lz2(Lz=2,Ly=2,Lx=48,dop=36,t=3,J=1,Jz= 1.2,dim=6000)
    df_14 = get_data_Lz2(Lz=2,Ly=2,Lx=48,dop=36,t=3,J=1,Jz= 1.4,dim=6000)
    #------------------------------
    fig = plt.figure(figsize=(4,12)) 
    ax1 = plt.subplot(1,1,1)
    plt.sca(ax1)  ##选择对ax1进行绘图
    ax1=plt.gca() #获得坐标轴的句柄
    # plot single layer
    '''
    slope = df_SL["slope"].values[0]
    intercept = df_SL["intercept"].values[0]
    label_data = "entropy"
    LSL, = ax1.plot(df_SL["logr"].values,df_SL["entropy"].values,
                   label=label_data,ls="-",lw=1.5,color="red",
             marker='o',alpha=1,markersize=8,markeredgewidth=1.5, markeredgecolor="red",
             markerfacecolor='None')
    label_fitdata = r"SL, Jz ={:.1f} c ={:.3f} g={:.3f}".format(0,slope,intercept)
    LSL_fit, = ax1.plot(df_SL["logr"].values,df_SL["fitentropy"].values,
                   label=label_fitdata,ls="--",lw=1.5,color="red",
             marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="red",
             markerfacecolor='w')
    '''
    
    # 取 single layer entropy 的 doublue
    slope = df_SL["slope"].values[0]
    intercept = df_SL["intercept"].values[0]
    label_data = "2*S (single layer)"
    LSL2, = ax1.plot(df_SL["logr"].values,2 * df_SL["entropy"].values,
                   label=label_data,ls="-",lw=1.5,color="red",
             marker='o',alpha=1,markersize=8,markeredgewidth=1.5, markeredgecolor="red",
             markerfacecolor='None')
    #label_fitdata = r"SL, Jz ={:.1f} c ={:.3f} g={:.3f}".format(0,slope,intercept)
    #LSL_fit, = ax1.plot(df_SL["logr"].values,df_SL["fitentropy"].values,
    #               label=label_fitdata,ls="--",lw=1.5,color="red",
    #         marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="red",
    #         markerfacecolor='w')
    
    '''
    # plot bilayer Jz = 0.01
    slope = df_001["slope"].values[0]
    intercept = df_001["intercept"].values[0]
    label_data = "entropy"
    L001, = ax1.plot(df_001["logr"].values,df_001["entropy"].values,
                   label=label_data,ls="-",lw=1.5,color="blue",
             marker='o',alpha=1,markersize=8,markeredgewidth=1.5, markeredgecolor="blue",
             markerfacecolor='None')
    label_fitdata = r"BL, Jz ={:.2f} c ={:.3f} g={:.3f}".format(0.01,slope,intercept)
    L001_fit, = ax1.plot(df_001["logr"].values,df_001["fitentropy"].values,
                   label=label_fitdata,ls="--",lw=1.5,color="blue",
             marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="blue",
             markerfacecolor='w')
    # plot bilayer Jz = 0.02
    
    slope = df_002["slope"].values[0]
    intercept = df_002["intercept"].values[0]
    label_data = "entropy"
    L002, = ax1.plot(df_002["logr"].values,df_002["entropy"].values,
                   label=label_data,ls="-",lw=1.5,color="green",
             marker='o',alpha=1,markersize=8,markeredgewidth=1.5, markeredgecolor="green",
             markerfacecolor='None')
    label_fitdata = r"BL, Jz ={:.2f} c ={:.3f} g={:.3f}".format(0.02,slope,intercept)
    L002_fit, = ax1.plot(df_002["logr"].values,df_002["fitentropy"].values,
                   label=label_fitdata,ls="--",lw=1.5,color="green",
             marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="green",
             markerfacecolor='w')
    '''
    # plot bilayer Jz = 0.1
    slope = df_01["slope"].values[0]
    intercept = df_01["intercept"].values[0]
    label_data = "entropy"
    L01, = ax1.plot(df_01["logr"].values,df_01["entropy"].values,
                   label=label_data,ls="-",lw=1.5,color="magenta",
             marker='o',alpha=1,markersize=8,markeredgewidth=1.5, markeredgecolor="magenta",
             markerfacecolor='None')
    label_fitdata = r"BL, Jz ={:.2f} c ={:.3f} g={:.3f}".format(0.1,slope,intercept)
    L01_fit, = ax1.plot(df_01["logr"].values,df_01["fitentropy"].values,
                   label=label_fitdata,ls="--",lw=1.5,color="magenta",
             marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="magenta",
             markerfacecolor='w')
    '''
    # plot bilayer Jz = 0.2
    slope = df_02["slope"].values[0]
    intercept = df_02["intercept"].values[0]
    label_data = "entropy"
    L02, = ax1.plot(df_02["logr"].values,df_02["entropy"].values,
                   label=label_data,ls="-",lw=1.5,color="cyan",
             marker='o',alpha=1,markersize=8,markeredgewidth=1.5, markeredgecolor="cyan",
             markerfacecolor='None')
    label_fitdata = r"BL, Jz ={:.2f} c ={:.3f} g={:.3f}".format(0.2,slope,intercept)
    L02_fit, = ax1.plot(df_02["logr"].values,df_02["fitentropy"].values,
                   label=label_fitdata,ls="--",lw=1.5,color="cyan",
             marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="cyan",
             markerfacecolor='w')
    # plot bilayer Jz = 0.4
    slope = df_04["slope"].values[0]
    intercept = df_04["intercept"].values[0]
    label_data = "entropy"
    L04, = ax1.plot(df_04["logr"].values,df_04["entropy"].values,
                   label=label_data,ls="-",lw=1.5,color="brown",
             marker='o',alpha=1,markersize=8,markeredgewidth=1.5, markeredgecolor="brown",
             markerfacecolor='None')
    label_fitdata = r"BL, Jz ={:.2f} c ={:.3f} g={:.3f}".format(0.4,slope,intercept)
    L04_fit, = ax1.plot(df_04["logr"].values,df_04["fitentropy"].values,
                   label=label_fitdata,ls="--",lw=1.5,color="brown",
             marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="brown",
             markerfacecolor='w')
    # plot bilayer Jz = 0.6
    slope = df_06["slope"].values[0]
    intercept = df_06["intercept"].values[0]
    label_data = "entropy"
    L06, = ax1.plot(df_06["logr"].values,df_06["entropy"].values,
                   label=label_data,ls="-",lw=1.5,color="violet",
             marker='o',alpha=1,markersize=8,markeredgewidth=1.5, markeredgecolor="violet",
             markerfacecolor='None')
    label_fitdata = r"BL, Jz ={:.2f} c ={:.3f} g={:.3f}".format(0.6,slope,intercept)
    L06_fit, = ax1.plot(df_06["logr"].values,df_06["fitentropy"].values,
                   label=label_fitdata,ls="--",lw=1.5,color="violet",
             marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="violet",
             markerfacecolor='w')
    # plot bilayer Jz = 0.8
    slope = df_08["slope"].values[0]
    intercept = df_08["intercept"].values[0]
    label_data = "entropy"
    L08, = ax1.plot(df_08["logr"].values,df_08["entropy"].values,
                   label=label_data,ls="-",lw=1.5,color="olive",
             marker='o',alpha=1,markersize=8,markeredgewidth=1.5, markeredgecolor="olive",
             markerfacecolor='None')
    label_fitdata = r"BL, Jz ={:.2f} c ={:.3f} g={:.3f}".format(0.8,slope,intercept)
    L08_fit, = ax1.plot(df_08["logr"].values,df_08["fitentropy"].values,
                   label=label_fitdata,ls="--",lw=1.5,color="olive",
             marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="olive",
             markerfacecolor='w')
    # plot bilayer Jz = 1.0
    slope = df_10["slope"].values[0]
    intercept = df_10["intercept"].values[0]
    label_data = "entropy"
    L10, = ax1.plot(df_10["logr"].values,df_10["entropy"].values,
                   label=label_data,ls="-",lw=1.5,color="steelblue",
             marker='o',alpha=1,markersize=8,markeredgewidth=1.5, markeredgecolor="steelblue",
             markerfacecolor='None')
    label_fitdata = r"BL, Jz ={:.2f} c ={:.3f} g={:.3f}".format(1.0,slope,intercept)
    L10_fit, = ax1.plot(df_10["logr"].values,df_10["fitentropy"].values,
                   label=label_fitdata,ls="--",lw=1.5,color="steelblue",
             marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="steelblue",
             markerfacecolor='w')
    # plot bilayer Jz = 1.2
    slope = df_12["slope"].values[0]
    intercept = df_12["intercept"].values[0]
    label_data = "entropy"
    L12, = ax1.plot(df_12["logr"].values,df_12["entropy"].values,
                   label=label_data,ls="-",lw=1.5,color="blue",
             marker='o',alpha=1,markersize=8,markeredgewidth=1.5, markeredgecolor="blue",
             markerfacecolor='None')
    label_fitdata = r"BL, Jz ={:.2f} c ={:.3f} g={:.3f}".format(1.2,slope,intercept)
    L12_fit, = ax1.plot(df_12["logr"].values,df_12["fitentropy"].values,
                   label=label_fitdata,ls="--",lw=1.5,color="blue",
             marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="blue",
             markerfacecolor='w')
    # plot bilayer Jz = 1.4
    slope = df_14["slope"].values[0]
    intercept = df_14["intercept"].values[0]
    label_data = "entropy"
    L14, = ax1.plot(df_14["logr"].values,df_14["entropy"].values,
                   label=label_data,ls="-",lw=1.5,color="green",
             marker='o',alpha=1,markersize=8,markeredgewidth=1.5, markeredgecolor="green",
             markerfacecolor='None')
    label_fitdata = r"BL, Jz ={:.2f} c ={:.3f} g={:.3f}".format(1.4,slope,intercept)
    L14_fit, = ax1.plot(df_14["logr"].values,df_14["fitentropy"].values,
                   label=label_fitdata,ls="--",lw=1.5,color="green",
             marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="green",
             markerfacecolor='w')
    # plot bilayer Jz = 2.6
    
    slope = df_26["slope"].values[0]
    intercept = df_26["intercept"].values[0]
    label_data = "entropy"
    L26, = ax1.plot(df_26["logr"].values,df_26["entropy"].values,
                   label=label_data,ls="-",lw=1.5,color="k",
             marker='o',alpha=1,markersize=8,markeredgewidth=1.5, markeredgecolor="k",
             markerfacecolor='None')
    label_fitdata = r"BL, Jz ={:.2f} c ={:.3f} g={:.3f}".format(2.6,slope,intercept)
    L26_fit, = ax1.plot(df_26["logr"].values,df_26["fitentropy"].values,
                   label=label_fitdata,ls="--",lw=1.5,color="k",
             marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k",
             markerfacecolor='w')
    '''
    #
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 26, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[LSL2, L01_fit,], 
                       loc = 4, bbox_to_anchor=(0.48, 0.0),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    # [L001_fit,L002_fit,L01_fit,L02_fit,L04_fit,L06_fit,L08_fit,L10_fit,]
    label_x = r"(1/6)log((Lx/$\pi$)sin(x$\pi$/Lx))"
    label_y = "Entropy"
    #plt.yscale("log")
    ax1.set_xlabel(label_x, size= 25)
    ax1.set_ylabel(label_y, size= 25)
    ax1.tick_params(labelsize = 25) # 设置坐标刻度对应数字的大小
    #ax1.set_xlim([0,8])
    #ax1.set_ylim([-0.1,1])
    #ax1.set_xticks([0,2,4,6,8])
    #ax1.set_yticks([-0.1,0,0.5,1])
    #ax1.text(0.3,-0.05, r'$\mathrm{(a)}$', fontsize=18)
    plt.title("$\delta$=%.4f"%(36/(4*48)), fontsize=25)
    plt.show()