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
    #
    Lz = 2
    Ly = 2
    Lx = 64
    dop = 48
    df_SL = get_data_Lz1(Lz=1,Ly=2,Lx=Lx,dop=24,t=3,J=1,Jz= 0,dim=6000) # !!! change dop
    #--------- Lz = 2, Jz = [0.01,2.6], dim = 6000
    df_01_6000 = get_data_Lz2(Lz=2,Ly=2,Lx=Lx,dop=dop,t=3,J=1,Jz= 0.1,dim=6000) # !!!
    df_01_8000 = get_data_Lz2(Lz=2,Ly=2,Lx=Lx,dop=dop,t=3,J=1,Jz= 0.1,dim=8000) # !!!
    df_01_10000 = get_data_Lz2(Lz=2,Ly=2,Lx=Lx,dop=dop,t=3,J=1,Jz= 0.1,dim=10000) # !!!
    #df_01_12000 = get_data_Lz2(Lz=2,Ly=2,Lx=48,dop=dop,t=3,J=1,Jz= 0.1,dim=12000) # !!!
    #------------------------------
    fig = plt.figure(figsize=(4,12)) 
    ax1 = plt.subplot(1,1,1)
    plt.sca(ax1)  ##选择对ax1进行绘图
    ax1=plt.gca() #获得坐标轴的句柄
    # plot single layer
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
    # plot bilayer Jz = 0.1  dim = 6000
    slope = df_01_6000["slope"].values[0]
    intercept = df_01_6000["intercept"].values[0]
    label_data = "entropy"
    L01_6000, = ax1.plot(df_01_6000["logr"].values,df_01_6000["entropy"].values,
                   label=label_data,ls="-",lw=1.5,color="magenta",
             marker='o',alpha=1,markersize=8,markeredgewidth=1.5, markeredgecolor="magenta",
             markerfacecolor='None')
    label_fitdata = r"BL, Jz ={:.2f} c ={:.3f} g={:.3f}, dim={:d}".format(0.1,slope,intercept,6000)
    L01_6000_fit, = ax1.plot(df_01_6000["logr"].values,df_01_6000["fitentropy"].values,
                   label=label_fitdata,ls="--",lw=1.5,color="magenta",
             marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="magenta",
             markerfacecolor='w')
    # plot bilayer Jz = 0.1  dim = 8000
    slope = df_01_8000["slope"].values[0]
    intercept = df_01_8000["intercept"].values[0]
    label_data = "entropy"
    L01_8000, = ax1.plot(df_01_8000["logr"].values,df_01_8000["entropy"].values,
                   label=label_data,ls="-",lw=1.5,color="blue",
             marker='o',alpha=1,markersize=8,markeredgewidth=1.5, markeredgecolor="blue",
             markerfacecolor='None')
    label_fitdata = r"BL, Jz ={:.2f} c ={:.3f} g={:.3f}, dim={:d}".format(0.1,slope,intercept,8000)
    L01_8000_fit, = ax1.plot(df_01_8000["logr"].values,df_01_8000["fitentropy"].values,
                   label=label_fitdata,ls="--",lw=1.5,color="blue",
             marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="blue",
             markerfacecolor='w')
    # plot bilayer Jz = 0.1  dim = 10000
    slope = df_01_10000["slope"].values[0]
    intercept = df_01_10000["intercept"].values[0]
    label_data = "entropy"
    L01_10000, = ax1.plot(df_01_10000["logr"].values,df_01_10000["entropy"].values,
                   label=label_data,ls="-",lw=1.5,color="cyan",
             marker='o',alpha=1,markersize=8,markeredgewidth=1.5, markeredgecolor="cyan",
             markerfacecolor='None')
    label_fitdata = r"BL, Jz ={:.2f} c ={:.3f} g={:.3f}, dim={:d}".format(0.1,slope,intercept,10000)
    L01_10000_fit, = ax1.plot(df_01_10000["logr"].values,df_01_10000["fitentropy"].values,
                   label=label_fitdata,ls="--",lw=1.5,color="cyan",
             marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="cyan",
             markerfacecolor='w')
    #
    # plot bilayer Jz = 0.1  dim = 12000
    '''
    slope = df_01_12000["slope"].values[0]
    intercept = df_01_12000["intercept"].values[0]
    label_data = "entropy"
    L01_12000, = ax1.plot(df_01_12000["logr"].values,df_01_12000["entropy"].values,
                   label=label_data,ls="-",lw=1.5,color="green",
             marker='o',alpha=1,markersize=8,markeredgewidth=1.5, markeredgecolor="green",
             markerfacecolor='None')
    label_fitdata = r"BL, Jz ={:.2f} c ={:.3f} g={:.3f}, dim={:d}".format(0.1,slope,intercept,12000)
    L01_12000_fit, = ax1.plot(df_01_12000["logr"].values,df_01_12000["fitentropy"].values,
                   label=label_fitdata,ls="--",lw=1.5,color="green",
             marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="green",
             markerfacecolor='w')
    '''
    #
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 26, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[LSL2, L01_6000_fit, L01_8000_fit, L01_10000_fit], 
                       loc = 4, bbox_to_anchor=(0.78, 0.0),
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
    plt.title("$\delta$=%.4f"%(dop/(Lz * Ly * Lx)), fontsize=25)
    plt.show()