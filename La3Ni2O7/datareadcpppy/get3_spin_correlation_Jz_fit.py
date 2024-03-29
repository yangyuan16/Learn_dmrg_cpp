# 对比不同 Jz 值下的自旋关联, 并同时画出fit曲线
import numpy as np
import pandas as pd
import matplotlib.pylab as plt

def get_data(Lz, Ly, Lx, dop, t, J, Jz, dim):
    # load data
    workpath = "E:\\WORK\\Work\\Project\\La3Ni2O7"
    filepath1 = "\\data_dmrgcpp\\Lz%d_Ly%d_Lx%d\\dop%g" % (Lz, Ly, Lx,dop)
    filepath2 = "\\t%d_J%d_Jz%g_dim%d" % (t, J, Jz, dim)
    filepath3 = "\\measurement_spin_correlation_ref74_fit.parquet"
    filename = workpath + filepath1 + filepath2 + filepath3
    print(filename)
    df = pd.read_parquet(filename)
    print(df.head())
    return df
#
if __name__ == "__main__":
    print()
    Lz = 2
    Ly = 3
    Lx = 48
    dop = 0.5
    t = 3
    J = 1
    dim = 6000 # dim cutoff
    #load data
    df05 = get_data(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=0.1,dim=dim)
    df1 = get_data(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=1,dim=dim)
    df2 = get_data(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=2,dim=dim)
    #
    r_05 = df05["r"].values
    corre_05 = df05["corre_abs"].values
    fitcorre_exp_05 = df05["fitcorre_exp"].values
    slope_exp_05 = df05["slope_exp"].values[0]
    #
    r_1 = df1["r"].values
    corre_1 = df1["corre_abs"].values
    fitcorre_exp_1 = df1["fitcorre_exp"].values
    slope_exp_1 = df1["slope_exp"].values[0]
    #
    r_2 = df2["r"].values
    corre_2 = df2["corre_abs"].values
    fitcorre_exp_2 = df2["fitcorre_exp"].values
    slope_exp_2 = df2["slope_exp"].values[0]
    #-----------------plot logr-r fig and logr-logr fig-------------------------
    fig = plt.figure(figsize=(6,10)) 
    ax1 = plt.subplot(1,1,1)
    plt.sca(ax1)  ##选择对ax2进行绘图
    ax1=plt.gca() #获得坐标轴的句柄
    
    L05, = ax1.plot(r_05,corre_05,label="Jz={}".format(0.5),ls="-",lw=1.5,color="r",
             marker='o',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="k",
             markerfacecolor='w')
    L05_exp, = ax1.plot(r_05,fitcorre_exp_05,label=r"$\xi$={:.4f}".format(-1/slope_exp_05),ls="--",lw=1.5,color="r",
             marker='o',alpha=1,markersize=0,markeredgewidth=1.5, markeredgecolor="k",
             markerfacecolor='w')
    #
    L1, = ax1.plot(r_1,corre_1,label="Jz={}".format(1),ls="-",lw=1.5,color="green",
             marker='s',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="k",
             markerfacecolor='w')
    L1_exp, = ax1.plot(r_1,fitcorre_exp_1,label=r"$\xi$={:.4f}".format(-1/slope_exp_1),ls="--",lw=1.5,color="green",
             marker='o',alpha=1,markersize=0,markeredgewidth=1.5, markeredgecolor="k",
             markerfacecolor='w')
    #
    L2, = ax1.plot(r_2,corre_2,label="Jz={}".format(2),ls="-",lw=1.5,color="blue",
             marker='^',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="k",
             markerfacecolor='w')
    L2_exp, = ax1.plot(r_2,fitcorre_exp_2,label=r"$\xi$={:.4f}".format(-1/slope_exp_2),ls="--",lw=1.5,color="blue",
             marker='o',alpha=1,markersize=0,markeredgewidth=1.5, markeredgecolor="k",
             markerfacecolor='w')
    ####图例设置
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 16, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[L05,L05_exp,L1,L1_exp,L2,L2_exp], loc = 4, bbox_to_anchor=(0.38, 0.18),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    
    label_x = r"|i-j|"
    label_y = "Spin Correlation"
    plt.yscale("log")
    #plt.xscale("log")      
    ax1.set_xlabel(label_x, size= 14)
    ax1.set_ylabel(label_y, size= 14)
    ax1.tick_params(labelsize = 15) # 设置坐标刻度对应数字的大小
    #ax1.set_xlim([0,8])
    #ax1.set_ylim([-0.1,1])
    #ax1.set_xticks([0,2,4,6,8])
    #ax1.set_yticks([-0.1,0,0.5,1])
    #ax1.text(0.3,-0.05, r'$\mathrm{(a)}$', fontsize=18)
    plt.show()
    