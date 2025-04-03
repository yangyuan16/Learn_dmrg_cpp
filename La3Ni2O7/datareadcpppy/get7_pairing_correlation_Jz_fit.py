# 对比不同 Jz 值下的density correlation function
import numpy as np
import pandas as pd
import matplotlib.pylab as plt
#
def get_data(Lz, Ly, Lx, dop, t, J, Jz, dim, bonds):
    workpath = "E:\\WORK\\Work\\Project\\La3Ni2O7"
    filepath1 = "\\data_dmrgcpp\\Lz%d_Ly%d_Lx%d\\dop%g" % (Lz, Ly, Lx,dop)
    filepath2 = "\\t%d_J%d_Jz%g_dim%d" % (t, J, Jz, dim)
    filepath3 = "\\measurement_pairing_"+bonds+"_ref74_fit.parquet"
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
    bonds = "zz" # zz or yy 方向的pairing
    #load data
    df01 = get_data(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=0.1,dim=dim,bonds=bonds)
    df025 = get_data(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=0.25,dim=dim,bonds=bonds)
    df05 = get_data(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=0.5,dim=dim,bonds=bonds)
    df075 = get_data(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=0.75,dim=dim,bonds=bonds)
    df1 = get_data(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=1,dim=dim,bonds=bonds)
    df2 = get_data(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=2,dim=dim,bonds=bonds)
    df4 = get_data(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=4,dim=dim,bonds=bonds)
    df6 = get_data(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=6,dim=dim,bonds=bonds)
    df8 = get_data(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=8,dim=dim,bonds=bonds)
    #
    #-----------------plot logr-r fig and logr-logr fig-------------------------
    fig = plt.figure(figsize=(6,10)) 
    ax1 = plt.subplot(1,1,1)
    plt.sca(ax1)  ##选择对ax2进行绘图
    ax1=plt.gca() #获得坐标轴的句柄
    #========= plot curve
    #L01, = ax1.plot(df01["r"],df01["corre_abs"],label="Jz={}".format(0.1),ls="-",lw=1.5,color="r",
    #                marker='o',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="k",markerfacecolor='w')
    #L025, = ax1.plot(df025["r"],df025["corre_abs"],label="Jz={}".format(0.25),ls="-",lw=1.5,color="blue",
    #                 marker='s',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="k",markerfacecolor='w')
    #L05, = ax1.plot(df05["r"],df05["corre_abs"],label="Jz={}".format(0.5),ls="-",lw=1.5,color="green",
    #                marker='^',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="k",markerfacecolor='w')
    #L075, = ax1.plot(df075["r"],df075["corre_abs"],label="Jz={}".format(0.75),ls="-",lw=1.5,color="orange",
    #                 marker='v',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="k",markerfacecolor='w')
    #L1, = ax1.plot(df1["r"],df1["corre_abs"],label="Jz={}".format(1),ls="-",lw=1.5,color="magenta",
    #               marker='v',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="k",markerfacecolor='w')
    
    #
    L2, = ax1.plot(df2["r"],df2["corre_abs"],label="Jz={}".format(2),ls="-",lw=1.5,color="cyan",
                   marker='^',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="k",markerfacecolor='w')
    L4, = ax1.plot(df4["r"],df4["corre_abs"],label="Jz={}".format(4),ls="-",lw=1.5,color="yellowgreen",
                   marker='*',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="k",markerfacecolor='w')
    L6, = ax1.plot(df6["r"],df6["corre_abs"],label="Jz={}".format(6),ls="-",lw=1.5,color="lime",
                   marker='h',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="k",markerfacecolor='w')
    L8, = ax1.plot(df8["r"],df8["corre_abs"],label="Jz={}".format(8),ls="-",lw=1.5,color="brown",
             marker='x',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="k",markerfacecolor='w')
    
    #============ plot logr-r curve
    #L01_exp, = ax1.plot(df01["r"],df01["fitcorre_exp"],label=r"$\xi$={:.4f}".format(-1/df01["slope_exp"].values[0]),ls="--",lw=1.5,color="r",
    #                    marker='o',alpha=1,markersize=0,markeredgewidth=1.5, markeredgecolor="k",markerfacecolor='w')
    #L025_exp, = ax1.plot(df025["r"],df025["fitcorre_exp"],label=r"$\xi$={:.4f}".format(-1/df025["slope_exp"].values[0]),ls="--",lw=1.5,color="blue",
    #                     marker='o',alpha=1,markersize=0,markeredgewidth=1.5, markeredgecolor="k",markerfacecolor='w')
    #L05_exp, = ax1.plot(df05["r"],df05["fitcorre_exp"],label=r"$\xi$={:.4f}".format(-1/df05["slope_exp"].values[0]),ls="--",lw=1.5,color="green",
    #                    marker='o',alpha=1,markersize=0,markeredgewidth=1.5, markeredgecolor="k",markerfacecolor='w')
    #L075_exp, = ax1.plot(df075["r"],df075["fitcorre_exp"],label=r"$\xi$={:.4f}".format(-1/df075["slope_exp"].values[0]),ls="--",lw=1.5,color="orange",
    #                     marker='o',alpha=1,markersize=0,markeredgewidth=1.5, markeredgecolor="k",markerfacecolor='w')
    #L1_exp, = ax1.plot(df1["r"],df1["fitcorre_exp"],label=r"$\xi$={:.4f}".format(-1/df1["slope_exp"].values[0]),ls="--",lw=1.5,color="magenta",
    #                   marker='o',alpha=1,markersize=0,markeredgewidth=1.5, markeredgecolor="k",markerfacecolor='w')
    
    #L2_exp, = ax1.plot(df2["r"],df2["fitcorre_exp"],label=r"$\xi$={:.4f}".format(-1/df2["slope_exp"].values[0]),ls="--",lw=1.5,color="cyan",
    #                   marker='o',alpha=1,markersize=0,markeredgewidth=1.5, markeredgecolor="k",markerfacecolor='w')
    #L4_exp, = ax1.plot(df4["r"],df4["fitcorre_exp"],label=r"$\xi$={:.4f}".format(-1/df4["slope_exp"].values[0]),ls="--",lw=1.5,color="yellowgreen",
    #                   marker='o',alpha=1,markersize=0,markeredgewidth=1.5, markeredgecolor="k", markerfacecolor='w')
    #L6_exp, = ax1.plot(df6["r"],df6["fitcorre_exp"],label=r"$\xi$={:.4f}".format(-1/df6["slope_exp"].values[0]),ls="--",lw=1.5,color="lime",
    #                   marker='o',alpha=1,markersize=0,markeredgewidth=1.5, markeredgecolor="k", markerfacecolor='w')
    #L8_exp, = ax1.plot(df8["r"],df8["fitcorre_exp"],label=r"$\xi$={:.4f}".format(-1/df8["slope_exp"].values[0]),ls="--",lw=1.5,color="brown",
    #                   marker='o',alpha=1,markersize=0,markeredgewidth=1.5, markeredgecolor="k", markerfacecolor='w')
    #============= plot logr-logr curve
    #
    #L01_pow, = ax1.plot(df01["r"],df01["fitcorre_pow"],label="$K$={:.4f}".format(-df01["slope_pow"].values[0]),ls="--",lw=1.5,color="r",
    #                    marker='o',alpha=1,markersize=0,markeredgewidth=1.5, markeredgecolor="k",markerfacecolor='w')
    #L025_pow, = ax1.plot(df025["r"],df025["fitcorre_pow"],label="$K$={:.4f}".format(-df025["slope_pow"].values[0]),ls="--",lw=1.5,color="blue",
    #                     marker='o',alpha=1,markersize=0,markeredgewidth=1.5, markeredgecolor="k",markerfacecolor='w')
    #L05_pow, = ax1.plot(df05["r"],df05["fitcorre_pow"],label="$K$={:.4f}".format(-df05["slope_pow"].values[0]),ls="--",lw=1.5,color="green",
    #                    marker='o',alpha=1,markersize=0,markeredgewidth=1.5, markeredgecolor="k",markerfacecolor='w')
    #L075_pow, = ax1.plot(df075["r"],df075["fitcorre_pow"],label="$K$={:.4f}".format(-df075["slope_pow"].values[0]),ls="--",lw=1.5,color="orange",
    #                     marker='o',alpha=1,markersize=0,markeredgewidth=1.5, markeredgecolor="k",markerfacecolor='w')
    #L1_pow, = ax1.plot(df1["r"],df1["fitcorre_pow"],label="$K$={:.4f}".format(-df1["slope_pow"].values[0]),ls="--",lw=1.5,color="magenta",
    #                   marker='o',alpha=1,markersize=0,markeredgewidth=1.5, markeredgecolor="k",markerfacecolor='w')    
    
    #-------------------
    L2_pow, = ax1.plot(df2["r"],df2["fitcorre_pow"],label=r"$K$={:.4f}".format(-df2["slope_pow"].values[0]),ls="--",lw=1.5,color="cyan",
                       marker='o',alpha=1,markersize=0,markeredgewidth=1.5, markeredgecolor="k", markerfacecolor='w')
    L4_pow, = ax1.plot(df4["r"],df4["fitcorre_pow"],label=r"$K$={:.4f}".format(-df4["slope_pow"].values[0]),ls="--",lw=1.5,color="yellowgreen",
                       marker='o',alpha=1,markersize=0,markeredgewidth=1.5, markeredgecolor="k", markerfacecolor='w')
    L6_pow, = ax1.plot(df6["r"],df6["fitcorre_pow"],label=r"$K$={:.4f}".format(-df6["slope_pow"].values[0]),ls="--",lw=1.5,color="lime",
                       marker='o',alpha=1,markersize=0,markeredgewidth=1.5, markeredgecolor="k", markerfacecolor='w')
    L8_pow, = ax1.plot(df8["r"],df8["fitcorre_pow"],label=r"$K$={:.4f}".format(-df8["slope_pow"].values[0]),ls="--",lw=1.5,color="brown",
                       marker='o',alpha=1,markersize=0,markeredgewidth=1.5, markeredgecolor="k", markerfacecolor='w')
    ####图例设置
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 16, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[L2,L2_pow,L4,L4_pow,L6,L6_pow,L8,L8_pow], loc = 4, bbox_to_anchor=(0.38, 0.02),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    # [L01,L01_exp,L025,L025_exp,L05,L05_exp,L075,L075_exp,L1,L1_exp]
    # [L2,L2_exp,L4,L4_exp,L6,L6_exp,L8,L8_exp]
    # [L01,L01_pow,L025,L025_pow,L05,L05_pow,L075,L075_pow,L1,L1_pow]
    # [L2,L2_pow,L4,L4_pow,L6,L6_pow,L8,L8_pow]
    # [L05,L05_exp,L1,L1_exp,L2,L2_exp]
    # [L05,L05_pow,L1,L1_pow,L2,L2_pow]

    label_x = r"|i-j|"
    label_y = bonds + " pairng corre."
    plt.yscale("log")
    plt.xscale("log")      
    ax1.set_xlabel(label_x, size= 14)
    ax1.set_ylabel(label_y, size= 14)
    ax1.tick_params(labelsize = 15) # 设置坐标刻度对应数字的大小
    #ax1.set_xlim([0,8])
    #ax1.set_ylim([-0.1,1])
    #ax1.set_xticks([0,2,4,6,8])
    #ax1.set_yticks([-0.1,0,0.5,1])
    #ax1.text(0.3,-0.05, r'$\mathrm{(a)}$', fontsize=18)
    plt.show()
    #