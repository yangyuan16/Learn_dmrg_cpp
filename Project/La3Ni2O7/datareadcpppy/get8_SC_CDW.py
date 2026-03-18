# 比较不同 Jz 下 YY ZZ pairing 以及 CDW powerlaw fitting 出来的曲线
import numpy as np
import pandas as pd
import matplotlib.pylab as plt
#
def get_data_pairng(Lz, Ly, Lx, dop, t, J, Jz, dim, bonds):
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
def get_data_density(Lz, Ly, Lx, dop, t, J, Jz, dim):
    # load data
    workpath = "E:\\WORK\\Work\\Project\\La3Ni2O7"
    filepath1 = "\\data_dmrgcpp\\Lz%d_Ly%d_Lx%d\\dop%g" % (Lz, Ly, Lx,dop)
    filepath2 = "\\t%d_J%d_Jz%g_dim%d" % (t, J, Jz, dim)
    filepath3 = "\\measurement_density_correlation_ref74_fit.parquet"
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
    Jz = 2
    dim = 6000 # dim cutoff
    bonds_zz = "zz" # zz or yy 方向的pairing
    bonds_yy = "yy" # yy 方向的pairing
    #
    df_cdw = get_data_density(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=Jz,dim=dim)
    df_yy_sc = get_data_pairng(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=Jz,dim=dim,bonds=bonds_yy)
    df_zz_sc = get_data_pairng(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=Jz,dim=dim,bonds=bonds_zz)
    #
    #-----------------plot logr-r fig and logr-logr fig-------------------------
    fig = plt.figure(figsize=(6,10)) 
    ax1 = plt.subplot(1,1,1)
    plt.sca(ax1)  ##选择对ax2进行绘图
    ax1=plt.gca() #获得坐标轴的句柄
    #========= plot curve
    Lcdw, = ax1.plot(df_cdw["r"],df_cdw["corre_abs"],label="$<n_in_j>$-$<n_i><n_j>$",ls="-",lw=1.5,color="r",
                    marker='o',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="k",markerfacecolor='w')
    Lsc_yy, = ax1.plot(df_yy_sc["r"],df_yy_sc["corre_abs"],label="$<\Delta_i^y\Delta_j^y>$",ls="-",lw=1.5,color="blue",
                     marker='s',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="k",markerfacecolor='w')
    Lsc_zz, = ax1.plot(df_zz_sc["r"],df_zz_sc["corre_abs"],label="$<\Delta_i^z\Delta_j^z>$",ls="-",lw=1.5,color="green",
                    marker='^',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="k",markerfacecolor='w')
    #============ plot logr-logr curve
    Lcdw_pow, = ax1.plot(df_cdw["r"],df_cdw["fitcorre_pow"],label="$K$={:.4f}".format(-df_cdw["slope_pow"].values[0]),ls="-",lw=1.5,color="r",
                    marker='o',alpha=1,markersize=0,markeredgewidth=1.5, markeredgecolor="k",markerfacecolor='w')
    Lsc_yy_pow, = ax1.plot(df_yy_sc["r"],df_yy_sc["fitcorre_pow"],label="$K$={:.4f}".format(-df_yy_sc["slope_pow"].values[0]),ls="-",lw=1.5,color="blue",
                     marker='s',alpha=1,markersize=0,markeredgewidth=1.5, markeredgecolor="k",markerfacecolor='w')
    Lsc_zz_pow, = ax1.plot(df_zz_sc["r"],df_zz_sc["fitcorre_pow"],label="$K$={:.4f}".format(-df_zz_sc["slope_pow"].values[0]),ls="-",lw=1.5,color="green",
                    marker='^',alpha=1,markersize=0,markeredgewidth=1.5, markeredgecolor="k",markerfacecolor='w')
    ####图例设置
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 16, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[Lcdw,Lcdw_pow,Lsc_yy,Lsc_yy_pow,Lsc_zz,Lsc_zz_pow], loc = 4, bbox_to_anchor=(0.38, 0.02),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    #
    label_x = r"|i-j|"
    label_y = "correlation"
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
    plt.title("Jz=%g"%Jz,fontsize=20)
    plt.show()
    #
