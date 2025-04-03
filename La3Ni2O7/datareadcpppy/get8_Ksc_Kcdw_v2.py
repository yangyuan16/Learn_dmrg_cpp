# 不同 Jz 下 Ksc * K cdw

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
    dim = 6000 # dim cutoff
    bonds_zz = "zz" # zz or yy 方向的pairing
    bonds_yy = "yy" # yy 方向的pairing
    #
    #load data
    Jz_list = [0.1,0.25,0.5,0.75,1,2,4,6,8]
    Kcdw_list = []
    Ksc_yy_list = []
    Ksc_zz_list = []
    for it in range(len(Jz_list)):
        Jz = Jz_list[it]
        df_cdw = get_data_density(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=Jz,dim=dim)
        Kcdw_list.append(-df_cdw["slope_pow"].values[0])
        df_yy_sc = get_data_pairng(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=Jz,dim=dim,bonds=bonds_yy)
        Ksc_yy_list.append(-df_yy_sc["slope_pow"].values[0])
        df_zz_sc = get_data_pairng(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=Jz,dim=dim,bonds=bonds_zz)
        Ksc_zz_list.append(-df_zz_sc["slope_pow"].values[0])
    #=======================================================
    KK_cdw_yy = np.array(Kcdw_list) * np.array(Ksc_yy_list) # 
    KK_cdw_zz = np.array(Kcdw_list) * np.array(Ksc_zz_list) #
    #
    fig = plt.figure(figsize=(6,10)) 
    ax1 = plt.subplot(1,1,1)
    plt.sca(ax1)  ##选择对ax2进行绘图
    ax1=plt.gca() #获得坐标轴的句柄
    #========= plot curve
    label = r"$K_{CDW} \cdot K_{sc}^{yy}$"  # !!! change this
    L1, = ax1.plot(Jz_list,KK_cdw_yy,label=label,ls="-",lw=1.5,color="red", #!!! change KK_cdw_yy
                   marker='v',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="k",markerfacecolor='w')
    L2, = ax1.plot([0.1,8],[1,1],label=r" ",ls="--",lw=1.5,color="k",
                   marker='s',alpha=1,markersize=0,markeredgewidth=1.5, markeredgecolor="k",markerfacecolor='w')
    ####图例设置
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 16, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[L1], loc = 4, bbox_to_anchor=(0.68, 0.42),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    label_x = r"Jz"
    label_y = label
    #plt.yscale("log")
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
    #
