# 对比不同 Jz 值下的green function
import numpy as np
import pandas as pd
import matplotlib.pylab as plt
#
def get_data_Lz1_Ly2(Lz=1,Ly=2,Lx=48,dop=18,t=3, J=1, Jz=0, dim=6000 ):
    # load data
    workpath = "E:\\WORK\\Work\\Project\\leg2-tj"
    filepath1 = "\\data_dmrgcpp\\Lz%d_Ly%d_Lx%d\\dop%g" % (Lz, Ly, Lx,dop)
    filepath2 = "\\t%d_J%d_Jz%d_dim%d" % (t, J, Jz, dim)
    filepath3 = "\\measurement_green_function.dat"
    filename = workpath + filepath1 + filepath2 + filepath3
    print(filename)
    df = pd.read_csv(filename, header=None, sep='\t',encoding='utf-8')
    df.rename(columns={0: "site1", 1: "site2", 2: "corre"},inplace=True)
    df = df[(df["site2"]-df["site1"]) % (Lz * Ly) == 0] 
    df.sort_values(["site2"],inplace=True)
    print(df.head())
    print(df.tail())
    print(len(df))
    #==============================================
    site1 = df["site1"].values
    site2 = df["site2"].values
    corre = np.abs(np.array(df["corre"].values))
    r = np.array((site2 - site1) / (Lz * Ly))
    print(r)
    return r, corre
#
def get_data(Lz, Ly, Lx, dop, t, J, Jz, dim):
    # load data
    workpath = "E:\\WORK\\Work\\Project\\La3Ni2O7"
    filepath1 = "\\data_dmrgcpp\\Lz%d_Ly%d_Lx%d\\dop%g" % (Lz, Ly, Lx,dop)
    filepath2 = "\\t%d_J%d_Jz%.2f_dim%d" % (t, J, Jz, dim)
    filepath3 = "\\measurement_green_function.dat"
    filename = workpath + filepath1 + filepath2 + filepath3
    print(filename)
    df = pd.read_csv(filename, header=None, sep='\t',encoding='utf-8')
    df.rename(columns={0: "site1", 1: "site2", 2: "corre"},inplace=True)
    df = df[(df["site2"]-df["site1"]) % (Lz * Ly) == 0] 
    df.sort_values(["site2"],inplace=True)
    print(df.head())
    print(df.tail())
    print(len(df))
    #==============================================
    site1 = df["site1"].values
    site2 = df["site2"].values
    corre = np.abs(np.array(df["corre"].values))
    r = np.array((site2 - site1) / (Lz * Ly))
    print(r)
    return r, corre
#
if __name__ == "__main__":
    print()
    Lz = 2
    Ly = 2
    Lx = 48
    dop = 72
    t = 3
    J = 1
    dim = 6000 # dim cutoff
    #load data
    r_Jz001, corre_Jz001 = get_data(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=0.01,dim=dim)
    r_Jz002, corre_Jz002 = get_data(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=0.02,dim=dim)
    r_Jz01, corre_Jz01 = get_data(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=0.1,dim=dim) 
    r_Jz02, corre_Jz02 = get_data(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=0.2,dim=dim)
    
    # load data of single layer leg-2 tj model
    r_SL, corre_SL = get_data_Lz1_Ly2(Lz=1,Ly=2,Lx=48,dop=36,t=3, J=1, Jz=0, dim=6000)
    #-----------------plot logr-r fig and logr-logr fig-------------------------
    fig = plt.figure(figsize=(6,10)) 
    ax1 = plt.subplot(1,1,1)
    plt.sca(ax1)  ##选择对ax2进行绘图
    ax1=plt.gca() #获得坐标轴的句柄
    LSL, = ax1.plot(r_SL,corre_SL,label="Jz={} (single layer)".format(0.0),ls="-",lw=1.5,color="r",
             marker='o',alpha=1,markersize=12,markeredgewidth=1.5, markeredgecolor="r",
             markerfacecolor='None')

    L001, = ax1.plot(r_Jz001,corre_Jz001,label="Jz={} (Bilayer)".format(0.01),ls="-",lw=1.5,color="magenta",
             marker='o',alpha=1,markersize=12,markeredgewidth=1.5, markeredgecolor="magenta",
             markerfacecolor='None')
    L002, = ax1.plot(r_Jz002,corre_Jz002,label="Jz={} (Bilayer)".format(0.02),ls="-",lw=1.5,color="cyan",
             marker='o',alpha=1,markersize=12,markeredgewidth=1.5, markeredgecolor="cyan",
             markerfacecolor='None')

    L01, = ax1.plot(r_Jz01,corre_Jz01,label="Jz={} (Bilayer)".format(0.1),ls="-",lw=1.5,color="blue",
             marker='o',alpha=1,markersize=12,markeredgewidth=1.5, markeredgecolor="blue",
             markerfacecolor='None')
    L02, = ax1.plot(r_Jz02,corre_Jz02,label="Jz={} (Bilayer)".format(0.2),ls="-",lw=1.5,color="g",
             marker='o',alpha=1,markersize=12,markeredgewidth=1.5, markeredgecolor="g",
             markerfacecolor='None')

    ####图例设置
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 20, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[LSL, L001, L002, L01, L02,], loc = 4, bbox_to_anchor=(0.38, 0.0),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    
    label_x = r"|i-j|"
    label_y = "Green function"
    plt.yscale("log")
    #plt.xscale("log")      
    ax1.set_xlabel(label_x, size= 25)
    ax1.set_ylabel(label_y, size= 25)
    ax1.tick_params(labelsize = 25) # 设置坐标刻度对应数字的大小
    #ax1.set_xlim([0,8])
    #ax1.set_ylim([-0.1,1])
    #ax1.set_xticks([0,2,4,6,8])
    #ax1.set_yticks([-0.1,0,0.5,1])
    #ax1.text(0.3,-0.05, r'$\mathrm{(a)}$', fontsize=18)
    plt.title("dop=%.4f, dim=%d "%((dop/(Lz*Ly*Lx)), dim),fontsize=25)
    plt.show()