# 对比不同 Jz 值下的density correlation function
import numpy as np
import pandas as pd
import matplotlib.pylab as plt
#
def get_data(Lz, Ly, Lx, dop, t, J, Jz, dim):
    #
    workpath = "E:\\WORK\\Work\\Project\\La3Ni2O7"
    filepath1 = "\\data_dmrgcpp\\Lz%d_Ly%d_Lx%d\\dop%g" % (Lz, Ly, Lx,dop)
    filepath2 = "\\t%d_J%d_Jz%.2f_dim%d" % (t, J, Jz, dim)
    filepath3 = "\\measurement_pairing.dat"
    filename = workpath + filepath1 + filepath2 + filepath3
    print(filename)
    # load the data
    df = pd.read_csv(filename, header=None, sep='\t',encoding='utf-8')
    df.rename(columns={0: "site1", 1: "site2", 2: "site3", 3: "site4", 4: "corre"}, inplace=True)
    df = df[(df["site3"]-df["site1"]) % (Lz * Ly) == 0] 
    df.sort_values(["site3"],inplace=True)
    print(df.head())
    print(df.tail())
    print(len(df))
    #
    #==============================================
    site1 = df["site1"].values
    site3 = df["site3"].values
    corre = np.abs(np.array(df["corre"].values))
    r = np.array((site3 - site1) / (Lz * Ly))
    print(r)
    return r, corre
#

if __name__ == "__main__":
    print()
    Lz = 2
    Ly = 3
    Lx = 48
    dop = 144
    t = 3
    J = 1
    dim = 8000 # dim cutoff
    #load data
    r_Jz06, corre_Jz06 = get_data(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=0.6,dim=dim)
    r_Jz10, corre_Jz10 = get_data(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=1.0,dim=dim)
    r_Jz20, corre_Jz20 = get_data(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=2.0,dim=dim)
    r_Jz30, corre_Jz30 = get_data(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=3.0,dim=dim)
    r_Jz40, corre_Jz40 = get_data(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=4.0,dim=dim)
    #-----------------plot logr-r fig and logr-logr fig-------------------------
    fig = plt.figure(figsize=(6,10)) 
    ax1 = plt.subplot(1,1,1)
    plt.sca(ax1)  ##选择对ax2进行绘图
    ax1=plt.gca() #获得坐标轴的句柄
    L06, = ax1.plot(r_Jz06,corre_Jz06,label="Jz={}".format(0.6),ls="-",lw=1.5,color="r",
             marker='o',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="r",
             markerfacecolor='None')
    L10, = ax1.plot(r_Jz10,corre_Jz10,label="Jz={}".format(1.0),ls="-",lw=1.5,color="blue",
             marker='s',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="blue",
             markerfacecolor='None')
    L20, = ax1.plot(r_Jz20,corre_Jz20,label="Jz={}".format(2.0),ls="-",lw=1.5,color="green",
             marker='^',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="green",
             markerfacecolor='None')
    L30, = ax1.plot(r_Jz30,corre_Jz30,label="Jz={}".format(3.0),ls="-",lw=1.5,color="cyan",
             marker='v',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="cyan",
             markerfacecolor='None')
    L40, = ax1.plot(r_Jz40,corre_Jz40,label="Jz={}".format(4.0),ls="-",lw=1.5,color="magenta",
             marker='h',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="magenta",
             markerfacecolor='None')  
    ####图例设置
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 20, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[L06,L10,L20,L30,L40], loc = 4, bbox_to_anchor=(0.38, 0.02),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    
    label_x = r"|i-j|"
    label_y = "pairing correlation"
    plt.yscale("log")
    plt.xscale("log")      
    ax1.set_xlabel(label_x, size= 25)
    ax1.set_ylabel(label_y, size= 25)
    ax1.tick_params(labelsize = 25) # 设置坐标刻度对应数字的大小
    #ax1.set_xlim([0,8])
    #ax1.set_ylim([-0.1,1])
    #ax1.set_xticks([0,2,4,6,8])
    #ax1.set_yticks([-0.1,0,0.5,1])
    #ax1.text(0.3,-0.05, r'$\mathrm{(a)}$', fontsize=18)
    plt.title("ZZ pairing",fontsize=25)
    plt.show()