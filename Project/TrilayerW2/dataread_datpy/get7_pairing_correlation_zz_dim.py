# 对比不同 Jz 值下的density correlation function
import numpy as np
import pandas as pd
import matplotlib.pylab as plt
#
def get_data(Lz, Ly, Lx, dop, t, J, Jz, dim):
    #
    workpath = "E:\\WORK\\Work\\Project\\TrilayerW2-tj-Jz-openz"
    filepath1 = "\\data_dmrgcpp\\Lz%d_Ly%d_Lx%d\\dop%g" % (Lz, Ly, Lx,dop)
    filepath2 = "\\t%.2f_J%.2f_Jz%.2f_dim%d" % (t, J, Jz, dim)
    filepath3 = "\\measurement_pairing_zz.dat"
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
    Lz = 3
    Ly = 2
    Lx = 48
    dop = 108
    t = 3
    J = 1
    Jz = 12.0
    #load data
    r_dim6000, corre_dim6000 = get_data(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=Jz,dim=6000)
    #r_Jz12_dim8000, corre_Jz12_dim8000 = get_data(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=Jz,dim=8000)
    r_dim10000, corre_dim10000 = get_data(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=Jz,dim=10000)
    r_dim12000, corre_dim12000 = get_data(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=Jz,dim=12000)
    r_dim14000, corre_dim14000 = get_data(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=Jz,dim=14000)
    #-----------------plot logr-r fig and logr-logr fig-------------------------
    fig = plt.figure(figsize=(6,10)) 
    ax1 = plt.subplot(1,1,1)
    plt.sca(ax1)  ##选择对ax2进行绘图
    ax1=plt.gca() #获得坐标轴的句柄
    L_6000, = ax1.plot(r_dim6000,corre_dim6000,label="Dim={}".format(6000),ls="-",lw=1.5,color="blue",
             marker='o',alpha=1,markersize=12,markeredgewidth=1.5, markeredgecolor="blue",
             markerfacecolor='None')
    #L12_8000, = ax1.plot(r_Jz12_dim8000,corre_Jz12_dim8000,label="Dim={}".format(8000),ls="-",lw=1.5,color="r",
    #         marker='s',alpha=1,markersize=12,markeredgewidth=1.5, markeredgecolor="r",
    #         markerfacecolor='None')
    L_10000, = ax1.plot(r_dim10000,corre_dim10000,label="Dim={}".format(10000),ls="-",lw=1.5,color="green",
             marker='^',alpha=1,markersize=12,markeredgewidth=1.5, markeredgecolor="green",
             markerfacecolor='None')
    #
    L_12000, = ax1.plot(r_dim12000,corre_dim12000,label="Dim={}".format(12000),ls="-",lw=1.5,color="magenta",
             marker='s',alpha=1,markersize=12,markeredgewidth=1.5, markeredgecolor="magenta",
             markerfacecolor='None')
    
    L_14000, = ax1.plot(r_dim14000,corre_dim14000,label="Dim={}".format(14000),ls="-",lw=1.5,color="cyan",
             marker='v',alpha=1,markersize=12,markeredgewidth=1.5, markeredgecolor="cyan",
             markerfacecolor='None')

    ####图例设置
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 20, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[L_6000,L_10000,L_12000,L_14000], loc = 4, bbox_to_anchor=(0.38, 0.0),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    
    label_x = r"|i-j|"
    label_y = "ZZ pairing correlation"
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
    plt.title("Lx=%d,dop=%.4f,Jz=%.2f,"%(Lx,dop/(Lz*Ly*Lx),Jz,), fontsize=25)
    plt.show()