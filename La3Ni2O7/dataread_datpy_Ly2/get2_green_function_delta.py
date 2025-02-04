# 对比不同 Jz 值下的green function
import numpy as np
import pandas as pd
import matplotlib.pylab as plt

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
    t = 3
    J = 1
    Jz = 1.0
    dim = 6000 # dim cutoff
    #load data
    r_01875, corre_01875 = get_data(Lz=Lz,Ly=Ly,Lx=Lx,dop=36,t=t,J=J,Jz=Jz,dim=dim) 
    r_025, corre_025 = get_data(Lz=Lz,Ly=Ly,Lx=Lx,dop=48,t=t,J=J,Jz=Jz,dim=dim)
    r_0375, corre_0375 = get_data(Lz=Lz,Ly=Ly,Lx=Lx,dop=72,t=t,J=J,Jz=Jz,dim=dim)
    #-----------------plot logr-r fig and logr-logr fig-------------------------
    fig = plt.figure(figsize=(6,10)) 
    ax1 = plt.subplot(1,1,1)
    plt.sca(ax1)  ##选择对ax2进行绘图
    ax1=plt.gca() #获得坐标轴的句柄
    L01875, = ax1.plot(r_01875,corre_01875,label="$\delta$={}".format(0.1875),ls="-",lw=1.5,color="blue",
             marker='o',alpha=1,markersize=12,markeredgewidth=1.5, markeredgecolor="blue",
             markerfacecolor='None')
    L025, = ax1.plot(r_025,corre_025,label="$\delta$={}".format(0.25),ls="-",lw=1.5,color="r",
             marker='o',alpha=1,markersize=12,markeredgewidth=1.5, markeredgecolor="r",
             markerfacecolor='None')
    L0375, = ax1.plot(r_0375,corre_0375,label="$\delta$={}".format(0.375),ls="-",lw=1.5,color="green",
             marker='o',alpha=1,markersize=12,markeredgewidth=1.5, markeredgecolor="green",
             markerfacecolor='None')
    ####图例设置
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 20, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[L01875,L025,L0375], loc = 4, bbox_to_anchor=(0.38, 0.0),
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
    plt.title("Jz = %.2f, dim=%d"%(Jz, dim),fontsize=25)
    plt.show()