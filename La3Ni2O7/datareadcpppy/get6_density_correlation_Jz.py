# 对比不同 Jz 值下的density correlation function
import numpy as np
import pandas as pd
import matplotlib.pylab as plt

def get_data(Lz, Ly, Lx, dop, t, J, Jz, dim):
    # load data
    workpath = "E:\\WORK\\Work\\Project\\La3Ni2O7"
    filepath1 = "\\data_dmrgcpp\\Lz%d_Ly%d_Lx%d\\dop%g" % (Lz, Ly, Lx,dop)
    filepath2 = "\\t%d_J%d_Jz%g_dim%d" % (t, J, Jz, dim)
    filepath3 = "\\measurement_density_correlation_ref74.dat"
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
    #------------------------------------
    print(" calculate <ninj> - <ni><nj> " )
    filepath4 = "\\measurement_electron_density.dat"
    filename2 = workpath + filepath1 + filepath2 + filepath4
    print(filename2)
    df2 = pd.read_csv(filename2, header=None, sep='\t',encoding='utf-8')
    df2.rename(columns={0: "site", 1: "density"},inplace=True)
    df2.sort_values(["site"],inplace=True)
    print(df2.head())
    print(df2.tail())
    print(len(df2))
    #
    site1_list = df['site1'].values
    site2_list = df['site2'].values
    corre_list = df['corre'].values
    corre2_list = []
    for it in range(len(df)):
        #print('it: ', it)
        site1 = site1_list[it]
        site2 = site2_list[it]
        ni = df2[df2['site']==site1]['density'].values[0]
        nj = df2[df2['site']==site2]['density'].values[0]
        corre2 = corre_list[it] - ni * nj  
        #print('site1:', site1, 'site2:', site2, 'ni:', ni, 'nj:', nj, 'corre2:', corre2)
        corre2_list.append(corre2)
    #
    df['corre2'] = corre2_list  
    corre2 = np.abs(np.array(df['corre2'].values))
    return r, corre2
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
    r_Jz01, corre_Jz01 = get_data(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=0.1,dim=dim)
    r_Jz025, corre_Jz025 = get_data(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=0.25,dim=dim)
    r_Jz05, corre_Jz05 = get_data(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=0.5,dim=dim)
    r_Jz075, corre_Jz075 = get_data(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=0.75,dim=dim)
    r_Jz1, corre_Jz1 = get_data(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=1,dim=dim)
    r_Jz2, corre_Jz2 = get_data(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=2,dim=dim)
    r_Jz4, corre_Jz4 = get_data(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=4,dim=dim)
    r_Jz6, corre_Jz6 = get_data(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=6,dim=dim)
    r_Jz8, corre_Jz8 = get_data(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=8,dim=dim)
    #-----------------plot logr-r fig and logr-logr fig-------------------------
    fig = plt.figure(figsize=(6,10)) 
    ax1 = plt.subplot(1,1,1)
    plt.sca(ax1)  ##选择对ax2进行绘图
    ax1=plt.gca() #获得坐标轴的句柄
    
    L01, = ax1.plot(r_Jz01,corre_Jz01,label="Jz={}".format(0.1),ls="-",lw=1.5,color="r",
             marker='o',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="k",
             markerfacecolor='w')
    L025, = ax1.plot(r_Jz025,corre_Jz025,label="Jz={}".format(0.25),ls="-",lw=1.5,color="r",
             marker='^',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="k",
             markerfacecolor='w')
    L05, = ax1.plot(r_Jz05,corre_Jz05,label="Jz={}".format(0.5),ls="-",lw=1.5,color="r",
             marker='v',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="k",
             markerfacecolor='w')
    L075, = ax1.plot(r_Jz075,corre_Jz075,label="Jz={}".format(0.75),ls="-",lw=1.5,color="r",
             marker='s',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="k",
             markerfacecolor='w')
    L1, = ax1.plot(r_Jz1,corre_Jz1,label="Jz={}".format(1),ls="-",lw=1.5,color="blue",
             marker='o',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="k",
             markerfacecolor='k')
    L2, = ax1.plot(r_Jz2,corre_Jz2,label="Jz={}".format(2),ls="-",lw=1.5,color="blue",
             marker='^',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="k",
             markerfacecolor='k')
    L4, = ax1.plot(r_Jz4,corre_Jz4,label="Jz={}".format(4),ls="-",lw=1.5,color="blue",
             marker='v',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="k",
             markerfacecolor='k')
    L6, = ax1.plot(r_Jz6,corre_Jz6,label="Jz={}".format(6),ls="-",lw=1.5,color="blue",
             marker='s',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="k",
             markerfacecolor='k')
    L8, = ax1.plot(r_Jz8,corre_Jz8,label="Jz={}".format(8),ls="-",lw=1.5,color="blue",
             marker='h',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="k",
             markerfacecolor='k')
    
    ####图例设置
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 16, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[L01,L025,L05,L075,L1,L2,L4,L6,L8], loc = 4, bbox_to_anchor=(0.38, 0.02),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    
    label_x = r"|i-j|"
    label_y = "Density correlation"
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