# 对比不同 Jz 值下的density correlation function
import numpy as np
import pandas as pd
import matplotlib.pylab as plt
#
def get_data(Lz, Ly, Lx, dop, t, J, Jz, dim, filetitle):
    #
    workpath = "E:\\WORK\\Work\\Project\\TrilayerW2-tj-Jz-openz"
    filepath1 = "\\data_dmrgcpp\\Lz%d_Ly%d_Lx%d\\dop%g" % (Lz, Ly, Lx,dop)
    filepath2 = "\\t%.2f_J%.2f_Jz%.2f_dim%d" % (t, J, Jz, dim)
    filepath3 = filetitle
    filename = workpath + filepath1 + filepath2 + filepath3
    print(filename)
    # load the data
    df = pd.read_csv(filename, header=None, sep='\t',encoding='utf-8')
    df.rename(columns={0: "site1", 1: "site2", 2: "ninj"},inplace=True)
    df = df[(df["site2"]-df["site1"]) % (Lz * Ly) == 0] 
    df.sort_values(["site2"],inplace=True)
    df.index = list(range(len(df)))
    print(df.head())
    print(df.tail())
    print(len(df))
    #
    #==============================================
    site1 = df["site1"].values
    site2 = df["site2"].values
    ninj = df["ninj"].values
    r = np.array((site2 - site1) / (Lz * Ly))
    print(r)
    #------------------------------------
    # calculate <ninj> - <ni><nj>
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
    ninj_list = df['ninj'].values
    corre_list = []
    for it in range(len(df)):
        #print('it: ', it)
        site1 = site1_list[it]
        site2 = site2_list[it]
        ni = df2[df2['site']==site1]['density'].values[0]
        nj = df2[df2['site']==site2]['density'].values[0]
        ncorre = ninj_list[it] - ni * nj  
        #print('site1:', site1, 'site2:', site2, 'ni:', ni, 'nj:', nj, 'corre2:', corre2)
        corre_list.append(ncorre)
    #
    df['corre'] = corre_list  
    corre = df['corre'].values
    corre_abs = np.abs(np.array(corre))
    #
    df["r"] = r
    df["corre_abs"] = corre_abs
    df["logr"] = np.log10(r)
    df["logcorre"] = np.log10(corre_abs)
    #
    print(df.head())
    print(df.tail())
    print(len(df))
    return df
#
if __name__ == "__main__":
    print()
    Lz = 3
    Ly = 2
    Lx = 48
    dop = 192
    t = 3
    J = 1
    Jz = 21.0
    Dim= 10000
    #load data
    file_up = "\\measurement_density_correlation_s1.dat" # up-layer
    file_in = "\\measurement_density_correlation.dat" # inner-layer
    file_dn = "\\measurement_density_correlation_s2.dat" # down-layer
    #
    df_up = get_data(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=Jz,dim=Dim,filetitle=file_up)
    df_in = get_data(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=Jz,dim=Dim,filetitle=file_in)
    df_dn = get_data(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=Jz,dim=Dim,filetitle=file_dn)
    #
    fig = plt.figure(figsize=(6,10)) 
    ax1 = plt.subplot(1,1,1)
    plt.sca(ax1)  ##选择对ax2进行绘图
    ax1=plt.gca() #获得坐标轴的句柄
    L_up, = ax1.plot(df_up["r"],df_up["corre_abs"],label="$<n_i\cdot n_j>-<n_i><n_j>$ (Up-layer)",ls="-",lw=1.5,color="red",
             marker='o',alpha=1,markersize=12,markeredgewidth=1.5, markeredgecolor="red",
             markerfacecolor='None')
    #
    L_in, = ax1.plot(df_in["r"],df_in["corre_abs"],label="$<n_i\cdot n_j>-<n_i><n_j>$ (In-layer)",ls="-",lw=1.5,color="blue",
             marker='o',alpha=1,markersize=12,markeredgewidth=1.5, markeredgecolor="blue",
             markerfacecolor='None')
    #
    L_dn, = ax1.plot(df_dn["r"],df_dn["corre_abs"],label="$<n_i\cdot n_j>-<n_i><n_j>$ (Dn-layer)",ls="-",lw=1.5,color="green",
             marker='o',alpha=1,markersize=12,markeredgewidth=1.5, markeredgecolor="green",
             markerfacecolor='None')
    #
    ####图例设置
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 20, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[L_up,L_in,L_dn], loc = 4, bbox_to_anchor=(0.68, 0.0),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    
    label_x = r"|i-j|"
    label_y = "Density Correlation"
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
    plt.title("Lx=%d,dop=%.4f,Jz=%.2f, dim=%d"%(Lx,dop/(Lz*Ly*Lx),Jz,Dim), fontsize=25)
    plt.show()