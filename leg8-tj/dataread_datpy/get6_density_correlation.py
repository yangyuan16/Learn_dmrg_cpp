import numpy as np
import pandas as pd
import matplotlib.pylab as plt
#
def plot_corre(r,corre,text, label_y):
    #-----------------plot <ninj>-------------------------
    fig = plt.figure(figsize=(6,10)) 
    ax1 = plt.subplot(3,1,1)
    ax2 = plt.subplot(3,1,2)
    ax3 = plt.subplot(3,1,3)
    plt.sca(ax1)  ##选择对ax1进行绘图
    ax1=plt.gca() #获得坐标轴的句柄
    text = text
    ax1.text(0.3, 0.009, text)
    ax1.plot(r,corre,label=r"G",ls="-",lw=1.5,color="red",
             marker='o',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="k",
             markerfacecolor='w')
    label_x = r"|i-j|"
    label_y = label_y
    ax1.set_xlabel(label_x, size= 14)
    ax1.set_ylabel(label_y, size= 14)
    ax1.tick_params(labelsize = 15) # 设置坐标刻度对应数字的大小
    #ax1.set_xlim([0,8])
    #ax1.set_ylim([-0.1,1])
    #ax1.set_xticks([0,2,4,6,8])
    #ax1.set_yticks([-0.1,0,0.5,1])
    #ax1.text(0.3,-0.05, r'$\mathrm{(a)}$', fontsize=18)
    #-----------------plot ax2--------------------------
    plt.sca(ax2)  ##选择对ax2进行绘图
    ax2=plt.gca() #获得坐标轴的句柄
    corre_abs = np.abs(np.array(corre))
    #corre1 = np.log10(corre_abs)
    ax2.plot(r,corre_abs,label=r"G",ls="-",lw=1.5,color="red",
             marker='o',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="k",
             markerfacecolor='w')
    label_x = r"|i-j|"
    plt.yscale("log")
    ax2.set_xlabel(label_x, size= 14)
    ax2.set_ylabel( label_y, size= 14)
    ax2.tick_params(labelsize = 15) # 设置坐标刻度对应数字的大小
    #ax2.set_xlim([0,8])
    #ax2.set_ylim([-0.1,1])
    #ax2.set_xticks([0,2,4,6,8])
    #ax2.set_yticks([-0.1,0,0.5,1])
    #ax2.text(0.3,-0.05, r'$\mathrm{(a)}$', fontsize=18)
    #----------------plot ax3---------------------------
    plt.sca(ax3)  ##选择对ax2进行绘图
    ax3=plt.gca() #获得坐标轴的句柄
    corre_abs = np.abs(np.array(corre))
    #corre1 = np.log10(corre_abs)
    #r1 = np.log(r)
    ax3.plot(r,corre_abs,label=r"G",ls="-",lw=1.5,color="red",
             marker='o',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="k",
             markerfacecolor='w')
    label_x = r"(|i-j|)"
    plt.xscale("log")
    plt.yscale("log")
    ax3.set_xlabel(label_x, size= 14)
    ax3.set_ylabel( label_y, size= 14)
    ax3.tick_params(labelsize = 15) # 设置坐标刻度对应数字的大小
    #ax3.set_xlim([0,8])
    #ax3.set_ylim([-0.1,1])
    #ax3.set_xticks([0,2,4,6,8])
    #ax3.set_yticks([-0.1,0,0.5,1])
    #ax3.text(0.3,-0.05, r'$\mathrm{(a)}$', fontsize=18)
    #
    plt.show()
    return

if __name__ == "__main__":
    print()
    Lz = 1
    Ly = 8
    Lx = 48
    dop = 48
    t = 3
    J = 1
    Jz = 0
    dim = 8000 # dim cutoff
    workpath = "E:\\WORK\\Work\\Project\\leg8-tj"
    filepath1 = "\\data_dmrgcpp\\Lz%d_Ly%d_Lx%d\\dop%g" % (Lz, Ly, Lx,dop)
    filepath2 = "\\t%d_J%d_Jz%g_dim%d" % (t, J, Jz, dim)
    filepath3 = "\\measurement_density_correlation.dat" 
    filename = workpath + filepath1 + filepath2 + filepath3
    print(filename)
    # load the data
    df = pd.read_csv(filename, header=None, sep='\t',encoding='utf-8')
    df.rename(columns={0: "site1", 1: "site2", 2: "corre"},inplace=True)
    df = df[(df["site2"]-df["site1"]) % (Lz * Ly) == 0] 
    df.sort_values(["site2"],inplace=True)
    print(df.head())
    print(df.tail())
    print(len(df))
    #
    #==============================================
    site1 = df["site1"].values
    site2 = df["site2"].values
    corre = df["corre"].values
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
    corre2 = df['corre2'].values 
    #--------------------------------------------------------------------
    text = "Lz%d_Ly%d_Lx%d_dop%g_" % (Lz, Ly, Lx,dop) + "t%d_J%d_Jz%g_dim%d" % (t, J, Jz, dim)
    label_y = '<ninj>'
    plot_corre(r=r,corre=corre,text=text,label_y=label_y)
    
    label_y = '<ninj>-<ni><nj>'
    plot_corre(r=r,corre=corre2,text=text,label_y=label_y)
    #--------------------------------------------------------------------
    # plot <ninj> - <ni><nj>
    







