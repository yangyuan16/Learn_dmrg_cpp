import numpy as np
import pandas as pd
import matplotlib.pylab as plt
#
if __name__ == "__main__":
    print()
    Lz = 2
    Ly = 2
    Lx = 48
    dop = 96
    t = 3
    J = 1
    Jz = 3.0  
    dim = 6000 # dim cutoff
    #============== ZZ pairing data ===========================
    workpath = "E:\\WORK\\Work\\Project\\La3Ni2O7"
    filepath1 = "\\data_dmrgcpp\\Lz%d_Ly%d_Lx%d\\dop%g" % (Lz, Ly, Lx,dop)
    filepath2 = "\\t%d_J%d_Jz%.2f_dim%d" % (t, J, Jz, dim)
    filepath3 = "\\measurement_pairing.dat" 
    filename = workpath + filepath1 + filepath2 + filepath3
    print(filename)
    # load the data
    df_zz = pd.read_csv(filename, header=None, sep='\t',encoding='utf-8')
    df_zz.rename(columns={0: "site1", 1: "site2", 2: "site3", 3: "site4", 4: "corre"}, inplace=True)
    print(df_zz.head())
    print(df_zz.tail())
    print(len(df_zz))
    df_zz = df_zz[(df_zz["site3"]-df_zz["site1"]) % (Lz * Ly) == 0] 
    df_zz.sort_values(["site3"],inplace=True)
    #print(df.head())
    #print(df.tail())
    #print(len(df))
    site1 = df_zz["site1"].values
    site3 = df_zz["site3"].values
    corre_zz = np.abs(np.array(df_zz["corre"].values))
    r_zz = np.array((site3 - site1) / (Lz * Ly))
    print(r_zz)
    #============= YY pairing data ===========================
    filepath3 = "\\measurement_pairing_yy.dat" 
    filename = workpath + filepath1 + filepath2 + filepath3
    print(filename)
    # load the data
    df_yy = pd.read_csv(filename, header=None, sep='\t',encoding='utf-8')
    df_yy.rename(columns={0: "site1", 1: "site2", 2: "site3", 3: "site4", 4: "corre"}, inplace=True)
    print(df_yy.head())
    print(df_yy.tail())
    print(len(df_yy))
    df_yy = df_yy[(df_yy["site3"]-df_yy["site1"]) % (Lz * Ly) == 0] 
    df_yy.sort_values(["site3"],inplace=True)
    site1 = df_yy["site1"].values
    site3 = df_yy["site3"].values
    corre_yy = np.abs(np.array(df_yy["corre"].values))
    r_yy = np.array((site3 - site1) / (Lz * Ly))
    print(r_yy)
    #============= XX pairing data ===========================
    filepath3 = "\\measurement_pairing_xx.dat" 
    filename = workpath + filepath1 + filepath2 + filepath3
    print(filename)
    # load the data
    df_xx = pd.read_csv(filename, header=None, sep='\t',encoding='utf-8')
    df_xx.rename(columns={0: "site1", 1: "site2", 2: "site3", 3: "site4", 4: "corre"}, inplace=True)
    print(df_xx.head())
    print(df_xx.tail())
    print(len(df_xx))
    df_xx = df_xx[(df_xx["site3"]-df_xx["site1"]) % (Lz * Ly) == 0] 
    df_xx.sort_values(["site3"],inplace=True)
    site1 = df_xx["site1"].values
    site3 = df_xx["site3"].values
    corre_xx = np.abs(np.array(df_xx["corre"].values))
    r_xx = np.array((site3 - site1) / (Lz * Ly))
    print(r_xx)
    #============= Green Function data ========================
    filepath3 = "\\measurement_green_function.dat"
    filename = workpath + filepath1 + filepath2 + filepath3
    print(filename)
    # load the data
    df_greenfun = pd.read_csv(filename, header=None, sep='\t',encoding='utf-8')
    df_greenfun.rename(columns={0: "site1", 1: "site2", 2: "corre"},inplace=True)
    df_greenfun = df_greenfun[(df_greenfun["site2"]-df_greenfun["site1"]) % (Lz * Ly) == 0] 
    df_greenfun.sort_values(["site2"],inplace=True)
    print(df_greenfun.head())
    print(df_greenfun.tail())
    print(len(df_greenfun))
    site1 = df_greenfun["site1"].values
    site2 = df_greenfun["site2"].values
    corre_greenfun2 = np.abs(np.array(df_greenfun["corre"].values)) * np.abs(np.array(df_greenfun["corre"].values)) / 4
    r_greenfun = np.array((site2 - site1) / (Lz * Ly))
    print(r_greenfun)
    #-----------------plot ax1-------------------------
    fig = plt.figure(figsize=(6,10)) 
    ax1 = plt.subplot(1,1,1)
    plt.sca(ax1)  ##选择对ax2进行绘图
    ax1=plt.gca() #获得坐标轴的句柄
    Lzz, = ax1.plot(r_zz,corre_zz,label="$<\Delta_i^{zz} \Delta_j^{zz}>$",ls="-",lw=1.5,color="red",
             marker='o',alpha=1,markersize=10,markeredgewidth=1.5, markeredgecolor="red",
             markerfacecolor='None')
    
    Lyy, = ax1.plot(r_yy,corre_yy,label="$<\Delta_i^{yy} \Delta_j^{yy}>$",ls="-",lw=1.5,color="green",
             marker='o',alpha=1,markersize=10,markeredgewidth=1.5, markeredgecolor="green",
             markerfacecolor='None')
    
    Lxx, = ax1.plot(r_xx,corre_xx,label="$<\Delta_i^{xx} \Delta_j^{xx}>$",ls="-",lw=1.5,color="brown",
             marker='o',alpha=1,markersize=10,markeredgewidth=1.5, markeredgecolor="brown",
             markerfacecolor='None')
    
    Lgreenfun, = ax1.plot(r_greenfun,corre_greenfun2,label="$G(r)^2/4$",ls="-",lw=1.5,color="blue",
             marker='o',alpha=1,markersize=10,markeredgewidth=1.5, markeredgecolor="blue",
             markerfacecolor='None')
    
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 22, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[Lzz,Lyy,Lxx,Lgreenfun], loc = 4, bbox_to_anchor=(0.48, 0.18),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    
    label_x = r"(|i-j|)"
    label_y = r"Correlation"
    plt.xscale("log") 
    plt.yscale("log") 
    ax1.set_xlabel(label_x, size= 25)
    ax1.set_ylabel(label_y, size= 25)
    ax1.tick_params(labelsize = 25) # 设置坐标刻度对应数字的大小
    #ax1.set_xlim([0,8])
    #ax1.set_ylim([-0.1,1])
    #ax1.set_xticks([0,2,4,6,8])
    #ax1.set_yticks([-0.1,0,0.5,1])
    #ax1.text(0.3,-0.05, r'$\mathrm{(a)}$', fontsize=18)
    #
    plt.title("Jz=%.2f"%Jz,fontsize=25)
    plt.show()