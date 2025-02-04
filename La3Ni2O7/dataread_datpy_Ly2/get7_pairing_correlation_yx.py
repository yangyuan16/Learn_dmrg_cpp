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
    workpath = "E:\\WORK\\Work\\Project\\La3Ni2O7"
    filepath1 = "\\data_dmrgcpp\\Lz%d_Ly%d_Lx%d\\dop%g" % (Lz, Ly, Lx,dop)
    filepath2 = "\\t%d_J%d_Jz%.2f_dim%d" % (t, J, Jz, dim)
    filepath3 = "\\measurement_pairing_yx.dat" 
    filename = workpath + filepath1 + filepath2 + filepath3
    print(filename)
    # load the data
    df = pd.read_csv(filename, header=None, sep='\t',encoding='utf-8')
    df.rename(columns={0: "site1", 1: "site2", 2: "site3", 3: "site4", 4: "corre"}, inplace=True)
    #print(df.head())
    #print(df.tail())
    #print(len(df))
    df = df[(df["site3"]-df["site1"]) % (Lz * Ly) == 0] 
    df.sort_values(["site3"],inplace=True)
    print(df.head())
    print(df.tail())
    print(len(df))
    #
    #==============================================
    site1 = df["site1"].values
    site3 = df["site3"].values
    corre = np.array(df["corre"].values)
    r = np.array((site3 - site1) / (Lz * Ly))
    print(r)
    sign_list = []
    for it in range(len(corre)):
        if corre[it] > 0:
            sign_list.append(1)
        else:
            sign_list.append(-1)
    #-----------------plot corre-------------------------
    fig = plt.figure(figsize=(6,10)) 
    ax2 = plt.subplot(1,1,1)
    plt.sca(ax2)  ##选择对ax2进行绘图
    ax2=plt.gca() #获得坐标轴的句柄
    #corre1 = np.log10(corre)
    ax2.plot(r,corre,label=r"G",ls="-",lw=1.5,color="red",
             marker='o',alpha=1,markersize=10,markeredgewidth=1.5, markeredgecolor="red",
             markerfacecolor='None')
    label_x = r"|i-j|"
    label_y = "YX Pairing Corre."
    #plt.yscale("log") 
    ax2.set_xlabel(label_x, size= 24)
    ax2.set_ylabel(label_y, size= 24)
    ax2.tick_params(labelsize = 25) # 设置坐标刻度对应数字的大小
    #ax2.set_xlim([0,8])
    #ax2.set_ylim([-0.1,1])
    #ax2.set_xticks([0,2,4,6,8])
    #ax2.set_yticks([-0.1,0,0.5,1])
    #ax2.text(0.3,-0.05, r'$\mathrm{(a)}$', fontsize=18)
    plt.title("Jz=%.2f"%Jz,fontsize=25)
    plt.show()
    #-----------------plot sign of corre-------------------------
    fig = plt.figure(figsize=(10,3)) 
    ax3 = plt.subplot(1,1,1)
    plt.sca(ax3)  ##选择对ax2进行绘图
    ax3=plt.gca() #获得坐标轴的句柄
    #corre1 = np.log10(corre)
    ax3.plot(r,sign_list,label=r"G",ls="--",lw=1.5,color="blue",
             marker='o',alpha=1,markersize=10,markeredgewidth=1.5, markeredgecolor="red",
             markerfacecolor='None')
    label_x = r"|i-j|"
    #label_y = "Sign of YY Pairing Corre." 
    ax3.set_xlabel(label_x, size= 25)
    #ax3.set_ylabel(label_y, size= 25)
    ax3.tick_params(labelsize = 25) # 设置坐标刻度对应数字的大小
    #ax3.set_xlim([0,8])
    #ax3.set_ylim([-0.1,1])
    #ax3.set_xticks([0,2,4,6,8])
    #ax3.set_yticks([-0.1,0,0.5,1])
    #ax3.text(0.3,-0.05, r'$\mathrm{(a)}$', fontsize=18)
    plt.title("Sign of YX Pairing Corre. Jz=%.2f"%Jz, fontsize=25)
    plt.show()