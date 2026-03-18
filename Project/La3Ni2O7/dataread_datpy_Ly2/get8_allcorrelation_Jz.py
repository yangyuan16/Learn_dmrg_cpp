#
import numpy as np
import pandas as pd
import matplotlib.pylab as plt
#
def get_data_sisj(Lz, Ly, Lx, dop, t, J, Jz, dim): # 自旋关联
    # load data
    workpath = "E:\\WORK\\Work\\Project\\La3Ni2O7"
    filepath1 = "\\data_dmrgcpp\\Lz%d_Ly%d_Lx%d\\dop%g" % (Lz, Ly, Lx,dop)
    filepath2 = "\\t%d_J%d_Jz%.2f_dim%d" % (t, J, Jz, dim)
    filepath3 = "\\measurement_spin_correlation.dat"
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
def get_data_cicj(Lz, Ly, Lx, dop, t, J, Jz, dim): # green function
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
def get_data_ninj(Lz, Ly, Lx, dop, t, J, Jz, dim):  # density correlation
    # load data
    workpath = "E:\\WORK\\Work\\Project\\La3Ni2O7"
    filepath1 = "\\data_dmrgcpp\\Lz%d_Ly%d_Lx%d\\dop%g" % (Lz, Ly, Lx,dop)
    filepath2 = "\\t%d_J%d_Jz%.2f_dim%d" % (t, J, Jz, dim)
    filepath3 = "\\measurement_density_correlation.dat"
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
def get_data_pairing_yy(Lz, Ly, Lx, dop, t, J, Jz, dim,): # yy pairing
    #
    workpath = "E:\\WORK\\Work\\Project\\La3Ni2O7"
    filepath1 = "\\data_dmrgcpp\\Lz%d_Ly%d_Lx%d\\dop%g" % (Lz, Ly, Lx,dop)
    filepath2 = "\\t%d_J%d_Jz%.2f_dim%d" % (t, J, Jz, dim)
    filepath3 = "\\measurement_pairing_yy.dat"
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
def get_data_pairing_zz(Lz, Ly, Lx, dop, t, J, Jz, dim): # zz pairing
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
    Ly = 2
    Lx = 48
    dop = 96
    t = 3
    J = 1
    Jz_list = [0.1,0.25,0.5,0.75,1.0,1.25,1.50,1.75,2.0,2.25,2.5,2.75,3.0,3.25,3.5,3.75,4.0]
    dim = 6000 # dim cutoff
    for it in range(len(Jz_list)):
        Jz = Jz_list[it]
        r_sisj, corre_sisj = get_data_sisj(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=Jz,dim=dim)
        r_cicj, corre_cicj = get_data_cicj(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=Jz,dim=dim)
        r_ninj, corre_ninj = get_data_ninj(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=Jz,dim=dim)
        r_pyy, corre_pyy = get_data_pairing_yy(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=Jz,dim=dim)
        r_pzz, corre_pzz = get_data_pairing_zz(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=Jz,dim=dim)
        #-----------------plot logr-r fig and logr-logr fig-------------------------
        fig = plt.figure(figsize=(6,10)) 
        ax1 = plt.subplot(1,1,1)
        plt.sca(ax1)  ##选择对ax2进行绘图
        ax1=plt.gca() #获得坐标轴的句柄
        Lsisj, = ax1.plot(r_sisj,corre_sisj,label="<SiSj>",ls="-",lw=1.5,color="r",
                 marker='o',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="k",
                 markerfacecolor='w')
        Lcicj, = ax1.plot(r_cicj,corre_cicj,label=r"$\langle c_i^+ c_j \rangle$",ls="-",lw=1.5,color="green",
                 marker='^',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="k",
                 markerfacecolor='w')
        Lninj, = ax1.plot(r_ninj,corre_ninj,label="$<n_in_j>$",ls="-",lw=1.5,color="blue",
                 marker='v',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="k",
                 markerfacecolor='w')
        Lpyy, = ax1.plot(r_pyy,corre_pyy,label="$<\Delta^{+,yy}_i \Delta_j>$",ls="-",lw=1.5,color="magenta",
                 marker='s',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="k",
                 markerfacecolor='w')
        Lpzz, = ax1.plot(r_pzz,corre_pzz,label="$<\Delta^{+,zz}_i \Delta_j>$",ls="-",lw=1.5,color="cyan",
                 marker='o',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="k",
                 markerfacecolor='k')
        ####图例设置
        legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 16, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
        legend1=plt.legend(handles=[Lsisj,Lcicj,Lninj,Lpyy,Lpzz], loc = 4, bbox_to_anchor=(0.38, 0.02),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    
        label_x = r"|i-j|"
        label_y = "correlation"
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
        plt.title("Jz=%.2f"%Jz)
        plt.show()