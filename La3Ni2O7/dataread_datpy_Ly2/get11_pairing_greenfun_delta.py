import numpy as np
import pandas as pd
import matplotlib.pylab as plt
#
def get_zzpairng(Lz, Ly, Lx, t, J, Jz, dim, dop):
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
    return r_zz, corre_zz
#
def get_yypairing(Lz, Ly, Lx, t, J, Jz, dim, dop):
    workpath = "E:\\WORK\\Work\\Project\\La3Ni2O7"
    filepath1 = "\\data_dmrgcpp\\Lz%d_Ly%d_Lx%d\\dop%g" % (Lz, Ly, Lx,dop)
    filepath2 = "\\t%d_J%d_Jz%.2f_dim%d" % (t, J, Jz, dim)
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
    return r_yy, corre_yy
#
def get_greenfun2(Lz, Ly, Lx, t, J, Jz, dim, dop):
    workpath = "E:\\WORK\\Work\\Project\\La3Ni2O7"
    filepath1 = "\\data_dmrgcpp\\Lz%d_Ly%d_Lx%d\\dop%g" % (Lz, Ly, Lx,dop)
    filepath2 = "\\t%d_J%d_Jz%.2f_dim%d" % (t, J, Jz, dim)
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
    return r_greenfun, corre_greenfun2

if __name__ == "__main__":
    print()
    Lz = 2
    Ly = 2
    Lx = 48
    t = 3
    J = 1
    Jz = 1.0  
    dim = 6000 # dim cutoff
    #============== dop = 36, delta=0.1875
    dop = 36
    r_zz_d36, corre_zz_d36 = get_zzpairng(Lz=Lz, Ly=Ly, Lx=Lx, t=t, J=J, Jz=Jz, dim=dim, dop=dop)
    r_yy_d36, corre_yy_d36 = get_yypairing(Lz=Lz, Ly=Ly, Lx=Lx, t=t, J=J, Jz=Jz, dim=dim, dop=dop)
    r_greenfun_d36, corre_greenfun2_d36 = get_greenfun2(Lz=Lz, Ly=Ly, Lx=Lx, t=t, J=J, Jz=Jz, dim=dim, dop=dop)
    #
    #============== dop = 72, delta = 0.375
    dop = 72
    r_zz_d72, corre_zz_d72 = get_zzpairng(Lz=Lz, Ly=Ly, Lx=Lx, t=t, J=J, Jz=Jz, dim=dim, dop=dop)
    r_yy_d72, corre_yy_d72 = get_yypairing(Lz=Lz, Ly=Ly, Lx=Lx, t=t, J=J, Jz=Jz, dim=dim, dop=dop)
    r_greenfun_d72, corre_greenfun2_d72 = get_greenfun2(Lz=Lz, Ly=Ly, Lx=Lx, t=t, J=J, Jz=Jz, dim=dim, dop=dop)
    #
    #-----------------plot ax1-------------------------
    fig = plt.figure(figsize=(6,10)) 
    ax1 = plt.subplot(1,1,1)
    plt.sca(ax1)  ##选择对ax2进行绘图
    ax1=plt.gca() #获得坐标轴的句柄
    Lzz36, = ax1.plot(r_zz_d36,corre_zz_d36,label="$\delta=%.4f,<\Delta_i^{zz} \Delta_j^{zz}>$"%(0.1875),ls="-",lw=1.5,color="red",
             marker='o',alpha=1,markersize=10,markeredgewidth=1.5, markeredgecolor="red",
             markerfacecolor='None')
    
    Lyy36, = ax1.plot(r_yy_d36,corre_yy_d36,label="$\delta=%.4f,<\Delta_i^{yy} \Delta_j^{yy}>$"%(0.1875),ls="-",lw=1.5,color="green",
             marker='o',alpha=1,markersize=10,markeredgewidth=1.5, markeredgecolor="green",
             markerfacecolor='None')
    
    Lgreenfun36, = ax1.plot(r_greenfun_d36,corre_greenfun2_d36,label="$\delta=%.4f,G(r)^2/4$"%(0.1875),ls="-",lw=1.5,color="blue",
             marker='o',alpha=1,markersize=10,markeredgewidth=1.5, markeredgecolor="blue",
             markerfacecolor='None')
    
    Lzz72, = ax1.plot(r_zz_d72,corre_zz_d72,label="$\delta=%.3f,<\Delta_i^{zz} \Delta_j^{zz}>$"%(0.375),ls="-",lw=1.5,color="red",
             marker='^',alpha=1,markersize=10,markeredgewidth=1.5, markeredgecolor="red",
             markerfacecolor='None')
    
    Lyy72, = ax1.plot(r_yy_d72,corre_yy_d72,label="$\delta=%.3f,<\Delta_i^{yy} \Delta_j^{yy}>$"%(0.375),ls="-",lw=1.5,color="green",
             marker='^',alpha=1,markersize=10,markeredgewidth=1.5, markeredgecolor="green",
             markerfacecolor='None')
    
    Lgreenfun72, = ax1.plot(r_greenfun_d72,corre_greenfun2_d72,label="$\delta=%.3f,G(r)^2/4$"%(0.375),ls="-",lw=1.5,color="blue",
             marker='^',alpha=1,markersize=10,markeredgewidth=1.5, markeredgecolor="blue",
             markerfacecolor='None')
    
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 22, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[Lzz36, Lyy36, Lgreenfun36, Lzz72, Lyy72, Lgreenfun72 ], loc = 4, bbox_to_anchor=(0.48, 0.18),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False) # [Lzz36, Lyy36, Lgreenfun36, Lzz72, Lyy72, Lgreenfun72]
    
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