import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pylab as plt
#
def get_data(Lz, Ly, Lx, dop, t, J, Jz, dim,):
    workpath = "E:\\WORK\\Work\\Project\\La3Ni2O7"
    filepath1 = "\\data_dmrgcpp\\Lz%d_Ly%d_Lx%d\\dop%g" % (Lz, Ly, Lx,dop)
    filepath2 = "\\t%d_J%d_Jz%.2f_dim%d" % (t, J, Jz, dim)
    filepath3 = "\\measurement_electron_density.dat"
    filename = workpath + filepath1 + filepath2 + filepath3
    print(filename)
    # load the data
    df = pd.read_csv(filename, header=None, sep='\t',encoding='utf-8')
    df.rename(columns={0: "site", 1: "density"},inplace=True)
    df.sort_values(['site'],inplace=True)
    print(df.head())
    print(len(df))
    #------- 沿着 x 方向的 electron density ----------------
    df_out = pd.DataFrame(columns = ["r0","dy0","r1","dy1","r2","dy2","r3","dy3","r4","dy4","r5","dy5","rmean","dymean"])
    #
    df_y0 = df[df['site'] % (Ly * Lz) ==0]
    df_out["r0"] = df_y0["site"].values
    df_out["dy0"] = df_y0["density"].values
    
    df_y1 = df[df['site'] % (Ly * Lz) ==1]
    df_out["r1"] = df_y1["site"].values
    df_out["dy1"] = df_y1["density"].values
    
    df_y2 = df[df['site'] % (Ly * Lz) ==2]
    df_out["r2"] = df_y2["site"].values
    df_out["dy2"] = df_y2["density"].values

    df_y3 = df[df['site'] % (Ly * Lz) ==3]
    df_out["r3"] = df_y3["site"].values
    df_out["dy3"] = df_y3["density"].values

    df_y4 = df[df['site'] % (Ly * Lz) ==4]
    df_out["r4"] = df_y4["site"].values
    df_out["dy4"] = df_y4["density"].values

    df_y5 = df[df['site'] % (Ly * Lz) ==5]
    df_out["r5"] = df_y5["site"].values
    df_out["dy5"] = df_y5["density"].values

    r_mean = range(len(df_y0))
    density_mean = (np.array(df_out["dy0"].values) + np.array(df_out["dy1"].values) + np.array(df_out["dy2"].values)  + 
                    np.array(df_out["dy3"].values) + np.array(df_out["dy4"].values) + np.array(df_out["dy5"].values)) / (Lz*Ly)
    #
    df_out["rmean"] = r_mean
    df_out["dymean"] = density_mean
    #
    return r_mean, density_mean 


if __name__ == "__main__":
    Lz = 2
    Ly = 3
    Lx = 48
    dop = 144
    t = 3
    J = 1
    Jz = 0.4 
    #
    r_dim4000, corre_dim4000 = get_data(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=Jz,dim=4000)
    #
    r_dim6000, corre_dim6000 = get_data(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=Jz,dim=6000)
    #
    r_dim8000, corre_dim8000 = get_data(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=Jz,dim=8000)
    
    r_dim10000, corre_dim10000 = get_data(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=Jz,dim=10000)

    r_dim12000, corre_dim12000 = get_data(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=Jz,dim=12000)
    #------------------------------------------------------------------------------------
    fig = plt.figure(figsize=(6,10)) 
    ax1 = plt.subplot(1,1,1)
    plt.sca(ax1)  ##选择对ax2进行绘图
    ax1=plt.gca() #获得坐标轴的句柄
    
    #L2000,= ax1.plot(r_dim2000,corre_dim2000,label="dim={}".format(2000),ls="-",lw=1.5,color="magenta",
    #         marker='h',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="magenta",
    #         markerfacecolor='None')  
    #
    #L3000,= ax1.plot(r_dim3000,corre_dim3000,label="dim={}".format(3000),ls="-",lw=1.5,color="cyan",
    #         marker='*',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="cyan",
    #         markerfacecolor='None')    
    
    L4000,= ax1.plot(r_dim4000,corre_dim4000,label="dim={}".format(4000),ls="-",lw=1.5,color="green",
             marker='+',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="green",
             markerfacecolor='None')  
    #
    L6000,= ax1.plot(r_dim6000,corre_dim6000,label="dim={}".format(6000),ls="-",lw=1.5,color="r",
             marker='o',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="r",
             markerfacecolor='None')
    
    L8000,= ax1.plot(r_dim8000,corre_dim8000,label="dim={}".format(8000),ls="-",lw=1.5,color="cyan",
             marker='o',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="cyan",
             markerfacecolor='None')    


    L10000,= ax1.plot(r_dim10000,corre_dim10000,label="dim={}".format(10000),ls="-",lw=1.5,color="b",
             marker='s',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="b",
             markerfacecolor='None')
    
    L12000,= ax1.plot(r_dim12000,corre_dim12000,label="dim={}".format(12000),ls="-",lw=1.5,color="magenta",
             marker='s',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="magenta",
             markerfacecolor='None')

    #
    ####图例设置
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 16, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[L4000, L6000, L8000, L10000, L12000], loc = 4, bbox_to_anchor=(0.38, 0.18),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    
    label_x = r"|i|"
    label_y = "charge density"    
    ax1.set_xlabel(label_x, size= 14)  
    ax1.set_ylabel(label_y, size= 14)
    ax1.tick_params(labelsize = 15) # 设置坐标刻度对应数字的大小
    plt.title("<ni> Jz=%g"%Jz,fontsize=25)
    #ax1.set_xlim([0,8])
    #ax1.set_ylim([-0.1,1])
    #ax1.set_xticks([0,2,4,6,8])
    #ax1.set_yticks([-0.1,0,0.5,1])
    #ax1.text(0.3,-0.05, r'$\mathrm{(a)}$', fontsize=18)
    ax1.set_xlabel(label_x, size= 20)
    ax1.set_ylabel(label_y, size= 20)
    ax1.tick_params(labelsize = 20) # 设置坐标刻度对应数字的大小
    plt.show()
