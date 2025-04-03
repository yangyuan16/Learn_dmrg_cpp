import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pylab as plt
#
def density_along_x_Ly3(df,Ly,Lz): # 沿着 x 方向的 electron density
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
    return df_out
#
def density_along_x_Ly2(df,Ly,Lz): # 沿着 x 方向的 electron density
    df_out = pd.DataFrame(columns = ["r0","dy0","r1","dy1","r2","dy2","r3","dy3","rmean","dymean"])
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

    r_mean = range(len(df_y0))
    density_mean = (np.array(df_out["dy0"].values) + np.array(df_out["dy1"].values) + np.array(df_out["dy2"].values)  + 
                    np.array(df_out["dy3"].values) ) / (Lz*Ly)
    #
    df_out["rmean"] = r_mean
    df_out["dymean"] = density_mean
    #
    return df_out
#
def get_density_amp(Lz, Ly, Lxdim,delta,Lmidbond):
    density_amp_list = []
    Lx_list = []
    for it in range(len(Lxdim)):
        Lx = Lxdim[it][0]
        dim = Lxdim[it][1]
        midbond = Lmidbond[it]
        dop = Lz * Ly * Lx * delta
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
        #print(df.head())
        #print(len(df))
        #------------------------------------------------------------
        #plot_loopJz_density(Jz=Jz)
        #plot_loopJz_density_curve(Jz=Jz)
        df_x = density_along_x_Ly2(df=df, Ly=Ly, Lz=Lz)
        dy0_amp = np.abs(df_x["dy0"][midbond[0]] - df_x["dy0"][midbond[1]])
        dy1_amp = np.abs(df_x["dy1"][midbond[0]] - df_x["dy1"][midbond[1]])
        dy2_amp = np.abs(df_x["dy2"][midbond[0]] - df_x["dy2"][midbond[1]])
        dy3_amp = np.abs(df_x["dy3"][midbond[0]] - df_x["dy3"][midbond[1]])
        density_amp = (dy0_amp+dy1_amp+dy2_amp+dy3_amp)/4
        density_amp_list.append(density_amp)
        Lx_list.append(Lx)
    return Lx_list, density_amp_list
#
if __name__ == "__main__":
    Lz = 2
    Ly = 2
    t = 3
    J = 1
    Jz = 0.1
    #dim = 6000 # dim cutoff
    #
    Lxdim_01875 = [(48,6000),(56,6000),(64,6000),(72,6000),(80,6000),(88,6000),(96,6000),(104,6000),(112,6000),(120,6000)]
    #midbond_01875=[(21,23), (25,27), (29,31), (33,35), (37,39), (41,43), (45,47), (49,51), (53,55), (58,60)]
    midbond_01875=[(21,23), (27,30), (29,31), (35,38), (37,39), (43,46), (45,47), (49,51), (53,55), (58,60)]
    Lxdim_025 = [(48,6000),(56,6000),(64,6000),(72,6000),(80,6000),(88,6000),(96,6000),(104,6000),(112,6000),(120,6000)]
    midbond_025=[(22,24), (26,28), (30,32), (34,36), (38,40), (42,44), (46,48), (50,52), (54,56), (58,60)]
    Lxdim_0375 = [(48,6000),(56,6000),(64,6000),(72,6000),(80,6000),(88,6000),(96,6000),(104,6000),(112,6000),(120,6000)]
    midbond_0375=[(21,22), (25,26), (29,30), (33,34), (37,38), (45,46), (45,46), (49,50), (53,54), (57,58)]
    # get density of delta 
    Lx_list_01875, density_01875 = get_density_amp(Lz=Lz, Ly=Ly, Lxdim=Lxdim_01875, delta=0.1875, Lmidbond=midbond_01875)
    Lx_list_025, density_025 = get_density_amp(Lz=Lz, Ly=Ly, Lxdim=Lxdim_025, delta=0.25, Lmidbond=midbond_025)
    Lx_list_0375, density_0375 = get_density_amp(Lz=Lz, Ly=Ly, Lxdim=Lxdim_0375, delta=0.375, Lmidbond=midbond_0375)
    #
    fig = plt.figure(figsize=(6,6)) 
    ax1 = plt.subplot(1,1,1)
    plt.sca(ax1)  ##选择对ax1进行绘图
    ax1=plt.gca() #获得坐标轴的句柄
    L01875, = ax1.plot(Lx_list_01875,density_01875,label=r"$\delta$=%.4f"%(0.1875),ls="-",lw=1.5,color="red", 
                       marker='o',alpha=1,markersize=12,markeredgewidth=1.5, markeredgecolor="red",
                       markerfacecolor='None')
    
    L025, = ax1.plot(Lx_list_025,density_025,label=r"$\delta$=%.4f"%(0.25),ls="-",lw=1.5,color="blue", 
                       marker='o',alpha=1,markersize=12,markeredgewidth=1.5, markeredgecolor="blue",
                       markerfacecolor='None')
    
    L0375, = ax1.plot(Lx_list_0375,density_0375,label=r"$\delta$=%.4f"%(0.375),ls="-",lw=1.5,color="green", 
                       marker='o',alpha=1,markersize=12,markeredgewidth=1.5, markeredgecolor="green",
                       markerfacecolor='None')
    #
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 22, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[L01875,L025,L0375], loc = 4, bbox_to_anchor=(0.78, 0.08),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    #
    label_x = r"Lx"
    label_y = r"$\Delta$ n"
    #plt.yscale("log")
    #plt.xscale("log")
    ax1.set_xlabel(label_x, size= 25)
    ax1.set_ylabel(label_y, size= 25)
    ax1.tick_params(labelsize = 25) # 设置坐标刻度对应数字的大小
    #ax1.set_xlim([0,8])
    #ax1.set_ylim([-0.1,1])
    #ax1.set_xticks([0,2,4,6,8])
    #ax1.set_yticks([-0.1,0,0.5,1])
    #ax1.text(0.3,-0.05, r'$\mathrm{(a)}$', fontsize=18)
    plt.title("Jz=%.2f"%(Jz), fontsize=25)
    plt.show()
    
