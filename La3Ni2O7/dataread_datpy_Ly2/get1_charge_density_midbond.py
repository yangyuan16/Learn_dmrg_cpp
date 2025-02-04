import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pylab as plt
#
def read_electron_density():
    return
#
def plot_density_curve(r,density, label_y, Jz):
    fig = plt.figure(figsize=(6,6)) 
    ax1 = plt.subplot(1,1,1)
    plt.sca(ax1)  ##选择对ax1进行绘图
    ax1=plt.gca() #获得坐标轴的句柄
    #ax1.text(0.3, 0.009, text)
    ax1.plot(r,density,label=r"G",ls="-",lw=1.5,color="red",
             marker='o',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="k",
             markerfacecolor='w')
    label_x = r"Lx"
    label_y = label_y
    ax1.set_xlabel(label_x, size= 25)
    ax1.set_ylabel(label_y, size= 25)
    #
    ax1.tick_params(labelsize = 25) # 设置坐标刻度对应数字的大小
    #ax1.set_xlim([0,8])
    #ax1.set_ylim([-0.1,1])
    #ax1.set_xticks([0,2,4,6,8])
    #ax1.set_yticks([-0.1,0,0.5,1])
    #ax1.text(0.3,-0.05, r'$\mathrm{(a)}$', fontsize=18)
    plt.title("dop = %.4f, Jz=%.2f, dim=%d, Lcut=%d"%(delta, Jz,dim, Lcut), fontsize=25)
    plt.show()
    return
#
def plot_density_curve_log(r,density, label_y, Jz):
    fig = plt.figure(figsize=(6,6)) 
    ax1 = plt.subplot(1,1,1)
    plt.sca(ax1)  ##选择对ax1进行绘图
    ax1=plt.gca() #获得坐标轴的句柄
    #ax1.text(0.3, 0.009, text)
    ax1.plot(r,density,label=r"G",ls="-",lw=1.5,color="red",
             marker='o',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="k",
             markerfacecolor='w')
    label_x = r"Lx"
    label_y = label_y
    plt.yscale("log")
    ax1.set_xlabel(label_x, size= 25)
    ax1.set_ylabel(label_y, size= 25)
    #
    ax1.tick_params(labelsize = 25) # 设置坐标刻度对应数字的大小
    #ax1.set_xlim([0,8])
    #ax1.set_ylim([-0.1,1])
    #ax1.set_xticks([0,2,4,6,8])
    #ax1.set_yticks([-0.1,0,0.5,1])
    #ax1.text(0.3,-0.05, r'$\mathrm{(a)}$', fontsize=18)
    plt.title("dop = %.4f, Jz=%.2f, dim=%d, Lcut=%d"%(delta, Jz,dim, Lcut), fontsize=25)
    plt.show()
    return
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
if __name__ == "__main__":
    Lz = 2
    Ly = 2
    Lx_list = [48,56,64,72,80,88,96,104,112,120]
    delta = 0.375
    Lcut = 6 # 从中间点向前向后取 Lcut 个点 
    t = 3
    J = 1
    Jz = 0.1
    dim = 6000 # dim cutoff
    density_amp_list = []
    for it in range(len(Lx_list)):
        Lx = Lx_list[it]
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
        dy0_amp = df_x["dy0"][round(Lx/2-Lcut):round(Lx/2+Lcut)].max() - df_x["dy0"][round(Lx/2-Lcut):round(Lx/2+Lcut)].min()
        dy1_amp = df_x["dy1"][round(Lx/2-Lcut):round(Lx/2+Lcut)].max() - df_x["dy1"][round(Lx/2-Lcut):round(Lx/2+Lcut)].min()
        dy2_amp = df_x["dy2"][round(Lx/2-Lcut):round(Lx/2+Lcut)].max() - df_x["dy2"][round(Lx/2-Lcut):round(Lx/2+Lcut)].min()
        dy3_amp = df_x["dy3"][round(Lx/2-Lcut):round(Lx/2+Lcut)].max() - df_x["dy3"][round(Lx/2-Lcut):round(Lx/2+Lcut)].min()
        density_amp = (dy0_amp+dy1_amp+dy2_amp+dy3_amp)/4
        density_amp_list.append(density_amp)
    print(density_amp_list)
    plot_density_curve(r=Lx_list,density=density_amp_list,label_y=r'$\Delta$ n',Jz=Jz)
    plot_density_curve_log(r=Lx_list,density=density_amp_list,label_y=r'$\Delta$ n',Jz=Jz)




