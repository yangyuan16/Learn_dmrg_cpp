import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pylab as plt
#
def read_electron_density():
    return
#
def plot_density_curve(r,density, label_y, Jz):
    fig = plt.figure(figsize=(10,3)) 
    ax1 = plt.subplot(1,1,1)
    plt.sca(ax1)  ##选择对ax1进行绘图
    ax1=plt.gca() #获得坐标轴的句柄
    text = "Lz%d_Ly%d_Lx%d_dop%g_" % (Lz, Ly, Lx,dop) + "t%d_J%d_Jz%.2f_dim%d" % (t, J, Jz, dim)
    ax1.text(0.3, 0.009, text)
    ax1.plot(r,density,label=r"G",ls="-",lw=1.5,color="red",
             marker='o',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="k",
             markerfacecolor='w')
    label_x = r"i"
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
    plt.title("Lx=%d,dop=%.4f,Jz=%.2f,Dim=%d"%(Lx,dop/(Lz*Ly*Lx),Jz,dim), fontsize=25)
    plt.show()
    return
#
def density_along_x_Lz3(df,Ly,Lz): # 沿着 x 方向的 electron density
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
def density_along_x_Lz2(df,Ly,Lz): # 沿着 x 方向的 electron density
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
    Lz = 3
    Ly = 2
    Lx = 48
    dop = 108
    t = 3.0
    J = 1.0
    Jz_list = [12.0]
    for it in range(len(Jz_list)):
        Jz = Jz_list[it]
        #============ data of dim = 6000
        dim = 6000 # dim cutoff
        workpath = "E:\\WORK\\Work\\Project\\TrilayerW2-tj-Jz-openz"
        filepath1 = "\\data_dmrgcpp\\Lz%d_Ly%d_Lx%d\\dop%g" % (Lz, Ly, Lx,dop)
        filepath2 = "\\t%.2f_J%.2f_Jz%.2f_dim%d" % (t, J, Jz, dim)
        filepath3 = "\\measurement_electron_density.dat"
        filename = workpath + filepath1 + filepath2 + filepath3
        print(filename)
        # load the data
        df_dim6000 = pd.read_csv(filename, header=None, sep='\t',encoding='utf-8')
        df_dim6000.rename(columns={0: "site", 1: "density"},inplace=True)
        df_dim6000.sort_values(['site'],inplace=True)
        df_x_dim6000 = density_along_x_Lz3(df=df_dim6000, Ly=Ly, Lz=Lz)
        #
        #=========== data of dim = 8000
        #dim = 8000 # dim cutoff
        #filepath2 = "\\t%d_J%d_Jz%.2f_dim%d" % (t, J, Jz, dim)
        #filename = workpath + filepath1 + filepath2 + filepath3
        #print(filename)
        # load the data
        #df_dim8000 = pd.read_csv(filename, header=None, sep='\t',encoding='utf-8')
        #df_dim8000.rename(columns={0: "site", 1: "density"},inplace=True)
        #df_dim8000.sort_values(['site'],inplace=True)
        #df_x_dim8000 = density_along_x_Ly2(df=df_dim8000, Ly=Ly, Lz=Lz)
        #
        #=========== data of dim = 10000
        dim = 10000 # dim cutoff
        filepath2 = "\\t%.2f_J%.2f_Jz%.2f_dim%d" % (t, J, Jz, dim)
        filename = workpath + filepath1 + filepath2 + filepath3
        print(filename)
        # load the data
        df_dim10000 = pd.read_csv(filename, header=None, sep='\t',encoding='utf-8')
        df_dim10000.rename(columns={0: "site", 1: "density"},inplace=True)
        df_dim10000.sort_values(['site'],inplace=True)
        df_x_dim10000 = density_along_x_Lz3(df=df_dim10000, Ly=Ly, Lz=Lz)
        #=========== data of dim = 12000
        dim = 12000 # dim cutoff
        filepath2 = "\\t%.2f_J%.2f_Jz%.2f_dim%d" % (t, J, Jz, dim)
        filename = workpath + filepath1 + filepath2 + filepath3
        print(filename)
        # load the data
        df_dim12000 = pd.read_csv(filename, header=None, sep='\t',encoding='utf-8')
        df_dim12000.rename(columns={0: "site", 1: "density"},inplace=True)
        df_dim12000.sort_values(['site'],inplace=True)
        df_x_dim12000 = density_along_x_Lz3(df=df_dim12000, Ly=Ly, Lz=Lz)
        #=========== data of dim = 14000
        dim = 14000 # dim cutoff
        filepath2 = "\\t%.2f_J%.2f_Jz%.2f_dim%d" % (t, J, Jz, dim)
        filename = workpath + filepath1 + filepath2 + filepath3
        print(filename)
        # load the data
        df_dim14000 = pd.read_csv(filename, header=None, sep='\t',encoding='utf-8')
        df_dim14000.rename(columns={0: "site", 1: "density"},inplace=True)
        df_dim14000.sort_values(['site'],inplace=True)
        df_x_dim14000 = density_along_x_Lz3(df=df_dim14000, Ly=Ly, Lz=Lz)
        #========== plot
        fig = plt.figure(figsize=(10,3))
        ax1 = plt.subplot(1,1,1)
        plt.sca(ax1)  ##选择对ax1进行绘图
        ax1=plt.gca() #获得坐标轴的句柄
        L6000, = ax1.plot(df_x_dim6000["rmean"],df_x_dim6000["dymean"],label=r"dim=6000",ls="-",lw=1.5,color="red",
             marker='o',alpha=1,markersize=10,markeredgewidth=1.5, markeredgecolor="red",
             markerfacecolor='None')
        #L8000, = ax1.plot(df_x_dim8000["rmean"],df_x_dim8000["dymean"],label=r"dim=8000",ls="-",lw=1.5,color="blue",
        #     marker='o',alpha=1,markersize=10,markeredgewidth=1.5, markeredgecolor="blue",
        #     markerfacecolor='None')
        L10000, = ax1.plot(df_x_dim10000["rmean"],df_x_dim10000["dymean"],label=r"dim=10000",ls="-",lw=1.5,color="blue",
             marker='o',alpha=1,markersize=10,markeredgewidth=1.5, markeredgecolor="blue",
             markerfacecolor='None')
        
        L12000, = ax1.plot(df_x_dim12000["rmean"],df_x_dim12000["dymean"],label=r"dim=12000",ls="-",lw=1.5,color="magenta",
             marker='o',alpha=1,markersize=10,markeredgewidth=1.5, markeredgecolor="magenta",
             markerfacecolor='None')
        
        L14000, = ax1.plot(df_x_dim14000["rmean"],df_x_dim14000["dymean"],label=r"dim=14000",ls="-",lw=1.5,color="cyan",
             marker='o',alpha=1,markersize=10,markeredgewidth=1.5, markeredgecolor="cyan",
             markerfacecolor='None')

        legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 22, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
        legend1=plt.legend(handles=[L6000,L10000,L12000,L14000], loc = 4, bbox_to_anchor=(0.58, 0.48),
                       ncol = 2,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
        label_x = r"i"
        label_y = "n(mean)"
        ax1.set_xlabel(label_x, size= 25)
        ax1.set_ylabel(label_y, size= 25)
        #
        ax1.tick_params(labelsize = 25) # 设置坐标刻度对应数字的大小
        #ax1.set_xlim([0,8])
        #ax1.set_ylim([-0.1,1])
        #ax1.set_xticks([0,2,4,6,8])
        #ax1.set_yticks([-0.1,0,0.5,1])
        #ax1.text(0.3,-0.05, r'$\mathrm{(a)}$', fontsize=18)
        plt.title("Lx=%d,dop=%.4f,Jz=%.2f"%(Lx,dop/(Lz*Ly*Lx),Jz,), fontsize=25)
        plt.show()

