import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pylab as plt

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
    ax1.set_xlabel(label_x, size= 14)
    ax1.set_ylabel(label_y, size= 14)
    #
    ax1.tick_params(labelsize = 15) # 设置坐标刻度对应数字的大小
    #ax1.set_xlim([0,8])
    #ax1.set_ylim([-0.1,1])
    #ax1.set_xticks([0,2,4,6,8])
    #ax1.set_yticks([-0.1,0,0.5,1])
    #ax1.text(0.3,-0.05, r'$\mathrm{(a)}$', fontsize=18)
    plt.title("Jz=%.2f"%Jz)
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
    t = 3
    J = 1
    Jz_list = [3.0]
    dim = 6000 # dim cutoff
    for it in range(len(Jz_list)):
        Jz = Jz_list[it]
        #============ data of Lx = 48
        Lx = 48
        dop = 88 
        workpath = "E:\\WORK\\Work\\Project\\La3Ni2O7"
        filepath1 = "\\data_dmrgcpp\\Lz%d_Ly%d_Lx%d\\dop%g" % (Lz, Ly, Lx,dop)
        filepath2 = "\\t%d_J%d_Jz%.2f_dim%d" % (t, J, Jz, dim)
        filepath3 = "\\measurement_electron_density.dat"
        filename = workpath + filepath1 + filepath2 + filepath3
        print(filename)
        # load the data
        df_Lx48 = pd.read_csv(filename, header=None, sep='\t',encoding='utf-8')
        df_Lx48.rename(columns={0: "site", 1: "density"},inplace=True)
        df_Lx48.sort_values(['site'],inplace=True)
        df_x_Lx48 = density_along_x_Ly2(df=df_Lx48, Ly=Ly, Lz=Lz)
        #plot_density_curve(r=df_x_Lx48["rmean"],density=df_x_Lx48["dymean"],label_y='<ni>(mean)',Jz=Jz)
        #=========== data of Lx = 64
        '''
        Lx = 64
        dop = 48 
        filepath1 = "\\data_dmrgcpp\\Lz%d_Ly%d_Lx%d\\dop%g" % (Lz, Ly, Lx,dop)
        filename = workpath + filepath1 + filepath2 + filepath3
        print(filename)
        # load the data
        df_Lx64 = pd.read_csv(filename, header=None, sep='\t',encoding='utf-8')
        df_Lx64.rename(columns={0: "site", 1: "density"},inplace=True)
        df_Lx64.sort_values(['site'],inplace=True)
        df_x_Lx64 = density_along_x_Ly2(df=df_Lx64, Ly=Ly, Lz=Lz)
        #plot_density_curve(r=df_x_Lx64["rmean"],density=df_x_Lx64["dymean"],label_y='<ni>(mean)',Jz=Jz)
        #=========== data of Lx = 128
        Lx = 128
        dop = 96 
        filepath1 = "\\data_dmrgcpp\\Lz%d_Ly%d_Lx%d\\dop%g" % (Lz, Ly, Lx,dop)
        filename = workpath + filepath1 + filepath2 + filepath3
        print(filename)
        # load the data
        df_Lx128 = pd.read_csv(filename, header=None, sep='\t',encoding='utf-8')
        df_Lx128.rename(columns={0: "site", 1: "density"},inplace=True)
        df_Lx128.sort_values(['site'],inplace=True)
        df_x_Lx128 = density_along_x_Ly2(df=df_Lx128, Ly=Ly, Lz=Lz)
        #plot_density_curve(r=df_x_Lx64["rmean"],density=df_x_Lx64["dymean"],label_y='<ni>(mean)',Jz=Jz)
        '''
        #=========== data of Lx = 120
        Lx = 120
        dop = 220
        filepath1 = "\\data_dmrgcpp\\Lz%d_Ly%d_Lx%d\\dop%g" % (Lz, Ly, Lx,dop)
        filename = workpath + filepath1 + filepath2 + filepath3
        print(filename)
        # load the data
        df_Lx120 = pd.read_csv(filename, header=None, sep='\t',encoding='utf-8')
        df_Lx120.rename(columns={0: "site", 1: "density"},inplace=True)
        df_Lx120.sort_values(['site'],inplace=True)
        df_x_Lx120 = density_along_x_Ly2(df=df_Lx120, Ly=Ly, Lz=Lz)
        #plot_density_curve(r=df_x_Lx64["rmean"],density=df_x_Lx64["dymean"],label_y='<ni>(mean)',Jz=Jz)


        #========== plot
        fig = plt.figure(figsize=(10,3)) 
        ax1 = plt.subplot(1,1,1)
        plt.sca(ax1)  ##选择对ax1进行绘图
        ax1=plt.gca() #获得坐标轴的句柄
        #L48, = ax1.plot(df_x_Lx48["rmean"] + 40,df_x_Lx48["dymean"],label=r"Lx=48",ls="-",lw=1.5,color="red",
        #     marker='o',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="red",
        #     markerfacecolor='None')
        #L64, = ax1.plot(df_x_Lx64["rmean"]+32,df_x_Lx64["dymean"],label=r"Lx=64",ls="-",lw=1.5,color="blue",
        #     marker='o',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="blue",
        #     markerfacecolor='None')
        #L128, = ax1.plot(df_x_Lx128["rmean"],df_x_Lx128["dymean"],label=r"Lx=128",ls="-",lw=1.5,color="green",
        #     marker='o',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="green",
        #     markerfacecolor='None')
        
        #L48, = ax1.plot(df_x_Lx48["r0"] * (8/3),df_x_Lx48["dy0"],label=r"Lx=48",ls="-",lw=1.5,color="red",
        #     marker='o',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="red",
        #     markerfacecolor='None')
        #L64, = ax1.plot(df_x_Lx64["r0"] * 2,df_x_Lx64["dy0"],label=r"Lx=64",ls="-",lw=1.5,color="blue",
        #     marker='o',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="blue",
        #     markerfacecolor='None')
        #L128, = ax1.plot(df_x_Lx128["r0"],df_x_Lx128["dy0"],label=r"Lx=128",ls="-",lw=1.5,color="green",
        #     marker='o',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="green",
        #     markerfacecolor='None')

        L48, = ax1.plot(df_x_Lx48["rmean"] + 36,df_x_Lx48["dymean"],label=r"Lx=48",ls="-",lw=1.5,color="red",
             marker='o',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="red",
             markerfacecolor='None')
        L120, = ax1.plot(df_x_Lx120["rmean"],df_x_Lx120["dymean"],label=r"Lx=128",ls="-",lw=1.5,color="green",
             marker='o',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="green",
             markerfacecolor='None')        


        legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 18, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
        legend1=plt.legend(handles=[L48,L120], loc = 4, bbox_to_anchor=(0.78, 0.08),
                       ncol = 3,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
        label_x = r"i"
        label_y = "n(i)"
        ax1.set_xlabel(label_x, size= 25)
        ax1.set_ylabel(label_y, size= 25)
        #
        ax1.tick_params(labelsize = 25) # 设置坐标刻度对应数字的大小
        #ax1.set_xlim([0,8])
        #ax1.set_ylim([-0.1,1])
        #ax1.set_xticks([0,2,4,6,8])
        #ax1.set_yticks([-0.1,0,0.5,1])
        #ax1.text(0.3,-0.05, r'$\mathrm{(a)}$', fontsize=18)
        plt.title("Jz=%.2f"%(Jz),fontsize=25)
        plt.show()