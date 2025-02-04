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
    Lx = 48
    dop = 96 
    t = 3
    J = 1
    dim = 6000
    #============ data of Jz = 0.1
    Jz = 0.1 # 
    workpath = "E:\\WORK\\Work\\Project\\La3Ni2O7"
    filepath1 = "\\data_dmrgcpp\\Lz%d_Ly%d_Lx%d\\dop%g" % (Lz, Ly, Lx,dop)
    filepath2 = "\\t%d_J%d_Jz%.2f_dim%d" % (t, J, Jz, dim)
    filepath3 = "\\measurement_electron_density.dat"
    filename = workpath + filepath1 + filepath2 + filepath3
    print(filename)
    # load the data
    df_01 = pd.read_csv(filename, header=None, sep='\t',encoding='utf-8')
    df_01.rename(columns={0: "site", 1: "density"},inplace=True)
    df_01.sort_values(['site'],inplace=True)
    df_x_01 = density_along_x_Ly2(df=df_01, Ly=Ly, Lz=Lz)
    #=========== data of Jz = 0.2
    Jz = 0.2 # 
    filepath2 = "\\t%d_J%d_Jz%.2f_dim%d" % (t, J, Jz, dim)
    filename = workpath + filepath1 + filepath2 + filepath3
    print(filename)
    # load the data
    df_02 = pd.read_csv(filename, header=None, sep='\t',encoding='utf-8')
    df_02.rename(columns={0: "site", 1: "density"},inplace=True)
    df_02.sort_values(['site'],inplace=True)
    df_x_02 = density_along_x_Ly2(df=df_02, Ly=Ly, Lz=Lz)
    #=========== data of Jz = 0.4
    Jz = 0.4 # 
    filepath2 = "\\t%d_J%d_Jz%.2f_dim%d" % (t, J, Jz, dim)
    filename = workpath + filepath1 + filepath2 + filepath3
    print(filename)
    # load the data
    df_04 = pd.read_csv(filename, header=None, sep='\t',encoding='utf-8')
    df_04.rename(columns={0: "site", 1: "density"},inplace=True)
    df_04.sort_values(['site'],inplace=True)
    df_x_04 = density_along_x_Ly2(df=df_04, Ly=Ly, Lz=Lz)
    #=========== data of Jz = 0.6
    Jz = 0.6 # 
    filepath2 = "\\t%d_J%d_Jz%.2f_dim%d" % (t, J, Jz, dim)
    filename = workpath + filepath1 + filepath2 + filepath3
    print(filename)
    # load the data
    df_06 = pd.read_csv(filename, header=None, sep='\t',encoding='utf-8')
    df_06.rename(columns={0: "site", 1: "density"},inplace=True)
    df_06.sort_values(['site'],inplace=True)
    df_x_06 = density_along_x_Ly2(df=df_06, Ly=Ly, Lz=Lz)
    #=========== data of Jz = 0.8
    Jz = 0.8 # 
    filepath2 = "\\t%d_J%d_Jz%.2f_dim%d" % (t, J, Jz, dim)
    filename = workpath + filepath1 + filepath2 + filepath3
    print(filename)
    # load the data
    df_08 = pd.read_csv(filename, header=None, sep='\t',encoding='utf-8')
    df_08.rename(columns={0: "site", 1: "density"},inplace=True)
    df_08.sort_values(['site'],inplace=True)
    df_x_08 = density_along_x_Ly2(df=df_08, Ly=Ly, Lz=Lz)
    #=========== data of Jz = 0.8
    '''
    Jz = 0.75 # 
    filepath2 = "\\t%d_J%d_Jz%.2f_dim%d" % (t, J, Jz, dim)
    filename = workpath + filepath1 + filepath2 + filepath3
    print(filename)
    # load the data
    df_075 = pd.read_csv(filename, header=None, sep='\t',encoding='utf-8')
    df_075.rename(columns={0: "site", 1: "density"},inplace=True)
    df_075.sort_values(['site'],inplace=True)
    df_x_075 = density_along_x_Ly2(df=df_075, Ly=Ly, Lz=Lz)
    '''
    #=========== data of Jz = 1.0
    Jz = 1.0 # Jz
    filepath2 = "\\t%d_J%d_Jz%.2f_dim%d" % (t, J, Jz, dim)
    filename = workpath + filepath1 + filepath2 + filepath3
    print(filename)
    # load the data
    df_10 = pd.read_csv(filename, header=None, sep='\t',encoding='utf-8')
    df_10.rename(columns={0: "site", 1: "density"},inplace=True)
    df_10.sort_values(['site'],inplace=True)
    df_x_10 = density_along_x_Ly2(df=df_10, Ly=Ly, Lz=Lz)
    #=========== data of Jz = 1.2
    '''
    Jz = 1.25 # Jz
    filepath2 = "\\t%d_J%d_Jz%.2f_dim%d" % (t, J, Jz, dim)
    filename = workpath + filepath1 + filepath2 + filepath3
    print(filename)
    # load the data
    df_125 = pd.read_csv(filename, header=None, sep='\t',encoding='utf-8')
    df_125.rename(columns={0: "site", 1: "density"},inplace=True)
    df_125.sort_values(['site'],inplace=True)
    df_x_125 = density_along_x_Ly2(df=df_125, Ly=Ly, Lz=Lz)
    #=========== data of Jz = 1.5
    Jz = 1.5 # Jz
    filepath2 = "\\t%d_J%d_Jz%.2f_dim%d" % (t, J, Jz, dim)
    filename = workpath + filepath1 + filepath2 + filepath3
    print(filename)
    # load the data
    df_15 = pd.read_csv(filename, header=None, sep='\t',encoding='utf-8')
    df_15.rename(columns={0: "site", 1: "density"},inplace=True)
    df_15.sort_values(['site'],inplace=True)
    df_x_15 = density_along_x_Ly2(df=df_15, Ly=Ly, Lz=Lz)
    #=========== data of Jz = 1.75
    Jz = 1.75 # Jz
    filepath2 = "\\t%d_J%d_Jz%.2f_dim%d" % (t, J, Jz, dim)
    filename = workpath + filepath1 + filepath2 + filepath3
    print(filename)
    # load the data
    df_175 = pd.read_csv(filename, header=None, sep='\t',encoding='utf-8')
    df_175.rename(columns={0: "site", 1: "density"},inplace=True)
    df_175.sort_values(['site'],inplace=True)
    df_x_175 = density_along_x_Ly2(df=df_175, Ly=Ly, Lz=Lz)
    '''
    #=========== data of Jz = 2
    Jz = 2.0 # Jz
    filepath2 = "\\t%d_J%d_Jz%.2f_dim%d" % (t, J, Jz, dim)
    filename = workpath + filepath1 + filepath2 + filepath3
    print(filename)
    # load the data
    df_20 = pd.read_csv(filename, header=None, sep='\t',encoding='utf-8')
    df_20.rename(columns={0: "site", 1: "density"},inplace=True)
    df_20.sort_values(['site'],inplace=True)
    df_x_20 = density_along_x_Ly2(df=df_20, Ly=Ly, Lz=Lz)
    #=========== data of Jz = 2.5
    '''
    Jz = 2.5 # Jz
    filepath2 = "\\t%d_J%d_Jz%.2f_dim%d" % (t, J, Jz, dim)
    filename = workpath + filepath1 + filepath2 + filepath3
    print(filename)
    # load the data
    df_25 = pd.read_csv(filename, header=None, sep='\t',encoding='utf-8')
    df_25.rename(columns={0: "site", 1: "density"},inplace=True)
    df_25.sort_values(['site'],inplace=True)
    df_x_25 = density_along_x_Ly2(df=df_25, Ly=Ly, Lz=Lz)
    '''
    #=========== data of Jz = 3.0
    Jz = 3.0 # Jz
    filepath2 = "\\t%d_J%d_Jz%.2f_dim%d" % (t, J, Jz, dim)
    filename = workpath + filepath1 + filepath2 + filepath3
    print(filename)
    # load the data
    df_30 = pd.read_csv(filename, header=None, sep='\t',encoding='utf-8')
    df_30.rename(columns={0: "site", 1: "density"},inplace=True)
    df_30.sort_values(['site'],inplace=True)
    df_x_30 = density_along_x_Ly2(df=df_30, Ly=Ly, Lz=Lz)
    #=========== data of Jz = 3.5
    '''
    Jz = 3.5 # Jz
    filepath2 = "\\t%d_J%d_Jz%.2f_dim%d" % (t, J, Jz, dim)
    filename = workpath + filepath1 + filepath2 + filepath3
    print(filename)
    # load the data
    df_35 = pd.read_csv(filename, header=None, sep='\t',encoding='utf-8')
    df_35.rename(columns={0: "site", 1: "density"},inplace=True)
    df_35.sort_values(['site'],inplace=True)
    df_x_35 = density_along_x_Ly2(df=df_35, Ly=Ly, Lz=Lz)
    '''
    #========== plot ==============
    fig = plt.figure(figsize=(10,6)) 
    ax1 = plt.subplot(1,1,1)
    plt.sca(ax1)  ##选择对ax1进行绘图
    ax1=plt.gca() #获得坐标轴的句柄
    #L01, = ax1.plot(df_x_01["rmean"],df_x_01["dymean"],label=r"Jz=0.1",ls="-",lw=1.5,color="red",
    #         marker='o',alpha=1,markersize=12,markeredgewidth=1.5, markeredgecolor="red",
    #         markerfacecolor='None')
    L02, = ax1.plot(df_x_02["rmean"],df_x_02["dymean"],label=r"Jz=0.2",ls="-",lw=1.5,color="blue",
             marker='o',alpha=1,markersize=12,markeredgewidth=1.5, markeredgecolor="blue",
             markerfacecolor='None')
    L04, = ax1.plot(df_x_04["rmean"],df_x_04["dymean"],label=r"Jz=0.4",ls="-",lw=1.5,color="green",
             marker='o',alpha=1,markersize=12,markeredgewidth=1.5, markeredgecolor="green",
             markerfacecolor='None')
    L06, = ax1.plot(df_x_06["rmean"],df_x_06["dymean"],label=r"Jz=0.6",ls="-",lw=1.5,color="magenta",
             marker='o',alpha=1,markersize=12,markeredgewidth=1.5, markeredgecolor="magenta",
             markerfacecolor='None')
    L08, = ax1.plot(df_x_08["rmean"],df_x_08["dymean"],label=r"Jz=0.8",ls="-",lw=1.5,color="cyan",
             marker='o',alpha=1,markersize=12,markeredgewidth=1.5, markeredgecolor="cyan",
             markerfacecolor='None')
    #L075, = ax1.plot(df_x_075["rmean"],df_x_075["dymean"],label=r"Jz=0.75",ls="-",lw=1.5,color="cyan",
    #         marker='o',alpha=1,markersize=12,markeredgewidth=1.5, markeredgecolor="cyan",
    #         markerfacecolor='None')
    L10, = ax1.plot(df_x_10["rmean"],df_x_10["dymean"],label=r"Jz=1.0",ls="-",lw=1.5,color="brown",
             marker='o',alpha=1,markersize=12,markeredgewidth=1.5, markeredgecolor="brown",
             markerfacecolor='None')
    #L125, = ax1.plot(df_x_125["rmean"],df_x_125["dymean"],label=r"Jz=1.25",ls="-",lw=1.5,color="k",
    #         marker='^',alpha=1,markersize=12,markeredgewidth=1.5, markeredgecolor="k",
    #         markerfacecolor='None')    

    #L15, = ax1.plot(df_x_15["rmean"],df_x_15["dymean"],label=r"Jz=1.5",ls="-",lw=1.5,color="red",
    #         marker='s',alpha=0.9,markersize=12,markeredgewidth=1.5, markeredgecolor="red",
    #         markerfacecolor='None')    
    
    #L175, = ax1.plot(df_x_175["rmean"],df_x_175["dymean"],label=r"Jz=1.75",ls="-",lw=1.5,color="blue",
    #         marker='s',alpha=0.9,markersize=12,markeredgewidth=1.5, markeredgecolor="blue",
    #         markerfacecolor='None')  

    L20, = ax1.plot(df_x_20["rmean"],df_x_20["dymean"],label=r"Jz=2.0",ls="-",lw=1.5,color="green",
             marker='s',alpha=0.9,markersize=12,markeredgewidth=1.5, markeredgecolor="green",
             markerfacecolor='None')

    #L25, = ax1.plot(df_x_25["rmean"],df_x_25["dymean"],label=r"Jz=2.5",ls="-",lw=1.5,color="magenta",
    #         marker='s',alpha=0.9,markersize=12,markeredgewidth=1.5, markeredgecolor="magenta",
    #         markerfacecolor='None')
    
    L30, = ax1.plot(df_x_30["rmean"],df_x_30["dymean"],label=r"Jz=3.0",ls="-",lw=1.5,color="cyan",
             marker='s',alpha=0.9,markersize=12,markeredgewidth=1.5, markeredgecolor="cyan",
             markerfacecolor='None')

    #L35, = ax1.plot(df_x_35["rmean"],df_x_35["dymean"],label=r"Jz=3.5",ls="-",lw=1.5,color="brown",
    #         marker='s',alpha=0.9,markersize=12,markeredgewidth=1.5, markeredgecolor="brown",
    #         markerfacecolor='None')
    #
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 22, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[L02,L04,L06,L08,L10,L20,L30], loc = 4, bbox_to_anchor=(0.78, 0.08),
                       ncol = 4,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False) #[L01,L025,L04,L06,L075,L10,L125,L15,L175,L20,L25,L30,L35]
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
    plt.title("dim=%d"%dim,fontsize=25)
    plt.show()