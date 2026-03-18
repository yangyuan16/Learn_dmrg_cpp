import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pylab as plt

#
def plot_density_curve(r,density, label_y):
    fig = plt.figure(figsize=(10,3)) 
    ax1 = plt.subplot(1,1,1)
    plt.sca(ax1)  ##选择对ax1进行绘图
    ax1=plt.gca() #获得坐标轴的句柄
    ax1.plot(r,density,label=r"G",ls="-",lw=1.5,color="red",
             marker='o',alpha=1,markersize=8,markeredgewidth=1.5, markeredgecolor="k",
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
    plt.title("Lx=%d,dop=%.4f,tnnn=%.2f,Jnnn=%.2f,Dim=%d"%(Lx,dop/(Lz*Ly*Lx),t1,J1,dim), fontsize=25)
    plt.show()
    return
#
def density_along_x_Ly4(df,Ly,Lz): # 沿着 x 方向的 electron density
    df_out = pd.DataFrame(columns = ["r0","dy0","r1","dy1","r2","dy2","rmean","dymean"])
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
    density_mean = (np.array(df_out["dy0"].values) + np.array(df_out["dy1"].values) + 
                    np.array(df_out["dy2"].values) + np.array(df_out["dy3"].values)) / (Lz*Ly)
    #
    df_out["rmean"] = r_mean
    df_out["dymean"] = density_mean
    #
    return df_out
#
def density_along_x_Ly3(df,Ly,Lz): # 沿着 x 方向的 electron density
    df_out = pd.DataFrame(columns = ["r0","dy0","r1","dy1","r2","dy2","rmean","dymean"])
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

    r_mean = range(len(df_y0))
    density_mean = (np.array(df_out["dy0"].values) + np.array(df_out["dy1"].values) + np.array(df_out["dy2"].values)) / (Lz*Ly)
    #
    df_out["rmean"] = r_mean
    df_out["dymean"] = density_mean
    #
    return df_out
#
def density_along_x_Ly2(df,Ly,Lz): # 沿着 x 方向的 electron density
    df_out = pd.DataFrame(columns = ["r0","dy0","r1","dy1","rmean","dymean"])
    #
    df_y0 = df[df['site'] % (Ly * Lz) ==0]
    df_out["r0"] = df_y0["site"].values
    df_out["dy0"] = df_y0["density"].values
    
    df_y1 = df[df['site'] % (Ly * Lz) ==1]
    df_out["r1"] = df_y1["site"].values
    df_out["dy1"] = df_y1["density"].values

    r_mean = range(len(df_y0))
    density_mean = (np.array(df_out["dy0"].values) + np.array(df_out["dy1"].values)) / (Lz*Ly)
    #
    df_out["rmean"] = r_mean
    df_out["dymean"] = density_mean
    #
    return df_out
#
if __name__ == "__main__":
    Lz = 1
    Ly = 4
    Lx = 48
    dop = 60
    t = 3
    J = 1
    t1 = -3
    J1 = 1.0
    dim = 6000 # dim cutoff
    workpath = "E:\\WORK\\Work\\Project\\SSLattice\\leg4"
    filepath1 = "\\data_dmrgcpp\\Lz%d_Ly%d_Lx%d\\dop%g" % (Lz, Ly, Lx,dop)
    filepath2 = "\\t%.2f_J%.2f_t%.2f_J%.2f_dim%d" % (t, J, t1, J1, dim)
    filepath3 = "\\measurement_electron_density.dat"
    filename = workpath + filepath1 + filepath2 + filepath3
    print(filename)
    # load the data
    df = pd.read_csv(filename, header=None, sep='\t',encoding='utf-8')
    df.rename(columns={0: "site", 1: "density"},inplace=True)
    df.sort_values(['site'],inplace=True)
    df_layer1 = df
    sites_layer1 = df_layer1["site"].values.reshape(-1,Ly).T
    density_layer1 = df_layer1["density"].values.reshape(-1, Ly).T
    print(df.head())
    print(len(df))
    #------------------------------------------------------------
    #
    if Ly == 2:
        df_x = density_along_x_Ly2(df=df, Ly=Ly, Lz=Lz)
    elif Ly == 3:
        df_x = density_along_x_Ly3(df=df, Ly=Ly, Lz=Lz)
    elif Ly == 4:
        df_x = density_along_x_Ly4(df=df, Ly=Ly, Lz=Lz)
    else:
        raise "wrong of Ly"
    #----plot density curve----------------
    if Ly == 2:
        plot_density_curve(r=df_x["r0"],density=df_x["dy0"],label_y='<ni>(y0)')
        plot_density_curve(r=df_x["r1"],density=df_x["dy1"],label_y='<ni>(y1)')
        #
        plot_density_curve(r=df_x["rmean"],density=df_x["dymean"],label_y='<ni>(mean)')
    elif Ly == 3:
        plot_density_curve(r=df_x["r0"],density=df_x["dy0"],label_y='<ni>(y0)')
        plot_density_curve(r=df_x["r1"],density=df_x["dy1"],label_y='<ni>(y1)')
        plot_density_curve(r=df_x["r2"],density=df_x["dy2"],label_y='<ni>(y2)')
        #
        plot_density_curve(r=df_x["rmean"],density=df_x["dymean"],label_y='<ni>(mean)')
    elif Ly == 4:
        plot_density_curve(r=df_x["r0"],density=df_x["dy0"],label_y='<ni>(y0)')
        plot_density_curve(r=df_x["r1"],density=df_x["dy1"],label_y='<ni>(y1)')
        plot_density_curve(r=df_x["r2"],density=df_x["dy2"],label_y='<ni>(y2)')
        plot_density_curve(r=df_x["r3"],density=df_x["dy3"],label_y='<ni>(y3)')
        #
        plot_density_curve(r=df_x["rmean"],density=df_x["dymean"],label_y='<ni>(mean)')
    else:
        raise Exception("wrong input of Ly")