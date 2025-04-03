import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pylab as plt
#
def read_electron_density():
    return
#
def density_fit(x,n0,A0,Kc,kf,phi,Lx):
    nx = []
    for it in x:
        Acdw = A0*(pow(it,-Kc/2)+pow((Lx + 1 - it),-Kc/2))
        ni = n0 + Acdw * np.cos(2*kf*it + phi)
        nx.append(ni)    
    return nx
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
    plt.title("Jz=%.2f"%Jz, fontsize=25)
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
def plot_loopJz_density_curve(Jz):
    if Ly == 2:
        df_x = density_along_x_Ly2(df=df, Ly=Ly, Lz=Lz)
    elif Ly == 3:
        df_x = density_along_x_Ly3(df=df, Ly=Ly, Lz=Lz)
    else:
        raise "wrong of Ly"

    #----plot density curve----------------
    if Ly == 2:
        plot_density_curve(r=df_x["r0"],density=df_x["dy0"],label_y='<ni>(y0)',Jz=Jz)
        plot_density_curve(r=df_x["r1"],density=df_x["dy1"],label_y='<ni>(y1)',Jz=Jz)
        plot_density_curve(r=df_x["r2"],density=df_x["dy2"],label_y='<ni>(y2)',Jz=Jz)
        plot_density_curve(r=df_x["r3"],density=df_x["dy3"],label_y='<ni>(y3)',Jz=Jz)  
        #
        plot_density_curve(r=df_x["rmean"],density=df_x["dymean"],label_y='<ni>(mean)',Jz=Jz)
    elif Ly == 3:
        #plot_density_curve(r=df_x["r0"],density=df_x["dy0"],label_y='<ni>(y0)',Jz=Jz)
        #plot_density_curve(r=df_x["r1"],density=df_x["dy1"],label_y='<ni>(y1)',Jz=Jz)
        #plot_density_curve(r=df_x["r2"],density=df_x["dy2"],label_y='<ni>(y2)',Jz=Jz)
        #plot_density_curve(r=df_x["r3"],density=df_x["dy3"],label_y='<ni>(y3)',Jz=Jz)
        #plot_density_curve(r=df_x["r4"],density=df_x["dy4"],label_y='<ni>(y4)',Jz=Jz)
        #plot_density_curve(r=df_x["r5"],density=df_x["dy5"],label_y='<ni>(y5)',Jz=Jz)
        #
        plot_density_curve(r=df_x["rmean"],density=df_x["dymean"],label_y='<ni>(mean)',Jz=Jz)
    else:
        raise Exception("wrong input of Ly")
    return
#
def plot_density_fit(r,density, label_y,x,nx,title):
    fig = plt.figure(figsize=(10,3)) 
    ax1 = plt.subplot(1,1,1)
    plt.sca(ax1)  ##选择对ax1进行绘图
    ax1=plt.gca() #获得坐标轴的句柄
    ax1.plot(r,density,label=r"G",ls="-",lw=1.5,color="red",marker='o',alpha=1,markersize=8,markeredgewidth=1.5, markeredgecolor="red",
             markerfacecolor='None')
    ax1.plot(x,nx,label=r"G",ls="--",lw=1.5,color="k",marker='o',alpha=1,markersize=0,markeredgewidth=1.5, markeredgecolor="k",
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
    plt.title(title,fontsize=18)
    plt.show()
    return
#
if __name__ == "__main__":
    Lz = 2
    Ly = 2
    Lx = 64
    dop = 48 
    t = 3
    J = 1
    Jz_list = [0.1,]
    dim = 6000 # dim cutoff
    for it in range(len(Jz_list)):
        Jz = Jz_list[it]
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
        #------------------------------------------------------------
        if Ly == 2:
            df_x = density_along_x_Ly2(df=df, Ly=Ly, Lz=Lz)
        elif Ly == 3:
            df_x = density_along_x_Ly3(df=df, Ly=Ly, Lz=Lz)
        else:
            raise "wrong of Ly"
        #----- fit dx["dymean"]----------------------
        print("=======================")
        print(df_x.head())
        print(df_x.tail())
        print(len(df_x))
        x_cut = list(df_x["rmean"].values)[1:32]
        x_cut = np.arange(x_cut[0],x_cut[-1],0.02)
        n0=0.81
        A0=0.05
        Kc=0.61
        kf=np.pi / 5.2
        phi= -0.4
        nx = density_fit(x=x_cut,n0=n0,A0=A0,Kc=Kc,kf=kf,phi=phi,Lx=Lx)
        print(nx)
        title = "n0=%.4f,A0=%.2f,Kc=%.2f,Kf=%.2f,phi=%.2f"%(n0,A0,Kc,kf,phi)
        plot_density_fit(r=df_x["rmean"],density=df_x["dymean"],label_y='<ni>(mean)',x=x_cut,nx=nx,title=title)
        #---------------save data --------------------------
        Issave = True
        df_fit = pd.DataFrame(columns=["x","nx","n0","A0","Kc","kf","phi"])
        df_fit['x'] = x_cut
        df_fit['nx'] = nx
        df_fit['n0'] = n0
        df_fit["A0"] = A0
        df_fit["Kc"] = Kc
        df_fit["kf"] = kf
        df_fit["phi"] = phi
        print("df_fit=\n",df_fit.head())
        if Issave:
            filepath4 = "\\measurement_density_fit.parquet"
            savepath = workpath + filepath1 + filepath2 + filepath4
            df_fit.to_parquet(savepath)