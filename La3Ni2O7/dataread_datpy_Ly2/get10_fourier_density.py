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
def plot_loopJz_density_curve(Jz):
    if Ly == 2:
        df_x = density_along_x_Ly2(df=df, Ly=Ly, Lz=Lz)
    elif Ly == 3:
        df_x = density_along_x_Ly3(df=df, Ly=Ly, Lz=Lz)
    else:
        raise "wrong of Ly"

    #----plot density curve----------------
    if Ly == 2:
        #plot_density_curve(r=df_x["r0"],density=df_x["dy0"],label_y='<ni>(y0)',Jz=Jz)
        #plot_density_curve(r=df_x["r1"],density=df_x["dy1"],label_y='<ni>(y1)',Jz=Jz)
        #plot_density_curve(r=df_x["r2"],density=df_x["dy2"],label_y='<ni>(y2)',Jz=Jz)
        #plot_density_curve(r=df_x["r3"],density=df_x["dy3"],label_y='<ni>(y3)',Jz=Jz)  
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
#---------------------------------------------------------------------------------------
def fourier_transform(k,density):
    N = len(density)
    dft = 0
    for n in range(N):
        dft += density[n] * np.exp(-1j * k * n ) / N
    return dft
def full_k_fourier_transform(k_list,density):
    dft_list = []
    for k in k_list:
        dft = fourier_transform(k=k,density=density)
        dft_list.append(dft)
    return dft_list
#
def plot_loopJz_density_dft_curve(k_list,dft_list,Jz):
    fig = plt.figure(figsize=(10,3)) 
    ax1 = plt.subplot(1,2,1)
    plt.sca(ax1)  ##选择对ax1进行绘图
    ax1=plt.gca() #获得坐标轴的句柄
    ax1.plot(np.array(k_list)/np.pi, np.real(dft_list), label=r"G",ls="-",lw=1.5,color="red",
             marker='o',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="k",
             markerfacecolor='w')
    label_x = r"k"
    label_y = "real part"
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
    #
    ax2 = plt.subplot(1,2,2)
    plt.sca(ax2)  ##选择对ax1进行绘图
    ax2=plt.gca() #获得坐标轴的句柄
    ax2.plot(np.array(k_list)/np.pi, np.abs(dft_list), label=r"G",ls="-",lw=1.5,color="blue",
             marker='o',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="k",
             markerfacecolor='w')
    label_x = r"k"
    label_y = "Magnitude"
    ax2.set_xlabel(label_x, size= 14)
    ax2.set_ylabel(label_y, size= 14)
    #
    ax2.tick_params(labelsize = 15) # 设置坐标刻度对应数字的大小
    #ax2.set_xlim([0,8])
    #ax2.set_ylim([-0.1,1])
    #ax2.set_xticks([0,2,4,6,8])
    #ax2.set_yticks([-0.1,0,0.5,1])
    #ax2.text(0.3,-0.05, r'$\mathrm{(a)}$', fontsize=18)
    plt.title("Jz=%.2f"%Jz)
    plt.show()
    return
#
if __name__ == "__main__":
    Lz = 2
    Ly = 2
    Lx = 48
    dop = 60
    t = 3
    J = 1
    Jz_list = [0.1,0.2,0.25,0.4,0.5,0.6,0.75,0.8,1.0,1.25,1.5,1.75,2.0,
               2.25,2.5,2.75]
    '''              0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,
               1.3,1.4,1.5,1.6,1.8,1.9,2.0,2.1,2.2,2.6,3.0,4.0'''
    k_list = list(np.arange(-np.pi,np.pi+np.pi/100,np.pi/100))
    print("k=", k_list) 
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
        df_x = density_along_x_Ly2(df=df, Ly=Ly, Lz=Lz)
        print(df_x.head())
        # calculate the full K spcace fourier transform results
        dft_list = full_k_fourier_transform(k_list=k_list,density=df_x["dymean"].values)
        print("dft_list:\n", dft_list)
        print(np.abs(dft_list))
        print(np.real(dft_list))
        #plot_loopJz_density_curve(Jz=Jz)
        plot_loopJz_density_dft_curve(k_list,dft_list,Jz)
        
        
        
