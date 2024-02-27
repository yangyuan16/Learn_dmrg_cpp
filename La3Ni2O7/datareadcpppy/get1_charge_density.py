import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pylab as plt

def read_electron_density():
    return
#
if __name__ == "__main__":
    Lz = 2
    Ly = 3
    Lx = 48
    dop = 0.25
    t = 3
    J = 1
    Jz = 2
    dim = 6000 # dim cutoff
    workpath = "E:\\WORK\\Work\\Project\\La3Ni2O7"
    filepath1 = "\\data_dmrgcpp\\Lz%d_Ly%d_Lx%d\\dop%g" % (Lz, Ly, Lx,dop)
    filepath2 = "\\t%d_J%d_Jz%g_dim%d" % (t, J, Jz, dim)
    filepath3 = "\\measurement_electron_density.dat"
    filename = workpath + filepath1 + filepath2 + filepath3
    print(filename)
    # load the data
    df = pd.read_csv(filename, header=None, sep='\t',encoding='utf-8')
    df.rename(columns={0: "site", 1: "density"},inplace=True)
    df_layer1 = df[df["site"] % 2 == 0] 
    df_layer2 = df[df["site"] % 2 == 1]
    sites_layer1 = df_layer1["site"].values.reshape(-1,Ly).T
    density_layer1 = df_layer1["density"].values.reshape(-1, Ly).T
    sites_layer2 = df_layer2["site"].values.reshape(-1, Ly).T
    density_layer2 = df_layer2["density"].values.reshape(-1, Ly).T
    #------------------------------------------------------------
    # plot the data
    fig = plt.figure(figsize=(12,2))
    ax1 = plt.subplot(2,1,1)
    ax2 = plt.subplot(2,1,2)
    #
    vmin1 = np.array(density_layer1).min()
    vmax1 = np.array(density_layer1).max()
    vmin2 = np.array(density_layer2).min()
    vmax2 = np.array(density_layer2).max()
    vmin = min(vmin1, vmin2)
    vmax = max(vmax1, vmax2)
    plt.sca(ax1)  ##选择对ax1进行绘图
    ax1=plt.gca() #获得坐标轴的句柄
    text = "Lz%d_Ly%d_Lx%d_dop%g_" % (Lz, Ly, Lx,dop) + "t%d_J%d_Jz%g_dim%d" % (t, J, Jz, dim)
    plt.text(0.3, -0.25, text)
    sns.heatmap(density_layer1, cmap="YlGnBu", annot=True, linewidths=0.5, vmax=vmax, vmin=vmin,
                annot_kws={'size':9,'weight':'bold', 'color':'white', "rotation":90})
    ax1.set_title("charge density, layer1")
    ax1.set_xlabel("")
    ax1.set_xticklabels([]) # 设置 x 轴图例为空置
    ax1.set_ylabel("")
    #---
    plt.sca(ax2)  ##选择对ax1进行绘图
    ax2=plt.gca() #获得坐标轴的句柄
    sns.heatmap(density_layer2, cmap="YlGnBu", annot=True, linewidths=0.5, vmax=vmax, vmin=vmin,
                annot_kws={'size':9,'weight':'bold', 'color':'white', "rotation":90})
    ax2.set_title("charge density, layer2")
    ax2.set_xlabel("")
    ax2.set_xticklabels([]) # 设置 x 轴图例为空置
    ax2.set_ylabel("")
    #f.savefig("charge_density.jpg", bbox_inches="tight")
    plt.show()