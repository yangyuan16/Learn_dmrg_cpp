import numpy as np
import pandas as pd
import matplotlib.pylab as plt
#
def entropy(Lx, x, c, g):
    S = []
    for it in x:
        s = (c/6)*np.log((Lx/np.pi)*np.sin(it*np.pi/Lx)) + g
        S.append(s)    
    return S
#
def df_along_Lx(df,Lz,Ly):
    print("---------df_res0---------")
    df_res0 = df[df["site"]%(Lz*Ly)==0]
    df_res0["r"] = range(1,len(df_res0)+1) 
    print(df_res0.head())
    print(df_res0.tail())
    print(len(df_res0))
    print("---------df_res1---------")
    df_res1 = df[df["site"]%(Lz*Ly)==1]
    df_res1["r"] = range(1,len(df_res1)+1) 
    print(df_res1.head())
    print(df_res1.tail())
    print(len(df_res1))
    return df_res0, df_res1
#
def plot_entropy(df0,df1,):
    #-----------------plot ax1---------------------------
    fig = plt.figure(figsize=(6,10)) 
    ax1 = plt.subplot(2,2,1)
    ax2 = plt.subplot(2,2,2)
    ax3 = plt.subplot(2,2,3)
    ax4 = plt.subplot(2,2,4)
    #----------------------------------------------------
    plt.sca(ax1)  ##选择对ax1进行绘图
    ax1=plt.gca() #获得坐标轴的句柄
    r = df0["r"].values
    L1, = ax1.plot(r,df0["entropy"],label="Jz=%.2f"%Jz,ls="-",lw=1.5,color="red",
             marker='o',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="red",
             markerfacecolor='None')
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 22, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[L1], loc = 4, bbox_to_anchor=(0.48, 0.18),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    label_x = r"x"
    label_y = "entropy"
    #plt.yscale("log")
    ax1.set_xlabel(label_x, size= 14)
    ax1.set_ylabel(label_y, size= 14)
    ax1.tick_params(labelsize = 15) # 设置坐标刻度对应数字的大小
    #ax1.set_xlim([0,8])
    #ax1.set_ylim([-0.1,1])
    #ax1.set_xticks([0,2,4,6,8])
    #ax1.set_yticks([-0.1,0,0.5,1])
    #ax1.text(0.3,-0.05, r'$\mathrm{(a)}$', fontsize=18)
    plt.title("case1",fontsize=16)
    #--------------------------------------------------    
    plt.sca(ax2)  ##选择对ax2进行绘图
    ax2=plt.gca() #获得坐标轴的句柄
    r = df1["r"].values
    L2, = ax2.plot(r,df1["entropy"],label="Jz=%.2f"%Jz,ls="-",lw=1.5,color="blue",
             marker='o',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="blue",
             markerfacecolor='None')
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 22, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[L2], loc = 4, bbox_to_anchor=(0.48, 0.18),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    label_x = r"x"
    label_y = "entropy"
    #plt.yscale("log")
    ax2.set_xlabel(label_x, size= 14)
    ax2.set_ylabel(label_y, size= 14)
    ax2.tick_params(labelsize = 15) # 设置坐标刻度对应数字的大小
    #ax2.set_xlim([0,8])
    #ax2.set_ylim([-0.1,1])
    #ax2.set_xticks([0,2,4,6,8])
    #ax2.set_yticks([-0.1,0,0.5,1])
    #ax2.text(0.3,-0.05, r'$\mathrm{(a)}$', fontsize=18)
    plt.title("case2",fontsize=16) 
    #--------------------------------------------------------
    plt.show()
    return
#
def plot_entropy_v2(df0,df1,):
    #-----------------plot ax1---------------------------
    fig = plt.figure(figsize=(6,10)) 
    ax1 = plt.subplot(1,1,1)
    #----------------------------------------------------
    plt.sca(ax1)  ##选择对ax1进行绘图
    ax1=plt.gca() #获得坐标轴的句柄
    r = df0["r"].values
    L1, = ax1.plot(r,df0["entropy"],label="case1",ls="-",lw=1.5,color="red",
             marker='o',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="red",
             markerfacecolor='None')
    #
    r = df1["r"].values
    L2, = ax1.plot(r,df1["entropy"],label="case2",ls="-",lw=1.5,color="blue",
             marker='o',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="blue",
             markerfacecolor='None')
    #
    #
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 22, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[L1,L2,], loc = 4, bbox_to_anchor=(0.48, 0.18),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    label_x = r"x"
    label_y = "entropy"
    #plt.yscale("log")
    ax1.set_xlabel(label_x, size= 14)
    ax1.set_ylabel(label_y, size= 14)
    ax1.tick_params(labelsize = 15) # 设置坐标刻度对应数字的大小
    #ax1.set_xlim([0,8])
    #ax1.set_ylim([-0.1,1])
    #ax1.set_xticks([0,2,4,6,8])
    #ax1.set_yticks([-0.1,0,0.5,1])
    #ax1.text(0.3,-0.05, r'$\mathrm{(a)}$', fontsize=18)
    plt.title("Jz=%.2f"%Jz,fontsize=16)
    plt.show()
    return
#
def plot_entropy_fit(df,r_ent,entropy_fit,c,g,color):
    #-----------------plot ax1---------------------------
    fig = plt.figure(figsize=(6,10)) 
    ax1 = plt.subplot(1,1,1)
    #----------------------------------------------------
    plt.sca(ax1)  ##选择对ax1进行绘图
    ax1=plt.gca() #获得坐标轴的句柄
    r = df["r"].values
    L1, = ax1.plot(r,df["entropy"],label=" ",ls="-",lw=1.5,color=color,
             marker='o',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor=color,
             markerfacecolor='None')
    #
    L2, = ax1.plot(r_ent,entropy_fit,label="c=%.2f,g=%.2f"%(c,g),ls="--",lw=1.5,color="k",
             marker='o',alpha=1,markersize=0,markeredgewidth=1.5, markeredgecolor="k",
             markerfacecolor='None')
    #
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 22, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[L1,L2,], loc = 4, bbox_to_anchor=(0.48, 0.18),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    label_x = r"x"
    label_y = "entropy"
    #plt.yscale("log")
    ax1.set_xlabel(label_x, size= 14)
    ax1.set_ylabel(label_y, size= 14)
    ax1.tick_params(labelsize = 15) # 设置坐标刻度对应数字的大小
    #ax1.set_xlim([0,8])
    #ax1.set_ylim([-0.1,1])
    #ax1.set_xticks([0,2,4,6,8])
    #ax1.set_yticks([-0.1,0,0.5,1])
    #ax1.text(0.3,-0.05, r'$\mathrm{(a)}$', fontsize=18)
    plt.title("Jz=%.2f"%Jz,fontsize=16)
    plt.show()
    return
#
if __name__ == "__main__":
    print()
    Lz = 1
    Ly = 2
    Lx = 128
    dop = 112
    t = 3
    J = 1
    Jz = 0
    dim = 1200 # dim cutoff
    workpath = "E:\\WORK\\Work\\Project\\La3Ni2O7"
    filepath1 = "\\data_xinlu\\Lz%d_Ly%d_Lx%d\\%g" % (Lz, Ly, Lx,dop)
    filepath2 = "\\t%d_J%d_Jz%d_dim%d_ent\\entanglement\\entropy" % (t, J, Jz, dim)
    filepath3 = "\\sys-ground_state"
    filename = workpath + filepath1 + filepath2 + filepath3
    print(filename)
    # load the data
    df = pd.read_csv(filename, header=None, sep='\t',encoding='utf-8')
    df.rename(columns={0: "site", 1: "entropy",}, inplace=True)
    L1 = Lz * Ly * Lx - 3
    L2 = len(df)
    df = df.iloc[L2-L1:L2]
    df.sort_values(["site"],inplace=True)
    print(df.head())
    print(df.tail())
    print(len(df))
    df_res0, df_res1, = df_along_Lx(df=df,Lz=Lz,Ly=Ly)
    plot_entropy(df0=df_res0,df1=df_res1,)
    plot_entropy_v2(df0=df_res0,df1=df_res1,)
    #=================================================
    c_list =[2,] 
    g_list = [0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4]
    color = "red"
    for c in c_list:
        for g in g_list:
            S_res0 = entropy(Lx=Lx,x=df_res0["r"].values,c=c,g=g)
            S_res1 = entropy(Lx=Lx,x=df_res1["r"].values,c=c,g=g)
            print("S_res0:\n", S_res0)
            print("S_res1:\n", S_res1)
            plot_entropy_fit(df=df_res0,r_ent = df_res0["r"].values,entropy_fit=S_res0,
                             c=c,g=g,color=color)


