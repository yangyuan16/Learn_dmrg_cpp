import numpy as np
import pandas as pd
import matplotlib.pylab as plt
# 得到最近邻格点之间的自旋关联
def pair_nn(df,bonds): # 近邻关联对
    return
#--------------------------------------------------------------------
def bonds_layer1_y(Lz, Ly, Lx):
    b_list = []
    #b = (Ly-1) * Lx + (Lx-1) * Ly
    for it1 in range(Lx):
        for it2 in range(Ly-1):
            b1 = it1 * Lz * Ly + it2 * Lz
            b2 = b1 + 2
            b_list.append([b1,b2])
    #print(b_list)
    return b_list
def bonds_layer1_x(Lz,Ly,Lx):
    b_list = []
    for it1 in range(Ly):
        for it2 in range(Lx-1):
            b1 = it2 * Lz * Ly + it1 * 2
            b2 = b1 + 6
            b_list.append([b1,b2])
    #print(b_list)    
    return b_list
def bonds_layer2_y(Lz, Ly, Lx):
    b_list = []
    #b = (Ly-1) * Lx + (Lx-1) * Ly
    for it1 in range(Lx):
        for it2 in range(Ly-1):
            b1 = it1 * Lz * Ly + it2 * Lz + 1
            b2 = b1 + 2 
            b_list.append([b1,b2])
    #print(b_list)
    return b_list
def bonds_layer2_x(Lz,Ly,Lx):
    b_list = []
    for it1 in range(Ly):
        for it2 in range(Lx-1):
            b1 = it2 * Lz * Ly + it1 * 2 + 1
            b2 = b1 + 6
            b_list.append([b1,b2])
    #print(b_list)    
    return b_list
def bonds_z(Ly,Lx):
    b_list = []
    for it in range(Ly*Lx):
        b1 = it * 2
        b2 = it  * 2 + 1
        b_list.append([b1,b2])
    #print(b_list)
    return b_list
def bonds_pattern():
    return
#------------------------------------------------------------------------
def get_bonds_corre(df,bonds):
    dff = pd.DataFrame(columns=['b_s1', 'b_s2', 'corre'])
    b_s1 = []
    b_s2 = []
    b_corre = []
    for it in range(len(bonds)):
        s1 = bonds[it][0]
        s2 = bonds[it][1]
        corre = df.loc[(df["site1"] == s1) & (df["site2"] == s2)]["corre"].values[0]
        b_s1.append(s1)
        b_s2.append(s2)
        b_corre.append(corre)
    dff['b_s1'] = b_s1
    dff['b_s2'] = b_s2
    dff['corre'] = b_corre
    #print(dff)
    return dff
#
def build_coordinates(Ly, Lx):
    df = pd.DataFrame(columns=["x","y","z","site"])
    x_list = []
    y_list = []
    z_list = []
    site_list = []
    for it1 in range(Ly):
        for it2 in range(Lx):
            x = it2
            y = it1
            z = 0
            site = x * 6 + 2 * y
            x_list.append(x)
            y_list.append(y)
            z_list.append(z)
            site_list.append(site)
    #
    for it1 in range(Ly):
        for it2 in range(Lx):
            x = it2
            y = it1
            z = 1
            site = x * 6 + 2 * y + 1
            x_list.append(x)
            y_list.append(y)
            z_list.append(z)
            site_list.append(site)
    df["x"] = x_list
    df["y"] = y_list
    df["z"] = z_list
    df["site"] = site_list
    return df
#
def plot_corre_layer1(df_bonds_corre, df_coord): # correlation on layer 1
    fig, ax = plt.subplots(figsize=(12,1))
    plt.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.95)
    plt.axis('off')
    df_bonds_corre['abs_corre'] = np.abs(df_bonds_corre['corre'])
    for it1 in range(len(df_bonds_corre)):
        s1 = df_bonds_corre['b_s1'][it1]
        s2 = df_bonds_corre['b_s2'][it1]
        s1_x = df_coord.loc[df_coord["site"]==s1]["x"].values[0]
        s1_y = df_coord.loc[df_coord["site"]==s1]["y"].values[0]
        s2_x = df_coord.loc[df_coord["site"]==s2]["x"].values[0]
        s2_y = df_coord.loc[df_coord["site"]==s2]["y"].values[0]
        lw = 6 * (df_bonds_corre['abs_corre'][it1] - df_bonds_corre['abs_corre'].min()) /  (
            df_bonds_corre['abs_corre'].max() - df_bonds_corre['abs_corre'].min())
        if df_bonds_corre['corre'][it1] > 0:
            ax.plot([s1_x, s2_x], [s1_y, s2_y], '-r', linewidth=lw) # >0, red
        else:
            ax.plot([s1_x, s2_x], [s1_y, s2_y], '-b', linewidth=lw) # <0, blue
        if s1_x == s2_x:
            ax.text(s1_x,s1_y+0.1, str(round(df_bonds_corre['corre'][it1],4)), rotation=90,
                    fontdict={'fontsize':8, 'color':'k'})
        else:
            ax.text(s1_x+0.1,s1_y, str(round(df_bonds_corre['corre'][it1],4)), rotation=0,
                    fontdict={'fontsize':8, 'color':'k'})
    plt.show()  
    return
#
def plot_corre_layer2(df_bonds_corre, df_coord): # correlation on layer 2
    fig, ax = plt.subplots(figsize=(12,1))
    plt.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.95)
    plt.axis('off')
    df_bonds_corre['abs_corre'] = np.abs(df_bonds_corre['corre'])
    for it1 in range(len(df_bonds_corre)):
        s1 = df_bonds_corre['b_s1'][it1]
        s2 = df_bonds_corre['b_s2'][it1]
        s1_x = df_coord.loc[df_coord["site"]==s1]["x"].values[0]
        s1_y = df_coord.loc[df_coord["site"]==s1]["y"].values[0]
        s2_x = df_coord.loc[df_coord["site"]==s2]["x"].values[0]
        s2_y = df_coord.loc[df_coord["site"]==s2]["y"].values[0]
        lw = 6 * (df_bonds_corre['abs_corre'][it1] - df_bonds_corre['abs_corre'].min()) /  (
            df_bonds_corre['abs_corre'].max() - df_bonds_corre['abs_corre'].min())
        if df_bonds_corre['corre'][it1] > 0:
            ax.plot([s1_x, s2_x], [s1_y, s2_y], '-r', linewidth=lw) # >0, red
        else:
            ax.plot([s1_x, s2_x], [s1_y, s2_y], '-b', linewidth=lw) # <0, blue
        if s1_x == s2_x:
            ax.text(s1_x,s1_y+0.1, str(round(df_bonds_corre['corre'][it1],4)), rotation=90,
                    fontdict={'fontsize':8, 'color':'k'})
        else:
            ax.text(s1_x+0.1,s1_y, str(round(df_bonds_corre['corre'][it1],4)), rotation=0,
                    fontdict={'fontsize':8, 'color':'k'})
    plt.show()  
    return
#
def plot_corre_z(df_bonds_corre, df_coord): # correlation between two layers
    fig, ax = plt.subplots(figsize=(12,1))
    plt.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.95)
    plt.axis('off')
    df_bonds_corre['abs_corre'] = np.abs(df_bonds_corre['corre'])
    for it1 in range(len(df_bonds_corre)):
        s1 = df_bonds_corre['b_s1'][it1]
        s1_x = df_coord.loc[df_coord['site']==s1]["x"].values[0]
        s1_y = df_coord.loc[df_coord['site']==s1]["y"].values[0]
        ms = 16 * (df_bonds_corre['abs_corre'][it1] - df_bonds_corre['abs_corre'].min()) /  (
            df_bonds_corre['abs_corre'].max() - df_bonds_corre['abs_corre'].min())
        if df_bonds_corre['corre'][it1] > 0:
            ax.plot(s1_x, s1_y, 'ro', markersize=ms)
        else:
            ax.plot(s1_x, s1_y, 'bo', markersize=ms)
        ax.text(s1_x+0.1,s1_y+0.1, str(round(df_bonds_corre['corre'][it1],2)), rotation=0,
                fontdict={'fontsize':6, 'color':'k'})
    plt.show()
    return    
#
if __name__ == "__main__":
    print()
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
    filepath3 = "\\measurement_sq_spin.dat"
    filename = workpath + filepath1 + filepath2 + filepath3
    print(filename)
    print('------------load the data------------------')
    df = pd.read_csv(filename, header=None, sep='\t',encoding='utf-8')
    df.rename(columns={0: "site1", 1: "site2", 2: "corre"},inplace=True)
    df.sort_values(["site1","site2"],inplace=True)
    #print(df.head())
    #print(df.tail())
    #print(len(df))
    df_coord = build_coordinates(Ly, Lx) # get the coordiantes of each sites
    #print(df_coord)
    print('---------creat bonds pattern----------------')
    b1y= bonds_layer1_y(Lz=Lz, Ly=Ly, Lx=Lx)
    b1x =bonds_layer1_x(Lz=Lz, Ly=Ly, Lx=Lx)
    bz = bonds_z(Ly=Ly, Lx=Lx)
    b2y = bonds_layer2_y(Lz=Lz, Ly=Ly, Lx=Lx)
    b2x = bonds_layer2_x(Lz=Lz, Ly=Ly, Lx=Lx)
    b1 = b1y + b1x
    df_bonds_corre1 = get_bonds_corre(df, bonds=b1) # layer1
    b2 = b2y + b2x
    df_bonds_corre2 = get_bonds_corre(df, bonds=b2) # layer2
    df_bonds_correz = get_bonds_corre(df, bonds=bz) # between layer1 and layer2    
    print('------------plot the data------------------')
    print('----df_bonds_coore1----')
    print(df_bonds_corre1)
    print('----df_bonds_coore2----')
    print(df_bonds_corre2)
    print('----df_bonds_coorez----')
    print(df_bonds_correz)
    print('--df_coord----------')
    print(df_coord)
    print('------------plot the data------------------')
    plot_corre_layer1(df_bonds_corre=df_bonds_corre1,df_coord=df_coord)
    plot_corre_layer2(df_bonds_corre=df_bonds_corre2,df_coord=df_coord)
    plot_corre_z(df_bonds_corre=df_bonds_correz, df_coord=df_coord)

