import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pylab as plt
#
def get_bond_corre_yy_1(Lz, Ly, Lx, df):
    # create yy_1 bonds
    bond_s1 = []
    bond_s2 = []
    bond_corre_list = []
    for it in range(Lx):
        s1 = it * Ly * Lz
        s2 = s1 + 2
        bond_s1.append(s1)
        bond_s2.append(s2)
        bond_corre_list.append(df[(df['site1']== s1) & (df['site2'] == s2)]['bond_corre'].values[0])

    df_bond = pd.DataFrame(columns=['site1','site2','bond_corre'])
    df_bond['site1'] = bond_s1
    df_bond['site2'] = bond_s2
    df_bond['bond_corre'] = bond_corre_list
    return df_bond
#
def get_bond_corre_yy_2(Lz, Ly, Lx, df):
    # create yy_1 bonds
    bond_s1 = []
    bond_s2 = []
    bond_corre_list = []
    for it in range(Lx):
        s1 = it * Ly * Lz + 1
        s2 = s1 + 2
        bond_s1.append(s1)
        bond_s2.append(s2)
        bond_corre_list.append(df[(df['site1']== s1) & (df['site2'] == s2)]['bond_corre'].values[0])

    df_bond = pd.DataFrame(columns=['site1','site2','bond_corre'])
    df_bond['site1'] = bond_s1
    df_bond['site2'] = bond_s2
    df_bond['bond_corre'] = bond_corre_list
    return df_bond
#
def get_bond_corre_xx_1(Lz, Ly, Lx, df):
    # create yy_1 bonds
    bond_s1 = []
    bond_s2 = []
    bond_corre_list = []
    for it in range(Lx-1):
        s1 = it * Ly * Lz 
        s2 = s1 + 4
        bond_s1.append(s1)
        bond_s2.append(s2)
        #print("it:",it)
        bond_corre_list.append(df[(df['site1']== s1) & (df['site2'] == s2)]['bond_corre'].values[0])

    df_bond = pd.DataFrame(columns=['site1','site2','bond_corre'])
    df_bond['site1'] = bond_s1
    df_bond['site2'] = bond_s2
    df_bond['bond_corre'] = bond_corre_list
    return df_bond
#
def get_bond_corre_xx_2(Lz, Ly, Lx, df):
    # create yy_1 bonds
    bond_s1 = []
    bond_s2 = []
    bond_corre_list = []
    for it in range(Lx-1):
        s1 = it * Ly * Lz + 1
        s2 = s1 + 4
        bond_s1.append(s1)
        bond_s2.append(s2)
        #print("it:",it)
        bond_corre_list.append(df[(df['site1']== s1) & (df['site2'] == s2)]['bond_corre'].values[0])

    df_bond = pd.DataFrame(columns=['site1','site2','bond_corre'])
    df_bond['site1'] = bond_s1
    df_bond['site2'] = bond_s2
    df_bond['bond_corre'] = bond_corre_list
    return df_bond
#
def get_bond_corre_xx_3(Lz, Ly, Lx, df):
    # create yy_1 bonds
    bond_s1 = []
    bond_s2 = []
    bond_corre_list = []
    for it in range(Lx-1):
        s1 = it * Ly * Lz + 2
        s2 = s1 + 4
        bond_s1.append(s1)
        bond_s2.append(s2)
        #print("it:",it)
        bond_corre_list.append(df[(df['site1']== s1) & (df['site2'] == s2)]['bond_corre'].values[0])

    df_bond = pd.DataFrame(columns=['site1','site2','bond_corre'])
    df_bond['site1'] = bond_s1
    df_bond['site2'] = bond_s2
    df_bond['bond_corre'] = bond_corre_list
    return df_bond
#
def get_bond_corre_xx_4(Lz, Ly, Lx, df):
    # create yy_1 bonds
    bond_s1 = []
    bond_s2 = []
    bond_corre_list = []
    for it in range(Lx-1):
        s1 = it * Ly * Lz + 3
        s2 = s1 + 4
        bond_s1.append(s1)
        bond_s2.append(s2)
        #print("it:",it)
        bond_corre_list.append(df[(df['site1']== s1) & (df['site2'] == s2)]['bond_corre'].values[0])

    df_bond = pd.DataFrame(columns=['site1','site2','bond_corre'])
    df_bond['site1'] = bond_s1
    df_bond['site2'] = bond_s2
    df_bond['bond_corre'] = bond_corre_list
    return df_bond
#
def get_bond_corre_zz_1(Lz, Ly, Lx, df):
    # create yy_1 bonds
    bond_s1 = []
    bond_s2 = []
    bond_corre_list = []
    for it in range(Lx):
        s1 = it * Ly * Lz 
        s2 = s1 + 1
        bond_s1.append(s1)
        bond_s2.append(s2)
        #print("it:",it)
        bond_corre_list.append(df[(df['site1']== s1) & (df['site2'] == s2)]['bond_corre'].values[0])

    df_bond = pd.DataFrame(columns=['site1','site2','bond_corre'])
    df_bond['site1'] = bond_s1
    df_bond['site2'] = bond_s2
    df_bond['bond_corre'] = bond_corre_list
    return df_bond
#
def get_bond_corre_zz_2(Lz, Ly, Lx, df):
    # create yy_1 bonds
    bond_s1 = []
    bond_s2 = []
    bond_corre_list = []
    for it in range(Lx):
        s1 = it * Ly * Lz + 2 
        s2 = s1 + 1
        bond_s1.append(s1)
        bond_s2.append(s2)
        #print("it:",it)
        bond_corre_list.append(df[(df['site1']== s1) & (df['site2'] == s2)]['bond_corre'].values[0])

    df_bond = pd.DataFrame(columns=['site1','site2','bond_corre'])
    df_bond['site1'] = bond_s1
    df_bond['site2'] = bond_s2
    df_bond['bond_corre'] = bond_corre_list
    return df_bond

def get_bond_corre(Lz,Ly,Lx, df ):
        df_bond_yy1 = get_bond_corre_yy_1(Lz=Lz,Ly=Ly,Lx=Lx,df=df)
        print('df_bond_yy1')
        print(df_bond_yy1.head())

        df_bond_yy2 = get_bond_corre_yy_2(Lz=Lz,Ly=Ly,Lx=Lx,df=df)
        print('df_bond_yy2')
        print(df_bond_yy2.head())

        df_bond_xx1 = get_bond_corre_xx_1(Lz=Lz,Ly=Ly,Lx=Lx,df=df)
        print('df_bond_xx1')
        print(df_bond_xx1.head())

        df_bond_xx2 = get_bond_corre_xx_2(Lz=Lz,Ly=Ly,Lx=Lx,df=df)
        print('df_bond_xx2')
        print(df_bond_xx2.head())

        df_bond_xx3 = get_bond_corre_xx_3(Lz=Lz,Ly=Ly,Lx=Lx,df=df)
        print('df_bond_xx3')
        print(df_bond_xx3.head())

        df_bond_xx4 = get_bond_corre_xx_4(Lz=Lz,Ly=Ly,Lx=Lx,df=df)
        print('df_bond_xx4')
        print(df_bond_xx4.head())

        df_bond_zz1 = get_bond_corre_zz_1(Lz=Lz,Ly=Ly,Lx=Lx,df=df)
        print('df_bond_zz1')
        print(df_bond_zz1.head())

        df_bond_zz2 = get_bond_corre_zz_2(Lz=Lz,Ly=Ly,Lx=Lx,df=df)
        print('df_bond_zz2')
        print(df_bond_zz2.head())
        #
        df_dict = {}
        df_dict['yy1'] = df_bond_yy1
        df_dict['yy2'] = df_bond_yy2
        df_dict['xx1'] = df_bond_xx1
        df_dict['xx2'] = df_bond_xx2
        df_dict['xx3'] = df_bond_xx3
        df_dict['xx4'] = df_bond_xx4
        df_dict['zz1'] = df_bond_zz1
        df_dict['zz2'] = df_bond_zz2

        return df_dict
    
# 
def plot_bond_density_curve(df_dict_density, Jz):
    #------- fig yy1
    fig = plt.figure(figsize=(10,3)) 
    ax1 = plt.subplot(1,1,1)
    plt.sca(ax1)  ##选择对ax1进行绘图
    ax1=plt.gca() #获得坐标轴的句柄
    r = range(len(df_dict_density['yy1']['bond_corre']))
    ax1.plot(r,df_dict_density['yy1']['bond_corre'],label=r"G",ls="-",lw=1.5,color="b",
             marker='o',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="b",
             markerfacecolor='None')
    label_x = r"i"
    label_y = "bond density YY 1"
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
    #------- fig yy2
    fig = plt.figure(figsize=(10,3)) 
    ax1 = plt.subplot(1,1,1)
    plt.sca(ax1)  ##选择对ax1进行绘图
    ax1=plt.gca() #获得坐标轴的句柄
    r = range(len(df_dict_density['yy2']['bond_corre']))
    ax1.plot(r,df_dict_density['yy2']['bond_corre'],label=r"G",ls="-",lw=1.5,color="b",
             marker='o',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="b",
             markerfacecolor='None')
    label_x = r"i"
    label_y = "bond density YY 2"
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
    #------- fig xx1
    fig = plt.figure(figsize=(10,3)) 
    ax1 = plt.subplot(1,1,1)
    plt.sca(ax1)  ##选择对ax1进行绘图
    ax1=plt.gca() #获得坐标轴的句柄
    r = range(len(df_dict_density['xx1']['bond_corre']))
    ax1.plot(r,df_dict_density['xx1']['bond_corre'],label=r"G",ls="-",lw=1.5,color="b",
             marker='o',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="b",
             markerfacecolor='None')
    label_x = r"i"
    label_y = "bond density XX 1"
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
    #------- fig xx2
    fig = plt.figure(figsize=(10,3)) 
    ax1 = plt.subplot(1,1,1)
    plt.sca(ax1)  ##选择对ax1进行绘图
    ax1=plt.gca() #获得坐标轴的句柄
    r = range(len(df_dict_density['xx2']['bond_corre']))
    ax1.plot(r,df_dict_density['xx2']['bond_corre'],label=r"G",ls="-",lw=1.5,color="b",
             marker='o',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="b",
             markerfacecolor='None')
    label_x = r"i"
    label_y = "bond density XX 2"
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
    #------- fig xx3
    fig = plt.figure(figsize=(10,3)) 
    ax1 = plt.subplot(1,1,1)
    plt.sca(ax1)  ##选择对ax1进行绘图
    ax1=plt.gca() #获得坐标轴的句柄
    r = range(len(df_dict_density['xx3']['bond_corre']))
    ax1.plot(r,df_dict_density['xx3']['bond_corre'],label=r"G",ls="-",lw=1.5,color="b",
             marker='o',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="b",
             markerfacecolor='None')
    label_x = r"i"
    label_y = "bond density XX 3"
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
    #------- fig xx4
    fig = plt.figure(figsize=(10,3)) 
    ax1 = plt.subplot(1,1,1)
    plt.sca(ax1)  ##选择对ax1进行绘图
    ax1=plt.gca() #获得坐标轴的句柄
    r = range(len(df_dict_density['xx4']['bond_corre']))
    ax1.plot(r,df_dict_density['xx4']['bond_corre'],label=r"G",ls="-",lw=1.5,color="b",
             marker='o',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="b",
             markerfacecolor='None')
    label_x = r"i"
    label_y = "bond density XX 4"
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
    #------- fig zz1
    fig = plt.figure(figsize=(10,3)) 
    ax1 = plt.subplot(1,1,1)
    plt.sca(ax1)  ##选择对ax1进行绘图
    ax1=plt.gca() #获得坐标轴的句柄
    r = range(len(df_dict_density['zz1']['bond_corre']))
    ax1.plot(r,df_dict_density['zz1']['bond_corre'],label=r"G",ls="-",lw=1.5,color="b",
             marker='o',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="b",
             markerfacecolor='None')
    label_x = r"i"
    label_y = "bond density ZZ 1"
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
    #------- fig zz2
    fig = plt.figure(figsize=(10,3)) 
    ax1 = plt.subplot(1,1,1)
    plt.sca(ax1)  ##选择对ax1进行绘图
    ax1=plt.gca() #获得坐标轴的句柄
    r = range(len(df_dict_density['zz2']['bond_corre']))
    ax1.plot(r,df_dict_density['zz2']['bond_corre'],label=r"G",ls="-",lw=1.5,color="b",
             marker='o',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="b",
             markerfacecolor='None')
    label_x = r"i"
    label_y = "bond density ZZ 2"
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
    #
    return
#
#
if __name__ == "__main__":
    Lz = 2
    Ly = 2
    Lx = 48
    dop = 96 
    t = 3
    J = 1
    Jz_list = [3.0,]
    dim = 6000 # dim cutoff
    for it in range(len(Jz_list)):
        Jz = Jz_list[it]
        workpath = "E:\\WORK\\Work\\Project\\La3Ni2O7"
        filepath1 = "\\data_dmrgcpp\\Lz%d_Ly%d_Lx%d\\dop%g" % (Lz, Ly, Lx,dop)
        filepath2 = "\\t%d_J%d_Jz%.2f_dim%d" % (t, J, Jz, dim)
        filepath3 = "\\measurement_bond_spin.dat"
        filename = workpath + filepath1 + filepath2 + filepath3
        print(filename)
        # load the bond spin data
        df_spin = pd.read_csv(filename, header=None, sep='\t',encoding='utf-8')
        print(df_spin.head())
        df_spin.rename(columns={0: "site1", 1: "site2", 2: "bond_corre"},inplace=True)
        print('-----bond spin data-----')
        df_dict_spin = get_bond_corre(Lz=Lz,Lx=Lx,Ly=Ly, df=df_spin)
        #
        #------load the bond density data
        filepath3 = "\\measurement_bond_density.dat"
        filename = workpath + filepath1 + filepath2 + filepath3
        print(filename)
        df_density = pd.read_csv(filename, header=None, sep='\t',encoding='utf-8')
        print(df_density.head())
        df_density.rename(columns={0: "site1", 1: "site2", 2: "bond_corre"},inplace=True)
        print('-----bond density data-----')
        df_dict_density = get_bond_corre(Lz=Lz,Lx=Lx,Ly=Ly, df=df_density)
        #
        #-------load the bond hopping data
        filepath3 = "\\measurement_bond_hopping.dat"
        filename = workpath + filepath1 + filepath2 + filepath3
        print(filename)
        df_hopping = pd.read_csv(filename, header=None, sep='\t',encoding='utf-8')
        print(df_hopping.head())
        df_hopping.rename(columns={0: "site1", 1: "site2", 2: "bond_corre"},inplace=True)
        print('-----bond density data-----')
        df_dict_hopping = get_bond_corre(Lz=Lz,Lx=Lx,Ly=Ly, df=df_hopping)
        #
        #
        #------ plot bond density
        plot_bond_density_curve(df_dict_density=df_dict_density,Jz=Jz)