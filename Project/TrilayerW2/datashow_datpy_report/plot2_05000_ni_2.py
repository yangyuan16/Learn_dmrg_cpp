import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

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

def get_data_ni(Lz, Ly, Lx, dop, t, J, Jz, dim,):
    workpath = "E:\\WORK\\Work\\Project\\TrilayerW2-tj-Jz-openz"
    filepath1 = "\\data_dmrgcpp\\Lz%d_Ly%d_Lx%d\\dop%g" % (Lz, Ly, Lx,dop)
    filepath2 = "\\t%.2f_J%.2f_Jz%.2f_dim%d" % (t, J, Jz, dim)
    filepath3 = "\\measurement_electron_density.dat"
    filename = workpath + filepath1 + filepath2 + filepath3
    print(filename)
    # load the data
    df = pd.read_csv(filename, header=None, sep='\t',encoding='utf-8')
    df.rename(columns={0: "site", 1: "density"},inplace=True)
    df.sort_values(['site'],inplace=True)
    print(df.head())
    print(len(df))
    df_out = density_along_x_Lz3( df=df,Ly=Ly,Lz=Lz)
    print(df_out.head())
    return df_out
#---------------------------------------------------
if __name__ =="__main__":
    Lz = 3
    Ly = 2
    Lx = 48
    t = 3
    J = 1
    #
    #---------------plot logr-r fig and logr-logr fig-----------------------------
    fig = plt.figure(figsize=(20,10))
    # plt.figure(facecolor='blue',edgecolor='black') # 设置画布的颜色
    params = {
        'axes.labelsize': '30',
        'xtick.labelsize':'20',
        'ytick.labelsize':'20',
        'ytick.direction':'in',
        'xtick.direction':'in',
        'lines.linewidth':10 ,
        'legend.fontsize': '27',
        # 'figure.figsize'   : '12, 9'    # set figure size
    }
    pylab.rcParams.update(params) # set figure parameter 更新绘图的参数
    #plt.rcParams['font.family'] = 'Times New Roman'  # 设置全局字体为 Times New Roman
    # 得到子图
    ax1 = plt.axes([0.1,0.55,0.24,0.4])
    ax2 = plt.axes([0.4,0.55,0.24,0.4])
    ax3 = plt.axes([0.7,0.55,0.24,0.4])
    ax4 = plt.axes([0.1,0.1,0.24,0.4])
    ax5 = plt.axes([0.4,0.1,0.24,0.4])
    ax6 = plt.axes([0.7,0.1,0.24,0.4])
    # plt.subplots_adjust(left=0.1, bottom=0.1, right=0.98, top=0.96, wspace=0.32, hspace=0.26)
    #-----------------------------------------------------------------------------------------------
    # 选择子图 ax1 进行绘图
    plt.sca(ax1) ## 选择对 ax1 进行绘图
    ax1 = plt.gca()
    df_ni_Jz90 = get_data_ni(Lz=Lz,Ly=Ly,Lx=Lx,dop=144,t=t,J=J,Jz=9.0,dim=10000)
    #
    L, =ax1.plot([0,50],[0.5,0.5],label=" ",ls="--",lw=1.5,color="k",
             marker='o',alpha=1,markersize=0,markeredgewidth=1.5, markeredgecolor="red", markerfacecolor='None')
    ni_layer1 = np.array(df_ni_Jz90['dy0'] /2 + df_ni_Jz90["dy1"] / 2)
    label = r"$\mu =1$"
    L90_mu1, =ax1.plot(df_ni_Jz90["rmean"].values,ni_layer1,label=label,ls="-",lw=1.5,color="red",
             marker='o',alpha=1,markersize=9,markeredgewidth=1.5, markeredgecolor="red", markerfacecolor='None')
    #
    ni_layer2 = np.array(df_ni_Jz90['dy2'] /2 + df_ni_Jz90["dy3"] / 2)
    label = r"$\mu =2$"
    L90_mu2, =ax1.plot(df_ni_Jz90["rmean"].values,ni_layer2,label=label,ls="-",lw=1.5,color="blue",
             marker='o',alpha=1,markersize=9,markeredgewidth=1.5, markeredgecolor="blue", markerfacecolor='None')
    #
    ni_layer3 = np.array(df_ni_Jz90['dy4'] /2 + df_ni_Jz90["dy5"] / 2)
    label = r"$\mu =3$"
    L90_mu3, =ax1.plot(df_ni_Jz90["rmean"].values,ni_layer3,label=label,ls="-",lw=1.5,color="magenta",
             marker='o',alpha=1,markersize=9,markeredgewidth=1.5, markeredgecolor="magenta", markerfacecolor='None')
    ####图例设置
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 16, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[L90_mu1,L90_mu2,L90_mu3], loc = 4, bbox_to_anchor=(0.56, 0.46),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)    
    #legend2=plt.legend(handles=[L04g,L10g,L20g], loc = 4, bbox_to_anchor=(0.64, 0.02),
    #                   ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    #plt.gca().add_artist(legend1)#####把图例legend1重新加载回来
    #
    label_x = r"x"
    label_y = "n(x)"
    #plt.yscale("log")
    #plt.xscale("log")
    ax1.set_xlabel(label_x, size= 16)
    ax1.set_ylabel(label_y, size= 16)
    ax1.tick_params(labelsize = 15) # 设置坐标刻度对应数字的大小
    ax1.set_xlim([0,48])
    #ax1.set_ylim([0.74,0.96])
    ax1.set_xticks([0,10,20,30,40,48])
    #ax1.set_yticks([0.75,0.8,0.85,0.9,0.95]) 
    #=========================================================
    ax1.text(2,0.42,"(a)",fontsize = 20, color='black', rotation = 0)
    ax1.text(10,0.42, r'$\mathrm{J_{\bot}} = 9.0$', fontsize = 18, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
    ax1.text(25,0.42, r'$\mathrm{\delta} = 0.5$', fontsize = 18, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
    ax1.text(32,0.52, r'D=10000', fontsize = 18, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
    #=======================================================================================
    # 坐标轴设置第一层
    labels = ax1.get_xticklabels() + ax1.get_yticklabels()
    #[label.set_fontname('Times New Roman') for label in labels]###设置ticket labled的字体格式
    ax1.xaxis.set_minor_locator(MultipleLocator(5))###设置次刻度的间隔
    #ax1.yaxis.set_minor_locator(MultipleLocator(10))###设置次刻度的间隔
    ax1.xaxis.set_major_formatter(FormatStrFormatter('%1.0f'))###设置X轴标签文本格式
    #ax1.yaxis.set_major_formatter(FormatStrFormatter('%1.1f'))###设置Y轴标签文本格式
    #
    #ax1.yaxis.get_major_locator().set_params(numticks=99)
    #ax1.yaxis.get_minor_locator().set_params(numticks=99, subs=[.2,.4,.6,.8]) # 将次要刻度显示出来 
    #----将次要刻度显示出来 
    #locmin = matplotlib.ticker.LogLocator(base=10.0, subs=(0.1,0.2,0.4,0.6,0.8,1,2,4,6,8,10 )) 
    #ax1.xaxis.set_minor_locator(locmin)
    #ax1.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    #
    #=====坐标轴的第二层： 坐标轴的设置
    ax1.spines['bottom'].set_linewidth(1.5) ###设置底部坐标轴的粗细
    ax1.spines['left'].set_linewidth(1.5)   ###设置左边坐标轴的粗细
    ax1.spines['right'].set_linewidth(1.5)  ###设置右边坐标轴的粗细
    ax1.spines['top'].set_linewidth(1.5)    ###设置上部坐标轴的粗细
    #ax1.spines['right'].set_color('none')# 将右边上边的两条边颜色设置为空 其实就相当于抹掉这两条边
    #ax1.spines['top'].set_color('none')
    #====坐标轴的第三层：  主刻度线的设置
    for line in ax1.xaxis.get_ticklines():
        #line is a Line2D instance
        #line.set_color('green')
        line.set_markersize(3)####设置刻度线的长度
        line.set_markeredgewidth(1.5)####设置刻度线的宽度
    for line in ax1.yaxis.get_ticklines():
        # line is a Line2D instance
        #line.set_color('green')
        line.set_markersize(3)####设置刻度线的长度
        line.set_markeredgewidth(1.5)####设置刻度线的宽度
    plt.tick_params(axis="x", which="minor", length=2.5, width=1.5, color="k")  ### 设置次要刻度 
    plt.tick_params(axis="y", which="minor", length=2.5, width=1.5, color="k")  ### 设置次要刻度  
    #-----------------------------------------------------------------------------------------------
    # 选择子图 ax2 进行绘图
    plt.sca(ax2) ## 选择对 ax2 进行绘图
    ax2 = plt.gca()
    df_ni_Jz120 = get_data_ni(Lz=Lz,Ly=Ly,Lx=Lx,dop=144,t=t,J=J,Jz=12.0,dim=14000)
    #
    L, =ax2.plot([0,50],[0.5,0.5],label=" ",ls="--",lw=1.5,color="k",
             marker='o',alpha=1,markersize=0,markeredgewidth=1.5, markeredgecolor="red", markerfacecolor='None')
    ni_layer1 = np.array(df_ni_Jz120['dy0'] /2 + df_ni_Jz120["dy1"] / 2)
    label = r"$\mu =1$"
    L120_mu1, =ax2.plot(df_ni_Jz120["rmean"].values,ni_layer1,label=label,ls="-",lw=1.5,color="red",
             marker='o',alpha=1,markersize=9,markeredgewidth=1.5, markeredgecolor="red", markerfacecolor='None')
    #
    ni_layer2 = np.array(df_ni_Jz120['dy2'] /2 + df_ni_Jz120["dy3"] / 2)
    label = r"$\mu =2$"
    L120_mu2, =ax2.plot(df_ni_Jz120["rmean"].values,ni_layer2,label=label,ls="-",lw=1.5,color="blue",
             marker='o',alpha=1,markersize=9,markeredgewidth=1.5, markeredgecolor="blue", markerfacecolor='None')
    #
    ni_layer3 = np.array(df_ni_Jz120['dy4'] /2 + df_ni_Jz120["dy5"] / 2)
    label = r"$\mu =3$"
    L120_mu3, =ax2.plot(df_ni_Jz120["rmean"].values,ni_layer3,label=label,ls="-",lw=1.5,color="magenta",
             marker='o',alpha=1,markersize=9,markeredgewidth=1.5, markeredgecolor="magenta", markerfacecolor='None')
    ####图例设置
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 16, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[L120_mu1,L120_mu2,L120_mu3], loc = 4, bbox_to_anchor=(0.56, 0.46),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)    
    #legend2=plt.legend(handles=[L04g,L10g,L20g], loc = 4, bbox_to_anchor=(0.64, 0.02),
    #                   ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    #plt.gca().add_artist(legend1)#####把图例legend1重新加载回来
    #
    label_x = r"x"
    label_y = "n(x)"
    #plt.yscale("log")
    #plt.xscale("log")
    ax2.set_xlabel(label_x, size= 16)
    ax2.set_ylabel(label_y, size= 16)
    ax2.tick_params(labelsize = 15) # 设置坐标刻度对应数字的大小
    ax2.set_xlim([0,48])
    #ax2.set_ylim([0.74,0.96])
    ax2.set_xticks([0,10,20,30,40,48])
    #ax2.set_yticks([0.75,0.8,0.85,0.9,0.95]) 
    #=========================================================
    ax2.text(2,0.41,"(b)",fontsize = 20, color='black', rotation = 0)
    ax2.text(10,0.41, r'$\mathrm{J_{\bot}} = 12.0$', fontsize = 18, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
    ax2.text(25,0.41, r'$\mathrm{\delta} = 0.5$', fontsize = 18, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
    ax2.text(35,0.52, r'D=14000', fontsize = 18, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
    #=======================================================================================
    # 坐标轴设置第一层
    labels = ax2.get_xticklabels() + ax2.get_yticklabels()
    #[label.set_fontname('Times New Roman') for label in labels]###设置ticket labled的字体格式
    ax2.xaxis.set_minor_locator(MultipleLocator(5))###设置次刻度的间隔
    #ax2.yaxis.set_minor_locator(MultipleLocator(10))###设置次刻度的间隔
    ax2.xaxis.set_major_formatter(FormatStrFormatter('%1.0f'))###设置X轴标签文本格式
    #ax2.yaxis.set_major_formatter(FormatStrFormatter('%1.1f'))###设置Y轴标签文本格式
    #
    #ax2.yaxis.get_major_locator().set_params(numticks=99)
    #ax2.yaxis.get_minor_locator().set_params(numticks=99, subs=[.2,.4,.6,.8]) # 将次要刻度显示出来 
    #----将次要刻度显示出来 
    #locmin = matplotlib.ticker.LogLocator(base=10.0, subs=(0.1,0.2,0.4,0.6,0.8,1,2,4,6,8,10 )) 
    #ax2.xaxis.set_minor_locator(locmin)
    #ax2.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    #
    #=====坐标轴的第二层： 坐标轴的设置
    ax2.spines['bottom'].set_linewidth(1.5) ###设置底部坐标轴的粗细
    ax2.spines['left'].set_linewidth(1.5)   ###设置左边坐标轴的粗细
    ax2.spines['right'].set_linewidth(1.5)  ###设置右边坐标轴的粗细
    ax2.spines['top'].set_linewidth(1.5)    ###设置上部坐标轴的粗细
    #ax2.spines['right'].set_color('none')# 将右边上边的两条边颜色设置为空 其实就相当于抹掉这两条边
    #ax2.spines['top'].set_color('none')
    #====坐标轴的第三层：  主刻度线的设置
    for line in ax2.xaxis.get_ticklines():
        #line is a Line2D instance
        #line.set_color('green')
        line.set_markersize(3)####设置刻度线的长度
        line.set_markeredgewidth(1.5)####设置刻度线的宽度
    for line in ax2.yaxis.get_ticklines():
        # line is a Line2D instance
        #line.set_color('green')
        line.set_markersize(3)####设置刻度线的长度
        line.set_markeredgewidth(1.5)####设置刻度线的宽度
    plt.tick_params(axis="x", which="minor", length=2.5, width=1.5, color="k")  ### 设置次要刻度 
    plt.tick_params(axis="y", which="minor", length=2.5, width=1.5, color="k")  ### 设置次要刻度     
    #-----------------------------------------------------------------------------------------------
    # 选择子图 ax3 进行绘图
    plt.sca(ax3) ## 选择对 ax3 进行绘图
    ax3 = plt.gca()
    df_ni_Jz150 = get_data_ni(Lz=Lz,Ly=Ly,Lx=Lx,dop=144,t=t,J=J,Jz=15.0,dim=10000)
    #
    L, =ax3.plot([0,50],[0.5,0.5],label=" ",ls="--",lw=1.5,color="k",
             marker='o',alpha=1,markersize=0,markeredgewidth=1.5, markeredgecolor="red", markerfacecolor='None')
    ni_layer1 = np.array(df_ni_Jz150['dy0'] /2 + df_ni_Jz150["dy1"] / 2)
    label = r"$\mu =1$"
    L150_mu1, =ax3.plot(df_ni_Jz150["rmean"].values,ni_layer1,label=label,ls="-",lw=1.5,color="red",
             marker='o',alpha=1,markersize=9,markeredgewidth=1.5, markeredgecolor="red", markerfacecolor='None')
    #
    ni_layer2 = np.array(df_ni_Jz150['dy2'] /2 + df_ni_Jz150["dy3"] / 2)
    label = r"$\mu =2$"
    L150_mu2, =ax3.plot(df_ni_Jz150["rmean"].values,ni_layer2,label=label,ls="-",lw=1.5,color="blue",
             marker='o',alpha=1,markersize=9,markeredgewidth=1.5, markeredgecolor="blue", markerfacecolor='None')
    #
    ni_layer3 = np.array(df_ni_Jz150['dy4'] /2 + df_ni_Jz150["dy5"] / 2)
    label = r"$\mu =3$"
    L150_mu3, =ax3.plot(df_ni_Jz150["rmean"].values,ni_layer3,label=label,ls="-",lw=1.5,color="magenta",
             marker='o',alpha=1,markersize=9,markeredgewidth=1.5, markeredgecolor="magenta", markerfacecolor='None')
    ####图例设置
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 16, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[L150_mu1,L150_mu2,L150_mu3], loc = 4, bbox_to_anchor=(0.56, 0.46),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)    
    #legend2=plt.legend(handles=[L04g,L10g,L20g], loc = 4, bbox_to_anchor=(0.64, 0.02),
    #                   ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    #plt.gca().add_artist(legend1)#####把图例legend1重新加载回来
    #
    label_x = r"x"
    label_y = "n(x)"
    #plt.yscale("log")
    #plt.xscale("log")
    ax3.set_xlabel(label_x, size= 16)
    ax3.set_ylabel(label_y, size= 16)
    ax3.tick_params(labelsize = 15) # 设置坐标刻度对应数字的大小
    ax3.set_xlim([0,48])
    #ax3.set_ylim([0.74,0.96])
    ax3.set_xticks([0,10,20,30,40,48])
    #ax3.set_yticks([0.75,0.8,0.85,0.9,0.95]) 
    #=========================================================
    ax3.text(2,0.38,"(c)",fontsize = 20, color='black', rotation = 0)
    ax3.text(10,0.38, r'$\mathrm{J_{\bot}} = 15.0$', fontsize = 18, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
    ax3.text(25,0.38, r'$\mathrm{\delta} = 0.5$', fontsize = 18, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
    ax3.text(35,0.52, r'D=10000', fontsize = 18, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
    #=======================================================================================
    # 坐标轴设置第一层
    labels = ax3.get_xticklabels() + ax3.get_yticklabels()
    #[label.set_fontname('Times New Roman') for label in labels]###设置ticket labled的字体格式
    ax3.xaxis.set_minor_locator(MultipleLocator(5))###设置次刻度的间隔
    #ax3.yaxis.set_minor_locator(MultipleLocator(10))###设置次刻度的间隔
    ax3.xaxis.set_major_formatter(FormatStrFormatter('%1.0f'))###设置X轴标签文本格式
    #ax3.yaxis.set_major_formatter(FormatStrFormatter('%1.1f'))###设置Y轴标签文本格式
    #
    #ax3.yaxis.get_major_locator().set_params(numticks=99)
    #ax3.yaxis.get_minor_locator().set_params(numticks=99, subs=[.2,.4,.6,.8]) # 将次要刻度显示出来 
    #----将次要刻度显示出来 
    #locmin = matplotlib.ticker.LogLocator(base=10.0, subs=(0.1,0.2,0.4,0.6,0.8,1,2,4,6,8,10 )) 
    #ax3.xaxis.set_minor_locator(locmin)
    #ax3.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    #
    #=====坐标轴的第二层： 坐标轴的设置
    ax3.spines['bottom'].set_linewidth(1.5) ###设置底部坐标轴的粗细
    ax3.spines['left'].set_linewidth(1.5)   ###设置左边坐标轴的粗细
    ax3.spines['right'].set_linewidth(1.5)  ###设置右边坐标轴的粗细
    ax3.spines['top'].set_linewidth(1.5)    ###设置上部坐标轴的粗细
    #ax3.spines['right'].set_color('none')# 将右边上边的两条边颜色设置为空 其实就相当于抹掉这两条边
    #ax3.spines['top'].set_color('none')
    #====坐标轴的第三层：  主刻度线的设置
    for line in ax3.xaxis.get_ticklines():
        #line is a Line2D instance
        #line.set_color('green')
        line.set_markersize(3)####设置刻度线的长度
        line.set_markeredgewidth(1.5)####设置刻度线的宽度
    for line in ax3.yaxis.get_ticklines():
        # line is a Line2D instance
        #line.set_color('green')
        line.set_markersize(3)####设置刻度线的长度
        line.set_markeredgewidth(1.5)####设置刻度线的宽度
    plt.tick_params(axis="x", which="minor", length=2.5, width=1.5, color="k")  ### 设置次要刻度 
    plt.tick_params(axis="y", which="minor", length=2.5, width=1.5, color="k")  ### 设置次要刻度  
    #------------------------------------------------------------------------------------------------
    # 选择子图 ax4 进行绘图
    plt.sca(ax4) ## 选择对 ax4 进行绘图
    ax4 = plt.gca()
    df_ni_Jz180 = get_data_ni(Lz=Lz,Ly=Ly,Lx=Lx,dop=144,t=t,J=J,Jz=18.0,dim=10000)
    #
    L, =ax4.plot([0,50],[0.5,0.5],label=" ",ls="--",lw=1.5,color="k",
             marker='o',alpha=1,markersize=0,markeredgewidth=1.5, markeredgecolor="red", markerfacecolor='None')
    ni_layer1 = np.array(df_ni_Jz180['dy0'] /2 + df_ni_Jz180["dy1"] / 2)
    label = r"$\mu =1$"
    L180_mu1, =ax4.plot(df_ni_Jz180["rmean"].values,ni_layer1,label=label,ls="-",lw=1.5,color="red",
             marker='o',alpha=1,markersize=9,markeredgewidth=1.5, markeredgecolor="red", markerfacecolor='None')
    #
    ni_layer2 = np.array(df_ni_Jz180['dy2'] /2 + df_ni_Jz180["dy3"] / 2)
    label = r"$\mu =2$"
    L180_mu2, =ax4.plot(df_ni_Jz180["rmean"].values,ni_layer2,label=label,ls="-",lw=1.5,color="blue",
             marker='o',alpha=1,markersize=9,markeredgewidth=1.5, markeredgecolor="blue", markerfacecolor='None')
    #
    ni_layer3 = np.array(df_ni_Jz180['dy4'] /2 + df_ni_Jz180["dy5"] / 2)
    label = r"$\mu =3$"
    L180_mu3, =ax4.plot(df_ni_Jz180["rmean"].values,ni_layer3,label=label,ls="-",lw=1.5,color="magenta",
             marker='o',alpha=1,markersize=9,markeredgewidth=1.5, markeredgecolor="magenta", markerfacecolor='None')
    ####图例设置
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 16, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[L180_mu1,L180_mu2,L180_mu3], loc = 4, bbox_to_anchor=(0.56, 0.46),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)    
    #legend2=plt.legend(handles=[L04g,L10g,L20g], loc = 4, bbox_to_anchor=(0.64, 0.02),
    #                   ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    #plt.gca().add_artist(legend1)#####把图例legend1重新加载回来
    #
    label_x = r"x"
    label_y = "n(x)"
    #plt.yscale("log")
    #plt.xscale("log")
    ax4.set_xlabel(label_x, size= 16)
    ax4.set_ylabel(label_y, size= 16)
    ax4.tick_params(labelsize = 15) # 设置坐标刻度对应数字的大小
    ax4.set_xlim([0,48])
    #ax4.set_ylim([0.74,0.96])
    ax4.set_xticks([0,10,20,30,40,48])
    #ax4.set_yticks([0.75,0.8,0.85,0.9,0.95]) 
    #=========================================================
    ax4.text(2,0.36,"(d)",fontsize = 20, color='black', rotation = 0)
    ax4.text(10,0.36, r'$\mathrm{J_{\bot}} = 18.0$', fontsize = 18, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
    ax4.text(25,0.36, r'$\mathrm{\delta} = 0.5$', fontsize = 18, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
    ax4.text(35,0.52, r'D=10000', fontsize = 18, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
    #=======================================================================================
    # 坐标轴设置第一层
    labels = ax4.get_xticklabels() + ax4.get_yticklabels()
    #[label.set_fontname('Times New Roman') for label in labels]###设置ticket labled的字体格式
    ax4.xaxis.set_minor_locator(MultipleLocator(5))###设置次刻度的间隔
    #ax4.yaxis.set_minor_locator(MultipleLocator(10))###设置次刻度的间隔
    ax4.xaxis.set_major_formatter(FormatStrFormatter('%1.0f'))###设置X轴标签文本格式
    #ax4.yaxis.set_major_formatter(FormatStrFormatter('%1.1f'))###设置Y轴标签文本格式
    #
    #ax4.yaxis.get_major_locator().set_params(numticks=99)
    #ax4.yaxis.get_minor_locator().set_params(numticks=99, subs=[.2,.4,.6,.8]) # 将次要刻度显示出来 
    #----将次要刻度显示出来 
    #locmin = matplotlib.ticker.LogLocator(base=10.0, subs=(0.1,0.2,0.4,0.6,0.8,1,2,4,6,8,10 )) 
    #ax4.xaxis.set_minor_locator(locmin)
    #ax4.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    #
    #=====坐标轴的第二层： 坐标轴的设置
    ax4.spines['bottom'].set_linewidth(1.5) ###设置底部坐标轴的粗细
    ax4.spines['left'].set_linewidth(1.5)   ###设置左边坐标轴的粗细
    ax4.spines['right'].set_linewidth(1.5)  ###设置右边坐标轴的粗细
    ax4.spines['top'].set_linewidth(1.5)    ###设置上部坐标轴的粗细
    #ax4.spines['right'].set_color('none')# 将右边上边的两条边颜色设置为空 其实就相当于抹掉这两条边
    #ax4.spines['top'].set_color('none')
    #====坐标轴的第三层：  主刻度线的设置
    for line in ax4.xaxis.get_ticklines():
        #line is a Line2D instance
        #line.set_color('green')
        line.set_markersize(3)####设置刻度线的长度
        line.set_markeredgewidth(1.5)####设置刻度线的宽度
    for line in ax4.yaxis.get_ticklines():
        # line is a Line2D instance
        #line.set_color('green')
        line.set_markersize(3)####设置刻度线的长度
        line.set_markeredgewidth(1.5)####设置刻度线的宽度
    plt.tick_params(axis="x", which="minor", length=2.5, width=1.5, color="k")  ### 设置次要刻度 
    plt.tick_params(axis="y", which="minor", length=2.5, width=1.5, color="k")  ### 设置次要刻度  
    #--------------------------------------------------------------------------------------------------
    # 选择子图 ax6 进行绘图
    plt.sca(ax6) ## 选择对 ax6 进行绘图
    ax6 = plt.gca()
    df_ni_Jz480 = get_data_ni(Lz=Lz,Ly=Ly,Lx=Lx,dop=144,t=t,J=J,Jz=24.0,dim=10000)
    #
    L, =ax6.plot([0,50],[0.5,0.5],label=" ",ls="--",lw=1.5,color="k",
             marker='o',alpha=1,markersize=0,markeredgewidth=1.5, markeredgecolor="red", markerfacecolor='None')
    ni_layer1 = np.array(df_ni_Jz480['dy0'] /2 + df_ni_Jz480["dy1"] / 2)
    label = r"$\mu =1$"
    L480_mu1, =ax6.plot(df_ni_Jz480["rmean"].values,ni_layer1,label=label,ls="-",lw=1.5,color="red",
             marker='o',alpha=1,markersize=9,markeredgewidth=1.5, markeredgecolor="red", markerfacecolor='None')
    #
    ni_layer2 = np.array(df_ni_Jz480['dy2'] /2 + df_ni_Jz480["dy3"] / 2)
    label = r"$\mu =2$"
    L480_mu2, =ax6.plot(df_ni_Jz480["rmean"].values,ni_layer2,label=label,ls="-",lw=1.5,color="blue",
             marker='o',alpha=1,markersize=9,markeredgewidth=1.5, markeredgecolor="blue", markerfacecolor='None')
    #
    ni_layer3 = np.array(df_ni_Jz480['dy4'] /2 + df_ni_Jz480["dy5"] / 2)
    label = r"$\mu =3$"
    L480_mu3, =ax6.plot(df_ni_Jz480["rmean"].values,ni_layer3,label=label,ls="-",lw=1.5,color="magenta",
             marker='o',alpha=1,markersize=9,markeredgewidth=1.5, markeredgecolor="magenta", markerfacecolor='None')
    ####图例设置
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 16, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[L480_mu1,L480_mu2,L480_mu3], loc = 4, bbox_to_anchor=(0.56, 0.4),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)    
    #legend2=plt.legend(handles=[L04g,L10g,L20g], loc = 4, bbox_to_anchor=(0.64, 0.02),
    #                   ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    #plt.gca().add_artist(legend1)#####把图例legend1重新加载回来
    #
    label_x = r"x"
    label_y = "n(x)"
    #plt.yscale("log")
    #plt.xscale("log")
    ax6.set_xlabel(label_x, size= 16)
    ax6.set_ylabel(label_y, size= 16)
    ax6.tick_params(labelsize = 15) # 设置坐标刻度对应数字的大小
    ax6.set_xlim([0,48])
    #ax6.set_ylim([0.74,0.96])
    ax6.set_xticks([0,10,20,30,40,48])
    #ax6.set_yticks([0.75,0.8,0.85,0.9,0.95]) 
    #=========================================================
    ax6.text(2,0.32,"(f)",fontsize = 20, color='black', rotation = 0)
    ax6.text(10,0.32, r'$\mathrm{J_{\bot}} = 24.0$', fontsize = 18, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
    ax6.text(25,0.32, r'$\mathrm{\delta} = 0.50$', fontsize = 18, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
    ax6.text(10,0.46, r'D=10000', fontsize = 18, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
    #=======================================================================================
    # 坐标轴设置第一层
    labels = ax6.get_xticklabels() + ax6.get_yticklabels()
    #[label.set_fontname('Times New Roman') for label in labels]###设置ticket labled的字体格式
    ax6.xaxis.set_minor_locator(MultipleLocator(5))###设置次刻度的间隔
    #ax6.yaxis.set_minor_locator(MultipleLocator(10))###设置次刻度的间隔
    ax6.xaxis.set_major_formatter(FormatStrFormatter('%1.0f'))###设置X轴标签文本格式
    #ax6.yaxis.set_major_formatter(FormatStrFormatter('%1.1f'))###设置Y轴标签文本格式
    #
    #ax6.yaxis.get_major_locator().set_params(numticks=99)
    #ax6.yaxis.get_minor_locator().set_params(numticks=99, subs=[.2,.4,.6,.8]) # 将次要刻度显示出来 
    #----将次要刻度显示出来 
    #locmin = matplotlib.ticker.LogLocator(base=10.0, subs=(0.1,0.2,0.4,0.6,0.8,1,2,4,6,8,10 )) 
    #ax6.xaxis.set_minor_locator(locmin)
    #ax6.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    #
    #=====坐标轴的第二层： 坐标轴的设置
    ax6.spines['bottom'].set_linewidth(1.5) ###设置底部坐标轴的粗细
    ax6.spines['left'].set_linewidth(1.5)   ###设置左边坐标轴的粗细
    ax6.spines['right'].set_linewidth(1.5)  ###设置右边坐标轴的粗细
    ax6.spines['top'].set_linewidth(1.5)    ###设置上部坐标轴的粗细
    #ax6.spines['right'].set_color('none')# 将右边上边的两条边颜色设置为空 其实就相当于抹掉这两条边
    #ax6.spines['top'].set_color('none')
    #====坐标轴的第三层：  主刻度线的设置
    for line in ax6.xaxis.get_ticklines():
        #line is a Line2D instance
        #line.set_color('green')
        line.set_markersize(3)####设置刻度线的长度
        line.set_markeredgewidth(1.5)####设置刻度线的宽度
    for line in ax6.yaxis.get_ticklines():
        # line is a Line2D instance
        #line.set_color('green')
        line.set_markersize(3)####设置刻度线的长度
        line.set_markeredgewidth(1.5)####设置刻度线的宽度
    plt.tick_params(axis="x", which="minor", length=2.5, width=1.5, color="k")  ### 设置次要刻度 
    plt.tick_params(axis="y", which="minor", length=2.5, width=1.5, color="k")  ### 设置次要刻度   
    #------------------------------------------
    #fig.tight_layout() # 自动调整 subplot 间的间隙参数
    #plt.savefig("E:\\WORK\\Work\\Project\\TrilayerW2-tj-Jz-openz\\datashow_datpy_report\\figs\\fig_05000_ni.eps",
    #            dpi=300, format='eps',bbox_inches='tight') # 白边紧凑型
    plt.show()