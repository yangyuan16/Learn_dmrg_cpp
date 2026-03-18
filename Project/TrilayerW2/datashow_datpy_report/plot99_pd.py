import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
#
if __name__ =="__main__":
    #---------------plot logr-logr fig-----------------------------
    fig = plt.figure(figsize=(8,6))
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
    #ax1 = plt.axes([0.12,0.57,0.8,0.4])
    ax2 = plt.axes([0.1,0.1,0.8,0.8])
    #ax3 = plt.axes([0.1,0.1,0.37,0.375])
    #ax4 = plt.axes([0.58,0.1,0.37,0.375])
    # plt.subplots_adjust(left=0.1, bottom=0.1, right=0.98, top=0.96, wspace=0.32, hspace=0.26)
    #-------------------------------------------------------------------------------------------
    # 选择子图 ax2 进行绘图
    plt.sca(ax2) ## 选择对 ax2 进行绘图
    ax2 = plt.gca()
    #------- J_\bot = 0.0
    ####图例设置
    #legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 15, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    #legend1=plt.legend(handles=[L04,L06,L08,], loc = 4, bbox_to_anchor=(0.99, 0.778),
    #                   ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    #legend2=plt.legend(handles=[Lyy10_36,Lyy10_48,Lyy10_72], loc = 4, bbox_to_anchor=(0.64, -0.01),
    #                   ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    #plt.gca().add_artist(legend1)#####把图例legend1重新加载回来
    label_x = r"$\delta$"
    label_y = r"$J_{\bot}$"
    #plt.yscale("log")
    #plt.xscale("log")
    #ax2.set_yscale("log",base=10,subs=[0.01,0.02,0.03])
    #ax2.set_xscale("log",base=10,subs=[0.01,0.02,0.03])       
    ax2.set_xlabel(label_x, size= 16)
    ax2.set_ylabel(label_y, size= 16)
    ax2.tick_params(labelsize = 14) # 设置坐标刻度对应数字的大小
    ax2.set_xlim([0.1,0.8])
    ax2.set_ylim([-0.02,24])
    ax2.set_xticks([0.125,0.25,0.375,0.5,0.6667,0.75])
    ax2.set_yticks([0,4,8,12,16,20,24])  
    #=========================================================
    #ax2.text(0.125,1.15,"(b)",fontsize = 20, color='black', rotation = 0)
    #ax2.text(2,0.5*1.0e-6, r'$\mathrm{\delta} = 0.375$', fontsize = 16, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
    #ax2.text(3,1.0e-9, r'$J_{\bot} = 1.0$', fontsize = 16, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
    #ax2.text(0.2,0.65, r'CDW', fontsize = 18, fontdict={'family' : 'Times New Roman'},color='olive', rotation = 0)
    #ax2.text(0.33,0.65, r'LEL', fontsize = 18, fontdict={'family' : 'Times New Roman'},color='olive', rotation = 0)
    #ax2.text(0.2,1.1, r'ZZ-SC', fontsize = 18, fontdict={'family' : 'Times New Roman'},color='olive', rotation = 0)
    #=========================================================
    # 坐标轴设置第一层
    labels = ax2.get_xticklabels() + ax2.get_yticklabels()
    #[label.set_fontname('Times New Roman') for label in labels]###设置ticket labled的字体格式
    ax2.xaxis.set_minor_locator(MultipleLocator(5))###设置次刻度的间隔
    #ax2.yaxis.set_minor_locator(MultipleLocator(10))###设置次刻度的间隔
    ax2.xaxis.set_major_formatter(FormatStrFormatter('%1.4f'))###设置X轴标签文本格式
    #ax2.yaxis.set_major_formatter(FormatStrFormatter('%1.1f'))###设置Y轴标签文本格式
    #
    #ax2.yaxis.get_major_locator().set_params(numticks=99)
    #ax2.yaxis.get_minor_locator().set_params(numticks=99, subs=[.2,.4,.6,.8]) # 将次要刻度显示出来 
    #----将次要刻度显示出来 
    #locmin = matplotlib.ticker.LogLocator(base=10.0, subs=(0.1,0.2,0.4,0.6,0.8,1,2,4,6,8,10 )) 
    #ax2.xaxis.set_minor_locator(locmin)
    #ax2.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
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
    #
    #fig.tight_layout() # 自动调整 subplot 间的间隙参数
    #plt.savefig("E:\\WORK\\Work\\Project\\TrilayerW2-tj-Jz-openz\\datashow_datpy_report\\figs\\fig_pd.eps",
    #            dpi=300, format='eps',bbox_inches='tight') # 白边紧凑型
    plt.show()
    