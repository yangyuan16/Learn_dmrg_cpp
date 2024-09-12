'''
    #-----------------------------------------------------------------------------------------------
    # 选择子图 ax2 进行绘图
    plt.sca(ax2) ## 选择对 ax1 进行绘图
    ax2 = plt.gca()
    #--------- J_\bot = 1.0
    slope = df_pyy_Jz10["slope_pow"].values[0]
    K_F = round(-slope,2)  
    label = r"$K_{SC}^{yy}$=%.2f"%(K_F)
    Lyy10, = ax2.plot(df_pyy_Jz10["r"],df_pyy_Jz10["corre_abs"],label=label,ls="-",lw=2,color="red",
             marker='s',alpha=1,markersize=8,markeredgewidth=1.5, markeredgecolor="red", markerfacecolor='None')
    Lyy10_, = ax2.plot(df_pyy_Jz10["r"],df_pyy_Jz10["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='None')
    #--------- J_\bot = 1.0
    slope = df_pzz_Jz10["slope_pow"].values[0]
    K_G = round(-slope,2)  
    label = r"$K_{SC}^{zz}$=%.2f"%(K_G)
    Lzz10, = ax2.plot(df_pzz_Jz10["r"],df_pzz_Jz10["corre_abs"],label=label,ls="-",lw=2,color="blue",
             marker='s',alpha=1,markersize=8,markeredgewidth=1.5, markeredgecolor="blue", markerfacecolor='None')
    Lzz10_, = ax2.plot(df_pzz_Jz10["r"],df_pzz_Jz10["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='None')
    #-------- J_\bot = 1.0
    slope = df_ninj_Jz10["slope_pow"].values[0]
    K_F = round(-slope,2)  
    label = r"$K_{CDW}$=%.2f"%(K_F)
    Lninj10, = ax2.plot(df_ninj_Jz10["r"][1:],df_ninj_Jz10["corre_abs"][1:],label=label,ls="-",lw=2,color="green",
             marker='s',alpha=1,markersize=8,markeredgewidth=1.5, markeredgecolor="green", markerfacecolor='None')
    Lninj10_, = ax2.plot(df_ninj_Jz10["r"][1:],df_ninj_Jz10["fitcorre_pow"][1:],label=label,ls="--",lw=1.5,color="k",
             marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='None')
    ####图例设置
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 16, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[Lyy10,Lzz10,Lninj10], loc = 4, bbox_to_anchor=(0.65, 0.02),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    #legend2=plt.legend(handles=[L10,L20,L40], loc = 4, bbox_to_anchor=(0.98, 0.70),
    #                   ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    #plt.gca().add_artist(legend1)
    label_x = r"r"
    label_y = "Correlation"
    #plt.yscale("log")
    #plt.xscale("log")
    ax2.set_xscale("log",base=10,subs=[0.01,0.02,0.03]) 
    ax2.set_yscale("log",base=10,subs=[0.01,0.02,0.03])      
    ax2.set_xlabel(label_x, size= 16)
    ax2.set_ylabel(label_y, size= 16)
    ax2.tick_params(labelsize = 14) # 设置坐标刻度对应数字的大小
    ax2.set_xlim([0,40])
    #ax2.set_xticks([0,10,20,30])
    #ax2.set_yticks([-1,-0.5,0,0.5,1]) 
    #
    #=========================================================
    ax2.text(25,0.015,"(b)",fontsize = 20, color='black', rotation = 0)
    ax2.text(2,1.5e-6, r'$\mathrm{\delta} = 0.375$', fontsize = 16, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
    ax2.text(2,1e-5, r'$J_{\bot}=1.0$', fontsize = 20, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
    #=========================================================
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
    plt.sca(ax3) ## 选择对 ax1 进行绘图
    ax3 = plt.gca()
    #--------- J_\bot = 2.0
    slope = df_pyy_Jz20["slope_pow"].values[0]
    K_F = round(-slope,2)  
    label = r"$K_{SC}^{yy}$=%.2f"%(K_F)
    Lyy20, = ax3.plot(df_pyy_Jz20["r"],df_pyy_Jz20["corre_abs"],label=label,ls="-",lw=2,color="red",
             marker='s',alpha=1,markersize=8,markeredgewidth=1.5, markeredgecolor="red", markerfacecolor='None')
    Lyy20_, = ax3.plot(df_pyy_Jz20["r"],df_pyy_Jz20["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='None')
    #--------- J_\bot = 2.0
    slope = df_pzz_Jz20["slope_pow"].values[0]
    K_G = round(-slope,2)  
    label = r"$K_{SC}^{zz}$=%.2f"%(K_G)
    Lzz20, = ax3.plot(df_pzz_Jz20["r"],df_pzz_Jz20["corre_abs"],label=label,ls="-",lw=2,color="blue",
             marker='s',alpha=1,markersize=8,markeredgewidth=1.5, markeredgecolor="blue", markerfacecolor='w')
    Lzz20_, = ax3.plot(df_pzz_Jz20["r"],df_pzz_Jz20["fitcorre_pow"],label=label,ls="--",lw=1.5,color="k",
             marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='None')
    #-------- J_\bot = 2.0
    slope = df_ninj_Jz20["slope_pow"].values[0]
    K_F = round(-slope,2)  
    label = r"$K_{CDW}$=%.2f"%(K_F)
    Lninj20, = ax3.plot(df_ninj_Jz20["r"][1:],df_ninj_Jz20["corre_abs"][1:],label=label,ls="-",lw=2,color="green",
             marker='s',alpha=1,markersize=8,markeredgewidth=1.5, markeredgecolor="green", markerfacecolor='None')
    Lninj20_, = ax3.plot(df_ninj_Jz20["r"][1:],df_ninj_Jz20["fitcorre_pow"][1:],label=label,ls="--",lw=1.5,color="k",
             marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='None')
    ####图例设置
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 16, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[Lyy20,Lzz20,Lninj20], loc = 4, bbox_to_anchor=(0.65, 0.02),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    #legend2=plt.legend(handles=[L10,L20,L40], loc = 4, bbox_to_anchor=(0.98, 0.70),
    #                   ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    #plt.gca().add_artist(legend1)
    label_x = r"r"
    label_y = "Correlation"
    #plt.yscale("log")
    #plt.xscale("log")
    ax3.set_xscale("log",base=10,subs=[0.01,0.02,0.03]) 
    ax3.set_yscale("log",base=10,subs=[0.01,0.02,0.03])      
    ax3.set_xlabel(label_x, size= 16)
    ax3.set_ylabel(label_y, size= 16)
    ax3.tick_params(labelsize = 14) # 设置坐标刻度对应数字的大小
    ax3.set_xlim([0,40])
    #ax3.set_xticks([0,10,20,30])
    #ax3.set_yticks([-1,-0.5,0,0.5,1]) 
    #
    #=========================================================
    ax3.text(25,0.01,"(c)",fontsize = 20, color='black', rotation = 0)
    ax3.text(2,0.15e-7, r'$\mathrm{\delta} = 0.375$', fontsize = 16, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
    ax3.text(2,0.3e-6, r'$J_{\bot}=2.0$', fontsize = 20, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
    #=========================================================
    # 坐标轴设置第一层
    labels = ax3.get_xticklabels() + ax3.get_yticklabels()
    #[label.set_fontname('Times New Roman') for label in labels]###设置ticket labled的字体格式
    ax3.xaxis.set_minor_locator(MultipleLocator(5))###设置次刻度的间隔
    #ax3.yaxis.set_minor_locator(MultipleLocator(10))###设置次刻度的间隔
    ax3.xaxis.set_major_formatter(FormatStrFormatter('%1.0f'))###设置X轴标签文本格式
    #ax3.yaxis.set_major_formatter(FormatStrFormatter('%1.1f'))###设置Y轴标签文本格式
    #
    #ax3.yaxis.get_major_locator().set_params(numticks=99)
    #ax3.yaxis.get_minor_locator().set_params(numticks=13, subs=[.2,.4,.6,.8]) # 将次要刻度显示出来 
    #----将次要刻度显示出来 
    #locmin = matplotlib.ticker.LogLocator(base=10.0, subs=(0.1,0.2,0.4,0.6,0.8,1,2,4,6,8,10 )) 
    #ax3.xaxis.set_minor_locator(locmin)
    #ax3.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
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
    #========================================================
    # 选择子图 ax4 进行绘图
    plt.sca(ax4) ## 选择对 ax1 进行绘图
    ax4 = plt.gca()

    LK2, = ax4.plot([0,2.5],[2,2],label=label,ls="--",lw=2,color="k",
             marker='s',alpha=1,markersize=0,markeredgewidth=1, markeredgecolor="green", markerfacecolor='None')
    label = "$K_{SC}^{yy}$"
    Lyy, = ax4.plot(Jz_list,K_sc_yy,label=label,ls="-",lw=2,color="red",
             marker='s',alpha=1,markersize=8,markeredgewidth=1.5, markeredgecolor="red", markerfacecolor='None')
    label = "$K_{SC}^{zz}$"
    Lzz, = ax4.plot(Jz_list,K_sc_zz,label=label,ls="-",lw=2,color="blue",
             marker='s',alpha=1,markersize=8,markeredgewidth=1.5, markeredgecolor="blue", markerfacecolor='None')
    label = "$K_{CDW}$"
    Lninj, = ax4.plot(Jz_list,K_cdw,label=label,ls="-",lw=2,color="green",
             marker='s',alpha=1,markersize=8,markeredgewidth=1.5, markeredgecolor="green", markerfacecolor='None')
    ####图例设置
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 16, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[Lyy,Lzz,Lninj], loc = 4, bbox_to_anchor=(0.95, 0.38),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    label_x = r"$J_{\bot}$"
    label_y = "K"
    #plt.yscale("log")
    #plt.xscale("log")
    #ax4.set_xscale("log",base=10,subs=[0.01,0.02,0.03]) 
    #ax4.set_yscale("log",base=10,subs=[0.01,0.02,0.03])      
    ax4.set_xlabel(label_x, size= 18)
    ax4.set_ylabel(label_y, size= 18)
    ax4.tick_params(labelsize = 16) # 设置坐标刻度对应数字的大小
    ax4.set_xlim([0,2.5])
    ax4.set_ylim([0,7])
    #ax4.set_xticks([0,10,20,30])
    #ax4.set_yticks([-1,-0.5,0,0.5,1]) 
    #
    #=========================================================
    ax4.text(2.1,6,"(d)",fontsize = 20, color='black', rotation = 0)
    ax4.text(0.8,6, r'$\mathrm{\delta} = 0.375$', fontsize = 16, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
    #=========================================================
    #======== 坐标轴设置第一层
    labels = ax4.get_xticklabels() + ax4.get_yticklabels()
    #[label.set_fontname('Times New Roman') for label in labels]###设置ticket labled的字体格式
    ax4.xaxis.set_minor_locator(MultipleLocator(0.25))###设置次刻度的间隔
    ax4.yaxis.set_minor_locator(MultipleLocator(0.5))###设置次刻度的间隔
    ax4.xaxis.set_major_formatter(FormatStrFormatter('%1.1f'))###设置X轴标签文本格式
    #ax4.yaxis.set_major_formatter(FormatStrFormatter('%1.1f'))###设置Y轴标签文本格式
    #
    #ax4.yaxis.get_major_locator().set_params(numticks=99)
    #ax4.yaxis.get_minor_locator().set_params(numticks=99, subs=[.2,.4,.6,.8]) # 将次要刻度显示出来 
    #----将次要刻度显示出来 
    #locmin = matplotlib.ticker.LogLocator(base=10.0, subs=(0.1,0.2,0.4,0.6,0.8,1,2,4,6,8,10 )) 
    #ax4.xaxis.set_minor_locator(locmin)
    #ax4.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
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
    #========================================================
    plt.sca(ax5) ## 选择对 ax1 进行绘图
    ax5 = plt.gca()
    #
    df_ni_Jz05 = get_data_ni(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=0.5,dim=dim)
    df_ni_Jz10 = get_data_ni(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=1.0,dim=dim)
    df_ni_Jz20 = get_data_ni(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=2.0,dim=dim)
    
    L, =ax5.plot([0,50],[0.625,0.625],label=" ",ls="--",lw=1.5,color="k",
             marker='o',alpha=1,markersize=0,markeredgewidth=1.5, markeredgecolor="red", markerfacecolor='None')
    label = r"$J_{\bot}$=0.5"
    L05, =ax5.plot(df_ni_Jz05["rmean"].values,df_ni_Jz05["dymean"].values,label=label,ls="--",lw=1.5,color="red",
             marker='o',alpha=1,markersize=9,markeredgewidth=1.5, markeredgecolor="red", markerfacecolor='None')
    label = r"$J_{\bot}$=1.0"
    L10, =ax5.plot(df_ni_Jz10["rmean"].values,df_ni_Jz10["dymean"].values,label=label,ls="--",lw=1.5,color="blue",
             marker='^',alpha=1,markersize=9,markeredgewidth=1.5, markeredgecolor="blue", markerfacecolor='None')
    label = r"$J_{\bot}$=2.0"
    L20, =ax5.plot(df_ni_Jz20["rmean"].values,df_ni_Jz20["dymean"].values,label=label,ls="--",lw=1.5,color="green",
             marker='v',alpha=1,markersize=9,markeredgewidth=1.5, markeredgecolor="green", markerfacecolor='None')
    ####图例设置
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 15, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[L05,L10,L20,], loc = 4, bbox_to_anchor=(0.88, 0.01),
                       ncol = 4,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)    
    #
    label_x = r"x"
    label_y = "n(x)"
    #plt.yscale("log")
    #plt.xscale("log")
    ax5.set_xlabel(label_x, size= 16)
    ax5.set_ylabel(label_y, size= 16)
    ax5.tick_params(labelsize = 15) # 设置坐标刻度对应数字的大小
    ax5.set_xlim([0,50])
    ax5.set_xticks([0,10,20,30,40,50])
    #ax5.set_yticks([-1,-0.5,0,0.5,1]) 
    #=========================================================
    ax5.text(4,0.675,"(e)",fontsize = 20, color='black', rotation = 0)
    ax5.text(20,0.675, r'$\mathrm{\delta} = 0.375$', fontsize = 16, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
    #=======================================================================================
    # 坐标轴设置第一层
    labels = ax5.get_xticklabels() + ax5.get_yticklabels()
    #[label.set_fontname('Times New Roman') for label in labels]###设置ticket labled的字体格式
    ax5.xaxis.set_minor_locator(MultipleLocator(5))###设置次刻度的间隔
    #ax5.yaxis.set_minor_locator(MultipleLocator(10))###设置次刻度的间隔
    ax5.xaxis.set_major_formatter(FormatStrFormatter('%1.0f'))###设置X轴标签文本格式
    #ax5.yaxis.set_major_formatter(FormatStrFormatter('%1.1f'))###设置Y轴标签文本格式
    #
    #ax5.yaxis.get_major_locator().set_params(numticks=99)
    #ax5.yaxis.get_minor_locator().set_params(numticks=99, subs=[.2,.4,.6,.8]) # 将次要刻度显示出来 
    #----将次要刻度显示出来 
    #locmin = matplotlib.ticker.LogLocator(base=10.0, subs=(0.1,0.2,0.4,0.6,0.8,1,2,4,6,8,10 )) 
    #ax5.xaxis.set_minor_locator(locmin)
    #ax5.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    #
    #=====坐标轴的第二层： 坐标轴的设置
    ax5.spines['bottom'].set_linewidth(1.5) ###设置底部坐标轴的粗细
    ax5.spines['left'].set_linewidth(1.5)   ###设置左边坐标轴的粗细
    ax5.spines['right'].set_linewidth(1.5)  ###设置右边坐标轴的粗细
    ax5.spines['top'].set_linewidth(1.5)    ###设置上部坐标轴的粗细
    #ax5.spines['right'].set_color('none')# 将右边上边的两条边颜色设置为空 其实就相当于抹掉这两条边
    #ax5.spines['top'].set_color('none')
    #====坐标轴的第三层：  主刻度线的设置
    for line in ax5.xaxis.get_ticklines():
        #line is a Line2D instance
        #line.set_color('green')
        line.set_markersize(3)####设置刻度线的长度
        line.set_markeredgewidth(1.5)####设置刻度线的宽度
    for line in ax5.yaxis.get_ticklines():
        # line is a Line2D instance
        #line.set_color('green')
        line.set_markersize(3)####设置刻度线的长度
        line.set_markeredgewidth(1.5)####设置刻度线的宽度
    plt.tick_params(axis="x", which="minor", length=2.5, width=1.5, color="k")  ### 设置次要刻度 
    plt.tick_params(axis="y", which="minor", length=2.5, width=1.5, color="k")  ### 设置次要刻度 
    #=========================================================================
    plt.sca(ax6) ## 选择对 ax1 进行绘图
    ax6 = plt.gca()
    #
    #
    df_sisj_Jz04 = get_data_sisj(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=0.4,dim=dim)
    df_sisj_Jz10 = get_data_sisj(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=1.0,dim=dim)
    df_sisj_Jz20 = get_data_sisj(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=2.0,dim=dim)
    # load green function data
    df_cicj_Jz04 = get_data_cicj(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=0.4,dim=dim)
    df_cicj_Jz10 = get_data_cicj(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=1.0,dim=dim)
    df_cicj_Jz20 = get_data_cicj(Lz=Lz,Ly=Ly,Lx=Lx,dop=dop,t=t,J=J,Jz=2.0,dim=dim)
    #
    #--------- J_\bot = 0.4
    slope = df_sisj_Jz04["slope_exp"].values[0]
    xi = round(-1/slope,2)  
    label = r"$J_{\bot}$=%.1f $\xi_F$=%.2f"%(0.4,xi)
    L04, = ax6.plot(df_sisj_Jz04["r"],df_sisj_Jz04["corre_abs"],label=label,ls="-",lw=1.5,color="blue",
             marker='s',alpha=1,markersize=6,markeredgewidth=1, markeredgecolor="blue", markerfacecolor='None')
    L04_, = ax6.plot(df_sisj_Jz04["r"],df_sisj_Jz04["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
             marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 1.0
    slope = df_sisj_Jz10["slope_exp"].values[0]
    xi = round(-1/slope,2)  
    label = r"$J_{\bot}$=%.1f $\xi_F$=%.2f"%(1.0,xi)
    L10, = ax6.plot(df_sisj_Jz10["r"],df_sisj_Jz10["corre_abs"],label=label,ls="-",lw=1.5,color="cyan",
             marker='+',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="cyan", markerfacecolor='None')
    L10_, = ax6.plot(df_sisj_Jz10["r"],df_sisj_Jz10["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
             marker='+',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 2.0
    slope = df_sisj_Jz20["slope_exp"].values[0]
    xi = round(-1/slope,2)  
    label = r"$J_{\bot}$=%.1f $\xi_F$=%.2f"%(2.0,xi)
    L20, = ax6.plot(df_sisj_Jz20["r"],df_sisj_Jz20["corre_abs"],label=label,ls="-",lw=1.5,color="pink",
             marker='*',alpha=1,markersize=8,markeredgewidth=1, markeredgecolor="pink", markerfacecolor='None')
    L20_, = ax6.plot(df_sisj_Jz20["r"],df_sisj_Jz20["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
             marker='*',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #
    #-------- Green Functions  <c_i c_j> ----------------- 
    #--------- J_\bot = 0.4
    slope = df_cicj_Jz04["slope_exp"].values[0]
    xi = round(-1/slope,2)  
    label = r"$J_{\bot}$=%.1f $\xi_G$=%.2f"%(0.4,xi)
    L04g, = ax6.plot(df_cicj_Jz04["r"],df_cicj_Jz04["corre_abs"],label=label,ls="-",lw=1.5,color="red",
             marker='s',alpha=1,markersize=6,markeredgewidth=1, markeredgecolor="red", markerfacecolor='None')
    L04g_, = ax6.plot(df_cicj_Jz04["r"],df_cicj_Jz04["fitcorre_exp"],label=label,ls="--",lw=1.5,color="k",
             marker='s',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 1.0
    slope = df_cicj_Jz10["slope_exp"].values[0]
    xi = round(-1/slope,2)  
    label = r"$J_{\bot}$=%.1f $\xi_G$=%.2f"%(1.0,xi)
    L10g, = ax6.plot(df_cicj_Jz10["r"],df_cicj_Jz10["corre_abs"],label=label,ls="-",lw=1.5,color="green",
             marker='+',alpha=1,markersize=6,markeredgewidth=1, markeredgecolor="green", markerfacecolor='None')
    L10g_, = ax6.plot(df_cicj_Jz10["r"],df_cicj_Jz10["fitcorre_exp"],label=label,ls="--",lw=1,color="k",
             marker='+',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #-------- J_\bot = 2.0
    slope = df_cicj_Jz20["slope_exp"].values[0]
    xi = round(-1/slope,2)  
    label = r"$J_{\bot}$=%.1f $\xi_G$=%.2f"%(2.0,xi)
    L20g, = ax6.plot(df_cicj_Jz20["r"],df_cicj_Jz20["corre_abs"],label=label,ls="-",lw=1.5,color="magenta",
             marker='*',alpha=1,markersize=6,markeredgewidth=1, markeredgecolor="magenta", markerfacecolor='None')
    L20g_, = ax6.plot(df_cicj_Jz20["r"],df_cicj_Jz20["fitcorre_exp"],label=label,ls="--",lw=1,color="k",
             marker='*',alpha=1,markersize=0,markeredgewidth=0, markeredgecolor="k", markerfacecolor='w')
    #
    #--------
    ####图例设置
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 12, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[L04,L10,L20], loc = 4, bbox_to_anchor=(0.99, 0.778),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    
    legend2=plt.legend(handles=[L04g,L10g,L20g], loc = 4, bbox_to_anchor=(0.64, 0.02),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    
    plt.gca().add_artist(legend1)#####把图例legend1重新加载回来

    label_x = r"r"
    label_y = "|F(r)|,|G(r)|"
    #plt.yscale("log")
    #plt.xscale("log")
    ax6.set_yscale("log",base=10,subs=[0.01,0.02,0.03])      
    ax6.set_xlabel(label_x, size= 14)
    ax6.set_ylabel(label_y, size= 14)
    ax6.tick_params(labelsize = 12) # 设置坐标刻度对应数字的大小
    ax6.set_xlim([0,40])
    #ax6.set_xticks([0,10,20,30])
    #ax6.set_yticks([-1,-0.5,0,0.5,1])  
    #=========================================================
    ax6.text(4,0.06,"(f)",fontsize = 20, color='black', rotation = 0)
    ax6.text(2,1.0e-8, r'$\mathrm{\delta} = 0.375$', fontsize = 16, fontdict={'family' : 'Times New Roman'},color='black', rotation = 0)
    #=========================================================
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
    #
    plt.sca(ax7)  ##选择对ax1进行绘图
    ax7=plt.gca() #获得坐标轴的句柄
    c_04 = 0.96
    g_04 = 3.35
    df_ent_04 = get_data_entropy(Lz=2,Ly=2,Lx=48,dop=72,t=3,J=1,Jz=0.4,dim=6000)
    S_res04 = entropy(Lx=Lx,x=df_ent_04["r"].values,c=c_04,g=g_04)
    #
    c_26 = 0.94
    g_26 = 2.50
    df_ent_26 = get_data_entropy(Lz=2,Ly=2,Lx=48,dop=72,t=3,J=1,Jz=2.6,dim=6000)
    S_res26 = entropy(Lx=Lx,x=df_ent_26["r"].values,c=c_26,g=g_26)
    #
    L04, = ax7.plot(df_ent_04['r'].values,df_ent_04["entropy"],label=r"$J_{\bot}$=%.1f"%(0.4),ls="-",lw=1.5,color='r',
             marker='o',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor='r',
             markerfacecolor='None')
    L04_fit, = ax7.plot(df_ent_04['r'].values,S_res04,label="c=%.2f,g=%.2f"%(c_04,g_04),ls="--",lw=1.5,color="k",
             marker='o',alpha=1,markersize=0,markeredgewidth=1.5, markeredgecolor="k",
             markerfacecolor='None')
    #
    L26, = ax7.plot(df_ent_26['r'].values,df_ent_26["entropy"],label=r"$J_{\bot}$=%.1f"%(2.6),ls="-",lw=1.5,color='blue',
             marker='o',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor='blue',
             markerfacecolor='None')
    L26_fit, = ax7.plot(df_ent_26['r'].values,S_res26,label="c=%.2f,g=%.2f"%(c_26,g_26),ls="--",lw=1.5,color="green",
             marker='o',alpha=1,markersize=0,markeredgewidth=1.5, markeredgecolor="k",
             markerfacecolor='None')
    #
    legfont = {'family' : 'Times New Roman','weight' : 'normal','size': 16, }###图例字体的大小###ncol 设置列的数量，使显示扁平化，当要表示的线段特别多的时候会有用
    legend1=plt.legend(handles=[L04,L04_fit,L26,L26_fit], loc = 4, bbox_to_anchor=(0.86, 0.30),
                       ncol = 1,prop=legfont,markerscale=1,fancybox=None,shadow=None,frameon=False)
    label_x = r"x"
    label_y = "entropy"
    #plt.yscale("log")
    ax7.set_xlabel(label_x, size= 16)
    ax7.set_ylabel(label_y, size= 16)
    ax7.tick_params(labelsize = 16) # 设置坐标刻度对应数字的大小
    #ax7.set_xlim([0,8])
    #ax7.set_ylim([-0.1,1])
    #ax7.set_xticks([0,2,4,6,8])
    #ax7.set_yticks([-0.1,0,0.5,1])
    #ax7.text(0.3,-0.05, r'$\mathrm{(a)}$', fontsize=18)
    ax7.text(2,3.9,"(g)",fontsize = 20, color='black', rotation = 0)
    #plt.title("Jz=%.2f"%Jz,fontsize=25)
    #
    #=======================================================================================
    # 坐标轴设置第一层
    labels = ax7.get_xticklabels() + ax7.get_yticklabels()
    #[label.set_fontname('Times New Roman') for label in labels]###设置ticket labled的字体格式
    ax7.xaxis.set_minor_locator(MultipleLocator(5))###设置次刻度的间隔
    #ax7.yaxis.set_minor_locator(MultipleLocator(10))###设置次刻度的间隔
    ax7.xaxis.set_major_formatter(FormatStrFormatter('%1.0f'))###设置X轴标签文本格式
    #ax7.yaxis.set_major_formatter(FormatStrFormatter('%1.1f'))###设置Y轴标签文本格式
    #
    #ax7.yaxis.get_major_locator().set_params(numticks=99)
    #ax7.yaxis.get_minor_locator().set_params(numticks=99, subs=[.2,.4,.6,.8]) # 将次要刻度显示出来 
    #----将次要刻度显示出来 
    #locmin = matplotlib.ticker.LogLocator(base=10.0, subs=(0.1,0.2,0.4,0.6,0.8,1,2,4,6,8,10 )) 
    #ax7.xaxis.set_minor_locator(locmin)
    #ax7.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    #
    #=====坐标轴的第二层： 坐标轴的设置
    ax7.spines['bottom'].set_linewidth(1.5) ###设置底部坐标轴的粗细
    ax7.spines['left'].set_linewidth(1.5)   ###设置左边坐标轴的粗细
    ax7.spines['right'].set_linewidth(1.5)  ###设置右边坐标轴的粗细
    ax7.spines['top'].set_linewidth(1.5)    ###设置上部坐标轴的粗细
    #ax7.spines['right'].set_color('none')# 将右边上边的两条边颜色设置为空 其实就相当于抹掉这两条边
    #ax7.spines['top'].set_color('none')
    #====坐标轴的第三层：  主刻度线的设置
    for line in ax7.xaxis.get_ticklines():
        #line is a Line2D instance
        #line.set_color('green')
        line.set_markersize(3)####设置刻度线的长度
        line.set_markeredgewidth(1.5)####设置刻度线的宽度
    for line in ax7.yaxis.get_ticklines():
        # line is a Line2D instance
        #line.set_color('green')
        line.set_markersize(3)####设置刻度线的长度
        line.set_markeredgewidth(1.5)####设置刻度线的宽度
    plt.tick_params(axis="x", which="minor", length=2.5, width=1.5, color="k")  ### 设置次要刻度 
    plt.tick_params(axis="y", which="minor", length=2.5, width=1.5, color="k")  ### 设置次要刻度 
    '''