def fcfoplot(Pth,LowerLimit,UpperLimit,log=True,colortable=['r','b','g','orange','magenta','cyan','black','purple','pink']):
   
    '''
    This function plots Fo-Fc graph from csv file.
    Input values are paths to csv file, lower and upper limits for plotting (the same value is used for x and y).
    You can supply one path as list for single plot, or several paths in list file, for multiple overlaid graphs.
    The csv is expected to be in the same format as one produced by Olex2.
        INPUT:
            Pth - one or more paths to csv files AS A LIST!
            LowerLimit - minimal value of Fo/Fc to be plotted,
            UpperLimit - maximal value of Fo/Fc to be plotted,
            log - by default set to True. Change to False to not plot logarithmic scale
            colortable - you can override colors used for plotting. By default table contains 9 colours:
                         red, blue, orange, magenta, cyan, black, purple, pink
        EXAMPLES:
            Single plot:
                FcFplot(['C:/data/SomeFCF.csv'],0.01,200)
            Multiple plots:
                        FcFplot(['C:/data/SomeFCF_1.csv','C:/data/SomeFCF_2.csv'],0.01,200)
            Changed color:
                        FcFplot([path],0.01,180,colortable=['green'])
    '''

    import pandas as pd
    import numpy as np
    from matplotlib.ticker import StrMethodFormatter
    import matplotlib.pyplot as plt
    plt.rcParams['figure.dpi'] = 150
    
    assert type(Pth) == list,"Please put path as a list! (Yes, even single item)"
    
    #TODO: makie it read Fc-Fo directly from shelx files
    
    clist = colortable
    f = plt.figure()
    ax=plt.gca()
    if log == True:
        ax.set_xscale('log')
        ax.set_yscale('log')
    f.set_figwidth(6)
    f.set_figheight(6)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.xlim(LowerLimit, UpperLimit)
    plt.ylim(LowerLimit, UpperLimit)
    for p in Pth:
        DF = pd.read_csv(p)
        plt.scatter(DF['x'],DF['y'],c=clist[0],s=(4*72./150)**2)
        clist = np.roll(clist,-1)
    plt.xlabel('Fc', fontsize=15)
    plt.ylabel('Fo', fontsize=15)
    plt.show()