# KrystalograPy - Small Scripts For Crystallography

This is a collection of functions I use in my crystallographic work. Most of them are used for plotting, but sometimes I add also different types of sctipts. Below I explain what each function do, and how to use it. Additionally, each function has docstring with explanation.

<b>Some functions are not finished yet. If you don't see a description here, it means that I'm still working on it.</b>

## Plotting functions


### Fc-Fo plotting
The function <b>fcfoplot</b> was made for comparing multiple Fc-Fo values. It can use two type of files as an imput: shelx FCF or CSV. The csv file should be in the same format as the one generated by Olex2 FcFo plot function. 
By default the logarithmic scale is used for both x and y to emphasize the low value region at which biggest disagreements between Fc and Fo are usually seen.

**How to use it?**

In this example I want to compare two plots. I assign paths to csv files to variables and then pass them to function. Also I set limits to 0.001 - 500000.

    from KrystalograPy import fcfoplot
    P = r"C:\structure\data1.fcf"
    P2 = r"C:\structure\data2.fcf"
    fcfoplot([P,P2],0.001,500000)
    
The result:
    
    
![image](https://user-images.githubusercontent.com/59794882/181138052-99c3d7f8-d7b6-4438-8826-871008399e4f.png)

You can also use this function for subplots, by passing axis object to the "axis" variable:

    from KrystalograPy import fcfoplot
    import matplotlib.pyplot as plt
    
    #A, B, C, D are paths to your data
    
    fig, axs = plt.subplots(nrows = 2, ncols = 2, figsize= (10, 6))
    fig.tight_layout(w_pad=2,h_pad=2)
    fcfoplot([A],0.001,500000,log=True,colortable='k',axis=axs[0][0])
    fcfoplot([B],0.001,500000,log=True,colortable='k',axis=axs[0][1])
    fcfoplot([C],0.001,500000,log=True,colortable='k',axis=axs[1][0])
    fcfoplot([D],0.001,500000,log=True,colortable='k',axis=axs[1][1])
    
The result:
![image](https://user-images.githubusercontent.com/59794882/181137969-ab07a3c8-c7eb-41da-b752-b0594a676c25.png)

### CIF data plotting

You can easily plot data from several CIF files with three functions presented below: <b>SearchInCif</b>, <b>MakeDF</b>, <b>plotSplit</b>.
The general idea is to search all CIF files in one location for specified data, then create DataFrame object from which the plotting is done. 
For plotting to be succesfull there are few things needed: (1) All CIF files should be in one directory, (2) All values must be a valid CIF keys*, (3) A index must be specified for x-axis.

*Unless you put your custom value in cif format, for all files

**How to use it?**
Let's say that we want to plot a, b, c and V in dependence of increasing temperature. First we need to perform search and create DataFrame object:

*!All data on plots comes from "faked" data i.e. I used standard CIF files but changed a, b, c, V and Temp. values to false ones, not associated with any real published data!*

    from KrystalograPy import SearchInCif, MakeDF, plotSplit
    
    root_path = r"C:\cifdata\"
    index = "_cell_measurement_temperature"
    
    #Notice that index and searched values are CIF keys
    #Remember, that index value also should be included in "searched" values!!!
    
    searched = ["_cell_length_a","_cell_length_b","_cell_length_c","_cell_volume","_cell_measurement_temperature"]
    
    results = SearchInCif(root_path,searched)
    DF = MakeDF(results,searched,index)

Now plotting can be performed. There are several ways of plotting within plotSplit function. Also there are some ways to customize the plot. We want to set labels for all axes, add units. First we want to see four separate plots in one figure. We do this by setting composition to "auto" (or "a"), or ignoring composition parameter (because "auto" is default behaviour).

    xlabel = "TEMP"
    ylabel = ['a','b','c','V']
    units_all  = ["[K]","[$\AA$]","[$\AA$]","[$\AA$]","[$\AA^3$]"]
    plotSplit(DF,searched,size=[10,6],x_label=xlabel,y_label=ylabel,units=units_all)

The result:
![image](https://user-images.githubusercontent.com/59794882/181144518-464f4401-8666-464f-9bdd-36a6a4c2c7de.png)

In some cases it is better to show such unit cell parameters change on one plot. In this case, different color will be assigned to each plot by us, and all data will be normalized to the initial one. **At the current version, labels and units must be specified for all subplots, despite x and y axis are shared.**

    colors    = ['r','g','b','k']
    labels    = ["Change of unit cell","Change of unit cell","Change of unit cell","Change of unit cell"]
    units_one = ['[K]','[%]','[%]','[%]','[%]']
    plotSplit(DF,searched,composition='s',size=[10,8],x_label=xlabel,y_label=labels,units=units_one,color=colors,legend=True)
    
The result:

![image](https://user-images.githubusercontent.com/59794882/181145662-08106bcc-27bb-4d95-a3aa-e88ffde82ac9.png)

