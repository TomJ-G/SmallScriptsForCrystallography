"""
Author: Tomasz Galica
Last updated: 01.07.2022
"""
#VandDerWaals radii dict.
#This is not finished, some atoms are missing.
VanDerWaals = {"H":[1.2,1], "Zn":[1.39,30], "He":[1.4,2], "Cu":[1.4,29], 
               "F":[1.47,9], "O":[1.52,8], "Ne":[1.54,10], "N":[1.55,7], 
               "Hg":[1.55,80], "Cd":[1.58,48], "Ni":[1.63,28], "Pd":[1.63,46],
               "Au":[1.66,79], "C":[1.7,6], "Ag":[1.72,47], "Mg":[1.73,12], 
               "Cl":[1.75,17], "Pt":[1.75,78], "P":[1.8,15], "S":[1.8,16], 
               "Li":[1.82,3], "As":[1.85,33], "Br":[1.85,35], "U":[1.86,92], 
               "Ga":[1.87,31], "Ar":[1.88,18], "Se":[1.9,34], "In":[1.93,49], 
               "Tl":[1.96,81], "I":[1.98,53], "Kr":[2.02,36], "Pb":[2.02,82], 
               "Te":[2.06,52], "Si":[2.10,14], "Xe":[2.16,54], "Sn":[2.17,50], 
               "Na":[2.27,11], "K":[2.75,19], }

def ArrayWIndex(Dataframe):
    """
    Parameters
    ----------
    Dataframe : pandas dataframe
        Dataframe to be converted to array.

    Returns
    -------
    new_array : numpy array
        Array with index value included.

    """
    import pandas as pd
    import numpy as np
    
    new_array = []
    
    Index = np.array(Dataframe.index)
    Values = Dataframe.values
    
    for i in range(0,len(Values)):
        new_array.append([Index[i],*Values[i]])
        
    new_array = np.array(new_array)
    
    return new_array

def popsmall(data):
    """
    Remove the smallest element from the list.
    :parameter data: list, list of numerical values
    :returns: list
    """
    from numpy import inf 
    small = inf
    n = 0
    index = 0
    
    for l in data:
        if l < small:
            small = l
            index = n
        n+=1
    data.pop(index)
    #data.sort()
    return(data,index)

def popouter(data,value):
    """
    Remove element the furthest of the supplied parameter.
    :parameter data: list, list of numerical values
    :parameter value: int|float, number with which list will be compared
    :returns: list
    """
    from numpy import inf, abs
    out = -inf
    n = 0
    index = 0
    
    for l in data:
        dif = abs(abs(l)-value)
        if dif > out:
            out = dif
            index = n
        n+=1
    data.pop(index)
    #data.sort()
    return(data,index)


def GetAtomNameVDW(n):
    """
    Based on WanDerWaals table, get atom type name from atomic number.
    :parameter n: int, Atomic number
    :return: string, Atom type name
    """
    for key,value in VanDerWaals.items():
        if value[1]==n:
            return(key)


def within(a,b,num,how='percent'):
    """
    Check if a is close to b. The acceptance range is defined by 'num' variable.
    There are three modes of comparison:
        percent  - range is definied by 'num'% of a. (for 10%, num = 0.1 etc)
        round    - a and b are rounded to 'num' decimal places and compared.
        constant - Difference between a and b must be less than 'num'.
    
    If a and b are list|tuple, each 1nth values will be compared (a[n] with b[n] etc.)
    
    Parameters
    ----------
    a : int|float
        1st number.
    b : int|float
        2nd number.
    num : int|float
        Comparing variable.
    how : str
        Comparing method: percent, round, constant.
    Returns
    -------
    Bool

    """
    iswithin = False
    
    assert type(a) in (int,float,tuple,list), "Type of A must be in int, float, tuple or list"
    assert type(b) in (int,float,tuple,list), "Type of B must be in int, float, tuple or list"
    assert how in ('percent','round','constant'), "How must be one of the follwing: percent, round, constant"
    
    if type(a) in (int,float) and type(b) in (int,float):
    
        if how == 'percent':
            if (b <= a*(1+num)) and (b >= a*(1-num)):
                iswithin = True
        elif how == 'round':
            if round(a,num) == round(b,num):
                iswithin = True
        elif how == 'constant':
            if abs(a-b) < num:
                iswithin = True
    
    elif type(a) in (list,tuple) and type(b) in (list,tuple):
        
        assert len(a) == len(b), "Size of the elements must be the same!"
        if how == 'percent':
            iswithin = all((b[n] <= a[n]*(1+num)) and (b[n] >= a[n]*(1-num)) for n in range(0,len(a)))
        elif how == 'round':
            iswithin = all(round(a[n],num) == round(b[n],num) for n in range(0,len(a)))
        elif how == 'constant':
            iswithin = all(abs(a[n]-b[n]) < num for n in range(0,len(a)))
    
    #Add functionality of mixed types
            
    return iswithin


def VDW_atom_search(atom,Cell,UnitCell):
    ###The purpose is to search for atoms which are within VDW radius.
    
    #1) Build unit cells surrounding this one.
    #2) Generate distances to each and every atom
    #3) Use WdV table to get distance.
    #4) Select only atoms which are below WvD
    
    #Cell == Cell_000 - this is the middle cell.
    #Cell_001 means that all atoms are +1 in 00L direction.
    
    
    return None


def BuildCell(atoms,codes_of_symmetry,acc=6):
    """
    Builds content of the unit cell.
    :param atoms: list, extracted atoms cartesian
    :param codes_of_symmetry: list, extracted from CIF
    :param acc: accuracy of rounding for duplicate check.
    :returns: list, all atoms in single unit cell
    """

    Cell = []
    
    for a in atoms:    
        temp = []
        x, y, z = a[2], a[3], a[4]
        
        #Transform for all symmetries and then compare positions.
        #If too close, then discard
        
        for s in codes_of_symmetry:   
            x, y, z = a[2], a[3], a[4]
            c1 = eval(codes_of_symmetry[s][0]) 
            c2 = eval(codes_of_symmetry[s][1])
            c3 = eval(codes_of_symmetry[s][2])
            
            cal = [(a[0]+"_"+str(s)),a[1],c1,c2,c3,s]
            
            #Temporary list of symmetrical atoms
            temp.append(cal)
        
        #Check if there are duplicates within accuracy
        pop = []
        for i in range(0,len(temp)):
                for j in range(i+1,len(temp)):
                    if within(temp[i][2:5],temp[j][2:5],acc,'round') == True:
                        pop.append(j)
        
        #List of indices to be removed from temp
        pop = list(set(pop))
        pop.sort(reverse=True)

        for i in pop:
            temp.pop(i)
        Cell += temp            
           
    return Cell


def LoadCif(Path,Q_peak_remove=True):
    #TODO
    #UPDATE DESCRIPTION
    """
    Reads CIF file and extracts unit cell info, symmetry operations, atom list,
    quality information, ect...
    :param Path: str, Path to CIF file
    :returns: symetry operators, atom list, parameter dictionary
    """
        
    from re import match, findall
        
    que_sym = r"\'(\s?\-?[xyzXYZ\+\-\d\/]+)\,\s?(\s?\-?[xyzXYZ\+\-\d\/]+)\,\s?(\s?\-?[xyzXYZ\+\-\d\/]+)\'"
    que_atm = r"(\w+\d?\w?)\s+(\d+)\s+(\d+\.?\d*)\s+(\d+\.?\d*)\s+(\d+\.?\d*)"
    
    pars = {"_cell_length_a":None,"_cell_length_b":None,
            "_cell_length_c":None,"_cell_angle_alpha":None,
            "_cell_angle_beta":None,"_cell_angle_gamma":None,
            "_cell_volume":None,"_cell_formula_units_Z":None,
            "_diffrn_reflns_av_R_equivalents":None,
            "_refine_ls_goodness_of_fit_ref":None,
            "_refine_ls_R_factor_all":None,
            "_refine_ls_R_factor_gt":None,"_refine_ls_wR_factor_gt":None,
            "_refine_ls_wR_factor_ref":None}
        
    #content = []
    
    symops = {}
    atoms = []
    symsearch, atmsearch, s = 0, 0, 0

    with open(Path,'r') as f:
        for line in f:

            if symsearch == 1:
                M = findall(que_sym,line)
                if len(M) != 0:
                    s+=1
                    symops['$'+str(s)] = ([M[0][0],M[0][1],M[0][2]])
                else:
                    symsearch = 0
            elif "_space_group_symop_operation_xyz" in str(line).lower():
                symsearch = 1
            elif atmsearch == 1:
                M = match(que_atm,line)
                if M:
                    if (Q_peak_remove == True) and ("Q" not in M[1]):
                        atoms.append([M[1],int(M[2]),float(M[3]),float(M[4]),
                                      float(M[5])])
                    elif Q_peak_remove == False:
                        atoms.append([M[1],int(M[2]),float(M[3]),float(M[4]),
                                      float(M[5])])
                elif "HKLF 4" in str(line).upper():
                    atmsearch == 0
            elif "FVAR" in str(line).upper():
                atmsearch = 1
            elif "_space_group_name_H-M_alt" in str(line):
                #This is not numeric value so different query is needed
                #If I need more non-numerics in future, I'll make separate par.
                L = findall(r"\s+\'(.+)\'",line)
            elif any(x in line for x in pars.keys()):
                M = findall(r"\s+([\d\.\(\)]+)",line)
                if len(M) !=0:
                    pars[str(line).split(" ")[0]] = M[0]

    pars["_space_group_name_H-M_alt"] = L[0]
            #This is for troubleshooting only
            #content.append(line)
    #return symops, atoms, par, content
    
    return symops,atoms, pars 


def Cartesian(coordinates,unit_cell):

    """
    Transforms crystal coordinates to cartesian.
    
    :param UnitCell: list|tuple|array, Unit cell parameters (a,b,c,alpha,beta,gamma)
    :param coordinates: list|tuple|array, coordinates to be transformed (x,y,z)
    
    :return: ndarray, Cartesian coordinates.
    """
    from numpy import cos, sin, sqrt, pi, array

    #If input is tuple, convert
    co = array(coordinates)
    
    #Unit cell
    a, b, c    = unit_cell[0], unit_cell[1], unit_cell[2]
    al,be,ga = unit_cell[3]*pi/180, unit_cell[4]*pi/180,unit_cell[5]*pi/180
    
    n = (cos(al) - cos(ga)*cos(be))/sin(ga)
    
    #Transform matrix
    M = array([[a,0,0],
         [b*cos(ga),b*sin(ga),0],
         [c*cos(be),c*n,c*sqrt(sin(be)**2-n**2)]])
    
    cartesian = co @ M
    
    return(cartesian)


def BondDist(atom1,atom2,UnitCell=None):
    """
    Calculates bond distance between pair of atoms.
    If coordinates are cartesian, only atom positions are needed.
    If coordinates are fractional (crystal system),
    also unit cell must be passed.
    
    :param atom1: list|tuple|array, Coordinates of 1st atom
    :param atom2: list|tuple|array, Coordinates of 2nd atom
    :param UnitCell: list|tuple|array, Optional parameter,
    If not supplied, it will be assumed that coordinates are alredy cartesian.
    Unit cell parameters to be passed: (a,b,c,alpha,beta,gamma)
    
    :return: float, Bond length
    """
    from numpy import sqrt
    
    if UnitCell == None:
        k1, k2 = atom1, atom2
        bond = sqrt((k1[0]-k2[0])**2+(k1[1]-k2[1])**2+(k1[2]-k2[2])**2)
    else:
        k1 = Cartesian(atom1,UnitCell)
        k2 = Cartesian(atom2,UnitCell)
        bond = sqrt((k1[0]-k2[0])**2+(k1[1]-k2[1])**2+(k1[2]-k2[2])**2)
    
    return bond


def Angle(atom1,atom2,atom3,UnitCell=None):
    """
    Calculates angle (in degrees), between 3 atoms.
    :param atom1: array, Coordinates of 1st atom
    :param atom2: array, Coordinates of 2nd atom
    :param atom3: array, Coordinates of 3rd atom
    :param UnitCell: list|tuple|array, Optional parameter,
    If not supplied, it will be assumed that coordinates are alredy cartesian.
    Unit cell parameters to be passed: (a,b,c,alpha,beta,gamma)
    
    :return: float, Angle in degrees
    """
    from numpy import arccos, array, dot, pi
    from numpy import linalg
    
    if UnitCell == None:
        k1, k2, k3 = array(atom1), array(atom2), array(atom3)
    else:
        k1 = Cartesian(atom1,UnitCell)
        k2 = Cartesian(atom2,UnitCell)
        k3 = Cartesian(atom3,UnitCell)
    
    #Calculate vectors
    v1 = k1 - k2
    v2 = k3 - k2

    coang = dot(v1, v2) / (linalg.norm(v1) * linalg.norm(v2))
    angle = arccos(coang) *180/pi
        
    return angle


def ExtractAtoms(Path):
    
    """
    Reads CIF, INS or RES file in search of atomic positions.
    :param Path: string, Full path to the CIF|INS|RES file
    :return: list, List of lists, with atom name, type and x, y, z positions.   
    """
    
    from re import match
    
    #Seek atoms
    content = []
    line_count = 0
    START = 0
    STOP = 0
    atoms = []

    #Put all content to variables and mark start-stop positions 
    with open(Path,'r') as f:
        for line in f:
            if "FVAR" in str(line).upper():
                START = line_count
            elif "HKLF 4" in str(line).upper() and STOP == 0:
                STOP = line_count
            content.append(line)
            line_count += 1

    #Extract atomic information to list
    #I do not care for displacement parameters at this moment
    query = r"(\w+\d?\w?)\s+(\d+)\s+(\d+\.?\d*)\s+(\d+\.?\d*)\s+(\d+\.?\d*)"
    for i in content[START+1:STOP]:
        M = match(query,i)
        if M:
            atoms.append([M[1],int(M[2]),float(M[3]),float(M[4]),float(M[5])])
        
    return atoms


def StructureCartesian(data,UnitCell):
    """
    Calculates structure to cartesian coordinates.
    :param data: string|list, Either path to CIF, INS or RES, or 
                              list of atoms in format:
                              [AtomName, AtomTypeNumber, X, Y, Z]
    
    :return: list, Atom list in cartesian system.
    """
    CartesianAtoms = []
    
    if type(data) == list:
        #This means that list of atoms is provided
        Atoms = data
        
    elif type(data) == str:
        #This means that atoms must be extracted
        Atoms = ExtractAtoms(data)
        
    for i in Atoms:
        a = Cartesian(i[2:],UnitCell)
        CartesianAtoms.append([i[0],i[1],*a])
   
    return CartesianAtoms


def ProcessCrystal(Path,Silent=True):

    #This is a full automatic processing of CIF file.
    #Processing steps:
        #1)Read CIF                           OK
        #2)Translate Atoms to Cartesian       OK
        #3)Build UnitCell class               OK
        #4)Construct full cell                OK
        #5)Prepare report of interactions     XX
    if Silent == False: print('Loading CIF file')
    Symm, Atom, Param = LoadCif(Path)
    UC = UnitCell(list(map(float,list(Param.values())[0:6])), list(Param.values())[-1])
    
    if Silent == False: print('Calculate to cartesian')
    AtomC = StructureCartesian(Atom, UC)
    
    if Silent == False: print('Build Unit Cell')
    Cell = BuildCell(AtomC, Symm)    

    return None


#=============================================================================#

def OMITME(maximal=999,minimal=0,r='high',error=10.0):
    
    """
    Prints OMIT command from the list of most disagreeable reflections.
    This function reads data from clipboard!
    Please copy " Most Disagreeable Reflections" to clip, then execute!
    Table must be of SHELXL format!
    
    :param maximal: float, Upper limit of fcf to be omitted.
    :param minimal: float, Lower limit of fcf to be omitted. Default is 0
    :param r: string, Rule by which omiting will be done;
                      default :'high'.
             
                Possible rules:
                    'high' - Omit only reflections with Fo^2 > Fc^2 |('h')
                    'low'  - Omit only reflections with Fo^2 < Fc^2 |('l')
                    'Error'- Omit reflections with Error higher than selected
                             value (check Error parameter). Fc/Fc(max) range 
                             will be ignored! |('e')
                             
    :param error: float, This variable will be considered only if r='Error'
                         (or r='e'). If reflection has error with higher value
                         than selected, it will be omited;
                         default: 10.0
    
    :examples:
        OMITME(0.014)          # will omit Fo^2 > Fc^2 below Fc/Fc(max) = 0.014
        OMITME(0.032,0.014)         # will omit Fo^2 > Fc^2 in Fc/Fc(max) range
        OMITME(0.032,0.014,rule='low')         # will omit Fo^2 < Fc^2 in range
        OMITME(0,r='Error',error=10)    # will omit reflections with error > 10
    """
    
    import pandas as pd
    DF = pd.read_clipboard()
    
    assert type(r) == str, "Rule must be a string!"
    r = r.lower()
    
    #If Fo^2 > Fc^2
    if r in ('high','h','Fo'):
        for n in range(0,len(DF)):
            if DF.iloc[n][6]<=0.031 and float(DF.iloc[n][3])>float(DF.iloc[n][4]):
                print("OMIT",int(DF.iloc[n][0]),int(DF.iloc[n][1]),int(DF.iloc[n][2]))
    #If Fo^2 < Fc^2
    elif r in ('low','l','Fc'):
        for n in range(0,len(DF)):
            if DF.iloc[n][6]<=0.031 and DF.iloc[n][3]<DF.iloc[n][4]:
                print("OMIT",int(DF.iloc[n][0]),int(DF.iloc[n][1]),int(DF.iloc[n][2]))
    elif r in ('error','e','error/esd'):
        for n in range(0,len(DF)):
            if DF.iloc[n][5]>=error:
                print("OMIT",int(DF.iloc[n][0]),int(DF.iloc[n][1]),int(DF.iloc[n][2]))
    else:
        print('Incorrect rule was selected!\n')
        print('Rule must be one of the following:\n')
        print('"high" (or "h", "Fo"); "low" (or "l", "Fc"); "error" (or "e", "error/esd")')

#Read file and extract content
def SearchInCif(path,columns):
    
    """
    Read all CIF files in single directory and seek for selected values.
    
    :param path: path, absolute path to folder in which CIF files are stored,
    :param columns: list, Values which should be plotted.
                    REMEMBER! INDEX value (x-axis) also should be included!
    :return: list, List of values from CIF. The same order as in "columns" 
    """
    import re
    from pathlib import Path

    
    Read_folder = Path(path).glob('*.cif')
    files = [x for x in Read_folder]
    
    Results = []
    for file in files:
        content = []
        with open(file,'r+') as InFile:
            for line in InFile:
                content.append(line)
        res = []
        err = []
        for query in columns:
            query="("+query+")"+"\s+(-?\d+\.?\d+)\(?(\d*)\)?"
            for line in content:
            #Create regex pattern from searched value
                M= re.match(query, line)
                if M:
                    res.append(M[2])
                    #Because of way errors in cif files are written extraction is needed
                    if M[3] !='':
                        digts = re.match(r'\d*.(\d*)',str(M[2]))
                        Ec = int(M[3])/(10**len(digts[1]))
                        err.append(Ec)
                    else:
                        #If no error is present for searched value, use zero instead
                        err.append(0)
        res=res+err
        Results.append(res)

    return Results


def MakeDF(data,columns,Index,Kelvin=False,SavePath=None):
    """
    Creates DataFrame from the values taken from CIF files.

    :param values: list, Results from ReadCif list of lists for searched values
    :param columns: list, Searched values. Used for headers of DataFrame
    :param index: string, Name of index column
    :param Kelvin: bool, If one of your values is temperature "TEMP",
                   if Kelvin=True, this will recalculate from C to K;
                   default: False
    :param SavePath: path, DataFrame will be stored as csv in selected path;
                     default: None
    :return: DataFrame
    """
    import pandas as pd
    
    search = columns.copy()
    #Add error columns
    for x in columns:
        if x != Index:
            search.append(x+"_error")
        else:
            n = (columns.index(x)*2)+1
            for i in data:
                del(i[n])
    #if found index val then delete "error" col from results
    
    #Create dataframe for graph    
    DF = pd.DataFrame(data, columns =search)

    #Convert Celcius to kelvin
    if Kelvin==True:
        if ("TEMP" in list(DF.columns)) == True:
            DF.TEMP = DF.TEMP.astype(float) + 273.15
        else:
            print("'TEMP' column not found. Check naming")
            
    #Modify Dataframe
    DF.set_index(Index,inplace=True)
    DF = DF.astype(float)
    DF.index = DF.index.astype(float)
    DF.sort_index(ascending=True, inplace=True)

    #Save file
    if SavePath != None:
        DF.to_csv(SavePath)
    
    return DF


def plotSplit(DF,columns,composition = "auto",color='b',size = [22,12],
              x_label = None,y_label=None,error_bars=False,
              units=None, y_limits=None,legend=False,font_size = 12):
    
    """
    Nicely composed plots, using data frame from 'MakeDF'.
    :param DF: DataFrame, result from MakeDF function 
    :param columns: list, Headers same as used for ReadCif and MakeDF functions
    :param composition: str, User can select 'auto', 'column' and 'single'; 
                      default: "auto"
                      
                      All compositions:
                      'auto' - Tries to compose plots in n-by-m matrix. ('a')
                      'column' - Each plot is below the previous one. (c')
                      'rows'   - Plots next to each other, in one row. ('r')
                      'single' - One plot with all the data.                       
                      1st value of each column normalized to 1. ('one','s')

    :param size: tuple|list, Plot size, the same as in figsize of pyplot
    :param x_label: string, Changes x-axis label, from the default one 
                    (the one you selected as index of DF).
    :param y_label: string, Changes y-axis names.
    :param error_bars: bool, If set to True, will plot error bars;
                       default: False
    :param units: string, adds units to labels (no matter if custom or default)
    :param y_limits: list(list(float))|tuple(tuple(float)), Limits y-axis 
                     You can limit any number of values by making list:
                     Ex.: If only second plot should have limits: 
                     [None,(1,5),None,...].
    :param color: string, Colors for markers. Acceots either 1 value for all 
                          or 1 for each column.
    :param legend: bool, enables or disables legend. ONLY IF composition=SINGLE
    :param font_size: int|float, change font size of ALL labels
        
    EXAMPLE:
        Step-by-step. Define path to CIF files, index column and searched vals.
        Then call reading cif, making dataframe and plotting:
            Path = r"C:/MyPath/CIFs"
            ColumnNames = ["_cell_length_a","_cell_length_b",
                           "_cell_measurement_temperature"] 
            res = SearchInCif(C:\MyPath\,ColumnNames)
            DF = MakeDF(res,ColumnNames,ColumnNames[2])
            plotSplit(DF,LookFor)
    """
    import pandas as pd
    from numpy import sqrt, ceil, roll
    import matplotlib.pyplot as plt
    
    def plotter(DF,columns,cmp,i,j,x_label,y_label,units,error_bars,color=color,y_limits=y_limits):
        """
        This is the core funtion of plotting. Variables the same as in main function
        """
        #Check how many colors in color table
        if type(color) == list or type(color) == tuple:
            clen = len(color)
            if clen > 0:
                clr = color[i]
        else:
            clr = color
        
        #Exclude index value
        if columns[i] != DF.index.name:
            plt.subplot(cmp[0], cmp[1], j+1)
            #Override labels values
            if x_label == None:
                plt.xlabel(DF.index.name+units[0])
            else:
                plt.xlabel(x_label+units[0])
            if y_label == None:
                plt.ylabel(columns[i]+units[i+1])
            else:
                plt.ylabel(y_label[i]+units[i+1])

            #Plot values
            #Check if limits exists
            if y_limits != None and (type(y_limits) in (list,tuple)):
                plt.scatter(DF.index,DF[columns[i]],c=clr,label=y_label)
                plt.ylim(y_limits[i])
            else:
                plt.scatter(DF.index,DF[columns[i]],c=clr,label=y_label)
            #Plot error bars if True
            if error_bars == True:
                plt.errorbar(DF.index,DF[columns[i]], DF[columns[i]+"_error"], fmt="o",ecolor=clr,color=clr)
        #=====================================================================================#

    
    #Check length of query
    length = len(columns) - 1
    
    #if there are label values, add space to each
    if units == None:
        unitsC = []
        for i in range(0,len(columns)):
            unitsC.append("")
    else:
        unitsC = units.copy()
        for i in range(0,len(columns)):
            unitsC[i] = " "+unitsC[i]
    
    #Check what composition should be used and execute plotter
    if composition == "auto" or composition == "a":    
        cmp = (int(round(sqrt(length))), int(ceil(sqrt(length))))
        plt.rc('font', size=font_size) 
        f = plt.figure(figsize=size)
        for i in range(0,length):
            plotter(DF,columns,cmp,i,i,x_label,y_label,unitsC,error_bars)
    elif composition == "column" or composition == "c":
        cmp = (length,1)
        plt.rc('font', size=font_size)
        f = plt.figure(figsize=size)
        for i in range(0,length):
            plotter(DF,columns,cmp,i,i,x_label,y_label,unitsC,error_bars)      
    elif composition == "rows" or composition == "r":
        cmp = (1,length)
        plt.rc('font', size=font_size)
        f = plt.figure(figsize=size)
        for i in range(0,length):
            plotter(DF,columns,cmp,i,i,x_label,y_label,unitsC,error_bars)
    elif composition == 'single' or composition == 'one' or composition == "s":
        lgnd = []
        cmp = (1,1)
        #Recalculates values - 1st value for each column is equal to 1.
        DF2 = DF.copy()
        plt.rc('font', size=font_size)
        f = plt.figure(figsize=size)
        for i in range(0,length):
            for n in columns:
                if n != DF.index.name:
                    start = DF[n].iloc[0]
                    temp = []
                    #CANNOT ADD LEGEND IF VALUE STARTS WITH "_"
                    lgnd.append(n.strip("_"))
                    for j in range(0,len(DF)):
                        temp.append(DF[n].iloc[j]/start)
                    DF2[n] = temp
            plotter(DF2,columns,cmp,i,0,x_label,y_label,unitsC,error_bars)
        if legend == True:
            plt.legend(lgnd)
    else:
        raise Exception("Incorrect composition! Please use one of the following: auto, column, rows, single.")
    
    return f
        
def TwoAxPlot(Data,colors=['r','b'],markers=["o","s"],select=[None,None,None],
              ylabels=None,xlabel=None,log=False,ax=None,Legend=None,box=None):
    """
    Two scatter plots on one figure, with shared x axis, but different y axes.
    By default it uses three first columns of data as x, y1 and y2.
    
    :param Data: array|list(list), Should consist of 3 columns (1st = index)
    :param colors: list(str)|tuple(str), Colors for markers of 1st, 2nd dataset
                   default: ['r','b']
    :param markers: list(str)|tuple(str), Markers shape,
                    default: ['o','s']
    :param select: list, column number which should be used for x, y1 and y2.
                    default: [0,1,2]
    :param ylabels: list(str)|tuple(str), Labels for 1st and 2nd y-axis
    :param xlabel: string, X-axis label
    :param log: bool, Set to True for logarithmic x-scale
                      default: False
    :param ax: axis, Axis on which plotting will be done
    :param Legend: list(int,int), If any value is used, it will be the location
                                  of the legend;
                   default: None               
    :param box: list(float,float), If any value is used, it will be used for
                                   legend box positioning.
                   default: None
                 
    :example:
        A = [[1,2,1],[2,4,4],[3,6,9],[4,8,16]]
        fig = plt.figure(figsize=(10,6))
        TwoAxPlot(A,ylabels=['2*x','x^2'],xlabel='x',ax=plt.axes(),
                  Legend="best",box=[1.01,1.12])

    """
    
    import matplotlib.pyplot as plt
    import numpy as np
    from pandas.core.frame import DataFrame
    
    if ax == None:
        ax=plt.gca()
    
    #Check data type. All except data frame translate to array.
    #DataFrame is exception, because index is ommited with translation!
    if type(Data) == DataFrame:
        Data = ArrayWIndex(Data)
    elif type(Data) != np.ndarray:
        Data = np.array(Data)
    
    #Data selection
    if select[0] == None:
        x_data = Data[:,0]
    else:
        x_data = Data[:,select[0]]
    
    if select[1] == None:
        y1_data = Data[:,1]
    else:
        y1_data = Data[:,select[1]]
    
    if select[2] == None:
        y2_data = Data[:,2]
    else:
        y2_data = Data[:,select[2]]
    
    if xlabel != None:
        ax.set_xlabel(xlabel)
    pl1 = ax.scatter(x = x_data, y=y1_data, c=colors[0], marker=markers[0])

    if ylabels != None:
        ax.set_ylabel(ylabels[0])
        ax.yaxis.label.set_color(colors[0])
    ax2 = ax.twinx()
    pl2 = ax2.scatter(x = x_data,y = y2_data,c=colors[1],marker=markers[1])

    if ylabels != None:
        ax2.set_ylabel(ylabels[1])
        ax2.yaxis.label.set_color(colors[1])
    if log == True:
        ax.set_xscale('log',base=10)
    if Legend !=None:
        if box != None:
            ax.legend([pl1, pl2], [ylabels[0], ylabels[1]],loc=Legend,bbox_to_anchor=box)
        else:
            ax.legend([pl1,pl2],ylabels,loc=Legend)
            
            
def fcfoplot(Pth,LowerLimit,UpperLimit,log=True,
             colortable=['r','b','g','orange','magenta','cyan','black','purple','pink'],
              axis=None):
   
    '''
    This function plots Fo-Fc graph from fcf or csv file.
    Input values are paths to csv or fcf files, lower and upper limits for plotting (the same value is used for x and y).
    You can supply one path as list for single plot, or several paths in list file, for multiple overlaid graphs.
    The csv is expected to be in the same format as one produced by Olex2 and fcf should be SHELX style.
    Paths to csv and fcf can be mixed!
        INPUT:
            Pth - one or more paths to fcf or csv files AS A LIST!
            LowerLimit - minimal value of Fo/Fc to be plotted,
            UpperLimit - maximal value of Fo/Fc to be plotted,
            log - by default set to True. Change to False to not plot logarithmic scale
            colortable - you can override colors used for plotting. By default table contains 9 colours:
                         red, blue, orange, magenta, cyan, black, purple, pink
        EXAMPLES:
            Single plot:
                FcFplot(['C:/data/SomeFCF.csv'],0.01,200)
            Multiple plots:
                        fcfoplot(['C:/data/SomeFCF_1.csv','C:/data/SomeFCF_2.fcf'],0.01,200)
            Changed color:
                        fcfoplot([paths],0.01,180,colortable=['green'])
            Different limits on each axis, logarithmic scale:
                        fcfoplot([paths],(0.5,0.5),40000,log=True)
    '''

    import re
    import pandas as pd
    import numpy as np
    from matplotlib.ticker import StrMethodFormatter

    import matplotlib.pyplot as plt
    
    plt.rcParams['figure.dpi'] = 150
    ext = r"(?:.+)(\.\w+)"
    fcfq = r"(\s*-?\d\s+-?\d\s+-?\d\s+)(\d+\.\d*\s+)(\d+\.\d*\s+)(\d+\.\d*\s+)"
    
    assert type(Pth) == list,"Please put path as a list! (Yes, even single item)"
    assert type(LowerLimit) in (int,float,list,tuple),"Unrecognized fromat detected for variable: LowerLimit."
    assert type(UpperLimit) in (int,float,list,tuple),"Unrecognized fromat detected for variable: UpperLimit."
    
    if type(LowerLimit) in (int,float):
        LLx = LowerLimit
        LLy = LowerLimit
    elif type(LowerLimit) in (list,tuple):
        LLx = LowerLimit[0]
        LLy = LowerLimit[1]
    
    if type(UpperLimit) in (int,float):
        ULx = UpperLimit
        ULy = UpperLimit
    elif type(UpperLimit) in (list,tuple):
        ULx = UpperLimit[0]
        ULy = UpperLimit[1]
    
    clist = colortable

    if axis == None:
        ax=plt.gca()
    else:
        ax = axis
        
    if log == True:
        ax.set_xscale('log')
        ax.set_yscale('log')
    #f.set_figwidth(6)
    #f.set_figheight(6)
    ax.set_aspect('auto', adjustable='box')
    ax.set_xlim(LLx, ULx)
    ax.set_ylim(LLy, ULy)
    for p in Pth:
        
        m = re.match(ext,p)
        #Check if fcf or csv
        if m[1]=='.fcf':
            content = []
            fcf = []
            with open(p,'r') as f:
                for i in f:
                    content.append(i)
                    M = re.match(fcfq,i)
                    if M:
                        temp = []
                        for val in (2,3,4):
                            t1 = str(M[val]).strip()
                            temp.append(float(t1))
                        fcf.append(temp)
            ax.scatter(np.array(fcf)[:,0],np.array(fcf)[:,1],c=clist[0],s=(5*72./150)**2) 
            clist = np.roll(clist,-1)
        elif m[1]=='.csv':
            DF = pd.read_csv(p)
            ax.scatter(DF['x'],DF['y'],c=clist[0],s=(4*72./150)**2)
            clist = np.roll(clist,-1)
        ax.set_xlabel('Fc', fontsize=15)
        ax.set_ylabel('Fo', fontsize=15)
    
    
#=============================================================================#
#=============================================================================#
#Clases
class UnitCell:
    """
    This class holds information about crystal UnitCell: 
    real and reciprocal parameters, volume, A',B' and G matrix.
    :param unit_cell: - list of UC parameters: [a,b,c,alpha,beta,gamma]
    Properties:
        a, b, c, al, be, ga, V - reals space UC
        ra, rb, rc, ral, rbe, rga - reciprocal UC
        rcosa, rcosb, rcosg - reciprocal cosines
        Ap - A' matrix
        Bp - B' matrix
        G  - metric matrix
    """
    def __init__(self,unit_cell,lattice) -> None:
        
        #Dependencies
        from numpy import sqrt, cos, sin, pi, array, arccos
        
        self.unit_cell = unit_cell
        self.lattice = lattice
        #Real-space unit cell parameters
        self.a  = unit_cell[0]
        self.b  = unit_cell[1]
        self.c  = unit_cell[2]
        self.al = unit_cell[3]
        self.be = unit_cell[4]
        self.ga = unit_cell[5]
        
        #Theese three are just to not repeat *pi/180 all the time.
        rda = self.al*pi/180
        rdb = self.be*pi/180
        rdg = self.ga*pi/180
        
        #Calculate the volume
        self.V  = self.a*self.b*self.c*sqrt(1-cos(rda)**2
                 -cos(rdb)**2-cos(rdg)**2
                 +2*(cos(rda)*cos(rdb)*cos(rdg)))
        
        #Reciprocal space unit cell parameters
        self.ra = self.b*self.c*sin(rda)/self.V
        self.rb = self.a*self.c*sin(rdb)/self.V
        self.rc = self.b*self.a*sin(rdg)/self.V
        #Cosine variables are left to minimalize re-calculation error.
        self.rcosa = (cos(rdb)*cos(rdg)-cos(rda))/(sin(rdb)*sin(rdg))
        self.rcosb = (cos(rda)*cos(rdg)-cos(rdb))/(sin(rda)*sin(rdg))
        self.rcosg = (cos(rda)*cos(rdb)-cos(rdg))/(sin(rda)*sin(rdb))
        #If i need to use other value than cosine
        self.ral = arccos(self.rcosa)*180/pi
        self.rbe = arccos(self.rcosb)*180/pi
        self.rga = arccos(self.rcosg)*180/pi

        #Matrixes
        self.Ap = array([[self.ra*sin(self.rbe),self.rb*(self.rcosg-self.rcosa*
                        self.rcosb)/sin(self.rbe),0],
                   [0,self.rbe*(((sin(self.ral))**2)-((self.rcosg-self.rcosa*
                    self.rcosb)**2)/sin(self.rbe)**2)**(0.5),0],
                   [self.ra*self.rcosb,self.rb*self.rcosa,self.rc]])
        self.Bp = array([[self.a*sin(self.ga),self.a*cos(self.ga),0],
                      [0,self.b,0],
                      [self.c*(cos(rdg)*cos(rda)-cos(rdb))/sin(rdg), self.c*
                       cos(rda),
                       self.c*(((sin(rda))**2)-((cos(rdg)*cos(rda)-
                      cos(rdb))**2)/sin(rdg)**2)**(0.5)]])        
        self.G  = array([[self.ra*self.ra,self.ra*self.rb*
                          self.rcosg,self.ra*self.rc*self.rcosb],
                       [self.ra*self.rb*self.rcosg,self.rb*
                        self.rb,self.rb*self.rc*self.rcosa],
                       [self.ra*self.rc*self.rcosb,self.rb*
                        self.rc*self.rcosa,self.rc*self.rc]])
        

    
    def __repr__(self):
        """Temporary representation function"""
        return(str([self.a,self.b,self.c,
             self.al,self.be,self.ga,self.V, self.lattice,
             self.ra,self.rb,self.rc,self.ral,self.rbe,self.rga,
             self.Ap,self.Bp,self.G]))
    
    def __getitem__(self,sliced):
        
        Structure =  [self.a,self.b,self.c,
                      self.al,self.be,self.ga,self.V,self.lattice,
                      self.ra,self.rb,self.rc,self.ral,self.rbe,self.rga,
                      self.Ap,self.Bp,self.G]
        
        return(Structure[sliced])
    