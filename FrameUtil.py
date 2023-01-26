#This is a module with utility functions for diffraction frame processing.
#Some updrades in the functions may be needed over time as some parameters are temporary fixed
#or some other issues are unresolved

#=================================== H5 MODULE =======================================#


from numba import njit

@njit
def assembly(I_Data,bad_px,Xv,Yv,Empty,x_size=8192,y_size=512):
    """
    Assemble SACLA H5 composite file into continous H5 which can be read as
    an array (useful for further processing) or loaded to programs such as CrysAlisPro.
    Numba was used to improve execution speed of this function.
    Because of this, empty array must be provided - cannot be created in function!
    Parameters:
        I_Data - array with pixel intensities
        bad_px - array with bad pixel positions
        Xv     - array with X positions (from H5 GEOM file)
        Yv     - array with Y positions (from H5 GEOM file)
        Empty  - prepared array of new H5 image
    """
    for x in range(0,x_size):
        for y in range(0,y_size):
            if bad_px[x][y] == 0:              
                xm = int(Xv[x][y])
                ym = int(Yv[x][y])
                Empty[xm,ym] = I_Data[x,y]
    return Empty

@njit
def assembly_ranged(I_Data,bad_px,Xv,Yv,Empty,pv,Lran,Uran,x_size=8192,y_size=512):
    """
    Assemble SACLA H5 composite file into continous H5 which can be read as
    an array (useful for further processing) or loaded to programs such as CrysAlisPro.
    The same as "assembly" function, except it filters only selected wavelengths.
    Numba was used to improve execution speed of this function.
    Because of this, empty array must be provided - cannot be created in function!
    Parameters:
        I_Data - array with pixel intensities
        bad_px - array with bad pixel positions
        Xv     - array with X positions (from H5 GEOM file)
        Yv     - array with Y positions (from H5 GEOM file)
        Empty  - prepared array of new H5 image
        pv     - wavelength of dataset,
        Lran, Uran - lower and upper limit for wavelength range,
    Output:
        Empty - assembled file
        Ec    - If wavelength out of range =0, if in range =0.
        """
    Ec = 0
    if (pv > Lran and pv < Uran):
        Ec = 1
        for x in range(0,x_size):
            for y in range(0,y_size):
                if bad_px[x][y] == 0:              
                    xm = int(Xv[x][y])
                    ym = int(Yv[x][y])
                    Empty[xm,ym] = I_Data[x,y]
    return Empty, Ec


def rename_blk(FPATH,Name,starti,digits,ext):
    """
    Rename ALL frame files in the location.
    Parameters:
        FPATH  - full path to folder, where files are stored
        Name   - Frame prefix
        starti - starting frame number
        digits - number of trailing digits
        ext    - Changed files extension
    Example:
        "Fr_", starti = 1, digits = 4 --> Fr_0001
    """
    import glob, os
     
    processed = glob.glob(FPATH+'/*.'+str(ext))
    for frame in processed:
        sname = FPATH+('/'+Name+str(starti).zfill(digits)+'.h5')
        os.rename(frame,sname)
        starti+=1


def Composite_powder(file_list,save_path,ranged=False,min=0.00,max=9.99):
    """
    Calculate composite powder image from H5 frames.
    By the default, all frames from the root folder are used.
    Instead of taking folder as input, list of frames is needed.
    
    Parameters:
        file_list  - list of H5 frames. Can be made with list_frames func.
        save_path  - path and name of the composite file.
        ranged     - False by default. User can select wavelength range for frames.
        min        - minimum wavelength.
        max        - maximum wavelength.
    Output
        composite   - array of composite image.
    """
    import h5py
    from numpy import full, maximum, array
    
    Composite = full([2400,2400],-1)
    
    #Iterate over all files in root
    for frame in file_list:
        h5f = h5py.File(frame, 'r')
        groups = list(h5f.keys())
        #Iterate over all frames in file
        for dset in groups :
            if dset not in 'metadata':
                #Get photon wavelength
                ds = h5f[dset]['photon_wavelength_A']
                pv = ds[()]
                #Only if ranged=True
                if (pv > min) and (pv < max):             
                    #Get array data                
                    ds = h5f[dset]['data']
                    Arr = array(ds)                
                    Sv = maximum(Composite,Arr)
                    Composite = Sv
            else:
                continue
            
        h5f.close()
    Save_hf = h5py.File(save_path, 'w')
    Save_hf.create_dataset('composite_file', data=Composite)
    Save_hf.close()
    
    return Composite



def average_background(frames,ftype):
    """
    Calculates average background from frame list.
    Uses FrameRead function (see for details).
    Parameters:
        frames - list of frames (adresses of frames, not frame content!)
        ftype - file type handled by FrameRead (TIF,IMG,H5)
    Output:
        avg - numpy array with averaged background (array is a single frame)
    """
    from numpy import ndarray, int32, add
    
    x,h, s = FrameRead(frames[0],ftype)
    N = len(frames)
    avg = ndarray(shape=s, dtype=int32)
    for f in frames:
        d,h, s = FrameRead(f,ftype)
        avg = add(avg,d)
    avg = avg/N
    return avg


def FindPeaks(Frame,Threshold,mdist,plot_on=False,maxlim=10):
    """
    Peaksearch basing on threshold value. Can also plot frame with shown peaks positions.
    Output is a list with [x,y] peak positions.
    Parameters:
        Frame - the frame data in a form of array,
        Threshold - minimum intensity of the peak,
        mdist - minimum distance between two peaks to count them as separate,
    Optional:
        plot_on - if True, the frame will be plotted, with peak positions
        maxlim - the colormap limit
    Output:
        coordinates - list of coordinates with elements [x,y]
    """
    from skimage.feature import peak_local_max
    import matplotlib.pyplot as plt
    
    coordinates = peak_local_max(Frame, min_distance=mdist,threshold_abs=Threshold)
        #Plot peaks
    if plot_on == True:
        f = plt.figure()
        ax=plt.gca()
        f.set_dpi(150)
        plt.imshow(Frame, cmap='inferno', vmin=0, vmax=maxlim)
        plt.plot(coordinates[:, 1], coordinates[:, 0], 'c.',markersize = 3)       
    return coordinates


def FrameRead(fname,ftype,header_size=4096,shape=(1043, 981)):
    """
    Reads single frame and outputs its data, header and shape. 
    It can read two formats: tif and h5
    Parameters:
        fname - frame address on the hard drive,
        ftype - frame extension. It can be 'tif','h5' and 'esperanto'
        header_size - header size when reading a tif file
        shape       - array shape
    Output:
        frame - array with intensity data,
        header - header of the image,
        fshape - size of the image
    
    """
    
    import h5py
    import numpy as np
  
    
    if ftype == 'tif':
        with open(fname, 'rb') as fr:
            header = fr.read(header_size).decode('unicode_escape')
            #SIZE IS HARDCODED!!! CHANGE FOR FUTURE WORK!!!
            frame = np.ndarray(shape=shape, dtype=np.int32, buffer=fr.read(4092732))
            frame = np.rot90(frame,k=-1)
            fshape = frame.shape
    elif ftype == 'h5':
        header = None
        h5f = h5py.File(fname, 'r')
        groups = list(h5f.keys())
        for dset in groups :
            if dset not in 'metadata':
                ds = h5f[dset]
                frame = np.array(ds)
                fshape = frame.shape
        h5f.close()
    return frame, header, fshape


def FramePlot(frame,maxlim,colormap='inferno',dpi=100,*kwargs):
    """
    Plot frame array on screen
    Parameters:
        frame    - Array to be ploted.
        maxlim   - Upper limit for the colormap.
        colormap - Color map to be used with plotting. The same values as for pyplot imshow are accepted.  
    Output:
        None. Image will be shown on screen.
    """
    import matplotlib.pyplot as plt
    f = plt.figure()
    ax=plt.gca()
    f.set_dpi(dpi)
    plt.imshow(frame,cmap=colormap, vmin=0, vmax=maxlim,*kwargs)
    plt.show()


def fract_hkl(coordinates,shape,detector_center,pixel_size,det_dist,wavelength,mat_B,mat_T):
    """
    Calculates fractional hkl (unknown hkl) as per W.Kabsch - J.Appl.Cryst. (1977) 10, 426-429.
    Parameters:
        coordinates     - list of spots, each element of a list is a tuple with [x,y] data.
        shape           - shape of the Frame array.
        detector_center - Beam center position in x,y format.
        pixel_size      - Size of pixel in mm.
        det_dist        - Distance to the detector in mm.
        wavelength      - wavelength in angstrom.
        mat_B           - B matrix.
        mat_T           - T matrix.
    Outputs:
        hkl             - List of fractional hkl. Each element of the list is in [h,k,l] format.
        r3d             - 3D space coordinates of the diffraction spots.
        theta           - List of theta angles for each spot.
    """
    ####Import all dependencies
    import numpy as np
    
    #Declare empty lists for spot positions storage
    r3d = []
    hkl = []
    theta = []
    #Calculate for each spot on the frame
    for peak in coordinates:
    #coordinates read from frame are in reversed order: y,x
        pos_x = peak[1]
        pos_y = shape[0] - peak[0]
    #Coordinates in detector coordinate system
        shft_x = (pos_x - detector_center[0])*pixel_size
        shft_y = (pos_y - detector_center[1])*pixel_size 
    #Angular information of reflection/unused but do not delete for now
        dcentr = (np.sqrt((detector_center[0]-pos_x)**2+(detector_center[1]-pos_y)**2))#*pixel_size
        tang = dcentr/det_dist
        th =  np.arctan(tang)
        #chi = np.arctan2(shft_y,shft_x) * 180 / np.pi
    #Reciprocal coordinates of reflections
        mianownik =wavelength*np.sqrt((shft_x**2)+(shft_y**2)+det_dist**2)
        rx = shft_x/(mianownik)
        ry = shft_y/(mianownik)
        rz = (det_dist/(mianownik)) - (1/wavelength)
    #save coordinates to list
        r3d.append([rx,ry,rz])
    #calculate xyz to fractional hkl
        hkl.append(np.dot(np.dot(mat_B,mat_T),np.array([[rx],[ry],[rz]])))
        theta.append(th)
    return hkl,r3d,theta


###This one requires KrystalograPy module
def HKL_check(coord,hkl,latt_type,rules,Tinv,Binv,xy,dd,wv,px):
    """
    Calculates integer hkl from fractional one, as per W.Kabsch - J.Appl.Cryst. (1977) 10, 426-429.
    Function calculates 8 close integers and checks which one is the closest.
    Parameters:
        coord     - list of peak coordinates
        hkl       - list of hkl coordinates (same order as coord)
        latt_type - type of lattice centering.
        rules     - zonal and serial rules for reflection to be observed (see zs_absence_check).
        Tinv - inverted T matrix
        Binv - inverted B matrix
        xy - primary beam position
        dd - detector distance
        wv - wavelength
        px - pixel size
    Output:
        recalculated_XY - list of XY spot positions, for integer hkl.
        dist_difference - list of distance differences between integer and fractional hkl.
        integer_hkl     - list of integer hkl closest to the fractional hkl.
    """    
    from numpy import floor, ceil, matmul, sqrt
    import KrystalograPy as kp
    
    recalculated_XY = []
    integer_hkl = []
    dist_difference = []
    for i in range(0,len(hkl)):
        ci = []
        num = hkl[i]
        #Get 8 integers closest to hkl
        ci.append([floor(float(num[0])),floor(float(num[1])),floor(float(num[2]))])
        ci.append([floor(float(num[0])),floor(float(num[1])),ceil(float(num[2]))])
        ci.append([floor(float(num[0])),ceil(float(num[1])),ceil(float(num[2]))])
        ci.append([ceil(float(num[0])),ceil(float(num[1])),ceil(float(num[2]))])
        ci.append([ceil(float(num[0])),ceil(float(num[1])),floor(float(num[2]))])
        ci.append([ceil(float(num[0])),floor(float(num[1])),floor(float(num[2]))])
        ci.append([ceil(float(num[0])),floor(float(num[1])),ceil(float(num[2]))])
        ci.append([floor(float(num[0])),ceil(float(num[1])),floor(float(num[2]))])
        pos = []
        #Exclude systematic absent reflections
        for c in ci:
            ec = kp.centering_absence_check(latt_type,*c)
            ec2 = kp.zs_absence_check(rules,*c)
            if ec == 1 and ec2 == 1:
                pos.append(c)
            elif ec == 99:
                print("PROGRAM INTERRUPTED!!!")
                break
        rcalc_coord = []
        for p in pos:
            #Take integer hkl and recalculate each one to x,y,z
            rcalc_xyz = matmul(Tinv,matmul(Binv,p))
            #Recalculate to X,Y (spot position on detector)
            rcalc_X = xy[0] - (rcalc_xyz[1]*dd)/(rcalc_xyz[2]+(1/(wv)))/px
            rcalc_Y = xy[1] + (rcalc_xyz[0]*dd)/(rcalc_xyz[2]+(1/(wv)))/px
            rcalc_coord.append([rcalc_X,rcalc_Y])
        #Initialize distance and index value
        dist_dfr = 99999
        #ival = 99
        #Obtain coordinate with smallest difference
        for o in range(0,len(pos)):
            rcc = rcalc_coord[o]
            clc_dfr = sqrt((rcc[0]-coord[i][0])**2 + ((rcc[1]-coord[i][1])**2))
            if clc_dfr < dist_dfr:
                dist_dfr = clc_dfr
                #rcalc_coordt = rcc
                #ival = o
        #Create objects for summary file
        recalculated_XY.append(rcalc_coord[o])
        integer_hkl.append(pos[o])
        dist_difference.append(dist_dfr)
    return recalculated_XY,dist_difference,integer_hkl  


def UnitCellExtraction(Amatrix, compare = False, ideal_values=[5,5,5,90,90,90], max_deviation=10):
    """
    This function extracts unit cell parameters from A* matrix calculated by Graph-search algorithm.
    Parameters:
        Amatrix      - 3x3 matrix from Graph-search (or similar) algorithm.
        compare      - if True, enables comparison with idealized cell. ONLY CELLS IN DEVIATION RANGE WILL BE RETURNED.
        ideal_values - the ideal values of unit cell, used for comparison. Used only if compare is True.
        max_ deviation    - the maximal % of deviation between calculated value and ideal one. Used only if compare is True.
    Output:
        cells        - If compare is False, returns all cells. If its True, returns list of unit cells in acceptable range.
        filtered     - Returned only if compare is True. List of matrixes in the defined range.
    """
    from numpy import sqrt, arccos, dot, array, pi, linalg
    def unit_vector(vector):
        return vector / linalg.norm(vector)
    
    ia,ib,ic = ideal_values[0], ideal_values[1], ideal_values[2]
    ial, ibe, iga = ideal_values[3], ideal_values[4], ideal_values[5]
    cells = []
    filtered = []
    for s in Amatrix:
        x1,x2,x3 = s[0][0], s[0][1], s[0][2]
        y1,y2,y3 = s[1][0], s[1][1], s[1][2]
        z1,z2,z3 = s[2][0], s[2][1], s[2][2]

        a = sqrt((x1**2)+(x2**2)+(x3**2))
        b = sqrt((y1**2)+(y2**2)+(y3**2))
        c = sqrt((z1**2)+(z2**2)+(z3**2))

        al = arccos(dot(array(unit_vector([y1,y2,y3])),array(unit_vector([z1,z2,z3])) )) * 180/pi
        be = arccos(dot(array(unit_vector([x1,x2,x3])),array(unit_vector([z1,z2,z3])) )) * 180/pi
        ga = arccos(dot(array(unit_vector([x1,x2,x3])),array(unit_vector([y1,y2,y3])) )) * 180/pi    

        pa  = abs(a-ia)/ia
        pb  = abs(b-ib)/ib
        pc  = abs(c-ic)/ic
        pal = abs(al-ial)/ial
        pbe = abs(be-ibe)/ibe
        pga = abs(ga-iga)/iga

        tempcell = [a,b,c,al,be,ga]
        
        d = max_deviation/100
        
        if compare == True:
            if (pa < d) and (pb < d) and (pc < d) and (pal < d) and (pbe < d) and (pga < d):
                cells.append(tempcell)
                filtered.append(s)
        else:
            cells.append(tempcell)
            
    if compare == True:
        return cells, filtered
    else :
        return cells


#Class to store some instrument parameters
class Instrument:

    def __init__(self,pixel_size,detector_distance,detector_center,wavelength,extype) -> None:
        
        self.pixel_size        = pixel_size
        self.detector_distance = detector_distance
        self.detector_center   = detector_center
        self.wavelength        = wavelength
        self.extype            = extype


#====================================__CCTBX_UTILITIES__===============================================#
#These functions were writen by me to work with some CCTBX files.
#Hence, use of these functions is very limited.

def extract_from_dials_show(path):
    """ 
    Generate DataFrame from dials.show output files.
    Parameter: path - full path to location with dials.show outputs
    Returns: Dataframe with peak parameters
    """
    
    import re
    import pandas as pd
    import numpy as np
    import KrystalograPy as kp
    
    data = kp.list_files(path,'txt',deep=False)
    DF = pd.DataFrame()

    for d in data:
        content = []
        with open(d,'r') as file:
            for line in file:
                content.append(line)

            source = content[5].strip().replace("reflections = ","")
            #Fix the header
            headers = re.split("\s+",content[42].strip())
            headers[0] = "h"
            headers.insert(1,"k")
            headers.insert(2,"l")
            headers[20] = 'xcal.mm'
            headers.insert(21,'ycal.mm')
            headers.insert(22,'zcal.mm')
            headers[23] = 'xcal.px'
            headers.insert(24,'ycal.px')
            headers.insert(25,'zcal.px')
            headers[27] = 'xobs.mm.value'
            headers.insert(28,'yobs.mm.value')
            headers.insert(29,'zobs.mm.value')
            headers[30] = 'xobs.mm.variance'
            headers.insert(31,'yobs.mm.variance')
            headers.insert(32,'zobs.mm.variance')
            headers[33] = 'xobs.px.value'
            headers.insert(34,'yobs.px.value')
            headers.insert(35,'zobs.px.value')
            headers[36] = 'xobs.px.variance'
            headers.insert(37,'yobs.px.variance')
            headers.insert(38,'zobs.px.variance')
            headers[39] = 's11'
            headers.insert(40,'s12')
            headers.insert(41,'s13')
            headers.insert(42,'source')
            frame = [re.split("\s+",c.strip().replace(",","")) for c in content[43:]]
            frame = [[*f,source[20:]] for f in frame]
            if len(DF) != 0:
                DF_temp = pd.DataFrame(frame,columns=headers)
                #DF = DF.append(DF_temp)
                DF = pd.concat([DF,DF_temp])
            else:
                DF = pd.DataFrame(frame,columns=headers)

    #Fix data frame index and set up optimal data types
    x = {h:(float if h!='source' else str) for h in headers }
    DF.reset_index(inplace=True,drop=True)
    DF = DF.astype(x)
    DF = DF.astype({'h': int,'k': int, 'l': int, 'd': float, 'id':int,'panel':int,'flags':int})

    DF['d.mm'] = np.sqrt( ((DF['xcal.mm']-DF['xobs.mm.value'])**2) + 
                          ((DF['ycal.mm']-DF['yobs.mm.value'])**2) + 
                          ((DF['zcal.mm']-DF['zobs.mm.value'])**2) )

    DF['d.px'] = np.sqrt( ((DF['xcal.px']-DF['xobs.px.value'])**2) + 
                          ((DF['ycal.px']-DF['yobs.px.value'])**2) + 
                          ((DF['zcal.px']-DF['zobs.px.value'])**2) )

    a, c = 4.77, 12.99
    DF['res.A'] = 1/(np.sqrt((4/3) * ((DF.h**2 + DF.h*DF.k + DF.k**2)/a**2) + (DF.l**2)/(c**2)))

    DF['bmean'] = DF['background.sum.value']/DF['num_pixels.background_used']
    DF['imean'] = DF['intensity.sum.value']/DF['num_pixels.foreground']
    return DF
  
def refl_seek(h,k,l,DFm):
    """
    Extract data for peak with given hkl index.
    Parameters: h, k, l - index of reflection
                DFm     - datframe generated by extract_from_dials_show
    """
    DF_x = DFm.loc[(DFm["h"]==h)&(DFm["k"]==k)&(DFm["l"]==l)][['h','k','l','background.mean', 'background.sum.value',
           'background.sum.variance', 'intensity.sum.value','intensity.sum.variance', 'num_pixels.background',
           'num_pixels.background_used', 'num_pixels.foreground', 'num_pixels.valid','source','bmean','imean',
                                                               'res.A','d.mm','d.px']]
    return(DF_x)

