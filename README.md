# Small Scripts For Crystallography

This is a collection of functions I use in my crystallographic work. Each function has docstring with explanation. At the moment there are two modules:  
KrystalograPy - mostly handles structural information and CIF files.   
  
FrameUtil     - made to read SACLA H5 images and show them in human redable form (not chopped, but assembled).  
  
Prepare.ps1   - Powershell script to extract data from SPring-8 files, to automatize partially import to CrysAlisPro. This version was written for output from BL02B1 beamline, when experimental data consisted of PHI.log, TIFF and INF files. If beamline changes the experimental output in the future, modifications to this script will be needed.  
  
Most functions require either CIF file or diffraction frames to work.

<b>Some functions are not finished yet. In such cases, I indicate this in docstring.</b>
