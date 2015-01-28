Calibration/verification model directory archive for USGS SIR titled:

  "Hydrologic Conditions in Urban Miami-Dade County, Florida and the Effect of 
   Groundwater Pumpage and Increased Sea Level on Canal Leakage and Regional 
   Groundwater Flow"
   
   By Joseph D. Hughes and Jeremy T. White
   
   07/28/2014
   08/11/2014 Updated based on Swain, E.D. review of model archive


Description: Input and output files for the calibration/verification model.
             Model results and post-processed results are in the Results\
             subdirectory. All other directories contain input data for the
             calibration/verification model or post-processing routines.
                 

System requirements:
  
  The models contained in this archive were run using the MD_mfnwt_x64.exe 
  executable in the bin\ directory in this model archive. MD_mfnwt_x64.exe
  was compiled using the Intel(R) Visual Fortran 64-bit compiler 
  (version 14.0.0092.11) for Windows and Visual C++ 2012.
  
  The models requires approximately 60 MB of available Random Access 
  Memory (RAM).
  
  The model has been run successfully on computers running the following
  Windows operating systems (OS):
  
    o 64-bit Windows 7 OS (Service Pack 1)
    o 64-bit Windows Server 2008 HPC OS (Service Pack 1)
    o 64-bit Windows Server 2012 R2 Standard
    
  The source code has also been successfully compiled using the gfortran 
  compiler and run on computers using the Mac OS X Mavericks (10.9) and 
  SUSE Linux operating systems.
                   

Running the model:

  The calibration/verification model can be run by typing the following command
  in a command line opened in this directory:
  
    ..\bin\MD_mfnwt_x64.exe UMD_fb.nam

  NOTE1: The model must be run using the MD_mfnwt_x64.exe executable
         in the bin\ directory in this model archive. The ZoneBud_x64.exe 
         executable must be used to used the zone budget post-processing 
         routines included in this archive.

  NOTE2: The files in the MET_DATA\ directory are required to run the models in 
	       the UMD.final directory. The MET_DATA\ directory must be located at the 
	       same level as the UMD.final directory so that relative path statements 
	       in the *.nam file are correct (..\MET_DATA\).
    

Subdirectories:

  .\
  
    Description: This directory contains the primary MODFLOW-NWT input files
                 for the calibration/verification model. The file UMD_fb.nam
                 is a catalog of files opened by the MODFLOW-NWT executable
                 (excluding files opened using an open/close statment). Example
                 of MODFLOW-NWT Package file extensions are (*.bas, *.dis, *.drn,
                 *.ets, *.gfb, *.ghb, *.hyd, *.lpf, *.nwt, *.oc, *.swi, *.swr, 
                 *.wel). This directory also includes the ZONEBUDGET input file
                 (UMD.zbg), ZONEBUDGET control file (zbg.in), and a python script 
                 for running ZONEBUDGET (runZoneBudget.py).
               
  
  BNDLIST\
  
    Description: ASCII list files containing daily data used by the MODFLOW-NWT
                 well package (in the BNDLIST\WEL\ subdirectory).
                 
               
  
  bndref\
  
    Description: ASCII and BINARY files containing daily data used by the MODFLOW-NWT
                 DRN, GFB, GHB, and SWR Packages/Processes.
    
    Sub-subdirectories:
    
      ag\
      
        Description: BINARY array files with daily estimated agricultural water
                     use included in the GFB file. Files include a standard MODFLOW
                     header and single-precision real data (NCOL by NROW).
                     
                     Daily agricultural water use data for every day between 
                     1/1/1996 and 12/31/2010 are included in this subdirectory.
                     
                       ag19960101.bin
                       ag19960102.bin
                       ag19960103.bin
                               .
                               .
                               .
                       ag20101229.bin
                       ag20101230.bin
                       ag20101231.bin
        
      
      drn_Obs\

        Description: BINARY list files with daily drain package data. Data
                     for each drain are written as single-precision real data 
                     (Layer Row Column Elevation Cond).
                     
                     Daily drain package data for every day between 1/1/1996 
                     and 12/31/2010 are included in this subdirectory.
                     
                       UMD_DRN_Obs_19960101.bin
                       UMD_DRN_Obs_19960102.bin
                       UMD_DRN_Obs_19960103.bin
                               .
                               .
                               .
                       UMD_DRN_Obs_20101229.bin
                       UMD_DRN_Obs_20101230.bin
                       UMD_DRN_Obs_20101231.bin
      
      flow\

        Description: ASCII SWR1 time series files used to define specified
                     flows for surface-water pumps included in the model.
                     See appendix 2 in doc\Hughes_TM6-A40.pdf for more information
                     on ASCII SWR1 time series files. Data for the following
                     time series are included in the model:

                        S151_C_SWR1.ref
                        S25B_P_SWR1.ref
                        S26_P_SWR1.ref
                        S331_P_SWR1.ref
                        S331_SIPH_SWR1.ref
                        S332B1_P_SWR1.ref
                        S332B2_P_SWR1.ref
                        S332C_P_SWR1.ref
                        S332DX1_C_SWR1.ref
                        S332D_P_SWR1.ref
                        S332_P_SWR1.ref
                      
                      Other files that are not used have been included for 
                      convenience to allow for possible future use.
      
      gate\

        Description: ASCII SWR1 time series files used to define effective gate
                     openings for operable movable-crest weirs and underflow
                     gates included in the model. See appendix 2 in 
                     doc\Hughes_TM6-A40.pdf for more information on ASCII SWR1 
                     time series files. See UMD_FB.nam for list of specific
                     SWR1 time series files that are being used. Other files that 
                     are not used have been included for convenience to allow for
                     possible future use.
      
      ghb_Obs\

        Description: BINARY list files with daily general head boundary (GHB) 
                     package data. Data for each drain are written as single-
                     precision real data (Layer Row Column Bhead Cond).
                     
                     Daily drain package data for every day between 1/1/1996 
                     and 12/31/2010 are included in this subdirectory.
                     
                       UMD_GHB_Obs_19960101.bin
                       UMD_GHB_Obs_19960102.bin
                       UMD_GHB_Obs_19960103.bin
                               .
                               .
                               .
                       UMD_GHB_Obs_20101229.bin
                       UMD_GHB_Obs_20101230.bin
                       UMD_GHB_Obs_20101231.bin
      
      rec\

        Description: BINARY array files with daily estimated recreational 
                     irrigation included in the GFB file. Files include a 
                     standard MODFLOW header and single-precision real data 
                     (NCOL by NROW).
                     
                     Daily recreational irrigation data for every day between 
                     1/1/1996 and 12/31/2010 are included in this subdirectory.
                     
                       rec19960101.bin
                       rec19960102.bin
                       rec19960103.bin
                               .
                               .
                               .
                       rec20101229.bin
                       rec20101230.bin
                       rec20101231.bin

      
      stage\
      
        Description: ASCII SWR1 time series files used to define specified
                     surface-water stages used in the model.
                     See appendix 2 in doc\Hughes_TM6-A40.pdf for more information
                     on ASCII SWR1 time series files. Data for the following
                     time series are included in the model:

                        VAKey_SWR1.ref
                      
                      Other files that are not used have been included for 
                      convenience to allow for possible future use.
      
      strcrit\
                 
        Description: ASCII SWR1 time series files used to define specified
                     surface-water criterion used to operate surface-water 
                     control structures in the model.
                     See appendix 2 in doc\Hughes_TM6-A40.pdf for more information
                     on ASCII SWR1 time series files. These data are not currently
                     used in the model but have been included for convenience to 
                     allow for possible future use.

  georef\
  
    Description: ASCII files cross-section geometry information for each 
                 unique SWR1 surface-water geometry item in the model. See 
                 appendix 2 in doc\Hughes_TM6-A40.pdf for more information
                 on ASCII SWR1 cross-section data. 
                 
                 The following cross-section data are included in this 
                 subdirectory.
                     
                   SWR_IGEONUM0001.ref
                   SWR_IGEONUM0002.ref
                   SWR_IGEONUM0003.ref
                           .
                           .
                           .
                   SWR_IGEONUM1007.ref
                   SWR_IGEONUM1008.ref
                   SWR_IGEONUM1009.ref

  GIS\
  
    Description: ESRI shapefiles that are used in several post-processing routines.
                 The shapefiles have been included for convenience to facilitate 
                 future use and have been given descriptive names. The ESRI 
                 shapefile technical description has been included in 
                 doc\shapefile.pdf for future users unfamiliar with the shapefile 
                 format (hello unknown future user). 


  obsref\
  
    Description: ASCII smp observed time series files used to compare with 
                 simulated model results. smp files have the following 
                 format:

                    DBL2_C 01/01/1996 12:00:00 528474.764093
                    DBL2_C 01/02/1996 12:00:00 528474.764093
                    DBL2_C 01/03/1996 12:00:00 538261.333798
                                       .
                                       .
                                       .
                    DBL2_C 12/30/2011 12:00:00 39146.2788217
                    DBL2_C 12/31/2011 12:00:00 29359.7091163
                    DBL2_C 01/01/2012 12:00:00 29359.7091163

                 
                 All smp files use model units of meters and cubic meters per day.
                 These units may be converted to other units by post-processing
                 routines.
    
    Sub-subdirectories:
    
      flow\
      
        Description: ASCII smp files with observed discharge results at surface-water
                     control structures in the active model domain. The data in this
                     subdirectory have units of cubic meters per day and are 
                     consistent with model length and time units.

      head\
      
        Description: ASCII smp files with observed groundwater heads at ground-
                     water monitoring locations used to calibrate/verify the model.
                     The data in this subdirectory have units of meters and are 
                     consistent with model length units.

      stage\
      
        Description: ASCII smp files with observed surface-water stages at stage
                     gages used to calibrate/verify the model. The data in this 
                     subdirectory have units of meters and are consistent with 
                     model length units.

                                  
  python\
  
    Description: python scripts used to post-process model results. These
                 scripts have been included for convenience to facilitate 
                 future use and have been given descriptive names. In order
                 to use the scripts Python 2.7 and the Numpy, Scipy, matplotlib,
                 ElementTree, and cluster_repo packages need to be installed.
                 The cluster_repo package is a custom package and has been included
                 in code\cluster_repo. 
                 
                 Instructions for installing cluster_repo in the python 2.7 
                 Lib\site-package\ subdirectory are given in README.txt in the
                 code\ directory.
                 
                 The post-processing python scripts are not required to run the 
                 UMD.final model or extract model results. Model results can be
                 extracted using standard MODFLOW or user-preferred post-
                 processing utilities.
                 
                 A brief description of the scripts are given below:
                 
                   Script                             Description
                   --------------------------------   ---------------------------------
                   UMDModel.py                        Dimensions and location of model
                   UMDUtils.py                        Common functionality used by one
                                                        or more post-processing scripts
                   ProcessUMDFinalStageHeadZeta.py    Extracts final stage, head, and
                                                       zeta surface from model. Use to
                                                       define initial condition for 
                                                       scenario simulations
                   ProcessUMDMFAnnualMET.py           Create maps showing simulated
                                                       meterological data (rainfall,
                                                       evapotranspiration, etc.).
                                                       Also creates summary tables of
                                                       simulated results.
                   ProcessUMDMFBOTM.py                Create maps showing bottom of 
                                                        model layers
                   ProcessUMDMFETSData.py             Create maps showing calibrated
                                                        ETS data file
                   ProcessUMDMFGeoCrosssections.py    Create cross-sections along
                                                        defined lines 
                                                        (defined in GIS\UMDCrossSections)
                   ProcessUMDMFGWProperties.py        Create maps showing hydraulic
                                                        conductivity and storage properties
                                                        used in the model
                   ProcessUMDMFHeads.py               Create maps of simulated water table
                                                        elevations
                   ProcessUMDMFObsHead.py             Create graphs showing simulated and
                                                        observed groundwater heads
                   ProcessUMDMFZeta.py                Create maps of simulated zeta 
                                                        surfaces
                   ProcessUMDMFZoneBudgetResults.py   Create tables with simulated 
                                                        ZONEBUDGET results
                   ProcessUMDSWRMannLeak.py           Create shapefiles with Manning's
                                                        n and canal bed leakance values
                   ProcessUMDSWRObsBudget.py          Create graphs showing simulated and
                                                        observed month net surface-water
                                                        subbasin discharge
                   ProcessUMDSWRObsFlow.py            Create graphs showing simulated and
                                                        observed surface-water discharge
                   ProcessUMDSWRObsStage.py           Create graphs showing simulated and
                                                        observed surface-water stage
                   ProcessUMDSWRUrbanBudget.py        Create table with summary of surface-
                                                        water budget for urban areas

  REF\
  
    Description: ASCII files containing arrays used to define two-dimensional data for
                 the model. These files are opened using open\close keywords in
                 MODFLOW Package files. Some of the files included in this directory
                 are not used directly by the model but have been been included for 
                 convenience to facilitate post-processing and to allow for possible 
                 future use. 
                                  
                                  
  Results\
  
    Description: Simulated results for the calibration/verification model.
    
    Sub-subdirectories:
    
      .\
      
        Description: Model results, ZONEBUDGET results, and post-processed tables created
                     by post-processing scripts in the python\ subdirectory. The post-
                     processed tables have been included for convenience. Model results 
                     can be extracted using standard MODFLOW or user-preferred post-
                     processing utilities.
      
      Figures\
      
        Description: Post-processed model results. These results are created using the
                     scripts in the python\ subdirectory and have been included for
                     convenience.

      GIS\
      
        Description: Shapefiles created by post-processing scripts in the python\ 
                     subdirectory. These shapefiles have been included for convenience.
                     The shapefile have been given descriptive names. The ESRI  
                     shapefile technical description has been included in 
                     doc\shapefile.pdf for future users unfamiliar with the shapefile 
                     format. 

  SWRREF\
  
    Description: ASCII files containing SWR1 Process data that are read by the SWR1 
                 Process. These files are either read using the EXTERNAL keyword and are 
                 defined in the MODFLOW-NWT name file (..\UMD_fb.nam) or using the
                 open\close keyword. Descriptive file names are used and define the
                 SWR1 Process dataset contained in the file. See appendix 2 in 
                 doc\Hughes_TM6-A40.pdf or doc\SWRProcessInputInstructions_v1.03.pdf
                 for  more information on the specific format of each of the SWR1 
                 data sets.

  xml\
  
    Description: xml files that are used by the post-processing python scripts in
                 the python\ subdirectory. These xml files are control files that
                 are used to define time series files used for particular
                 observations, time periods to evaluate, time sampling period, etc.
                 Descriptive tags have been added to make the xml files human
                 readable. The xml files have been included for convenience to 
                 facilitate future use.
                                  