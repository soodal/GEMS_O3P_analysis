; IDL Version 8.4.1 (linux x86_64 m64)
; Journal File for soodal@srsl115
; Working directory: /home/soodal/works/GEMS_O3P_analysis
; Date: Fri Oct  8 15:50:49 2021
 
merra2gemso3p.latitude[50, 450]
; % Expression must be a structure in this context: MERRA2GEMSO3P.
merra2gemsl2o3p.latitude[50, 450]
;     -0.35180873
merra2gemsl2o3p.effectivecloudfractionuv[50, 450]
;      0.50000000000000000
merra2gemsl2o3p.effectivecloudfractionuv[51, 450]
;      0.50000000000000000
merra2gemsl2o3p.effectivecloudfractionuv[53, 450]
;      0.50000000000000000
merra2gemsl2o3p.effectivecloudfractionuv[63, 450]
;      0.50000000000000000
merra2gemsl2o3p.effectivecloudfractionuv[74, 450]
;      0.50000000000000000
merra2gemsl2o3p.effectivecloudfractionuv[100, 450]
;      0.50000000000000000
ecf = merra2gemsl2o3p.effectivecloudfractionuv[100, 450]
ecf = merra2gemsl2o3p.effectivecloudfractionuv
showme, ecf
;================================================================================
;Double-precision floating: 174 x 512, 2 dimension.
;--------------------------------------------------------------------------------
;Max   value      0.50000000
;Mean  value  -5.2835400e+28 STDEV value    2.2370602e+29
;Min   value  -1.0000000e+30
;    Finite number count        89088
;Not Finite number count            0
;Finite number Ratio       1.00000
;Not a Number Ratio        0.00000
;First 6 elements and Last 6 elements are
;                                                                                
;  -1.0000000e+30  -1.0000000e+30                     . . .                   -1.0000000e+30  -1.0000000e+30
;                                                                 0.50000000  -1.0000000e+30
;                                     . . .                                      
;  -1.0000000e+30      0.30000001
;  -1.0000000e+30  -1.0000000e+30                     . . .                   -1.0000000e+30  -1.0000000e+30
;--------------------------------------------------------------------------------
.reset
.r test0077_test_max_finder
; % Error opening file. File: test0077_test_max_finder
.r test0077_test_local_max_finder
; % Not a legal system variable: !MAGIC.
