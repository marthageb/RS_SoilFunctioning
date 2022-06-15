#!/usr/bin/python3
#python analyze_reflect.py --infile_prefix Z:\SRER\Martha\hyperspec\flightlines\output\noNDVI_thick_sparse\NEON --outfile FoliarN\results.h5

doc = """
Primer on mapInfo and UTM...  Assuming SRER/Reflectance/Metadata/Coordinate_System/Map_Info is the following:

UTM,  1.000,  1.000,  505834.000,  3515363.000,  1.0000000000e+000,  1.0000000000e+000,  12,  North,  WGS-84,  units=Meters, 0

We assign these text fields into the mapInfo dictionary as follows:

            mapInfo = {
                'projection_name'       : mapInfo_split[0].strip(),             # usually 'UTM'
                'reference_pixel_x'     : float(mapInfo_split[1].strip()),      # not sure of use
                'reference_pixel_y'     : float(mapInfo_split[2].strip()),      # not sure of use
                'pixel_easting'         : float(mapInfo_split[3].strip()),      # xMin - westernmost boundary (array is x-ascending order)
                'pixel_northing'        : float(mapInfo_split[4].strip()),      # yMax - northernmost boundary (array is y-descending order)
                'x_pixel_size'          : float(mapInfo_split[5].strip()),      # x resolution
                'y_pixel_size'          : float(mapInfo_split[6].strip()),      # y resolution
                'utm_zone'              : int(mapInfo_split[7].strip()),        # for SRER will always be 12
                'utm_ns'                : mapInfo_split[8].strip(),             # 'North' or 'South'... for northern-hemisphere, always 'North'
                'datum'                 : mapInfo_split[9].strip(),             # usually 'WGS-84' - DOD GPS standard
                'units'                 : mapInfo_split[10].strip(),            # usually 'units=Meters'
                'rotation_angle'        : float(mapInfo_split[11].strip()),     # usually 0
                }
To display this during run, specify -v parameter

The "projection_name" mapInfo key value says that we are using UTM (Universal Transverse Mercator) coordinate system which depends on:
    1) Which UTM zone you are in (per above, we are in mapInfo keys "utm_zone-utm_ns" or 12-North).
    2) Which reference datum you are you using in this zone (per above, we are using mapInfo key "datum" value of WGS-84 which is the DOD reference system for GPS).
       This defines the absolute southwest point that all "northing" and "easting" displacements are measured from.  North America also, less-commonly, uses NAD-27.
    3) The first, or X number of a UTM coordinate is the "easting" which is usually the number of meters east of your datum
       (per above, we have mapInfo key "pixel_easting" value 505834.000).  IMPORTANT... For Neon this is the "minimum" or west-most X boundary.
    4) The second, or Y number of a UTM coordinate is the "northing" which is usually the number of meters north of your datum
       (per above, we have mapInfo key "pixel_northing" value 3515363.000).  IMPORTANT... For Neon, this is the "maximum" or north-most Y boundary.

Neon uses Y:X:Z arrangement of data in reflectance_data_ds as follows... 
Y axis is in Y-descending sequence starting at the mapInfo key "pixel_northing" value 3515363.000
X axis is in X-ascending sequence starting at the mapInfo key "pixel_easting" value 505834.000
Z axis is the 400+ wavelength reflectances measured at this Y:X coordinate
"""

import csv
import argparse
import sys
##### Parse arguments
parser = argparse.ArgumentParser(description='Reflectance Analyzer', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
arg_group = parser.add_mutually_exclusive_group(required=True)
arg_group.add_argument('--infile',          type=str, dest='h5_input_file',                             help='Input h5 filename')
arg_group.add_argument('--infile_prefix',   type=str, dest='h5_input_file_prefix',                      help='Input h5 filename prefix')
arg_group.add_argument('--doc',             action='store_true',                                        help='Print documentation on UTM and Neon MapInfo')

parser.add_argument('--outfile',        type=str,       dest='h5_output_file', default='results.h5',    help='Output h5 filename')
parser.add_argument('--numprocess',     type=int,                                                       help='Overrides the calculated number of subprocesses to run')
parser.add_argument('--verbose', '-v',  action='count', default=0,                                      help='Verbosity of output "-v" or "-vv" for even more')
parser.add_argument('--status_interval',type=int,       default=5,                                      help='The interval (in minutes) when each subtask prints status message - 0=no messages')
parser.add_argument('--csv',            action='store_true',                                            help='Generate csv in addition to h5 output')

parser.add_argument('--y',              type=int,                                                       help='"y" value if only wanting detail for single row')
parser.add_argument('--x',              type=int,                                                       help='"x" value if only wanting detail for single row')
numba_group = parser.add_mutually_exclusive_group(required=False)
numba_group.add_argument('--no_use_numba',   action='store_true',                                       help='Do not use numba jit for matrix operations')
numba_group.add_argument('--use_numba_cuda',   action='store_true',                                     help='Use numba cuda jit for matrix operations - requires nvidia gpu')
parser.add_argument('--timings',        action='store_true',                                            help='collect array operation timings')
args = parser.parse_args()
if args.doc:
    print(doc)
    sys.exit()

from multiprocessing import Pool, cpu_count
import h5py
from scipy.signal import savgol_filter
import pprint
import time
import pdb
import numpy as np
import re
import os
import datetime
import hashlib
from collections import namedtuple
from termcolor import colored, cprint
import colorama
if args.use_numba_cuda:
    cuda_base = "C:\\Program Files\\NVIDIA GPU Computing Toolkit\\CUDA\\v10.0"
    cuda_env = {
        "NUMBAPRO_CUDA_DRIVER"  : "C:\\windows\\system32\\nvcuda.dll",
        "NUMBAPRO_NVVM"         : cuda_base + "\\nvvm\\bin\\nvvm64_33_0.dll",
        "NUMBAPRO_LIBDEVICE"    : cuda_base + "\\nvvm\\libdevice",
        "CUDA_VISIBLE_DEVICES"  : "0",
        }
    for key, val in cuda_env.items():
        pass
        os.environ[key] = val
        print('set {}={}'.format(key,val))


reflectance_ds_name = 'SRER/Reflectance'
reflectance_data_ds_name = 'BRDF/Correction'
script_start_time = datetime.datetime.now()

def timestamp():
    return time.strftime('%H:%M:%S', time.localtime(time.time()))

###### Define function to convert index->coordinate
coord_tuple = namedtuple('coord_tuple',['x','y'])
def calc_coord(x_sub, y_sub, mapInfo):
    ##################################################################################################
    # IMPORTANT:  See doc block at beginning of this script detailing mapInfo and UTM usage
    ##################################################################################################
    # The Reflectance_Data dataset is organized as a Y,X,Z array
    # Y axis is in Y-descending sequence starting at Reflectance/Metadata/Coordinate_System/Map_Info[4] which is in mapInfo['pixel_northing']
    # X axis is in X-ascending sequence starting at Reflectance/Metadata/Coordinate_System/Map_Info[3] which is in mapInfo['pixel_easting']
    ##################################################################################################
    x = float(mapInfo['pixel_easting']  + (x_sub * mapInfo['x_pixel_size']))        # gives x coordinate via positive displacement from pixel_easting/xMin
    y = float(mapInfo['pixel_northing'] - (y_sub * mapInfo['y_pixel_size']))        # gives y coordinate via negative displacement from pixel_northing/yMax
    return coord_tuple(x,y)

#reg_coeff_list = [0,0,0,0,0,0,0,0,0,0,0,-1.477479419,-11.3331926,-10.94994106,-6.019146564,-8.974330121,-4.093858536,-0.618652025,-9.535672496,-9.308052279,-2.130628387,-3.928565114,-5.209352102,-10.18442702,-11.2324472,-8.617415704,-10.91546154,-10.40709097,-9.809190561,-10.91820924,-11.71126793,-6.387920853,-2.286901227,-5.474094385,-12.46710644,-7.357762347,-3.409187407,-6.031183434,-6.074014977,-4.094906774,-1.615777779,-9.221004827,-8.991123754,-1.168331962,-6.675592366,-5.994222272,-0.65754865,-4.076150719,-3.283082796,-4.195463506,-4.872447313,1.04338934,-3.038682292,-2.648941731,4.377918289,6.67289889,4.358633525,1.072362628,-4.00852391,-2.23714642,-3.435186005,-20.90929233,-17.69223883,-10.75708428,3.72155717,21.04452923,16.47756072,15.08713967,13.10842931,4.665190617,0.018229015,4.528115626,2.36804487,-0.248117201,-8.331787737,-27.10883044,-2.56248429,18.77844928,9.084333317,1.03038007,-0.96015422,-2.167448143,-2.646525225,2.238229485,1.551759008,-2.00260465,-3.231356352,4.300054568,10.85752738,0.144466058,-6.9847028,4.208438879,4.51211189,0.493640828,1.128306516,3.878899084,0.910527757,3.57393949,-1.019929203,-7.239410176,4.080577147,5.641055727,1.535759837,-5.039926317,-6.660639206,0.911094603,-1.439061696,1.190939497,3.720695349,7.123774049,2.591752838,0.324242537,14.68977708,-4.362039584,8.744848384,11.10163101,-9.072260726,3.090588005,4.763364336,0.809996893,0.195507738,-0.184544076,1.026635458,3.818750622,1.180457516,-1.938614404,-5.764106698,-1.003528727,2.106511166,0.994390897,3.324617753,-4.454964988,3.318434656,-3.884239266,-4.591643496,2.413991886,-8.575777334,0.044257335,5.526675004,4.226329159,-3.551654887,-4.146653226,5.530095629,-1.894747024,-0.275512724,4.711381843,3.07982569,-4.797508271,-1.138831815,3.797505363,7.323540394,-3.137185689,-0.767506827,12.0659501,8.065403574,0.594218104,-12.4220461,1.821625892,-1.030757271,-3.71576381,-4.303670994,-4.52212954,-6.66329654,-6.039823179,2.414257427,4.613732355,6.541022159,6.91435898,4.16264615,1.354855715,8.151667207,3.097864755,-3.185056531,-5.560316571,-6.399130015,-1.554572762,-4.40779125,-8.702575552,4.045565957,15.59372692,6.339992952,1.763584048,2.985177934,-0.656808797,0.790244574,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-4.950771724,1.675186734,-2.339328085,2.046425703,2.624254717,-5.371734517,-10.54150924,-0.863305907,-3.092781621,-6.850107081,-1.379557869,-3.781795573,-6.396216757,-3.881520867,-0.508191529,1.55581903,-1.632762372,-4.214811535,-10.41678676,-5.222579272,2.830372326,-3.932380648,-3.625189784,-3.836854084,-7.734786727,-14.96204576,-10.80221067,3.256123098,-4.278970343,-7.350017948,-3.070076853,-5.627008215,-7.720439707,-5.896064548,-5.13201971,-7.991148485,-6.779263119,-13.22300406,-11.29562179,-8.570943265,-8.608725655,-5.802619993,-6.358729601,-4.137327515,-14.32490117,-16.25119247,-8.003990018,-8.788486502,-4.408897842,-4.244809114,0.217211024,8.313376773,17.86285789,6.542428754,2.83332887,15.4690469,7.801009263,15.45092982,24.17336104,10.15141906,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5.47593095,7.273733506,1.132103776,8.95037083,16.22369201,5.058510212,9.157905401,12.55595607,4.512254264,3.319766551,10.86862457,6.137755286,-11.00064807,-7.162366739,0.747702215,-13.00571408,-11.02857464,-1.942384144,-8.773086065,-8.801975076,-15.47338017,-21.40510617,-4.774451699,-3.21431088,-8.110271485,-12.69240837,-1.361546824,1.572822367,-10.48595801,-2.633190107,3.975354404,-0.172520615,-9.1032043,3.412345842,9.869080098,-8.170413914,-4.971034351,16.40185834,12.79604258,2.282207617,7.881345793,12.39481018,-0.573512645,-8.899912373,0.148833923,5.111740301,-5.525913243,-8.9018296,2.348820558,-6.968629192,-13.06697038,-2.482166313,0.170783119,3.644279003,13.32213092,24.27672864,25.96841046,16.00160489,-1.815189982,-2.2341157,-2.101349714,-9.107594853,6.370853807,11.64859499,5.030429017,2.872926123,-4.397563912,-4.342109229,1.127145962,7.274365624,10.58714391,0.931769533,-17.8725291,-8.186274328,14.78881542,20.77281135,8.668399429,1.340477397,21.47232456,17.88014474,11.23953084,-5.360210629,0,0,0,0,0,0,0,0,0,0,0,0,0]
reg_coeff_list = []

with open("PLSRcoef.csv", newline="") as csvfile:
    reader = csv.reader(csvfile)
    for rownum, row in enumerate(reader):
        if rownum == 0:
            if row[1] != "PLSRcoef":
                raise Exception("PLSRcoef not found in column 2")
        else:
            reg_coeff_list.append(float(row[1]))

reg_coeff_ary = np.array(reg_coeff_list, np.float32)
timings = {
    'count'             : 0,
    'cell_transform'    : 0,
    'savgol'            : 0,
    'calc'              : 0,
    }
    
if not args.no_use_numba:
    import numba
    if args.use_numba_cuda:
        from numba import cuda
        @cuda.jit(device=True, fastmath=False)
        def setit(x_cell_val, constants):
            if x_cell_val == constants[0]:
                return np.nan
            else:
                return x_cell_val / constants[1]
                
        @cuda.jit()
        def massage(x_data_ary, z_ary, constants):
            i = cuda.grid(1)
            z_ary[i] = setit(x_data_ary[i],constants)
            constants[8] = float(constants[4] + (constants[2] * constants[6]))        # gives x coordinate via positive displacement from pixel_easting/xMin
            constants[9] = float(constants[5] - (constants[3] * constants[7]))        # gives y coordinate via negative displacement from pixel_northing/yMax
            """
            for i in range(x_data_ary.shape[0]):
                if x_data_ary[i] == constants[0]:
                    z_ary[i] = np.nan
                else:
                    z_ary[i] = x_data_ary[i] / constants[1]
            """
            return
        stream = cuda.stream()
    else:
        from numba import njit
        @njit(cache=True)
        def massage(x_data_ary, z_ary, noDataValue, scaleFactor, x_sub, y_sub, pixel_easting, pixel_northing, x_pixel_size, y_pixel_size, coord_ary):
            for i in range(x_data_ary.shape[0]):
                if x_data_ary[i] == noDataValue:
                    z_ary[i] = np.nan
                else:
                    z_ary[i] = x_data_ary[i] / scaleFactor                   
            coord_ary[0] = float(pixel_easting  + (x_sub * x_pixel_size))        # gives x coordinate via positive displacement from pixel_easting/xMin
            coord_ary[1] = float(pixel_northing - (y_sub * y_pixel_size))        # gives y coordinate via negative displacement from pixel_northing/yMax
            return
            
###### Define function to calculate each row "result"  
def calc_val(x_data_ary, y_sub, x_sub, scaleFactor, noDataValue, mapInfo, diag=False):
    if args.timings:
        timings['count'] += 1
        t0 = time.time()
    if args.no_use_numba:
        z_ary = []
        for cell in x_data_ary:
            if cell == noDataValue:
                z_ary.append(np.nan)
            else:
                z_ary.append(cell / scaleFactor)
        coords = calc_coord(x_sub, y_sub, mapInfo)
        x_coord = coords.x
        y_coord = coords.y
    elif args.use_numba_cuda:
        with stream.auto_synchronize():
            asize           = len(x_data_ary)
            d_x_data_ary    = cuda.to_device(x_data_ary, stream)
            z_ary           = np.empty(asize)
            d_z_ary         = cuda.to_device(z_ary, stream)
            constants       = [noDataValue,scaleFactor, x_sub, y_sub, mapInfo['pixel_easting'], mapInfo['pixel_northing'], mapInfo['x_pixel_size'], mapInfo['y_pixel_size'], float(0), float(0)]
            d_constants     = cuda.to_device(constants, stream)
            threadsperblock = 512
            blockspergrid   = ((asize * 2) + 2 + (threadsperblock - 1)) // threadsperblock
            massage[blockspergrid, threadsperblock, stream](d_x_data_ary, d_z_ary, d_constants)
            x_coord = constants[8]
            y_coord = constants[9]
    else:
        z_ary = np.empty(len(x_data_ary))
        coord_ary = np.empty(2)
        massage(np.float32(x_data_ary), z_ary, noDataValue, scaleFactor, x_sub, y_sub, mapInfo['pixel_easting'], mapInfo['pixel_northing'], mapInfo['x_pixel_size'], mapInfo['y_pixel_size'], coord_ary)
        x_coord, y_coord = coord_ary
        
    if args.timings:
        t1 = time.time()
        timings['cell_transform'] += t1 - t0
        t0 = t1
    #savgol_result = savgol_filter(x=z_ary, window_length=3, polyorder=1, deriv = 1)         # obtain savgol coefficients
    sqrtss = np.sqrt(np.sum(np.power(z_ary[~np.isnan(z_ary)],2))) # to get rid of wavelengths with nan vals
    #sqrtss = np.sqrt(np.sum(np.power(z_ary,2)))

    savgol_result = z_ary/sqrtss

    if args.timings:
        t1 = time.time()
        timings['savgol'] += t1 - t0
        t0=t1
    # multiply each savgol coefficient by above fixed coefficients, add all of them together, then add fixed value
    result = 3.17275573709091 + np.sum(np.multiply(np.nan_to_num(savgol_result),reg_coeff_ary))
    if args.timings:
        t1 = time.time()
        timings['calc'] += t1 - t0
    if diag:
        np.set_printoptions(suppress=True)
        print("*********************** subscript:[{}:{}] ***********************************".format(y_sub,x_sub))
        print("x_data_ary:{}\n".format(x_data_ary))
        print("z_ary:{}\n".format(z_ary))
        print("sqrtss:{}\n".format(sqrtss))
        print("savgol_result:{}\n".format(savgol_result))
        print("result:{}\n".format(result))
    return [y_coord, x_coord, result]
    
###### Define function to iterate through all y:x within y_start:y_end range
def iter_y(process_num, h5_input_file_full, y_start, y_end, scaleFactor, noDataValue, mapInfo):
    colorama.init()
    pid = os.getpid()
    with h5py.File(h5_input_file_full,'r') as f:                                        # must be re-opened within each subprocess
        reflectance_data_ds = f[reflectance_data_ds_name]
        num_y, num_x, num_z = reflectance_data_ds.shape
        num_y = y_end-y_start                                                           # number of y values
        cprint("{} Subprocess:{}/{} subsetting y[{}:{},:,:]".format(timestamp(), process_num, pid, y_start, y_end),'green',attrs=['bold'])
        subset = reflectance_data_ds[y_start:y_end,:,:]                                 # subset ds into memory with only specified y-range - highly i/o and memory intensive
    totrows = num_y * num_x                                                             # total rows in this y-range
    nan_compare_ary = np.full(subset.shape[2],noDataValue)                              # used later for comparison of all-nan
    results = []
    cprint("{} iter_y subprocess {}/{} now iterating y index {}:{} ({} columns)".format(timestamp(), process_num, pid, y_start, y_end-1, num_y),'green',attrs=['bold'])
    stime = time.time()                                                                 # initial iteration time used for status calculations below
    last_status_interval_processed = 0
    status_interval_seconds = args.status_interval * 60                                 # want status message every this-many seconds
    rows_processed_last_interval = 0
    for y_sub_orig, y_data_ary in enumerate(subset):                                    # iterate across the subset y array
        y_sub = y_sub_orig + y_start + 1                                                # this is the y-subscript adjusted to the entire dataset
        processed_y = y_sub_orig * num_x                                                # number of processed rows
        for x_sub, x_data_ary in enumerate(y_data_ary):                                 # now iterate across the second-level y_data x array
            if np.array_equal(x_data_ary, nan_compare_ary):                             # if entire x_data_ary is nan
                coords = calc_coord(x_sub, y_sub, mapInfo)
                results.append([coords.y,coords.x,np.nan])                                     # then result is nan
            else:                                                                       # otherwise we need to calculate
                results.append(calc_val(x_data_ary, y_sub, x_sub, scaleFactor, noDataValue, mapInfo))   # calc_val will return "result list"
            if args.status_interval:                                                    # only do stats if --status_interval is non-zero
                interval = time.time() - stime                                          # seconds since this subprocess started
                if not(int(interval % status_interval_seconds)):                        # do we need stats on this interval?
                    this_interval = int(interval // status_interval_seconds)            # calculate this interval
                    if last_status_interval_processed != this_interval:                 # have we already presented message for this interval?
                        last_status_interval_processed = this_interval                  # no then reset this_interval and present the status message
                        processed_rows = processed_y + x_sub                            # total rows processed
                        lper = int((processed_rows/totrows)*100)                        # percentage of all total rows that I have to process
                        rrate = int((processed_rows - rows_processed_last_interval) / status_interval_seconds) # this interval run-rate (rows/sec)
                        rows_processed_last_interval = processed_rows                   # for next time around
                        #rrate=int(processed_rows / interval)                            # run-rate in rows per second
                        remain_secs=int((totrows - processed_rows) / rrate)             # using past history project the remaining seconds
                        etc_text = time.strftime('%H:%M:%S', time.localtime(time.time() + remain_secs)) # get formatted wall-clock time of estimated completion
                        remain_secs_fmt = str(datetime.timedelta(seconds = remain_secs))# formatted hh:mm:ss of remaining seconds
                        cprint("{} Subprocess:{}/{} rrate:{} coordinates/sec ({:,}/{:,} or {}%) expected complete in {} at {}".format(timestamp(), process_num, pid, rrate, processed_rows, totrows, lper, remain_secs_fmt, etc_text),'green',attrs=['bold'])
                        if args.timings:
                            cprint("Timings: iters={:,.0f}, tot/avg times: cell_transform {:,.2f}/{:.6f}, savgol {:,.2f}/{:.6f}, calc {:,.2f}/{:.6f}".format(timings['count'], timings['cell_transform'], timings['cell_transform']/timings['count'], timings['savgol'], timings['savgol']/timings['count'], timings['calc'], timings['calc']/timings['count']),'green',attrs=['bold'])
    cprint("{} Subprocess {}/{} completed iteration in {}".format(timestamp(),process_num, pid, str(datetime.timedelta(seconds = time.time()-stime))),'green',attrs=['bold'])
    return results

###### Define function to break up y into separate subprocesses and invoke iter_y for each
def iter_reflect(h5_input_file_full, filepref, scaleFactor, noDataValue, mapInfo):
    stime = time.time()
    totrows = num_x_cols*num_y_cols                                                 # total rows in the dataset
    if args.numprocess:
        num_subprocesses = args.numprocess
    else:
        num_subprocesses = cpu_count()                                              # set to all available processors
        if cpu_count() == 2: num_subprocesses -= 1                                  # if only 2 processors then only use one of them
        elif 3 <= cpu_count() <=4: num_subprocesses -= 1                            # 3-4 processors:  leave one unused
        elif 5 <= cpu_count() <=7: num_subprocesses -= 1                            # 5-7 processors:  leave one unused 
        elif 8 <= cpu_count() <=15: num_subprocesses -= 1                           # 8-15 processors: leave one unused
        else: num_subprocesses += 1                                                 # > 15 processors: overload by 1 
        # above calculate more than number of processors we have assuming begin and end processes will go away quickly due to nan's
    cols_per_process = int(np.ceil(num_y_cols / num_subprocesses))                  # evenly divide y subscripts across all subprocesses and round up
    cprint("{} This box has {} processors.  Function iter_reflect spawning {} iter_y subprocesses with {:,} y cols per process and {:,} total y columns".format(timestamp(), cpu_count(), num_subprocesses, cols_per_process, num_y_cols),'green',attrs=['bold'])
    iter_parms = []
    for i in range(num_subprocesses):                                               # build parameter list for each subprocess
        y_start = int(i * cols_per_process)                                         # the only parameters that vary are y_start and y_end
        if i == num_subprocesses-1:                                                 # if this is last subprocess
            y_end = num_y_cols                                                      # then null y_end to subset remainder of dataset
        else:                                                                       # for all except the last subprocess
            y_end = y_start + cols_per_process                                      # set y_end to one below next subset
        parm_tuple = (i, h5_input_file_full, y_start, y_end, scaleFactor, noDataValue, mapInfo)
        iter_parms.append(parm_tuple)                                               # append parmlist for this subprocess
    with Pool(num_subprocesses) as p:                                               # instantiate num_subprocesses subprocesses
        results = p.starmap(iter_y, iter_parms)                                     # and spawn iter_y into each - will block until all processes complete
    cprint('{} Processed {:,} coordinates in {}. Now aggregating and sorting full result list...'.format(timestamp(),totrows, str(datetime.timedelta(seconds = time.time()-stime))),'green',attrs=['bold'])
    full_result_list = []
    for subprocess_results in results:                                              # iterate through each subprocess result list in results
        for result in subprocess_results:                                           # now iterate through each result from this subprocess
            full_result_list.append(result)                                         # and add from subprocess list
    full_result_list = np.array(full_result_list)                                   # convert to numpy for lexsort
    # Since reflectance data is in Y-desc, X-asc sequence, we re-sort to Y-asc, X-asc for readability and later analyze_overlap merging
    full_result_list = full_result_list[np.lexsort((full_result_list[:,1],full_result_list[:,0]))]
    cprint('{} Aggregation and sort of full result list complete.  Proceeding with hdf5 dataset creation to {}...'.format(timestamp(), args.h5_output_file),'green',attrs=['bold'])
    
    col_labels = [
        'y_coordinate',                                                             # the y absolute coordinate
        'x_coordinate',                                                             # the x absolute coordinate
        'result',                                                                   # the result
        ]
        
    if args.csv:                                                                    # if you want a csv
        import csv
        csv_filename = "{}.csv".format(filepref)
        with open(csv_filename,'w+',newline='') as csvfile:                         # create the new csvfile
            csvwriter = csv.writer(csvfile)                                         # and establish the csv writer object
            csvwriter.writerow(col_labels)                                          # then, output the header row
            for result in full_result_list:                                         # now iterate through each result from all subprocesses
                csvwriter.writerow(result)                                          # and write it to csv
        cprint("{} Created {}".format(timestamp(),csv_filename),'green',attrs=['bold'])

    def get_dir_name(dsn):
        return int(hashlib.sha1(dsn.encode('ascii','ignore')).hexdigest(), 16) % 10 # return mod10 remainder of crc of dsn
        
    with h5py.File(args.h5_output_file, 'a', libver='latest') as f:                 # open existing h5 results file - create if non-existant
        subdir = get_dir_name(filepref)                                             # subdir will be one character digit 0-9
        dsname = '/results/{}/{}'.format(subdir,filepref)                           # build the fq hdf5 dataset name
        if dsname in f:                                                             # if this dataset is already in the file...
            del f[dsname]                                                           # then it must be deleted before recreation
            cprint("{} Deleted existing dataset {} from hdf5 file {}".format(timestamp(), dsname, args.h5_output_file), 'green', attrs=['bold'])
        ds = f.create_dataset(dsname, data=full_result_list, compression="gzip", compression_opts=9) # create the result dataset from result_list
        ds.attrs['Creation_Date'] = time.strftime('%Y-%m-%dT%H:%M:%S+00:00', time.gmtime()).encode('ascii','ignore') # add creation date attribute
        ds.attrs['Array_Column_Labels'] = np.string_(col_labels)                    # add column labels attribute
        ds.attrs['Input_File'] = h5_input_file.encode('ascii','ignore')             # add input file attribute
        ds.attrs['Sort'] = 'y_coordinate ascending then x_coordinate ascending'.encode('ascii','ignore')              # add sort attribute
        ds.attrs['Created_By_Script']   = __file__.encode('ascii','ignore')         # add this script name attribute
        ds.attrs['Map_Info'] = "UTM Zone {}-{} with datum {} and {:.0f}/{:.0f} meter x/y resolution".format(mapInfo['utm_zone'], mapInfo['utm_ns'], mapInfo['datum'], mapInfo['x_pixel_size'], mapInfo['y_pixel_size']).encode('ascii','ignore')
        cprint("{} Added dataset {} with shape {} to hdf5 file {}".format(timestamp(), dsname, ds.shape, args.h5_output_file), 'green', attrs=['bold'])
    full_result_list = None                                                         # reclaim this memory space


#######################   begin mainline ############################
if __name__ == '__main__':
    colorama.init()
    if args.h5_input_file:
        h5_input_file_dir, h5_input_file_prefix = os.path.split(args.h5_input_file)
        if h5_input_file_prefix:
            h5_input_file_prefix = '^' + h5_input_file_prefix + '$'
        else:
            h5_input_file_prefix = '^' + h5_input_file_dir + '$'
            h5_input_file_dir = os.getcwd()
    elif args.h5_input_file_prefix:
        h5_input_file_dir, h5_input_file_prefix = os.path.split(args.h5_input_file_prefix)
        if h5_input_file_prefix:
            h5_input_file_prefix = '^' + h5_input_file_prefix + '.+?(?:h5|hdf5)$'
        else:
            h5_input_file_prefix = '^' + h5_input_file_dir + '.+?(?:h5|hdf5)$'
            h5_input_file_dir = os.getcwd()
    else:
        raise RuntimeError('neither --infile nor --infile_prefix found')
    h5_input_file_regex = re.compile(h5_input_file_prefix,flags=re.IGNORECASE)
    filesuff = '_reflectance.*?\\.(?:h5|hdf5)$'
    cprint('Selecting all files in directory "{}" matching pattern "{}"'.format(h5_input_file_dir, h5_input_file_prefix),'green',attrs=['bold'])
    for h5_input_file in next(os.walk(h5_input_file_dir))[2]:
        if not h5_input_file_regex.search(h5_input_file):
            continue    
        filepref = re.search('(.+){}'.format(filesuff),h5_input_file).group(1)
        ##### open/instantiate the h5 file and then print key/value detail if running in "-vv" mode
        h5_input_file_full = os.path.join(h5_input_file_dir, h5_input_file)
        cprint('Now processing input file "{}"'.format(h5_input_file_full),'green',attrs=['bold'])
        with h5py.File(h5_input_file_full, 'r') as f:
            ##### Obtain objects for Reflectance, Reflectance_Data datasets
            reflectance_ds = f[reflectance_ds_name]
            reflectance_data_ds = f[reflectance_data_ds_name]
            
            ##### Obtain scale_factor, data_ignore
            #scaleFactor = reflectance_data_ds.attrs['Scale_Factor']
            scaleFactor = 10000
            #noDataValue = reflectance_data_ds.attrs['Data_Ignore_Value']
            noDataValue = -10000


            ##### Extract Map_Info data
            mapInfo_ds = reflectance_ds['Metadata/Coordinate_System/Map_Info']
            # UTM,  1.000,  1.000,  505834.000,  3515363.000,  1.0000000000e+000,  1.0000000000e+000,  12,  North,  WGS-84,  units=Meters, 0
            mapInfo_split = mapInfo_ds.value.decode("utf-8").split(",")   # split the string (looks like above example) into substrings using the separator ","
            mapInfo = {
                'projection_name'       : mapInfo_split[0].strip(),             # usually 'UTM'
                'reference_pixel_x'     : float(mapInfo_split[1].strip()),      # not sure of use
                'reference_pixel_y'     : float(mapInfo_split[2].strip()),      # not sure of use
                'pixel_easting'         : float(mapInfo_split[3].strip()),      # xMin - westernmost boundary (array is x-ascending order)
                'pixel_northing'        : float(mapInfo_split[4].strip()),      # yMax - northernmost boundary (array is y-descending order)
                'x_pixel_size'          : float(mapInfo_split[5].strip()),      # x resolution
                'y_pixel_size'          : float(mapInfo_split[6].strip()),      # y resolution
                'utm_zone'              : int(mapInfo_split[7].strip()),        # for SRER will always be 12
                'utm_ns'                : mapInfo_split[8].strip(),             # 'North' or 'South'... for northern-hemisphere, always 'North'
                'datum'                 : mapInfo_split[9].strip(),             # usually 'WGS-84' - DOD GPS standard
                'units'                 : mapInfo_split[10].strip(),            # usually 'units=Meters'
                'rotation_angle'        : float(mapInfo_split[11].strip()),     # usually 0
                }
            ##################################################################################################
            ## IMPORTANT:  See doc block at beginning of this script detailing mapInfo and UTM usage
            ################################################################################################## 
            
            if mapInfo['x_pixel_size'] != 1 or mapInfo['y_pixel_size'] != 1:
                cprint('*** WARNING *** This script has not been tested at resolutions other than 1.  Output may not be valid.','yellow',attrs=['bold'])
                cprint('Resolution:{}/{}'.format(mapInfo['x_pixel_size'], mapInfo['y_pixel_size']),'yellow',attrs=['bold'])             

            num_y_cols, num_x_cols, num_z_cols = reflectance_data_ds.shape
            cprint('Reflectance Data Dimensions: y:{}, x:{} wavelengths:{} yielding {:,} discrete wavelength collections'.format(num_y_cols, num_x_cols, num_z_cols, (num_y_cols*num_x_cols)),'green',attrs=['bold'])
            
            if args.verbose:
                if args.verbose > 1:
                    def print_attrs(name, obj):
                        cprint("    Group: {}".format(name),'cyan',attrs=['bold'])
                        for key, val in obj.attrs.items():
                            cprint("        {}: {}".format(key, val),'cyan',attrs=['bold'])
                    cprint("\n**** H5 file {} Group Item Attribute Data ****".format(h5_input_file_full),'cyan',attrs=['bold'])
                    f.visititems(print_attrs)
                    cprint("**** End Group Item Attribute Data ****\n",'cyan',attrs=['bold'])
                    
                cprint('Resolution: {}/{} (x_pixel_size/y_pixel_size)'.format(mapInfo['x_pixel_size'], mapInfo['y_pixel_size']),'cyan',attrs=['bold'])    
                wavelength_ds = reflectance_ds['Metadata/Spectral_Data/Wavelength']
                cprint('min wavelength:{}nm  bandwidth:{}nm'.format(np.amin(wavelength_ds),(wavelength_ds.value[1]-wavelength_ds.value[0])),'cyan',attrs=['bold'])
                cprint('max wavelength:{}nm  bandwidth:{}nm'.format(np.amax(wavelength_ds),(wavelength_ds.value[-1]-wavelength_ds.value[-2])),'cyan',attrs=['bold'])
                cprint('Map Info: {}'.format(mapInfo),'cyan',attrs=['bold'])
                dataset_coords_dict = {
                    'xMin': mapInfo['pixel_easting'],   # from above example, this will be 505834.000 which is the Western-most x position which is the x minimum
                    'xMax': mapInfo['pixel_easting']  + (num_x_cols * mapInfo['x_pixel_size']),     # xMin + (# of columns * resolution) Eastern-most extent
                    'yMin': mapInfo['pixel_northing'] - (num_y_cols * mapInfo['y_pixel_size']),     # yMax - (# of columns * resolution) Southern-most extent
                    'yMax': mapInfo['pixel_northing'],  # from above example, this will be 3515363.000 which is the Northern-most y position which is the y maximum
                    }
                cprint('Global Extent Coordinates of all data in dataset:{}'.format(dataset_coords_dict),'cyan',attrs=['bold'])
                cprint('Extent Definitions: xMin=West, xMax=East, yMin=South, yMax=North','cyan',attrs=['bold'])
                cprint('Scale Factor:{}'.format(scaleFactor),'cyan',attrs=['bold'])
                cprint('Data Ignore Value:{}'.format(noDataValue),'cyan',attrs=['bold'])
                sq_km = ((dataset_coords_dict['xMax'] - dataset_coords_dict['xMin']) * (dataset_coords_dict['yMax'] - dataset_coords_dict['yMin'])) / 1000**2
                spatial_extent = reflectance_data_ds.attrs['Spatial_Extent_meters']
                cprint('Spatial Extent (meters):{} ({} sq-km)'.format(spatial_extent, sq_km),'cyan',attrs=['bold'])

            if args.x or args.y:
                if args.x and args.y:
                    calc_val(reflectance_data_ds[args.y,args.x], args.y, args.x, scaleFactor, noDataValue, diag=True)
                else:
                    cprint('Must specify both x and y value','red',attrs=['bold'])

        if not (args.x or args.y):
            iter_reflect(h5_input_file_full, filepref, scaleFactor, noDataValue, mapInfo)