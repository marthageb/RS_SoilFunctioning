#!/usr/bin/python3
from pprint import pprint, pformat
import time
import pdb
import numpy as np
import csv
import argparse
import re
import os
import sys
from datetime import datetime, timedelta
import copy
import warnings
from termcolor import colored, cprint
import colorama
import gc
import psutil

def timestamp():
    return time.strftime('%H:%M:%S', time.localtime(time.time()))
    
def mprint(text, color='green'):
    cprint('{} {}'.format(timestamp(),text),color,attrs=['bold'])
    
def convert_bytes(num):
    step_unit = 1000.0 #1024 bad the size
    for x in ['bytes', 'KB', 'MB', 'GB', 'TB']:
        if num < step_unit:
            return "%3.1f %s" % (num, x)
        num /= step_unit
        
def dict_key(begend,dsn):
    return '{:08.0f}_{}'.format(begend,dsn)    
    
def add_dict(dsn, ds, dspath, type, starty, endy, startx, endx, numrows):
    mdict = {                                           # this structure will be used in both ydict and xdict
                'type'      : type,                     # either 'start' or 'end'
                'dsn'       : dsn,                      # the original dsn that this came from
                'ds'        : ds,                       # the h5py object or internal array object representing this numpy dataset
                'dspath'    : dspath,
                'starty'    : starty,                   # the starting 'y' value in this ds
                'endy'      : endy,                     # the ending 'y' value in this ds
                'yrowcount' : int(endy)-int(starty)+1,  # the range of endy-starty+1
                'startx'    : startx,                   # the starting 'x' value in this ds
                'endx'      : endx,                     # the ending 'x' value in this ds
                'xrowcount' : int(endx)-int(startx)+1,  # the range of endx-startx+1
                'numrows'   : numrows,                  # the number of rows in this dataset
                'index'     : None,                     # used within xdict as the placeholder for where we are at in the iteration
                }
    if type == 'end':
        ydict.setdefault(dict_key(endy,dsn),mdict)
    else:
        ydict.setdefault(dict_key(starty,dsn),mdict)
    xdict.setdefault(dsn,mdict)
    
def dump_mosaic_to_dict():
    global dump_mosaic_iter, f_mosaic_tmp, mosaic
    opath = "{:06.0f}".format(dump_mosaic_iter)                                         # build the output hdf5 dataset path
    dump_mosaic_iter += 1                                                               # bump the iter for next time through
    f_mosaic_tmp.create_dataset(opath, data=mosaic)                                     # dump mosaic to hdf5 dataset - no compression
    if args.verbose >1: mprint('dumped mosaic with {:,.0f} rows to {}'.format(len(mosaic), opath),'cyan')
    mosaic = []                                                                         # reset mosaic to empty     

def fmt_dsns(xdict,ylined):
    dstxt = ''
    for ds in ylined['dsns']:     
        if dstxt:
            dstxt += ','
        numx = xdict[ds]['endx']-xdict[ds]['startx']+1
        numx_rows = numx * (ylined['end']-ylined['start']+1)
        dstxt += "\n         '{}' x:{:.0f}-{:.0f} (numx:{:,.0f},numxrows:{:,.0f})".format(ds, xdict[ds]['startx'], xdict[ds]['endx'], numx, numx_rows )
    return dstxt

script_start_time = datetime.now()
colorama.init()
##### Parse arguments
parser = argparse.ArgumentParser(description='Reflectance Results Overlap Analyzer', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--infile', dest='h5_input_file', default='results.h5', type=str, help='The input h5 filename')
parser.add_argument('--outfile', dest='h5_output_file', default='mosaic.h5', type=str, help='The output h5 filename')
parser.add_argument('--verbose', '-v', action='count', default=0, help='Verbosity of output "-v" or "-vv" for even more')

arg_group = parser.add_mutually_exclusive_group()
arg_group.add_argument('--use_h5py_mem', dest='use_option', action='store_const', const='h5py_mem', default='h5py_mem', help='Uses h5py hdf5 file access with memory preload... Close to deepdish performance but significantly more memory-conservative.')
arg_group.add_argument('--use_deepdish', dest='use_option', action='store_const', const='deepdish', help='Uses deepdish hdf5 file access.  Memory intensive but marginally faster than h5py_mem.')
arg_group.add_argument('--use_h5py_file', dest='use_option',action='store_const', const='h5py_file', help='Uses h5py direct hdf5 file access.  Very memory-conservative but > 4x slower than h5py_mem.')

parser.add_argument('--stats_interval', default=5, type=int, help='Statistics interval in minutes')
parser.add_argument('--skipiters', action='store_true', help='Only perform first iteration - testing only')
parser.add_argument('--debug', action='store_true', help='Perform debug sections - testing only')
args = parser.parse_args()
mprint('Script {} starts at {} using input file "{}" and hdf5 option {}'.format( __file__, datetime.now(), args.h5_input_file, args.use_option))

if args.use_option == 'deepdish':                                                   # if using deepdish
    statinfo = os.stat(args.h5_input_file)                                          # obtain hdf5 file stats
    size_limit = .5 * 10**9                                                         # .5G - my arbitrary size limit for deepdish
    if statinfo.st_size > size_limit:                                               # does the file size exceed this limit
        mprint('The file size of {} is {:,.0f}, exceeding {}, which may lead to memory problems with deepdish.  Switching to "--use_h5py_mem" mode'.format(args.h5_input_file, statinfo.st_size, convert_bytes(size_limit)),'yellow')
        args.use_option = 'h5py_mem'                                                # yes - change to h5py_mem mode
    else:
        mprint('Using deepdish for hdf5 file access')
        import deepdish as dd
        f = dd.io.load(args.h5_input_file)                                          # deepdish loads everything (all datasets) from this h5 into the "f" dict
        
if args.use_option in ['h5py_mem', 'h5py_file']:                                    # if using h5py
    if args.use_option == 'h5py_mem':
        mprint('Using h5py "copy array to memory" mode for hdf5 file access')
    else:
        mprint('Using h5py "direct file access" mode for hdf5 file access')
    import h5py
    f = h5py.File(args.h5_input_file, 'r')                                          # h5py instantiates an instance in the "f" object - does NOT load datasets
    mosaic_tmp_filename = '{}.tmp'.format(args.h5_output_file)
    f_mosaic_tmp = h5py.File(mosaic_tmp_filename, 'w')                              # open tmp file (with truncation) to hold intermediate mosaic dumps

ydict = {}
xdict = {}
ary = {}                                                  
totrows = 0
yline = []
inprog = {}

#### This section builds the xdict and ydict dictionaries containing all of the result datasets demographics
for dir in sorted(f['results']):
    for dsn in sorted(f['results'][dir]):
        ds = f['results'][dir][dsn]
        dspath = 'results/{}/{}'.format(dir,dsn)
        if args.use_option == 'h5py_mem':
            useds = None
        else:
            useds = ds
        add_dict(dsn, useds, dspath, 'start', ds[0,0], ds[-1,0], ds[0,1], ds[-1,1], ds.shape[0])
        add_dict(dsn, useds, dspath, 'end',   ds[0,0], ds[-1,0], ds[0,1], ds[-1,1], ds.shape[0])
        totrows += xdict[dsn]['numrows']

#### This section builds the yline list (from the above ydict dictionary) which is extensively used to determine when to start or stop using a particular dataset        
if args.verbose > 1: cprint("\n*** ydict values ***\n",'cyan',attrs=['bold'])
for ykey in sorted(ydict):
                                             # this will give us ykey sorted in y-descending order
    if args.verbose > 1: cprint('{:5s} y:{:.0f}-{:.0f} x:{:.0f}-{:.0f} {}'.format(ydict[ykey]['type'],ydict[ykey]['starty'],ydict[ykey]['endy'],ydict[ykey]['startx'],ydict[ykey]['endx'],ydict[ykey]['dsn']),'cyan',attrs=['bold'])

    if ydict[ykey]['type'] == 'start':                                              # if a new dataset is starting
        inprog.setdefault(ykey,ydict[ykey])                                         # add it to the in_progress dict
    else:                                                                           # a dataset is ending
        del inprog[dict_key(ydict[ykey]['starty'],ydict[ykey]['dsn'])]              # remove it from in_progress dict
    ylined = {                                                                      # init basic yline dict
        'start'     : None,
        'end'       : None,
        'xrowcount' : 0,
        'dsns'      : [],
        }
    if ydict[ykey]['type'] == 'start':                                              # if dataset is starting
        ylined['start'] = ydict[ykey]['starty']                                     # then set ylined start value
        try:
            yline[-1]['end'] = ydict[ykey]['starty'] - 1                            # if there is a prior ylined then set its end to one less than this starty
        except:
            pass                                                                    # otherwise ignore
    else:                                                                           # dataset is ending
        ylined['start'] = ydict[ykey]['endy'] + 1                                   # set ylined start to one more than this endy
        try:
            yline[-1]['end'] = ydict[ykey]['endy']                                  # if there is a prior ylined then set its end to this endy
        except:
            pass                                                                    # otherwise ignore
    # iterate through inprog to add all dsn's in progress
    for inprog_key in sorted(inprog):                                               # iterate through all dsn's thet are in_progress
        dsn = inprog[inprog_key]['dsn']
        ylined['xrowcount'] += xdict[dsn]['xrowcount']                              # increment ylined xrowcount by the xdict xrowcount
        ylined['dsns'].append(dsn)                                                  # append dsn to ylined dsns
    yline.append(ylined)                                                            # append ylined to yline
del yline[-1]                                                                       # trash last one
yd_rowcount = 0                                                                     # for following error check
for ylined in yline:                                                                # iterate over all yline's
    ylined['rowcount'] = ylined['xrowcount'] * ((int(ylined['end']) - int(ylined['start'])) + 1)    # sed ylined rowcount to real rowcount (x*y)
    yd_rowcount += ylined['rowcount']                                               # for following error check
assert yd_rowcount == totrows,"ydict rowcount {:,d} does not equal rowcount {:,d} of all datasets".format(yd_rowcount,totrows)

if args.verbose:
    cprint("\n*** yline values ***\n",'cyan',attrs=['bold'])
    for ylined in yline:
        dstxt = fmt_dsns(xdict,ylined)
        cprint('y:{:.0f}-{:.0f} (numy={:5.0f}, numrows={:,d}) [{}]'.format(ylined['start'],ylined['end'] ,ylined['end']-ylined['start']+1, ylined['rowcount'], dstxt),'cyan',attrs=['bold'])
    print("")

#### This section does all the work.  Using the yline list from above, start iterating through each dataset    
mprint('Processing {:,} rows'.format(totrows))
mosaic = []                                                                         # this is the final output list containing combined datasets                                                                       
numylined = 0
stime = time.time()                                                                 # initial iteration time used for status calculations below
last_status_interval_processed = 0
status_interval_minutes = args.stats_interval                                       # want status message every this-many minutes
status_interval_seconds = status_interval_minutes * 60                              # want status message every this-many seconds
processed_rows = 0                                                                  # incremented for each row accessed in each input ds  
last_interval_processed_rows = 0                                                    # for stats                          
expected_cum_rows = 0                                                               # post-incremented in each yline iteration 
expected_rows = 0                                                                   # used to perform above post-increment
found_nan = False
dump_mosaic_interval = 10**7                                                        # dump the mosaic dict whenever it gets this big - 10**7 is about 1-1.5G memory
dump_mosaic_iter = 0

for ylined in yline:                                                                # iterate over each yline overlap
    numylined += 1                                                                  # bump for stats message
    msgstr = 'Processing Y subset {} : {:.0f}-{:.0f} ({:.0f}) with {:,} rows. There have been {:,} rows ({:.2f}%) processed thus far.'.format(numylined, ylined['start'], ylined['end'], ylined['end']-ylined['start']+1, ylined['rowcount'], processed_rows, processed_rows/totrows * 100)
    if args.verbose:
        mprint('{} {}'.format(msgstr, fmt_dsns(xdict,ylined)))
    else:
        mprint(msgstr)
    expected_cum_rows += expected_rows                                              # expected_frows from prior iteration
    expected_rows = ylined['rowcount']                                              # for next time around
    if expected_cum_rows != processed_rows:                                         # if there is a mismatch between expected:actual
        mprint('WARNING: Expected cumulative rows {} and processed rows {} are not equal'.format(expected_cum_rows, processed_rows),'red')
    lowestx = 99999999                                                              # set high to find lower
    highestx = 0                                                                    # set low to find higher
    
    if args.use_option == 'h5py_mem':                                               # if we are using h5py mem then need to delete any unused ary members
        for dsn in list(ary.keys()):                                                # must convert iterator to list to allow del in iteration
            if dsn not in ylined['dsns']:                                           # if this dataset is no longer needed
                del ary[dsn]                                                        # then remove it from dict
                mprint('Deleted ary[{}] from memory.  There are now {} datasets loaded'.format(dsn, len(ary)))
                xdict[dsn]['ds'] = None                                             # and set xdict that it is gone also
        # note... at this point, some ary members may have been deleted but the memory usage MAY not decrease as python sometimes retains the memory for future use
        # that memory will be reused for future allocations, including instantiation of new ary members
        gc.collect()                                                                # perform garbage collection anyhow... it won't hurt and may help
            
    for dsn in ylined['dsns']:                                                      # iterate all dsns in this ylined
        if args.use_option == 'h5py_mem' and dsn not in ary:                        # if this dataset is not in memory then...
            ary[dsn] = f[xdict[dsn]['dspath']][:]                                   # read entire dataset into memory
            mprint('Added ary[{}] to memory.  There are now {} datasets loaded'.format(dsn, len(ary)))
            xdict[dsn]['ds'] = ary[dsn]                                             # and set dataset object in xdict
        if not xdict[dsn]['index']:                                                 # if we are not already processing this ylined 
            xdict[dsn]['index'] = np.searchsorted(xdict[dsn]['ds'][:,0], ylined['start'])   # in this ds, the index of the first [0] that matches this y start 
        if int(xdict[dsn]['endx']) > highestx:                                      # of all ds's is this the highest x?
            highestx = int(xdict[dsn]['endx'])                                      # yes, set new highest
        if int(xdict[dsn]['startx']) < lowestx:                                     # of all ds's is this the lowest x?
            lowestx = int(xdict[dsn]['startx'])                                     # yes, set new lowest
            
    break_new_ylined = False                                                        # initialize no-break
    while not break_new_ylined:                                                     # while we don't need to break to new ylined, process new y
        if (args.use_option in ['h5py_file', 'h5py_mem']) and len(mosaic) >= dump_mosaic_interval:   # if using h5py and mosaic is longer than interval
            dump_mosaic_to_dict()                                                   # then dump/clear mosaic
        for x in range(lowestx,highestx+1):                                         # process all x for this y
            mean_list = []                                                          # initialize it for this x
            y_exceeded = 0
            for dsn in ylined['dsns']:                                              # for this x, y then process all dsns in this ylined
                xd = xdict[dsn]                                                     # for easier reference below
                if xd['index'] >= xd['numrows']:                                    # if reached the end of the array then...
                    y_exceeded += 1                                                 # signal that y is exceeded for this dataset
                    continue                                                        # and iter to next dataset
                curr_y, curr_x, curr_val = xdict[dsn]['ds'][xd['index']]            # assign scalar names to list members
                if curr_y > ylined['end']:                                          # if we are at the end of this ylined subset then...
                    y_exceeded += 1                                                 # signal that y is exceeded for this dataset
                    continue                                                        # and iter to next dataset
                if int(curr_x) != x:                                                # found non-consecutive x in input dataset
                    if xd['startx'] <= x <= xd['endx']:                             # if there SHOULD HAVE been an x:
                        raise Exception('RuntimeError','xdict[{}]["ds"][{}]={} prev mo:{} expecting x={} but detected non-consecutive x'.format(dsn, xd['index'], xdict[dsn]['ds'][xd['index']], mo, x))
                    else:
                        continue                                                    # if we've stepped into new y then iterate next dsn
                if np.isnan(curr_val):                                              # if curr_val is NaN
                    found_nan = True                                                # then signal at least one NaN was found
                else:                                                               # otherwise, curr_val should be a float
                    mean_list.append(curr_val)                                      # append this value to overlap mean_list
                gath_y = curr_y                                                     # set gath_y to the curr_y that this x was collected at
                xd['index'] += 1                                                    # bump index for next iteration for this ds
                processed_rows += 1                                                 # bump for statistics
                # iterate next dsn in ylined['dsns'] for same x
            if y_exceeded >= len(ylined['dsns']):                                   # if past ylined['end'] on all datasets
                break_new_ylined = True                                             # signal break to outer loops
                break                                                               # and break out of for x: iter
            # now mean_list contains all values for this y:x in all overlapping datasets
            if len(mean_list):                                                      # if we are in-x-range on at least one dataset
                mo = [gath_y, float(x), sum(mean_list) / float(len(mean_list))]     # set mean of mean_list
                mosaic.append(mo)                                                   # and add mosaic y/x averaged value
            elif found_nan:                                                         # if there was no mean_list and at least one NaN was found
                mo = [gath_y, float(x), np.NaN]                                     # mo is mosaic list to be appended... use NaN as value
                mosaic.append(mo)                                                   # and add y/x nan value
            found_nan = False                                                       # reset flag
            interval = time.time() - stime                                          # seconds since this process started
            if not(int(interval % status_interval_seconds)):                        # do we need stats on this interval?
                this_interval = int(interval // status_interval_seconds)            # calculate this interval
                if last_status_interval_processed != this_interval:                 # have we already presented message for this interval?
                    last_status_interval_processed = this_interval                  # no then reset this_interval and present the status message
                    interval_processed_rows = processed_rows - last_interval_processed_rows # gives us rows processed in this interval
                    last_interval_processed_rows = processed_rows                   # for next interval calculation
                    lper = processed_rows/totrows * 100                             # percentage of all total y rows that I have to process
                    rrate = interval_processed_rows / status_interval_seconds       # run-rate in rows per second
                    remain_secs = (totrows - processed_rows) / rrate                # using past history project the remaining seconds
                    etc_text = time.strftime('%H:%M:%S', time.localtime(time.time() + remain_secs)) # get formatted wall-clock time of estimated completion
                    remain_secs_fmt = str(timedelta(seconds = remain_secs))         # formatted hh:mm:ss of remaining seconds
                    mprint("Processing Y:{:.0f} rrate:{:.0f} ds rows/sec {:,.0f} rows this interval ({:,.0f}/{:,.0f} or {:.2f}%) expected complete in {} at {} memory:{}".format(curr_y, rrate, interval_processed_rows, processed_rows, totrows, lper, remain_secs_fmt, etc_text, convert_bytes(psutil.Process().memory_info().vms)),'cyan')
            # iterate "for x in range"
        # iterate "while not break_new_ylined:"
        if args.skipiters: break
    # iterate "for ylined in yline:"
    if args.skipiters: break
# mosaic loading is now complete... take any follow-up action required before writing it to --outfile
f.close()                                                                           # close results input file

ods = "mosaic"                                                                      # output dsname    
if args.use_option == 'deepdish':                                                   # deepdish?
    mosaic = np.array(mosaic)                                                       # convert list to numpy
    mosaic_rows = np.shape(mosaic)[0]                                               # obtain # rows in array
    mprint('Creation of mosaic array consisting of the averages of {:,d} input rows is complete.  Now writing dataset "{}" with {:,.0f} rows to hdf5 file "{}" (may take several minutes for large rowcount)...'.format(processed_rows, ods, mosaic_rows, args.h5_output_file))
    dd.io.save(args.h5_output_file, {ods: mosaic})                                  # and dump mosaic numpy to the hdf5 dataset
else:                                                                               # using h5py
    mprint('Creation of mosaic array to tmp datasets consisting of the averages of {:,d} input rows is complete.  Now rebuilding in-memory mosaic from tmp datasets...'.format(processed_rows))
    dump_mosaic_to_dict()                                                           # dump and clear mosaic
    mosaic = np.empty([0,3])                                                        # create 0 x 3 empty numpy array                     
    for tmp in sorted(f_mosaic_tmp):                                                # iter over all sorted dsnames in mosaic_tmp
        mosaic = np.append(mosaic, f_mosaic_tmp[tmp], axis=0)                       # append tmp to mosaic
        if args.verbose >1:
            mrows = np.shape(f_mosaic_tmp[tmp])[0]                                  # rows in tmp
            trows = np.shape(mosaic)[0]                                             # rows in mosaic
            mprint('appended {} to mosaic with {:,.0f} rows which now has {:,.0f} rows'.format(tmp, mrows, trows),'cyan')
    f_mosaic_tmp.close()                                                            # close mosaic tmp file
    try:
        f_mosaic_tmp.close()                                                        # close it again (for some reason, twice is required)
    except:
        pass
    os.remove(mosaic_tmp_filename)                                                  # remove the tmp file
    mosaic_rows = np.shape(mosaic)[0]                                               # obtain # rows in array
    mprint('mosaic numpy recreation complete... now writing dataset "{}" with {:,.0f} rows to file "{}" (may take several minutes for large rowcount)...'.format(ods, mosaic_rows, args.h5_output_file))
    with h5py.File(args.h5_output_file, 'w') as f_mosaic:                           # open the output file (with truncation)
        # note... using szip compression on create_dataset because, for 343,146,362 x 3, gzip compression took > 45 minutes and szip took 2 minutes and compressed slightly better
        ds = f_mosaic.create_dataset(ods, data=mosaic, compression="szip")          # dump mosaic numpy to the hdf5 dataset
        ds.attrs['Creation_Date_UTC']   = time.strftime('%Y-%m-%dT%H:%M:%S+00:00', time.gmtime()).encode('ascii','ignore')  # add creation date attribute
        col_labels = [
            'y_coordinate',                                                                                                 # the y absolute coordinate
            'x_coordinate',                                                                                                 # the x absolute coordinate
            'result',                                                                                                       # the result
            ]
        ds.attrs['Array_Column_Labels'] = np.string_(col_labels)                                                            # add column labels attribute
        ds.attrs['Input_File']          = args.h5_input_file.encode('ascii','ignore')                                       # add input file attribute
        ds.attrs['Sort']                = 'y_coordinate ascending then x_coordinate ascending'.encode('ascii','ignore')     # add sort sequence attribute       
        ds.attrs['Created_By_Script']   = __file__.encode('ascii','ignore')                                                 # add this script name attribute

mprint('Created dataset "{}" in file "{}" with {:,d} rows'.format(ods, args.h5_output_file, mosaic_rows))
mprint('{:.1f}% of the input rows were eliminated/averaged due to overlap'.format((processed_rows - mosaic_rows) / processed_rows * 100))
mprint('Script {} ends at {} in {}'.format(__file__, datetime.now(), datetime.now() - script_start_time))