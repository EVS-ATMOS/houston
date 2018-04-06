from matplotlib import use
use('agg')
import pyart
import os
import sys
from glob import glob
import platform
import copy
import netCDF4
import datetime
import numpy as np
import time
import gc
import dask.bag as db
import xarray
import matplotlib.dates as mdates

from matplotlib import pyplot as plt
from distributed import Client, wait
from scipy import interpolate
from netCDF4 import num2date

def retrieve_qvp(radar, hts, flds = None):
    if flds == None:
        flds = ['differential_phase',
            'cross_correlation_ratio',
            'spectrum_width',
            'reflectivity', 'differential_reflectivity']
    desired_angle = 20.0
    index = abs(radar.fixed_angle['data'] - desired_angle).argmin()
    ss = radar.sweep_start_ray_index['data'][index]
    se = radar.sweep_end_ray_index['data'][index]
    mid = int((ss + se)/2)
    z = radar.gate_altitude['data'][mid, :]
    texture = pyart.retrieve.calculate_velocity_texture(radar)    
    qvp = {}
    for fld in flds:
        this_fld = radar.get_field(index, fld)
        shapes = this_fld.shape
        this_fld = np.ma.masked_where(texture['data'][:shapes[0], :shapes[1]] > 3, this_fld)
        this_fld = np.ma.mean(this_fld, axis=0) 
        intery = interpolate.interp1d(z,
                this_fld, bounds_error=False, fill_value='extrapolate')
        ithis = intery(hts)
        qvp.update({fld:ithis})
    qvp.update({'time': radar.time})
    qvp.update({'height': hts})
    return qvp


def get_file_tree(start_dir, pattern):
    """
    Make a list of all files matching pattern
    above start_dir
    Parameters
    ----------
    start_dir : string
        base_directory
    pattern : string
        pattern to match. Use * for wildcard
    Returns
    -------
    files : list
        list of strings
    """

    files = []

    for dir, _, _ in os.walk(start_dir):
        files.extend(glob(os.path.join(dir, pattern)))
    return files

                
def get_qvp_from_file(file_name, hts, flds = None):
    try:
       radar = pyart.io.read(file_name)
       qvp = retrieve_qvp(radar, hts, flds)
       del radar
       return qvp
    except KeyError:
       return []

if __name__=='__main__':
    #hello_world()
    top = ['/lcrc/group/earthscience/radar/houston/data/' + sys.argv[1] + '/']

    odir_data = '/lcrc/group/earthscience/rjackson/houston_radar_set/grid_files/'
    odir_images = '/lcrc/group/earthscience/radar/rjackson/houston_radar_set/plots/'

    if len(sys.argv) > 1:
        patt = sys.argv[1]
    else:
        patt = None

    all_files = []
    for paths in top:
        all_files = all_files + glob(paths + '*sur*.nc')
    print(all_files[1])
    print('found ', len(all_files), ' files')
    good_files = []
    for ffile in all_files:
        statinfo = os.stat(ffile)
        if statinfo.st_size > 1*1e6:
            good_files.append(ffile)

    if patt is None:
        really_good = good_files
    else:
        really_good = []
        for ffile in good_files:
            if patt in ffile:
                really_good.append(ffile)

    print(len(good_files))
    print(len(really_good))
    
    packing = []
    file_list = good_files
    
    client = Client(scheduler_file='/home/rjackson/scheduler.json')
    start_time = time.time()
    the_bag = db.from_sequence(good_files)
    n_cores = sum(client.ncores().values())
    heights = np.arange(0.0, 20000.0, 250.0)
    print("Dask cluster has " + str(n_cores) + " cores")
    get_qvp = lambda file_name: get_qvp_from_file(file_name, heights)
    qvps = the_bag.map(get_qvp).compute()
    #create an empty dictionary to be filled
    time_series = {}
    #work out how many profiles
    try:
        qvps.remove([])
    except:
        print('Yay!')
    print(qvps)
    n_times = len(qvps)
    print(qvps[0].keys())
   
    zz =  np.zeros([n_times, len(heights)])

    #Loop over all fields 
    for fld in ['differential_phase', 
	        'cross_correlation_ratio',
		'spectrum_width',
		'reflectivity', 'differential_reflectivity']:
        #Create empty time series
        this_fld = np.zeros([n_times, len(qvps[0][fld])])
    
        #create time array
        times = np.zeros(n_times)

        #Loop over all times
        for i in range(n_times):

	    #Create a datetime object
            dateobj = num2date(qvps[i]['time']['data'][0], qvps[i]['time']['units'])

            #Append a numerical date using Matplotlib's time reference
            times[i] = mdates.date2num(dateobj)

	    #Update the profile to the time series
            this_fld[i,:] = qvps[i][fld]
    
            #append the time series to the dictionary
        time_series.update({fld: (['time', 'z'], this_fld[times.argsort(), :])})
    
    for i in range(n_times):
	#Update the profile to the time series
        zz[i,:] = qvps[i]['height']

	#append the time series to the dictionary
        zz = zz[times.argsort(), :]
    times.sort()
        # Write QVP to netCDF file using xarray
    
    print("Time to complete: " + str(time.time() - start_time))
    ds = xarray.Dataset(time_series, coords={'time': times, 'z': heights,
                                             'reference_time': min(times)},
                        attrs={'units': ['deg', ' ', 'm/s', 'dBZ', 'dB']})
    ds = ds.sortby('time', ascending=True)
    ds.to_netcdf(path=('QVP_Houston' + sys.argv[1] + '.cdf'), mode='w')
    print('## Processing done in ' + str(time.time() - bt) + ' s')
    Client.shutdown(timeout=10)





