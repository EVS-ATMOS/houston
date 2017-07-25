from netCDF4 import Dataset
from dask import delayed, compute
from dask.diagnostics import Profiler, ResourceProfiler, visualize
import dask.multiprocessing
from distributed import Client, LocalCluster
import numpy as np
from glob import glob
import time
import sys
from datetime import datetime, timedelta
import xarray
import gc

# Set path to grid files. Use the highest level directory as
# Glob can search for files recursively
houston_path = '/lcrc/group/earthscience/radar/houston/data/'

ref_thresholds = np.arange(0, 40, 5)
n_levels = 31
fields_not_used = ['velocity', 'cross correlation ratio', 'spectrum_width',
                   'specific_differential_phase',
                   'unfolded_differential_phase', 'differential_reflectivity',
                   'cylp_processed_phase']


def calculate_scp(file):
    grid = Dataset(file)
    # Get reflectivity and height levels
    try:
        refl = grid.variables['reflectivity'][:]
    except:
        SCP_array = np.nan*np.ones((len(ref_thresholds), n_levels))
        grid_time = datetime(1900, 1, 1, 1, 1, 1)
        z_levels = np.nan*np.ones(n_levels)
        return SCP_array, grid_time, z_levels

    z_levels = grid.variables['z'][:]
    refl = refl[0]
    # Convert grid time to datetime object
    time_units = grid.variables['time'].units
    time_data = grid.variables['time'][:]

    del grid
    grid_time = datetime.strptime(time_units,
                                  "seconds since %Y-%m-%dT%H:%M:%SZ")
    grid_time += timedelta(seconds=time_data[0])
    grid_shape = refl.shape
    print('...Loading time ' + str(grid_time))

    # Initalize Statistical Coverage Product array
    SCP_array = np.zeros((len(ref_thresholds), grid_shape[0]))
    total_points = grid_shape[1]*grid_shape[2]

    for i in range(0, len(ref_thresholds)):
        for j in range(0, len(z_levels)):
            masked_refl = np.where(refl[j] > ref_thresholds[i])
            count = len(masked_refl[0])
            SCP_array[i, j] = float(count)/float(total_points)*100

    del refl
    return SCP_array, grid_time, z_levels

# Get list of grid files
print('## Getting list of grid files... ##')
bt = time.time()
grid_list = glob((houston_path + sys.argv[1] + '**/*grid*.nc'),
                 recursive=True, )
ignore_list = open('ignore_file', 'r')
for files in ignore_list:
    try:
        grid_list.remove(files[:-1])
    except ValueError:
        pass

# Corrupt file
print('## Done! Time = ' + str(time.time() - bt) + ' s ##')
print('## Starting dask cluster ##')
print('## About to process ' + str(len(grid_list)) + ' grids')

# Start dask cluster
client = Client('127.0.0.1:8786')

# Map workers to file
scp = delayed(calculate_scp)
chunks = np.ceil(len(grid_list)/(256))
chunk_size = int(len(grid_list)/chunks)
print('Separating into ' + str(chunk_size) + ' size chunks')
SCPs = []
for i in range(0, len(grid_list), chunk_size):
    time_chunk = time.time()
    SCP_chunk = [scp(file) for file in grid_list[i:i+chunk_size]]
    SCPs.append(compute(*SCP_chunk, get=client.get))
    print(str(i) + 'th chunk done in ' + str(time.time() - time_chunk) + ' s')


SCPs = np.concatenate(SCPs)
SCP_array = [a for (a, b, c) in SCPs]
times = [b for (a, b, c) in SCPs]
z_levels = [c for (a, b, c) in SCPs]
z_levels = z_levels[0]
SCP_array = np.stack(SCP_array, axis=0)

# Write to xarray file
ds = xarray.Dataset({'SCP': (['time', 'reflectivity', 'z'], SCP_array)},
                    coords={'time': times, 'z': z_levels,
                            'reference_time': min(times),
                            'reflectivity': ref_thresholds},
                    attrs={'units': '%', 'long_name': ('May and Lane (2009)' +
                                                       ' Statistical' +
                                                       ' Coverage Product')})
ds = ds.sortby('time', ascending=True)
ds.to_netcdf(path=('SCP_Houston' + sys.argv[1] + '.cdf'), mode='w')
print('## Processing done in ' + str(time.time() - bt) + ' s')
Client.shutdown(timeout=10)
