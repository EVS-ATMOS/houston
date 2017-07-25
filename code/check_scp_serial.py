from netCDF4 import Dataset
from dask import delayed, compute
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
n_workers = int(sys.argv[1])
ref_thresholds = np.arange(0, 40, 5)
n_levels = 31
fields_not_used = ['velocity', 'cross correlation ratio', 'spectrum_width',
                   'specific_differential_phase',
                   'unfolded_differential_phase', 'differential_reflectivity',
                   'cylp_processed_phase']


def calculate_scp(file):
    print(file)
    grid = Dataset(file)
    # Get reflectivity and height levels
    try:
        refl = grid.variables['reflectivity'][:]
    except:
        SCP_array = np.nan*np.ones((len(ref_thresholds), n_levels))
        grid_time = np.nan
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
ignore_file = open('ignore_file', 'r')
for files in ignore_file:
    try:
        print(grid_list[0])
        grid_list.remove(files[:-1])
        print('Removed ' + files[:-1])
    except ValueError:
        print(("Not removing file " + files[:-1]))
print('## Done! Time = ' + str(time.time() - bt) + ' s ##')
print('## About to process ' + str(len(grid_list)) + ' grids')

# Map workers to file
chunks = np.ceil(len(grid_list)/(n_workers*8))
chunk_size = int(len(grid_list)/chunks)
print('Separating into ' + str(chunk_size) + ' size chunks')
SCPs = []
for i in range(0, len(grid_list), chunk_size):
    time_chunk = time.time()
    SCP_chunk = [calculate_scp(file) for file in grid_list[i:i+chunk_size]]
    print(str(i) + 'th chunk done in ' + str(time.time() - time_chunk) + ' s')
