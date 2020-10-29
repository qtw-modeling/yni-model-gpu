import subprocess
import timeit
import numpy as np
import sys
import numpy.random as rnd
#import resource
import os


mode_ = sys.argv[1]
variants = [mode_] # ['gpu_local']

timings = {key: [] for key in variants}
exe_files = {key: ("exe_yni_%s.exe" % (str(key)[0:3])) for key in variants} # for Windows; [0:3] == 'cpu' or 'gpu', as they consist of 3 symbols, each
#exe_files = {key: "./exe_yni_%s" % (str(key)) for key in variants} # for Linux
print('exe-files are: ', exe_files, '\n')
#input()

nums_cells_x_dim = [128, 256, 512, 1024] # [(int(2**n) - 1) for n in range(7, 12)] # (8, 15) --- original] # number of segments in x-dimension
# minus 1 for every dim size: dim size must be ODD
print('Cell maps dims: ', nums_cells_x_dim)

T = int(100.) # endtime in ms.

output_timing_dir = './timings'

# loop for Nxes = 32, 64, 128, ..., values
for num_cells_x_dim in nums_cells_x_dim:
    print ('Num cells in x-direction: %d; calculations begin\n' % num_cells_x_dim) # TODO: find out how to continue printing on the same string
    #num_segments_x = Nx - 1

    for mode in variants:
        subprocess.call( [exe_files[mode], str(num_cells_x_dim), str(T), mode]) # , str(tmp_output_file)] )#, str(serie_of_launches_num)] ) # Ny == Nx, for now

