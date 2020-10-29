import subprocess
import timeit
import numpy as np
import sys
import numpy.random as rnd
#import resource
import os


mode_ = sys.argv[1]
variants = [mode_] # ['cpu', 'gpu_local', 'gpu_cloud']

timings = {key: [] for key in variants}
speedups = {key: [] for key in variants} # init with an empty list for now

nums_cells_x_dim = [128, 256, 512, 1024] # [(int(2**n) - 1) for n in range(7, 12)] # (8, 15) --- original] # number of segments in x-dimension
# minus 1 for every dim size: dim size must be ODD
print('Cell maps dims: ', nums_cells_x_dim)

T = int(100.) # endtime in ms.


output_timing_dir = './timings'

# after starting all series of launches: 
# go through the directory with timings and extract data into arrays;
# then, write a csv-file with all timings
for num_cells_x_dim in nums_cells_x_dim:
    for mode in variants:
        out_file_timing_name = mode + '_timing_dim_%dx%d_cells_T_%d_ms.txt' % (num_cells_x_dim, num_cells_x_dim, T)
        out_file_timing_full_name = output_timing_dir +  '/' + out_file_timing_name
        
        out_file_timing = open(out_file_timing_full_name, 'r')
        elapsed_time = out_file_timing.read()

        # old version of reading the file
        # tmp_file = open(tmp_output_file, 'r') 
        # elapsed_time = tmp_file.read() 
        
        timings[mode].append(float(elapsed_time))


# final actions
for mode in variants:
    timings[mode] = np.array(timings[mode])

# calculating speedups
for mode in variants:
    speedups[mode] = timings['cpu']/timings[mode]

nums_cells_x_dim_str = [(str(elem) + 'x' + str(elem)) for elem in nums_cells_x_dim]

# writing the data to a csv-file
#output_2d = np.vstack(nums_cells_x_dim_str, timings['cpu_local'], timings['gpu_local'], timings['gpu_cloud'], speedups)).transpose()
output_2d = np.vstack(nums_cells_x_dim_str, speedups['cpu_local'], speedups['gpu_local'], speedups['gpu_cloud']).transpose()
current_dir_name = os.path.basename(os.path.normpath(os.getcwd())) # displays only the last 'part' of the current working dir
np.savetxt('speedups_%s.txt' % current_dir_name, output_2d, header='Число клеток,Время расчета CPU,Время расчета GPU,Ускорение', fmt='%.7f', delimiter=',') # 7 --- for detailed outputsl
########


###### old actions
'''
#if len(variants) == 1: # only 'gpu' vers. is present
#    print('Elapsed times: ', timings[variants[0]], ' s.')

if len(variants) == 2: # if 'cpu' and 'gpu' modes are present
    speedups = timings['cpu']/timings['gpu']
    print ('Speedups:', speedups)


    # writing the data to a file
    output_2d = np.vstack((np.array(Nxes)**2, timings['cpu'], timings['gpu'], speedups)).transpose()
    current_dir_name = os.path.basename(os.path.normpath(os.getcwd())) # displays only the last 'part' of the current working dir
    np.savetxt('speedups_%s.txt' % current_dir_name, output_2d, header='Число клеток,Время расчета CPU,Время расчета GPU,Ускорение', fmt='%.7f', delimiter=',') # 7 --- for detailed outputsl
########
'''
