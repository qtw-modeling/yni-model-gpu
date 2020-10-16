import subprocess
import timeit
import numpy as np
import sys
import numpy.random as rnd
#import resource
import os


variants = ['gpu']

timings = {key: [] for key in variants}
#print(timings)
#input()

# old vers.
#timings = {variants[0]: []} # , variants[1]: []} 
#exe_files = {variants[0]: './ap_model_cpu', variants[1]: './ap_model_gpu'}
#exe_files = {variants[0]: "exe_yni_gpu.exe"} #{variants[0]: "exe_yni_cpu.exe", variants[1]: "exe_yni_gpu.exe"}

exe_files = {key: "exe_yni_%s.exe" % (str(key)) for key in variants}
#print(exe_files)
#input()


#num_launches = 1
nums_cells_x_dim = [(int(2**n) - 1) for n in range(7, 12)] # (8, 15) --- original] # number of segments in x-dimension
# minus 1 for every dim size: dim size must be ODD
print('Cell maps dims: ', nums_cells_x_dim)

T = int(600.) # endtime in ms.

#tmp_output_file = "tmp_output.txt"

#rnd.seed() # with no arguments the seed is set to current wall time;
# for explicit setting use this: timeit.default_timer())
#serie_of_launches_num = rnd.randint(0, 100)


output_timing_dir = './timings'

# loop for Nxes = 32, 64, 128, ..., values
#Nxes = Nxes[0:-3]  # uncomment this if testing
for num_cells_x_dim in nums_cells_x_dim:
    print ('Num cells in x-direction: %d; calculations begin...' % num_cells_x_dim) # TODO: find out how to continue printing on the same string
    #num_segments_x = Nx - 1

    for mode in variants:
        subprocess.call( [exe_files[mode], str(num_cells_x_dim), str(T)]) # , str(tmp_output_file)] )#, str(serie_of_launches_num)] ) # Ny == Nx, for now
    	
        out_file_timing_name = 'timing_dim_%dx%d_cells_T_%d_ms.txt' % (num_cells_x_dim, num_cells_x_dim, T)
        out_file_timing_full_name = output_timing_dir +  '/' + out_file_timing_name

        #info = resource.getrusage(resource.RUSAGE_CHILDREN)
    	#elapsed_time = timeit.default_timer() - start
    	#elapsed_time = info[0]
        
        out_file_timing = open(out_file_timing_full_name, 'r')
        elapsed_time = out_file_timing.read()

        # old version of reading the file
        # tmp_file = open(tmp_output_file, 'r') 
        # elapsed_time = tmp_file.read() 
        
        timings[mode].append(float(elapsed_time))
	
    if len(variants) == 2:
        speedup_local = timings['cpu'][-1]/timings['gpu'][-1]
        print ('Speedup local: %.2f' % speedup_local) # 7 --- for detailed output


# final actions
for mode in variants:
    timings[mode] = np.array(timings[mode])

if len(variants) == 1: # only 'gpu' vers. is present
    print('Elapsed times: ', timings[variants[0]], ' s.')

if len(variants) == 2: # if 'cpu' and 'gpu' modes are present
    speedups = timings['cpu']/timings['gpu']
    print ('Speedups:', speedups)


    # writing the data to a file
    output_2d = np.vstack((np.array(Nxes)**2, timings['cpu'], timings['gpu'], speedups)).transpose()
    current_dir_name = os.path.basename(os.path.normpath(os.getcwd())) # displays only the last 'part' of the current working dir
    np.savetxt('speedups_%s.txt' % current_dir_name, output_2d, header='Число клеток,Время расчета CPU,Время расчета GPU,Ускорение', fmt='%.7f', delimiter=',') # 7 --- for detailed outputsl