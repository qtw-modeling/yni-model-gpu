compile_cpu:
	g++ -std=c++11 -g main_yni.cpp CurrentSlow.cpp CurrentNa.cpp CurrentK.cpp CurrentLeak.cpp CurrentHyperpolar.cpp -o yni_exec

#compile_gpu:
#	pgc++ main_yni.cpp -acc -Minfo=accel -ta=nvidia -o yni_exec

run:
	./yni_exec

clear_output:
	rm output/*.vtk