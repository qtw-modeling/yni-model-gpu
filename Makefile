compile_cpu:
	g++ main_yni.cpp CurrentNa.cpp CurrentK.cpp -o yni_exec

#compile_gpu:
#	pgc++ main_yni.cpp -acc -Minfo=accel -ta=nvidia -o yni_exec

run:
	./yni_exec
