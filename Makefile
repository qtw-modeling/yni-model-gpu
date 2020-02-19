compile_for_cpu:
	g++ main_yni.cpp -o yni_exec

compile_for_gpu:
	pgc++ main_yni.cpp -acc -Minfo=accel -ta=nvidia -o yni_exec

run:
	./yni_exec
