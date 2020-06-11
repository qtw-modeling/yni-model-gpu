all: compile_cpu compile_gpu

compile_cpu:
	pgcc main_yni.c CurrentSlow.c CurrentNa.c CurrentK.c CurrentLeak.c CurrentHyperpolar.c extra.c -o exec_yni_cpu

compile_gpu:
	pgcc -acc -fast -ta=nvidia -Minfo=accel main_yni.c CurrentSlow.c CurrentNa.c CurrentK.c CurrentLeak.c CurrentHyperpolar.c extra.c -o exec_yni_gpu

clear_output:
	rm output/*.vtk

phase:
	pgcc WritePhaseDataIntoFile.c CurrentSlow.c CurrentNa.c CurrentK.c CurrentLeak.c CurrentHyperpolar.c -o exec_yni_phase