all: cpu gpu

cpu:
	pgcc main_yni.c CurrentSlow.c CurrentNa.c CurrentK.c CurrentLeak.c CurrentHyperpolar.c extra.c -o exe_yni_cpu

gpu:
	pgcc -acc -fast -ta=tesla -Minfo=accel main_yni.c CurrentSlow.c CurrentNa.c CurrentK.c CurrentLeak.c CurrentHyperpolar.c extra.c -o exe_yni_gpu

clear_out:
	rm output/*.vtk

phase:
	pgcc WritePhaseDataIntoFile.c CurrentSlow.c CurrentNa.c CurrentK.c CurrentLeak.c CurrentHyperpolar.c -o exe_yni_phase