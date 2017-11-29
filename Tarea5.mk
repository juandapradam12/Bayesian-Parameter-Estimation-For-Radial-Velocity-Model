#### TAREA 5 

# compilo con make -f Tarea5.mk

*.png : *.dat Plots.py
	python Plots.py

*.dat : CurvaRotacion.x 
	./CurvaRotacion.x 

CurvaRotacion.x : CurvaRotacion.c
	cc CurvaRotacion.c -o CurvaRotacion.x -lm

clean :
	rm *.dat  *.x *.png *.log *.aux 
