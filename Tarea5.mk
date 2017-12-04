#### TAREA 5 

# compilo con make -f Tarea5.mk

Resultados_hw5.pdf : Resultados_hw5.tex *.png
	pdflatex $<

*.png : *.dat Plots.py
	python Plots.py

*.dat : CurvaRotacion.x 
	./CurvaRotacion.x 

CurvaRotacion.x : CurvaRotacion.c
	cc CurvaRotacion.c -o CurvaRotacion.x -lm

clean :
	rm *.dat  *.x *.png *.log *.aux 
