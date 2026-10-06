# Legacy C pipeline (educational reference)
# make -f Tarea5.mk

.PHONY: all clean

all: DatosyModelo.png

CurvaRotacion.x: CurvaRotacion.c
	cc -std=c99 -O2 CurvaRotacion.c -o CurvaRotacion.x -lm

Datos.dat Ajuste.dat chain.dat: CurvaRotacion.x ../data/RadialVelocities.dat
	./CurvaRotacion.x

DatosyModelo.png: Datos.dat Ajuste.dat Plots.py
	python3 Plots.py

clean:
	rm -f *.dat *.x *.png *.log *.aux chain.dat
