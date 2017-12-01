#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* COMPILO CON cc CurvaRotacion.c -o CurvaRotacion.x -lm Y ./CurvaRotacion.x */

#define MAG_CAR 32
#define FILAS 302
#define COLS 2 
/* CTES DIMENSIONALES EN KPC */
#define B_b 0.2497
#define B_d 5.16
#define A_d 0.3105
#define A_h 64.3

void Lector_Datos(void); // LEE LOS DATOS, LOS GUARDA EN MATRICES Y GENERA ARCHIVO
double V_Rad(double M_b, double M_d, double M_h); // FUNCION V_rad(R)
double Error(double v_rand);

double *R;
double *V;

int main(){

	Lector_Datos();

	return 0; 

}

void Lector_Datos(void){

	FILE *datos = fopen("RadialVelocities.dat","r");
	
	double **M_Datos;
    	M_Datos = malloc(sizeof(double *)*FILAS*COLS);

	char linea[MAG_CAR];
	char *split_linea;
		
	int i = 0, j = 0;
    	i = 0;
	while(fgets(linea, MAG_CAR, datos) != NULL){
		j = 0; 
		M_Datos[i] = malloc(sizeof(double)*MAG_CAR);
		split_linea = strtok(linea, " ");
		
		while (split_linea != NULL){
			M_Datos[i][j] = atof(split_linea);
            		split_linea = strtok(NULL, " ");
            		j += 1;
        	}
		i +=1;
	}
	
	FILE *Datos = fopen("Datos.dat", "w");
	
	for(int i = 0; i < FILAS-1; i++){
            	fprintf(Datos, "%lf" " " "%lf" ,  M_Datos[i][0], M_Datos[i][1]);
        	fprintf(Datos, "\n");
    	}
	fclose(Datos);

	R = malloc(sizeof(double *)*FILAS);
	
	for(int i = 0; i < FILAS; i++){
            	R[i] = M_Datos[i][0];
    	}
	
	V = malloc(sizeof(double *)*FILAS);
	
	for(int i = 0; i < FILAS; i++){
            	R[i] = M_Datos[i][1];
    	}
}

double V_Rad(double M_b, double M_d, double M_h){
	double R;
	return (sqrt(M_b)*R/pow(pow(R, R)+pow(B_b, B_b), 3/4))+(sqrt(M_d)*R/pow(pow(R, R)+pow((B_b+A_d), (B_b+A_d)), 3/4))+(sqrt(M_h)/pow(pow(R, R)+pow(A_h, A_h), 1/4));
}

double Error(double v_rand){
	double V; 
	double suma = 0.0; 
	for(int i = 0; i < FILAS; i++){
            	suma += (pow(V-v_rand, 2));
    	}
	return suma;
}




