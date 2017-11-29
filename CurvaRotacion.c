#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAG_CAR 32
#define FILAS 301
#define COLS 2 

#define b_b 0.2497
#define b_d 5.16
#define a_d 0.3105
#define a_h 64.3

int main(){

	/* LEE LOS DATOS */	
	
	// LEE DATOS 

	FILE *datos; // datos es el archivo 
	datos = fopen("RadialVelocities.dat","r"); //cond_ini_cuerda

	// DEFINO ARRAY EN MEMORIA VARIABLE
	
	double **M_Datos;
    	M_Datos = malloc(sizeof(double *)*FILAS*COLS);

	// SEPARA LOS DATOS EN FILAS Y DESPUES EN COLUMNAS DE ACUERDO A EL DELIMITADOR

	char linea[MAG_CAR];
	char *split_linea;
		
	int i = 0, j = 0;
    	i = 0;
	while(fgets(linea, MAG_CAR, datos) != NULL){
		j = 0; 
		M_Datos[i] = malloc(sizeof(double)*MAG_CAR);
		split_linea = strtok(linea, " ");
		
		while (split_linea != NULL)
        	{
			M_Datos[i][j] = atof(split_linea);
            		split_linea = strtok(NULL, " ");
            		j += 1;
        	}
		i +=1;
	}

	// IMPRIME LOS DATOS	
	/*
	for(int j = 0; j < FILAS; j++){
        	for (int i = 0; i < COLS; i++){
			if( j != 0 ){
				printf("Datos[%d][%d]=%lf", j, i, M_Datos[j][i] );
        			printf("\n");
			}
		}
	}
	*/

	/* GENERA ARCHIVOS DE DATOS */ 

	FILE *Datos = fopen("Datos.dat", "w");
	
	for(int i = 0; i < FILAS; i++)
	{
            	fprintf(Datos, "%lf" " " "%lf" ,  M_Datos[i][0], M_Datos[i][1]);
        	fprintf(Datos, "\n");
    	}
	fclose(Datos);
	
	return 0; 

}
