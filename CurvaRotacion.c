#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* COMPILO CON cc CurvaRotacion.c -o CurvaRotacion.x -lm Y ./CurvaRotacion.x */

#define MAG_CAR 32
#define FILAS 301
#define COLS 2 
/* CTES DIMENSIONALES EN KPC */
#define B_b 0.2497
#define B_d 5.16
#define A_d 0.3105
#define A_h 64.3

void Lector_Datos(void); // LEE LOS DATOS, LOS GUARDA EN MATRICES Y GENERA ARCHIVO
double V_Rad(double R, double M_b, double M_d, double M_h); // FUNCION V_rad(R)
double Probabilidad(double D_exp, double D_rad); // Likelihood observados, modelo 
double Aleatorio(void);
double Prom(double *X);
void Imprime_en_Cons(void);
void M_H(int N); 
void Ajuste(void);

double *R_exp;
double *V_exp;
double *P;
double *M_b;
double *M_d;
double *M_h;
double *V_rad_ajuste;
double M_b_prom;
double M_d_prom;
double M_h_prom;
int N = 1000;

int main(){

	Lector_Datos();
	M_H(N);
	
	Imprime_en_Cons();
	
	Ajuste();

	return 0; 

}


void Lector_Datos(void){

	FILE *datos = fopen("RadialVelocities.dat","r");
	
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


	/* GENERA ARCHIVOS DE DATOS */ 

	FILE *Datos = fopen("Datos.dat", "w");
	
	for(int i = 0; i < FILAS; i++){
            	fprintf(Datos, "%lf" " " "%lf" ,  M_Datos[i][0], M_Datos[i][1]);
        	fprintf(Datos, "\n");
    	}
	fclose(Datos);

	/* MATRICES DE DATOS EXPERIMENTALES */ 

	R_exp = malloc(sizeof(double *)*FILAS);
	
	for(int i = 0; i < FILAS; i++){
            	R_exp[i] = M_Datos[i][0];
		//printf("R[%d]=%lf", i, R[i] );
        	//printf("\n");
    	}

	V_exp = malloc(sizeof(double *)*FILAS);

	for(int i = 0; i < FILAS; i++){
            	V_exp[i] = M_Datos[i][1];
		//printf("V[%d]=%lf", i, V[i]);
        	//printf("\n");
    	}
}


double Probabilidad(double D_exp, double D_rad){
	double chi_2 = pow((D_exp - D_rad)/0.1, 2); 
	return exp(-(1.0/2.0)*chi_2);
}


double V_Rad(double R_exp, double M_b, double M_d, double M_h){
	return (sqrt(M_b)*R_exp/pow(pow(R_exp, R_exp)+pow(B_b, B_b), 3/4))+(sqrt(M_d)*R_exp/pow(pow(R_exp, R_exp)+pow((B_b+A_d), (B_b+A_d)), 3/4))+(sqrt(M_h)/pow(pow(R_exp, R_exp)+pow(A_h, 2), 1/4));
}


double Aleatorio(void){
	double r = (double) rand()/RAND_MAX; // Lo arreglo para los bethas de N-H
	if(r !=1 ){
		return r; // LO NORMALIZO PORQUE NO QUIERO NUMEROS TAN GRANDES
	}
}


void M_H(int N){
	M_b = malloc(sizeof(double *)*N);
	M_d = malloc(sizeof(double *)*N);
	M_h = malloc(sizeof(double *)*N);
	P = malloc(sizeof(double *)*N); 

	double *Mb_prima = malloc(sizeof(double *)*N);
	double *Md_prima = malloc(sizeof(double *)*N);
	double *Mh_prima = malloc(sizeof(double *)*N);
	double *P_prima = malloc(sizeof(double *)*N);
	double *P_temp = malloc(sizeof(double *)*N);
		
	M_b[0] = Aleatorio();
	M_d[0] = Aleatorio();
	M_h[0] = Aleatorio();
	P[0] = Probabilidad(R_exp[0], V_Rad(R_exp[0], M_b[0], M_d[0], M_h[0]));

	/*printf("%lf", M_b[0]);
	printf("\n");
	printf("%lf", M_d[0]);
	printf("\n");
	printf("%lf", M_h[0]);
	printf("\n");
	printf("%lf", P[0]);
	printf("\n");*/
			

	for(int i = 0; i < N; i++){

        	Mb_prima[i] = Aleatorio();
		Md_prima[i] = Aleatorio();
		Mh_prima[i] = Aleatorio();
		P_prima[i] = Probabilidad(R_exp[i], V_Rad(R_exp[i], Mb_prima[i], Md_prima[i], Mh_prima[i]));
		
		P_temp[i] = Probabilidad(R_exp[i], V_Rad(R_exp[i], M_b[i], M_d[i], M_h[i]));
	
		double alpha = (double)abs(P_prima[i-1]/P_temp[i-1]);
	
		//printf("[%d]%lf", i, alpha);
		//printf("\n");
		
		if( alpha >=1.0 ){

			M_b[i] = Mb_prima[i-1];
			M_d[i] = Md_prima[i-1];
			M_h[i] = Mh_prima[i-1];
			P[i] = P_prima[i-1];

		}else{

			double betha = 1.0/Aleatorio();
	
			if( betha <= alpha ){

				M_b[i] = Mb_prima[i-1];
				M_d[i] = Md_prima[i-1];
				M_h[i] = Mh_prima[i-1];
				P[i] = P_prima[i-1];

			}else{

				M_b[i] = Mb_prima[i-1];
				M_d[i] = Md_prima[i-1];
				M_h[i] = Mh_prima[i-1];
				P[i] = P_temp[i-1];

			}
		}	
		//printf("%lf" " " "%lf", P_prima[i], P[i]);
		//printf("\n");
	
		//printf("[%d]%lf", i, M_b[i]);
		//printf("\n");
	}

}

double Prom(double *X){
	float suma = 0.0;
	for(int i = 0; i < N; i++){
            	suma += X[i];
    	}
	float prom = suma/N;

	return prom;
}

void Imprime_en_Cons(void){
	M_b_prom = Prom(M_b);
	M_d_prom = Prom(M_d);
	M_h_prom = Prom(M_h);

	printf("M_b=%lf", M_b_prom);
	printf("\n");
	printf("M_d=%lf", M_d_prom);
	printf("\n");
	printf("M_h=%lf", M_h_prom);
	printf("\n");
}

void Ajuste(void){
	M_b_prom = Prom(M_b);
	M_d_prom = Prom(M_d);
	M_h_prom = Prom(M_h);
		
	V_rad_ajuste = malloc(sizeof(double *)*N);

	for(int i = 0; i < N; i++){
            	V_rad_ajuste[i] = (sqrt(M_b_prom)*R_exp[i]/pow(pow(R_exp[i], R_exp[i])+pow(B_b, B_b), 3/4))+(sqrt(M_d_prom)*R_exp[i]/pow(pow(R_exp[i], R_exp[i])+pow((B_b+A_d), (B_b+A_d)), 3/4))+(sqrt(M_h_prom)/pow(pow(R_exp[i], R_exp[i])+pow(A_h, 2), 1/4));

	//printf("[%d]%lf", i, V_rad_ajuste[i]);
	//printf("\n");

    	}
	
	
	
	FILE *Ajuste = fopen("Ajuste.dat", "w");
	
	for(int i = 0; i < FILAS; i++){
            	fprintf(Ajuste, "%lf", V_rad_ajuste[i]);
        	fprintf(Ajuste, "\n");
    	}
	fclose(Ajuste);
	
}

