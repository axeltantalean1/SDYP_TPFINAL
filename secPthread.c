// Compilar:
//		gcc -o n_body_simple_NOGL n_body_simple.c -lm
// Ejecutar:
//		./n_body_simple_NOGL <nro de cuerpos> <DT> <Pasos>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include <pthread.h>

//
// Para tiempo de ejecucion
//

double dwalltime()
{
	double sec;
	struct timeval tv;

	gettimeofday(&tv,NULL);
	sec = tv.tv_sec + tv.tv_usec/1000000.0;
	return sec;
}

double tIni, tFin, tTotal;

//
// Constantes para Algoritmo de gravitacion
//
#define PI (3.141592653589793)
#define G 6.673e-11
#define ESTRELLA 0
#define POLVO 1
#define H2 2 //Hidrogeno molecular

// ===============
// ===== CPU =====
// ===============

//
// Estructuras y variables para Algoritmo de gravitacion
//
typedef struct cuerpo cuerpo_t;
struct cuerpo{
	float masa;
	float px;
	float py;
	float pz;
	float vx;
	float vy;
	float vz;
	float r;
	float g;
	float b;
	int cuerpo;
};

float *fuerza_totalX,*fuerza_totalY, *fuerza_totalZ, *fuerzas_parcialesX,*fuerzas_parcialesY,*fuerzas_parcialesZ;
float toroide_alfa;
float toroide_theta;
float toroide_incremento;
float toroide_lado;
float toroide_r;
float toroide_R;

cuerpo_t *cuerpos;
cuerpo_t *cuerposSec;
int dt = 1.0f; //Intervalo de tiempo, longitud de un paso
int pasos;
int N;
int CP;

pthread_barrier_t barrier;

//
// Funciones para Algoritmo de gravitacion
//

void calcularFuerzas(int id){
    int cuerpo1, cuerpo2,idN;
    float dif_X, dif_Y, dif_Z;
    float distancia;
    float F;
    idN=id*N;
	for(cuerpo1 = id; cuerpo1<N-1 ; cuerpo1+=CP){
		for(cuerpo2 = cuerpo1 + 1; cuerpo2<N ; cuerpo2++){
			if ( (cuerpos[cuerpo1].px == cuerpos[cuerpo2].px) && (cuerpos[cuerpo1].py == cuerpos[cuerpo2].py) && (cuerpos[cuerpo1].pz == cuerpos[cuerpo2].pz))
                continue;

            dif_X = cuerpos[cuerpo2].px - cuerpos[cuerpo1].px;
            dif_Y = cuerpos[cuerpo2].py - cuerpos[cuerpo1].py;
            dif_Z = cuerpos[cuerpo2].pz - cuerpos[cuerpo1].pz;
            distancia = sqrt(dif_X*dif_X + dif_Y*dif_Y + dif_Z*dif_Z);

            F = (G*cuerpos[cuerpo1].masa*cuerpos[cuerpo2].masa)/(distancia*distancia);
            dif_X *= F;
            dif_Y *= F;
            dif_Z *= F;


            //Se podria realizar una matriz triangular y ahorrar computo
            fuerzas_parcialesX[idN+cuerpo1] += dif_X;
            fuerzas_parcialesY[idN+cuerpo1] += dif_Y;
            fuerzas_parcialesZ[idN+cuerpo1] += dif_Z;

            fuerzas_parcialesX[idN+cuerpo2] -= dif_X;
            fuerzas_parcialesY[idN+cuerpo2] -= dif_Y;
            fuerzas_parcialesZ[idN+cuerpo2] -= dif_Z;

		}
	}
}

void moverCuerpos(cuerpo_t *cuerpos, int N, int dt,int id){
    int cuerpo;
    int ini=id*N/CP;
    int fin=N*(id+1)/CP;
	for(cuerpo=ini;cuerpo<fin;cuerpo++){

        fuerza_totalX[cuerpo] *= 1/cuerpos[cuerpo].masa;
        fuerza_totalY[cuerpo] *= 1/cuerpos[cuerpo].masa;
        //fuerza_totalZ[cuerpo] *= 1/cuerpos[cuerpo].masa;

        cuerpos[cuerpo].vx += fuerza_totalX[cuerpo]*dt;
        cuerpos[cuerpo].vy += fuerza_totalY[cuerpo]*dt;
        //cuerpos[cuerpo].vz += fuerza_totalZ[cuerpo]*dt;

        cuerpos[cuerpo].px += cuerpos[cuerpo].vx *dt;
        cuerpos[cuerpo].py += cuerpos[cuerpo].vy *dt;
        //cuerpos[cuerpo].pz += cuerpos[cuerpo].vz *dt;

        fuerza_totalX[cuerpo] = 0.0;
		fuerza_totalY[cuerpo] = 0.0;
		fuerza_totalZ[cuerpo] = 0.0;

    	}
}

void sumarFuerzasParciales(int id){
    int colIni=N*id/CP;
    int colFin=N*(id+1)/CP;
    int i,j;
    for(i=colIni;i<colFin;i++){
        for(j=0;j<CP;j++){
            fuerza_totalX[i]+=fuerzas_parcialesX[j*N+i];
            fuerza_totalY[i]+=fuerzas_parcialesY[j*N+i];
            fuerza_totalZ[i]+=fuerzas_parcialesZ[j*N+i];
			fuerzas_parcialesX[j*N+i] = 0;
			fuerzas_parcialesY[j*N+i] = 0;
			fuerzas_parcialesZ[j*N+i] = 0;
        }
    }
}

void gravitacionCPU(int id){

    calcularFuerzas(id);

    pthread_barrier_wait(&barrier);
	
    sumarFuerzasParciales(id);
moverCuerpos(cuerpos,N,dt,id);
    pthread_barrier_wait(&barrier);

}

void *funcionThread(void *arg){
    int id = (int)arg;
	int paso;
	for(paso=0; paso<pasos; paso++){
		gravitacionCPU(id);
	}    
}

void inicializarEstrella(cuerpo_t *cuerpo,int i,double n){

    cuerpo->masa = 0.001*8; 

        if ((toroide_alfa + toroide_incremento) >=2*M_PI){
            toroide_alfa = 0;
            toroide_theta += toroide_incremento;
        }else{
            toroide_alfa+=toroide_incremento;
        }

	cuerpo->px = (toroide_R + toroide_r*cos(toroide_alfa))*cos(toroide_theta);
	cuerpo->py = (toroide_R + toroide_r*cos(toroide_alfa))*sin(toroide_theta);
	cuerpo->pz = toroide_r*sin(toroide_alfa);

    	cuerpo->vx = 0.0;
	cuerpo->vy = 0.0;
	cuerpo->vz = 0.0;

		cuerpo->r = 1.0; //(double )rand()/(RAND_MAX+1.0);
		cuerpo->g = 1.0; //(double )rand()/(RAND_MAX+1.0);
		cuerpo->b = 1.0; //(double )rand()/(RAND_MAX+1.0);
}

void inicializarPolvo(cuerpo_t *cuerpo,int i,double n){

    cuerpo->masa = 0.001*4;
	
        if ((toroide_alfa + toroide_incremento) >=2*M_PI){
            toroide_alfa = 0;
            toroide_theta += toroide_incremento;
        }else{
            toroide_alfa+=toroide_incremento;
        }

	cuerpo->px = (toroide_R + toroide_r*cos(toroide_alfa))*cos(toroide_theta);
	cuerpo->py = (toroide_R + toroide_r*cos(toroide_alfa))*sin(toroide_theta);
	cuerpo->pz = toroide_r*sin(toroide_alfa);
	
	cuerpo->vx = 0.0;
	cuerpo->vy = 0.0;
	cuerpo->vz = 0.0;
    
	cuerpo->r = 1.0; //(double )rand()/(RAND_MAX+1.0);
	cuerpo->g = 0.0; //(double )rand()/(RAND_MAX+1.0);
	cuerpo->b = 0.0; //(double )rand()/(RAND_MAX+1.0);
}

void inicializarH2(cuerpo_t *cuerpo,int i,double n){

    cuerpo->masa = 0.001;

	if ((toroide_alfa + toroide_incremento) >=2*M_PI){
            toroide_alfa = 0;
            toroide_theta += toroide_incremento;
	}else{
            toroide_alfa+=toroide_incremento;
	}

	cuerpo->px = (toroide_R + toroide_r*cos(toroide_alfa))*cos(toroide_theta);
	cuerpo->py = (toroide_R + toroide_r*cos(toroide_alfa))*sin(toroide_theta);
	cuerpo->pz = toroide_r*sin(toroide_alfa);

	cuerpo->vx = 0.0;
	cuerpo->vy = 0.0;
	cuerpo->vz = 0.0;

	cuerpo->r = 1.0; //(double )rand()/(RAND_MAX+1.0);
	cuerpo->g = 1.0; //(double )rand()/(RAND_MAX+1.0);
	cuerpo->b = 1.0; //(double )rand()/(RAND_MAX+1.0);
}

void inicializarCuerpos(cuerpo_t *cuerpos,int N){
 int cuerpo;
 double n = N;

	toroide_alfa = 0.0;
	toroide_theta = 0.0;
	toroide_lado = sqrt(N);
	toroide_incremento = 2*M_PI / toroide_lado;
	toroide_r = 1.0;
	toroide_R = 2*toroide_r;
	
	srand(time(NULL));

	for(cuerpo = 0; cuerpo < N; cuerpo++){
	
        fuerza_totalX[cuerpo] = 0.0;
		fuerza_totalY[cuerpo] = 0.0;
		fuerza_totalZ[cuerpo] = 0.0;

		cuerpos[cuerpo].cuerpo = 0;//(rand() %3);

		if (cuerpos[cuerpo].cuerpo == ESTRELLA){
			inicializarEstrella(&cuerpos[cuerpo],cuerpo,n);
		}else if (cuerpos[cuerpo].cuerpo == POLVO){
			inicializarPolvo(&cuerpos[cuerpo],cuerpo,n);
		}else if (cuerpos[cuerpo].cuerpo == H2){
			inicializarH2(&cuerpos[cuerpo],cuerpo,n);
		}

	}

		//cuerpos[0].masa = 2.0e2;
/*	        cuerpos[0].px = 0.0;
		cuerpos[0].py = 0.0;
		cuerpos[0].pz = 0.0;
		cuerpos[0].vx = -0.000001;
		cuerpos[0].vy = -0.000001;
		cuerpos[0].vz = 0.0;

		cuerpos[1].masa = 1.0e1;
	        cuerpos[1].px = -1.0;
		cuerpos[1].py = 0.0;
		cuerpos[1].pz = 0.0;
		cuerpos[1].vx = 0.0;
		cuerpos[1].vy = 0.0001;
		cuerpos[1].vz = 0.0;
*/
}


void calcularFuerzasSec(cuerpo_t *cuerpos, int N, int dt){
int cuerpo1, cuerpo2;
float dif_X, dif_Y, dif_Z;
float distancia;
float F;

	for(cuerpo1 = 0; cuerpo1<N-1 ; cuerpo1++){
		for(cuerpo2 = cuerpo1 + 1; cuerpo2<N ; cuerpo2++){
			if ( (cuerpos[cuerpo1].px == cuerpos[cuerpo2].px) && (cuerpos[cuerpo1].py == cuerpos[cuerpo2].py) && (cuerpos[cuerpo1].pz == cuerpos[cuerpo2].pz))
                continue;

	            	dif_X = cuerpos[cuerpo2].px - cuerpos[cuerpo1].px;
			dif_Y = cuerpos[cuerpo2].py - cuerpos[cuerpo1].py;
			dif_Z = cuerpos[cuerpo2].pz - cuerpos[cuerpo1].pz;
                
			distancia = sqrt(dif_X*dif_X + dif_Y*dif_Y + dif_Z*dif_Z);

	                F = (G*cuerpos[cuerpo1].masa*cuerpos[cuerpo2].masa)/(distancia*distancia);

	                dif_X *= F;
			dif_Y *= F;
			dif_Z *= F;

	                fuerza_totalX[cuerpo1] += dif_X;
	                fuerza_totalY[cuerpo1] += dif_Y;
	                fuerza_totalZ[cuerpo1] += dif_Z;

	                fuerza_totalX[cuerpo2] -= dif_X;
	                fuerza_totalY[cuerpo2] -= dif_Y;
	                fuerza_totalZ[cuerpo2] -= dif_Z;
		}
	}
}

void moverCuerposSec(cuerpo_t *cuerpos, int N, int dt){
 int cuerpo;
	for(cuerpo = 0; cuerpo<N ; cuerpo++){

        	fuerza_totalX[cuerpo] *= 1/cuerpos[cuerpo].masa;
        	fuerza_totalY[cuerpo] *= 1/cuerpos[cuerpo].masa;
        	//fuerza_totalZ[cuerpo] *= 1/cuerpos[cuerpo].masa;

        	cuerpos[cuerpo].vx += fuerza_totalX[cuerpo]*dt;
        	cuerpos[cuerpo].vy += fuerza_totalY[cuerpo]*dt;
        	//cuerpos[cuerpo].vz += fuerza_totalZ[cuerpo]*dt;

        	cuerpos[cuerpo].px += cuerpos[cuerpo].vx *dt;
        	cuerpos[cuerpo].py += cuerpos[cuerpo].vy *dt;
        	//cuerpos[cuerpo].pz += cuerpos[cuerpo].vz *dt;

        	fuerza_totalX[cuerpo] = 0.0;
		fuerza_totalY[cuerpo] = 0.0;
		fuerza_totalZ[cuerpo] = 0.0;

    	}
}

void gravitacionCPUSec(cuerpo_t *cuerpos, int N, int dt){
	calcularFuerzasSec(cuerpos,N,dt);
	moverCuerposSec(cuerpos,N,dt);
}

void finalizar(void){
	free(cuerpos);
	free(fuerza_totalX);
	free(fuerza_totalY);
	free(fuerza_totalZ);
    free(fuerzas_parcialesX);
    free(fuerzas_parcialesY);
    free(fuerzas_parcialesZ);
    pthread_barrier_destroy(&barrier);
}

void copiarCuerpos(cuerpo_t *cuerpos, cuerpo_t *cuerposSec, int N){
    for(int i=0;i<N;i++){
        cuerposSec[i].masa = cuerpos[i].masa;
        cuerposSec[i].px = cuerpos[i].px;
        cuerposSec[i].py = cuerpos[i].py;
        cuerposSec[i].pz = cuerpos[i].pz;
        cuerposSec[i].vx = cuerpos[i].vx;
        cuerposSec[i].vy = cuerpos[i].vy;
        cuerposSec[i].vz = cuerpos[i].vz;
        cuerposSec[i].r = cuerpos[i].r;
        cuerposSec[i].g = cuerpos[i].g;
        cuerposSec[i].b = cuerpos[i].b;
        cuerposSec[i].cuerpo = cuerpos[i].cuerpo;
    }
}

int main(int argc, char * argv[]) {

	if (argc < 5){
		printf("Ejecutar: %s <nro. de cuerpos> <DT> <pasos> <cant procesos>\n",argv[0]);
		return -1;
	}

	N = atoi(argv[1]);
	dt = atof(argv[2]);
	pasos = atoi(argv[3]);
    CP = atoi(argv[4]);

	pthread_t misThreads[CP];
    pthread_barrier_init(&barrier,NULL,CP);

    cuerpos = (cuerpo_t*)malloc(sizeof(cuerpo_t)*N);
    cuerposSec = (cuerpo_t*)malloc(sizeof(cuerpo_t)*N); // Agregado malloc para cuerposSec
    fuerza_totalX = (float*)calloc(N,sizeof(float));
    fuerza_totalY = (float*)calloc(N,sizeof(float));
    fuerza_totalZ = (float*)calloc(N,sizeof(float));
    fuerzas_parcialesX = (float*)calloc(N*CP,sizeof(float));
    fuerzas_parcialesY = (float*)calloc(N*CP,sizeof(float));
    fuerzas_parcialesZ = (float*)calloc(N*CP,sizeof(float));

	inicializarCuerpos(cuerpos,N);
    copiarCuerpos(cuerpos,cuerposSec,N);
	
	tIni = dwalltime(); 

    for(int id=0;id<CP;id++){

		pthread_create(&misThreads[id],NULL,&funcionThread,(void*)id);
	}
	
	for(int id=0; id<CP;id++){
		pthread_join(misThreads[id],NULL);
	}

	tFin =	dwalltime();
	tTotal = tFin - tIni;
	
	printf("Tiempo en segundos paralelo: %f\n",tTotal);

	tIni = dwalltime(); 
    for(int i=0; i < pasos; i++){
        gravitacionCPUSec(cuerposSec,N,dt);
    }

	tFin =	dwalltime();
	tTotal = tFin - tIni;
	
	printf("Tiempo en segundos secuencial: %f\n",tTotal);

	// Comparar resultados entre versión paralela y secuencial
	int errores = 0;
	for(int i = 0; i < N; i++) {
					printf("Error en cuerpo %d:\n", i);
			printf("Paralelo  - px: %f, py: %f, pz: %f\n", 
				cuerpos[i].px, cuerpos[i].py, cuerpos[i].pz);
			printf("Secuencial- px: %f, py: %f, pz: %f\n",
				cuerposSec[i].px, cuerposSec[i].py, cuerposSec[i].pz);
		
	}
	
	if(errores == 0) {
		printf("Los resultados paralelos y secuenciales coinciden\n");
	} else {
		printf("Se encontraron %d diferencias entre resultados\n", errores);
	}

	finalizar();
    return(0);

}