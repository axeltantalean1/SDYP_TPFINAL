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

//
// Para tiempo de ejecucion
//

double dwalltime(){
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

cuerpo_t *cuerpos, *cuerpos_recibidos;
int dt = 1.0f; //Intervalo de tiempo, longitud de un paso
int pasos;
int N;
int idW;
int subN;
int CP;
int T;
pthread_barrier_t barrier;

//
// Funciones para Algoritmo de gravitacion
//

//rango va a ser de 0 a subN 
// N/T * (idW + 1) a N


void calcularFuerzas(int id,int length, cuerpos_t *cuerpos, int isLocal){
    int cuerpo1, cuerpo2,idN;
    float dif_X, dif_Y, dif_Z;
    float distancia;
    float F;
	int inicio;
	int idN = id*N;
	inicio = subN * (idW + 1);
	for(cuerpo1 = id; cuerpo1< subN ; cuerpo1+CP){
		if(isLocal){
			inicio = cuerpo1 + 1;
		}
		for(cuerpo2 = inicio; cuerpo2 < length ; cuerpo2++){
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

            
            fuerzas_parcialesX[idN + cuerpo1] += dif_X;
            fuerzas_parcialesY[idN + cuerpo1] += dif_Y;
            fuerzas_parcialesZ[idN + cuerpo1] += dif_Z;

			fuerzas_parcialesX[idN + cuerpo2] -= dif_X;
			fuerzas_parcialesY[idN + cuerpo2] -= dif_Y;
			fuerzas_parcialesZ[idN + cuerpo2] -= dif_Z;
		}
	}
}

void moverCuerpos(cuerpo_t *cuerpos, int length, int dt){
    int cuerpo;
    int ini=id*length/CP;
    int fin=length*(id+1)/CP;
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

// inicio global subN * idW
// fin global N


// N - subN * (idW)

void sumarFuerzasParciales(int id){
	int total = subN * (T - idW);
    int colIni=total*id/CP;
    int colFin=total*(id+1)/CP;
    int i,j;
    for(i=colIni;i<colFin;i++){
        for(j=0;j<CP;j++){
            fuerza_totalX[i]+=fuerzas_parcialesX[j*N+i];
            fuerza_totalY[i]+=fuerzas_parcialesY[j*N+i];
            fuerza_totalZ[i]+=fuerzas_parcialesZ[j*N+i];
        }
    }
}

void gravitacionCPU(int id){
	//Barrera para empezar
    calcularFuerzas(id, subN, cuerpos, 1);
    pthread_barrier_wait(&barrier); //barrera para esperar el resto de posiciones
	calcularFuerzas(id, N, cuerpos, 0); // barrera para esperar que todos terminen 
	pthread_barrier_wait(&barrier); // barrera para esperar que todos terminen sus fuerzas
	sumarFuerzasParciales(id);
	pthread_barrier_wait(&barrier); // Le avisamos que tenemos todas las fuerzas al "root"
	pthread_barrier_wait(&barrier); //barrera de comunicacion
	//sumarFuerzasParciales(id);
	//pthread_barrier_wait(&barrier); // barrera para esperar que todos terminen sus fuersza
	moverCuerpos(id, subN, dt);
	pthread_barrier_wait(&barrier); // avisar que todos los cuerpos se han movido

}

void *funcionThread(void *arg){
    int id = (int)id;
	for(i=0;i<pasos;i++){
		barrera
		calcularFuerzas
		
	}
    
}



void funcionProcesoA(){
	

	for (int i=0 ;i<idW;i++){
		MPI_Isend(cuerpos,subN*sizeof(cuerpo_t),MPI_BYTE,i,0,MPI_COMM_WORLD,MPI_Request *req); //SOLUCIONAR MPI_REQUEST
	}

	//levantar barrera de pthreads para procesamiento local

	for(int j=id + 1;j<T;j++){
		MPI_Recv(cuerpos + subN*(j-idW + 1),subN*sizeof(cuerpo_t),MPI_BYTE,j,0,MPI_COMM_WORLD,MPI_Status *status); //SOLUCIONAR MPI_REQUEST
	}

	//levanto barrera de pthreads para procesamiento total

	//SOLUCION A

	for(int j=idW + 1;j<T;j++){
		MPI_Isend(fuerza_totalX + subN*(j-idW + 1),subN*sizeof(float),MPI_FLOAT,j,0,MPI_COMM_WORLD);
		MPI_Isend(fuerza_totalY + subN*(j-idW + 1),subN*sizeof(float),MPI_FLOAT,j,0,MPI_COMM_WORLD);
		MPI_Isend(fuerza_totalZ + subN*(j-idW + 1),subN*sizeof(float),MPI_FLOAT,j,0,MPI_COMM_WORLD);
	}

	//HABRIA QUE CREAR UNA MATRIZ DE FUERZAS
	for(int j=0 ;j<idW;j++){
		MPI_Recv(fuerza_totalX + subN*j,subN*sizeof(float),MPI_FLOAT,j,0,MPI_COMM_WORLD,MPI_Status *status); //Ver el tema de status
		MPI_Recv(fuerza_totalY + subN*j,subN*sizeof(float),MPI_FLOAT,j,0,MPI_COMM_WORLD,MPI_Status *status);
		MPI_Recv(fuerza_totalZ + subN*j,subN*sizeof(float),MPI_FLOAT,j,0,MPI_COMM_WORLD,MPI_Status *status);
	}
	
	//Levantar barrera para hacer la suma de fuerzas

	//SOLUCION B

	for(int j = idW;j<T;j++){
		MPI_Reduce(fuerza_totalX + subN*(j-idW), fuerza_totalX + subN*(j-idW), subN, MPI_FLOAT, MPI_SUM, j, MPI_COMM_WORLD);
		MPI_Reduce(fuerza_totalY + subN*(j-idW), fuerza_totalY + subN*(j-idW), subN, MPI_FLOAT, MPI_SUM, j, MPI_COMM_WORLD);
		MPI_Reduce(fuerza_totalZ + subN*(j-idW), fuerza_totalZ + subN*(j-idW), subN, MPI_FLOAT, MPI_SUM, j, MPI_COMM_WORLD);
	}


	//Levantar barrera para mover cuerpos
	

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

		cuerpos[cuerpo].cuerpo = (rand() %3);

		if (cuerpos[cuerpo].cuerpo == ESTRELLA){
			inicializarEstrella(&cuerpos[cuerpo],cuerpo,n);
		}else if (cuerpos[cuerpo].cuerpo == POLVO){
			inicializarPolvo(&cuerpos[cuerpo],cuerpo,n);
		}else if (cuerpos[cuerpo].cuerpo == H2){
			inicializarH2(&cuerpos[cuerpo],cuerpo,n);
		}

	}

		cuerpos[0].masa = 2.0e2;
	        cuerpos[0].px = 0.0;
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
}
void inicializarMemoria(){
    int indice = subN * (T - idW);
    cuerpos=(cuerpo_t*)malloc(sizeof(cuerpo_t)*indice);
    fuerza_totalX = (float*)malloc(sizeof(float)*indice);
    fuerza_totalY = (float*)malloc(sizeof(float)*indice);
    fuerza_totalZ = (float*)malloc(sizeof(float)*indice);
    fuerzas_parcialesX = (float*)malloc(sizeof(float)*indice*CP);
    fuerzas_parcialesY = (float*)malloc(sizeof(float)*indice*CP);
    fuerzas_parcialesZ = (float*)malloc(sizeof(float)*indice*CP);
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

void procesoA(){
    int i;
    for(i=0;i<pasos;i++){
        funcionProcesoA();
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
	subN = N/T;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&idW);
    MPI_Comm_size(MPI_COMM_WORLD,&T);
    inicializarMemoria();
    if(idW==0){
        inicializarCuerpos(cuerpos,N);        
    }
    else{

    }

	pthread_t misThreads[CP];
    pthread_barrier_init(&barrier,NULL,CP + 1);

	inicializarCuerpos(cuerpos,N);
	
	tIni = dwalltime(); 
	
    for(int id=0;id<CP;id++){
		pthread_create(&misThreads[id],NULL,&funcionThread,(void*)id);
	}
	
	for(int id=0; id<T;id++){
		pthread_join(misThreads[id],NULL);
	}


	int paso;
	for(paso=0; paso<pasos; paso++){
		gravitacionCPU(cuerpos,N,delta_tiempo);
	}

	tFin =	dwalltime();
	tTotal = tFin - tIni;
	
	printf("Tiempo en segundos: %f\n",tTotal);

	finalizar();
    MPI_Finalize();
    return(0);

}