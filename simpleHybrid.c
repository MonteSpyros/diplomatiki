#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>
#include <pthread.h>
#include <omp.h>

#define NUM_THREADS 4

pthread_mutex_t Alock = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t Block = PTHREAD_MUTEX_INITIALIZER;

double x[10][64000];
int Ap = 0;
int Bp = 0;
double parallel_sum = 0.0;
double parallel_standardDeviation = 0.0;
double global_result = 0.0;
double previewsR = 0.0;
struct thread_data args;
int part = 0;

struct thread_data {
 int N;
 int m;
 double r;
}; 

void error(const char * s)
{
	puts(s);
	exit(-1);
}

void initPart(){
	part = 0;
}
void initSum(){
	parallel_sum = 0.0;
}

void initStandardDeviation(){
	parallel_standardDeviation = 0.0;
}

void initResult(){
	global_result = 0.0;
}

void initR(){
	args.r = previewsR;
}

void initAB(){
	Ap = 0;
	Bp = 0;
}

double parallel_result(){
	 if (Ap*Bp==0)
         return log((args.N-args.m)*(args.N-args.m-1));
    else
        return -log10(1.0*Ap/Bp);
}


void parallel_Computations(double *x, int N, int m, double r){
		int Nm = N-m;
	    int i,j,k;
	    initSum();
	    initStandardDeviation();

	    #pragma omp parallel for schedule(dynamic) num_threads(4) private(i) reduction(+: parallel_sum)
	    for(i=0; i<N; i++){
	        parallel_sum += x[i];
	    }

	    double mean = parallel_sum/N;        
	  
	    #pragma omp parallel for schedule(dynamic) num_threads(4) private(i) reduction(+: parallel_standardDeviation)
	    for (i=0; i<N; i++){
	        parallel_standardDeviation += (x[i]-mean)*(x[i]-mean);
	    }

	    parallel_standardDeviation = sqrt(parallel_standardDeviation/N);
	    args.r=args.r*parallel_standardDeviation;	
} 

void *thrfunc(void *arg){
	int i,j,k;
	int line = (int) arg;
	int thread_part = part++;
	int slash = 0.0;
	int A = 0;
	int B = 0;
	int counter1=0;
	int counter2=0;

	if(thread_part == 3){
		slash = (args.N/NUM_THREADS) + (args.N%NUM_THREADS);
		for(i = thread_part * (args.N/NUM_THREADS) ; i < (thread_part + 1) * (args.N / NUM_THREADS) + (args.N%NUM_THREADS); i++){
			for(j=i+1; j<args.N-args.m; j++){
				for(k=0; k<args.m;k++){
					if (fabs(x[line][i+k]-x[line][j+k])>args.r){
						break;
					}
				}
				if (k==args.m){
					B++;
					if (fabs(x[line][i+args.m]-x[line][j+args.m])<=args.r){
						A++;
					}	
				}
			}
		}
	}
	else{
		slash = (args.N/NUM_THREADS);
 		for(i = thread_part * ((args.N-args.m)/NUM_THREADS) ; i < (thread_part + 1) * ((args.N-args.m) / NUM_THREADS) ; i++){
			for(j=i+1; j<args.N-args.m; j++){
				for(k=0; k<args.m; k++){
					if (fabs(x[line][i+k]-x[line][j+k])>args.r){
						
						break;
					}
				}
				if(k == args.m){
					B++;
					if (fabs(x[line][i+args.m]-x[line][j+args.m])<=args.r){
						A++;
					}
				}
			}
		}
	}
	pthread_mutex_lock(&Block);
	Ap = Ap + A;
	Bp = Bp + B;
 	pthread_mutex_unlock(&Block);
}

double parallelResult(){
	// return Sampen
    if (Ap*Bp==0){
    	printf("IF parallel\n");
    	return log((args.N-args.m)*(args.N-args.m-1));
    }
    else
        return -log10(1.0*Ap/Bp); 
}


double serial_simple(double *x, int N, int m, double r) 
{
    int Nm = N-m;
    int i,j,k;

    double sum=0.0;

    for(i=0; i<N; i++){
        sum += x[i];
    }

    double mean = sum/N;        
    double standardDeviation=0.0;

    for (i=0; i<N; i++){
        standardDeviation += (x[i]-mean)*(x[i]-mean);
    }

    standardDeviation = sqrt(standardDeviation/N);
    r=r*standardDeviation;
   
//compute A,B
    int c = 0;
	int A=0; int B=0;
	for (i=0;i<Nm;i++)
	{
		for (j=i+1;j<Nm;j++)
		{
			for (k=0;k<m;k++)
				if (fabs(x[i+k]-x[j+k])>r){
					
					break;
				}
			if (k==m)
			{
				B++;
				if (fabs(x[i+m]-x[j+m])<=r)
					A++;
			}
		}
	}
printf("TO A %d kai B %d \n", A,B);
// return Sampen
    if (A*B==0){
    	printf("IF\n");
         return log((N-m)*(N-m-1));
    }
    else
        return -log10(1.0*A/B);      
}   

int main(){

	int N; printf("\nN: "); if(scanf("%d", &N)!=1) error("invalid input in scanf");
	int m; printf("m: "); if(scanf("%d", &m)!=1) error("invalid input in scanf");
	double r; printf("r: "); if(scanf("%lf", &r)!=1)  error("invalid input in scanf");
	puts("");
	int i,k;
	int thrCounter;
	pthread_t tid[NUM_THREADS];

	args.N = N;
	args.m = m;
	args.r = r;
	previewsR = r;
	
	// file read
    FILE *fp;
	char s[100];
	for (k=2001;k<=2010;k++)
	{
		sprintf(s,"/Users/spirosmontesantos/Desktop/Διπλωματική/entropysource/nsr/nsr%d.rr",k);
		
		fp=fopen(s,"r");

		if (fp==NULL){
			error("error in opening input file");
		} 
		
		//skip the first line
		fscanf(fp, "%*[^\n]\n", NULL);

		for(i=0;i<N;i++){  //pare ta N prwta shmata kai valta ston pinaka x
			if (fscanf(fp,"%lf",&x[k-2001][i])!=1){ 
				error("error in reading file");
			} 
			else{
				//printf("h timh pou phra einai %lf\n" , x[k-2001][i]);
			}
		}
		
	}

	FILE * fpout;
	sprintf(s,"/Users/spirosmontesantos/Desktop/Διπλωματική/entropysource/m%d_r%d_N%d",m,(int)(r*100),N);
	fpout=fopen(s,"w");
	if (fpout==NULL) error("error in opening output file");


	double elapsed_serial_time = 0;
	double elapsed_parallel_time = 0;
	double new_spent = 0;
	double p_spent = 0;
	struct timeval start,half, parallel_beg, end;
	double result;
	double parallel_result;

	fprintf(stderr,"\nN=%d, m=%d, r=%.2f\n\n",N,m,r);
	fprintf(fpout,"\nN=%d, m=%d, r=%.2f\n\n",N,m,r);

	//SERIAL

	for (i=0;i<1;i++){
		gettimeofday(&start, NULL);	
	 	result=serial_simple(x[i],N,m,r);
		gettimeofday(&half, NULL);
		new_spent = half.tv_sec - start.tv_sec;
		new_spent = new_spent + (half.tv_usec - start.tv_usec)/1000000.0;
		elapsed_serial_time = elapsed_serial_time + new_spent;
	}

		printf("Serial calculation time is : %lf secs\n" , elapsed_serial_time);
		printf("Serial result is %lf\n" , result);
	// fprintf(stderr,"simple: \t %5.2f msec\n",time_spent*2);
	// fprintf(fpout,"simple: \t %5.2f msec\n",time_spent*2);


	//PARALLEL
	
	for(i=0;i<1;i++){
		gettimeofday(&parallel_beg, NULL);

		parallel_Computations(x[i],N,m,r);

		for(thrCounter=0; thrCounter<NUM_THREADS; thrCounter++){
			pthread_create(&(tid[thrCounter]),NULL,thrfunc,(void*)(i));
		}
		
		for(thrCounter = 0; thrCounter<NUM_THREADS; thrCounter++){
		    pthread_join(tid[thrCounter], NULL);
		}
		printf("Ap %d kai Bp %d\n", Ap, Bp);
		gettimeofday(&end, NULL);
		global_result = parallelResult();
		printf("To result einai %lf\n" , global_result);
		initPart();
		initR();
		initAB();
		initResult();
		p_spent = end.tv_sec - parallel_beg.tv_sec;
		p_spent = p_spent + (end.tv_usec - parallel_beg.tv_usec)/1000000.0;
		elapsed_parallel_time = elapsed_parallel_time + p_spent;
	}

	printf(" Parallel calculation time is : %lf secs\n" , elapsed_parallel_time);
	// fprintf(stderr,"simple: \t %5.2f msec\n",time_spent*2);
	// fprintf(fpout,"simple: \t %5.2f msec\n",time_spent*2);
	// fprintf(stderr,"parallel: \t %5.2f msec\n",parallel_time_spent*2);
	// fprintf(fpout,"parallel: \t %5.2f msec\n",parallel_time_spent*2);

	fclose(fp);
	fclose(fpout);



	return 0;
}