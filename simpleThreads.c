#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <pthread.h>

#define NUM_THREADS 4
#define K 10 //DIASTHMATA

void *parallel_simple(void *arg);

pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t stand_lock = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t sumLock = PTHREAD_MUTEX_INITIALIZER;

double x[10][64000];
double global_sum = 0;
double temp_Global_sum = 0;
double global_result;
double standardDeviation = 0;
clock_t parallel_end;
struct thread_data args;
int t = 0;
int NTASK;
int A = 0;
int B = 0;
double test_Stand;


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

void calculateFirst(){
	parallel_simple(t);
}

void *thrfunc(void *arg){
	while(1){
		pthread_mutex_lock(&lock);
		t = t + NTASK;
		pthread_mutex_unlock(&lock);
		//printf("To taski einai: %d kai to NTASK einai %d\n" , t, NTASK);
		if(t >= args.N){
			break;
		}
		//printf("To taski einai: %d kai to NTASK einai %d\n" , t, NTASK);
		parallel_simple(t);
	}
}

void *parallel_simple(void *arg){

 	int Nm = (args.N)-(args.m);
    int i,j,k, me = (int *) arg;
   
    double sum = 0.0;

    for(i=me; i<(me+NTASK); i++){
        sum += x[0][i];
    }

    pthread_mutex_lock(&sumLock);
    global_sum = global_sum + sum;
    pthread_mutex_unlock(&sumLock);


    //edw vazw to barrier kai perimenoun ta alla nhmata! 
   	double mean = temp_Global_sum/args.N;        
    double temp_standardDeviation = 0.0;

    //printf("To mean einai: %lf\n" , mean);

    
   for (i=me; i<(me+NTASK); i++){
   		 pthread_mutex_lock(&stand_lock);
         standardDeviation += (x[0][i]-mean)*(x[0][i]-mean);
         pthread_mutex_unlock(&stand_lock);
    }

    temp_standardDeviation = sqrt(standardDeviation/args.N);
    //printf("to standardDeviation einai %lf\n" , standardDeviation);
	//deytero barrier  
 //     standardDeviation = sqrt(standardDeviation/(i-me));
 //     pthread_mutex_lock(&stand_lock);
 //     test_Stand = test_Stand + standardDeviation;
 //     pthread_mutex_unlock(&stand_lock);
 //     printf("to standardDeviation einai %lf\n" , standardDeviation);
 //     double temp_r = args.r * standardDeviation;
 //     int mem = me - args.m;
    

 //     //compute A,B
	// for (i=0;i<mem;i++)
	// {
	// 	for (j=i+1;j<mem;j++)
	// 	{
	// 		for (k=0;k<args.m;k++)
	// 			if (fabs(x[0][i+k]-x[0][j+k]) > temp_r)
	// 				break;
	// 		if (k == args.m)
	// 		{
	// 			pthread_mutex_lock(&lock);
	// 			B++;
	// 			pthread_mutex_unlock(&lock);
	// 			if (fabs(x[0][i+args.m]-x[0][j+args.m]) <= temp_r)
	// 				pthread_mutex_lock(&lock);
	// 				A++;
	// 				pthread_mutex_unlock(&lock);
	// 		}
	// 	}
	// }

// //return Sampen
// 	double exit_value = 0.0;
//     if (A*B==0){
//     	pthread_mutex_lock(&lock);
//     	exit_value = log(((data->N)-(data->m)*((data->N)-(data->m)-1)));
//     	pthread_mutex_unlock(&lock);
//     	//printf("To result ths diergasias %d einai %lf\n" , temp_Counter, exit_value);
//     	if(temp_Counter == NUM_THREADS-1){
//     		global_result = exit_value;
//     		parallel_end = clock();
//     	}
//     }
//     else{
// 		exit_value = -log10(1.0*A/B);  
// 		//printf("To result ths diergasias %d einai %lf\n" , temp_Counter, exit_value);
// 		if(temp_Counter == NUM_THREADS-1){
//     		global_result = exit_value;
//     		parallel_end = clock();
//     	}   
//     }
 

 	
} 

double serial_simple(double *x, int N, int m, double r) 
{
    int Nm = N-m;
    int i,j,k;

    double sum=0.0;

    for(i=0; i<N; i++){
        sum += x[i];
    }
 
    temp_Global_sum = sum;
    double mean = sum/N;        
    double standardDeviation=0.0;

    //printf("To serial mean einai %lf\n" , mean);

    for (i=0; i<N; i++){
        standardDeviation += (x[i]-mean)*(x[i]-mean);
    }

    standardDeviation = sqrt(standardDeviation/N);
   // printf("to standardDeviation SERIAL einai %lf\n" , standardDeviation);
    r=r*standardDeviation;
   
//compute A,B
	int As=0; int Bs=0;
	for (i=0;i<Nm;i++)
	{
		for (j=i+1;j<Nm;j++)
		{
			for (k=0;k<m;k++)
				if (fabs(x[i+k]-x[j+k])>r)
					break;
			if (k==m)
			{
				B++;
				if (fabs(x[i+m]-x[j+m])<=r)
					A++;
			}
		}
	}

	printf("To A kai to B einai: %d %d\n" , As, Bs);

// return Sampen
    if (A*B==0)
         return log((N-m)*(N-m-1));
    else
        return -log10(1.0*A/B);      
}   

int main(){

	int N; printf("\nN: "); if(scanf("%d", &N)!=1) error("invalid input in scanf");
	int m; printf("m: "); if(scanf("%d", &m)!=1) error("invalid input in scanf");
	double r; printf("r: "); if(scanf("%lf", &r)!=1)  error("invalid input in scanf");
	puts("");

	int i,k;
	pthread_t tid[NUM_THREADS];
	//struct thread_data args;

	args.N = N;
	args.m = m;
	args.r = r;
	NTASK = N/K;

	
	printf("data of struct %d %d %lf:\n " , args.N , args.m , args.r);
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


	clock_t begin,end;
	clock_t parallel_begin;
	double time_spent = 0;
	double parallel_time_spent = 0;
	
	double result;
	
	fprintf(stderr,"\nN=%d, m=%d, r=%.2f\n\n",N,m,r);
	fprintf(fpout,"\nN=%d, m=%d, r=%.2f\n\n",N,m,r);

	//for(k=0;k<50;k++){ //giati ginetai auto?!?
	begin = clock();
		for (i=0;i<1;i++)	
	 		result=serial_simple(x[i],N,m,r);
		end = clock();
		time_spent += (double)(end - begin) / CLOCKS_PER_SEC;
	//}

	fprintf(stderr,"simple: \t %5.2f msec\n",time_spent*2);
	fprintf(fpout,"simple: \t %5.2f msec\n",time_spent*2);

	//parallel 
	parallel_begin = clock();
	calculateFirst();
	for(i=0;i<4;i++){
		pthread_create(&(tid[i]),NULL,thrfunc,NULL); 
	}

	for(i = 0; i<NUM_THREADS; i++){
		pthread_join(tid[i], NULL);
	}

	printf("To sum einai: %lf\n" , global_sum);
	parallel_time_spent += (double)(parallel_end - parallel_begin) / CLOCKS_PER_SEC;




	printf("To result einai %lf\n" , result);
	printf("To parallel result einai %lf\n" , global_result);
	printf("To A kai to B einai: %d %d\n" , A, B);
	printf("to standardDeviation parall einai %lf\n" , test_Stand);
	// fprintf(stderr,"simple: \t %5.2f msec\n",time_spent*2);
	// fprintf(fpout,"simple: \t %5.2f msec\n",time_spent*2);
	fprintf(stderr,"parallel: \t %5.2f msec\n",parallel_time_spent*2);
	fprintf(fpout,"parallel: \t %5.2f msec\n",parallel_time_spent*2);

	fclose(fp);
	fclose(fpout);



	return 0;
}