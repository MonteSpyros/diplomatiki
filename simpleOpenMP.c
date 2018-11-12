#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>
#include <omp.h>

void error(const char * s){
	puts(s);
	exit(-1);
}

double parallel_simple(double *x, int N, int m, double r){
		int Nm = N-m;
	    int i,j,k;

	    double sum=0.0;

	    #pragma omp parallel for schedule(dynamic) num_threads(4) private(i) reduction(+: sum)
	    for(i=0; i<N; i++){
	        sum += x[i];
	    }

	    double mean = sum/N;        
	    double standardDeviation=0.0;

	    #pragma omp parallel for schedule(dynamic) num_threads(4) private(i) reduction(+: standardDeviation)
	    for (i=0; i<N; i++){
	        standardDeviation += (x[i]-mean)*(x[i]-mean);
	    }

	    standardDeviation = sqrt(standardDeviation/N);
	    r=r*standardDeviation;

	    //barriers sta 2 panw kai dynamika alla grammh sta threads! 
	    int A=0; int B=0;
	#pragma omp parallel for schedule(dynamic) num_threads(4) private(i,j,k) reduction(+: A) reduction(+: B)
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

// return Sampen
    if (A*B==0)
         return log((N-m)*(N-m-1));
    else
        return -log10(1.0*A/B);     
   	
} 

double serial_simple(double *x, int N, int m, double r) {
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
	int A=0; int B=0;
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
printf("TO a einai %d\n", A);
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
	double x[10][64000];
	int i,k;
	
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

	for (i=0;i<10;i++){
		gettimeofday(&start, NULL);	
	 	result=serial_simple(x[i],N,m,r);
		gettimeofday(&half, NULL);
		new_spent = half.tv_sec - start.tv_sec;
		new_spent = new_spent + (half.tv_usec - start.tv_usec)/1000000.0;
		elapsed_serial_time = elapsed_serial_time + new_spent;
	}

		printf("Serial calculation time is : %lf secs\n" , elapsed_serial_time);
	// fprintf(stderr,"simple: \t %5.2f msec\n",time_spent*2);
	// fprintf(fpout,"simple: \t %5.2f msec\n",time_spent*2);

	//parallel 
	
	for(i=0;i<10;i++){
		gettimeofday(&parallel_beg, NULL);
		parallel_result = parallel_simple(x[i],N,m,r);
		gettimeofday(&end, NULL);
		p_spent = end.tv_sec - parallel_beg.tv_sec;
		p_spent = p_spent + (end.tv_usec - parallel_beg.tv_usec)/1000000.0;
		elapsed_parallel_time = elapsed_parallel_time + p_spent;
	}

	printf(" Parallel calculation time is : %lf secs\n" , elapsed_parallel_time);
	printf("To result einai %lf\n" , result);
	printf("To parallel result einai %lf\n" , parallel_result);
	// fprintf(stderr,"simple: \t %5.2f msec\n",time_spent*2);
	// fprintf(fpout,"simple: \t %5.2f msec\n",time_spent*2);
	// fprintf(stderr,"parallel: \t %5.2f msec\n",parallel_time_spent*2);
	// fprintf(fpout,"parallel: \t %5.2f msec\n",parallel_time_spent*2);

	fclose(fp);
	fclose(fpout);



	return 0;
}