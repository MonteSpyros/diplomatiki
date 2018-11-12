

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>

#define EMPTY -1

void error(const char * s)
{
	puts(s);
	exit(-1);
} 

typedef struct {
	double Dato;
	unsigned int p_ini;
} d_o;

int sort_function( const void *a, const void *b );

double serial_Lightweight(double *x, int N, int m, double r){
	int Nm = N-m;
    int i,j,k,ii,jj;
    int A=0,B=0;
	int r_sup, i_inf, i_sup, i_mez;
	double *X;
	double D_piu_r;
	d_o  *_D_ordinati;
	int *originalPositions;



////////////////////

    //value of r
    double sum=0.0;
    for(i=0; i<N; i++)
        sum += x[i];
    double mean = sum/N;        
    double standardDeviation=0.0;
    for (i=0; i<N; i++)
        standardDeviation += (x[i]-mean)*(x[i]-mean);
    standardDeviation = sqrt(standardDeviation/N);
    r=r*standardDeviation;
    


	
	X = (double *)malloc(Nm*sizeof(double));
	for (i=0;i<Nm;i++) 
	{
		X[i]=x[i];
		for (j=1;j<m;j++)
			X[i]+=x[i+j];
	}

	for(i=0; i<Nm; i++){
		printf("Simple: %lf\n" , X[i]);
	}
	
	_D_ordinati = (d_o *)malloc(Nm*sizeof(d_o));
	for(i=0; i<Nm; i++)
	{
		_D_ordinati[i].Dato=X[i];
		_D_ordinati[i].p_ini=i;
	}
	
	qsort((void *)_D_ordinati, (size_t)Nm, sizeof(d_o), sort_function);

	originalPositions = (int *)malloc(Nm*sizeof(int));
	
	for(i=0; i<Nm; i++)
	{
		X[i]=_D_ordinati[i].Dato;
		originalPositions[i]=_D_ordinati[i].p_ini;
	}

	free(_D_ordinati);
	
	for(i=0; i<Nm; i++)
	{
		D_piu_r=X[i]+m*r;
		if(D_piu_r >= X[Nm-1])
			r_sup=Nm-1;
		else 
		{
			i_inf=i;
			i_sup=Nm-1;
			while(i_sup-i_inf>1) 
			{
				i_mez=(i_inf+i_sup)>>1;
				if( X[i_mez] > D_piu_r )
					i_sup=i_mez;
				else
					i_inf=i_mez;
			}
			r_sup=i_inf;
		}
		ii=originalPositions[i];
		for(j=i+1; j<=r_sup; j++) 
		{
			jj=originalPositions[j];
			
			for (k=0;k<m;k++)
				if (fabs(x[ii+k]-x[jj+k])>r)
					break;
			if (k==m)
			{
				B++;
				//printf("found %d,%d\n",i,j);
				if (fabs(x[ii+m]-x[jj+m])<=r)
					A++;
			}
		}    
	}

//////////////////////////////////////////
// return Sampen
	//printf("A=%d B=%d \n",A,B);
    if (A*B==0)
         return log((N-m)*(N-m-1));
    else
        return -log10(1.0*A/B);    
}



double parallel_Lightweight(double *x, int N, int m, double r){
 
    int Nm = N-m;
    int i,j,k,ii,jj;
    int A=0,B=0;
	int r_sup, i_inf, i_sup, i_mez;
	double *X;
	double D_piu_r;
	d_o  *_D_ordinati;
	int *originalPositions;

    //value of r
    double sum=0.0;
    #pragma omp parallel for schedule(dynamic) num_threads(4) private(i) reduction(+: sum)
    for(i=0; i<N; i++)
        sum += x[i];


    double mean = sum/N;        
    double standardDeviation=0.0;


    #pragma omp parallel for schedule(dynamic) num_threads(4) private(i) reduction(+: standardDeviation)
    for (i=0; i<N; i++)
        standardDeviation += (x[i]-mean)*(x[i]-mean);


    standardDeviation = sqrt(standardDeviation/N);
    r=r*standardDeviation;
    

	X = (double *)malloc(Nm*sizeof(double));

	//#pragma omp parallel for schedule(static) num_threads(4) private(i,j) reduction(+: X[i])
	for (i=0;i<Nm;i++) 
	{
		X[i]=x[i];
		for (j=1;j<m;j++)
			X[i]+=x[i+j];
	}

	
	_D_ordinati = (d_o *)malloc(Nm*sizeof(d_o));
	for(i=0; i<Nm; i++)
	{
		_D_ordinati[i].Dato=X[i];
		_D_ordinati[i].p_ini=i;
	}
	
	qsort((void *)_D_ordinati, (size_t)Nm, sizeof(d_o), sort_function);

	originalPositions = (int *)malloc(Nm*sizeof(int));
	
	for(i=0; i<Nm; i++)
	{
		X[i]=_D_ordinati[i].Dato;
		originalPositions[i]=_D_ordinati[i].p_ini;
	}

	free(_D_ordinati);
	
	for(i=0; i<Nm; i++)
	{
		D_piu_r=X[i]+m*r;
		if(D_piu_r >= X[Nm-1])
			r_sup=Nm-1;
		else 
		{
			i_inf=i;
			i_sup=Nm-1;
			while(i_sup-i_inf>1) 
			{
				i_mez=(i_inf+i_sup)>>1;
				if( X[i_mez] > D_piu_r )
					i_sup=i_mez;
				else
					i_inf=i_mez;
			}
			r_sup=i_inf;
		}
		ii=originalPositions[i];
		for(j=i+1; j<=r_sup; j++) 
		{
			jj=originalPositions[j];
			
			for (k=0;k<m;k++)
				if (fabs(x[ii+k]-x[jj+k])>r)
					break;
			if (k==m)
			{
				B++;
				//printf("found %d,%d\n",i,j);
				if (fabs(x[ii+m]-x[jj+m])<=r)
					A++;
			}
		}    
	}

//////////////////////////////////////////
// return Sampen
	//printf("A=%d B=%d \n",A,B);
    if (A*B==0)
         return log((N-m)*(N-m-1));
    else
        return -log10(1.0*A/B);    
////////////////////////////////////////  

}


int sort_function( const void *a, const void *b)
{
	return ( ((d_o *)a)->Dato > ((d_o *)b)->Dato ) ? 1 : -1;
}



int main(){
	
	int N; printf("\nN: "); if(scanf("%d", &N)!=1) error("invalid input in scanf"); //N MEGALO
	int m; printf("m: "); if(scanf("%d", &m)!=1) error("invalid input in scanf"); // m = 2
	double r; printf("r: "); if(scanf("%lf", &r)!=1)  error("invalid input in scanf"); // r = 0.2
	puts("");

	int i,k;
	double x[10][64000];
	
////////////////////

    // file read
    FILE *fp;   // pointer to input file
	char s[100];
	for (k=2001;k<=2010;k++)
	{
	
		sprintf(s,"/Users/spirosmontesantos/Desktop/Διπλωματική/entropysource/nsr/nsr%d.rr",k);
		
		fp=fopen(s,"r");

		if (fp==NULL){
			error("error in opening input file");
		} 

		fscanf(fp, "%*[^\n]\n", NULL);

		for(i=0;i<N;i++) // to N einai megalo 
			if (fscanf(fp,"%lf",&x[k-2001][i])!=1){
				error("error in reading file");
			} 
			else{
				//printf("h timh pou phra einai %lf\n" , x[k-2001][i]);
			}
	}
	
	FILE * fpout;
	sprintf(s,"/Users/spirosmontesantos/Desktop/Διπλωματική/entropysource/m%d_r%d_N%d",m,(int)(r*100),N);
	fpout=fopen(s,"w");
	if (fpout==NULL) error("error in opening output file");
		

	
	//clock_t begin,end;
	//double time_spent = 0;
	
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
			result=serial_Lightweight(x[i],N,m,r);
			gettimeofday(&half, NULL);
			new_spent = half.tv_sec - start.tv_sec;
			new_spent = new_spent + (half.tv_usec - start.tv_usec)/1000000.0;
			elapsed_serial_time = elapsed_serial_time + new_spent;
		}
		
		printf("Serial calculation time is : %lf secs\n" , elapsed_serial_time);

		for(i=0;i<10;i++){
			result = parallel_Lightweight(x[i],N,m,r);
		}
				
	
	//fprintf(stderr,"lw: \t\t %5.2f msec\n",time_spent*2);
	//fprintf(fpout,"lw: \t\t %5.2f msec\n",time_spent*2);
				
				
	fclose(fp);
	fclose(fpout);
}
    
    
    
