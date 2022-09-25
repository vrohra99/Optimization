/*Problem_3 :
Following methods are used in this code:
Bounding phase Method
Bisection Method
Conjugate Gradient Method
Bracket operator penalty method*/

#include<stdio.h>
#include<conio.h>
#include<math.h>
#include<stdlib.h>

double original_function(double* x,int n);                                          /*Main objective function*/
double objective_function(double* x,int n,double R);								/*penalty function*/
double function(double alpha,double* xx,double* s,int n,double R);					/*Single variable function for unidirectional search*/
double derivative(double x,double* xx,double* s,int n,double R);						/*Derivative of single variable function for unidirectional search*/
double* main_derivative(double* xx,int n,double R);									/*To calculate -grad(f(x))*/
void bounding_phase(double x0,double d,double* xfinal,double* yfinal,double* xx,double* s,int n,double R);		/*function for bounding phase method*/
double bisection_method(double x1,double x2,double* xx,double* s,int n,double R);	/*function for bisection method*/
double* unidirectional_search(double* xx,double* s,int n,double R);					/*function for unidirectional search*/
double* fr_equation(double* x0,double* x1,double* s0,int n,double R);				/*function for Fletcher & Reeves equation*/
double norm(double* vector,int n);											/*To calculate norm of any vector*/
double condition_1(double* s0,double* s1,int n);							/*linear independency check between two search directions*/
double condition_2(double* x1,double* x2,int n);							/*To check ||x2-x1||/||x1|| <= e2*/													
double condition_3(double* x,int n,double R);										/*To check ||grad(f(x))|| <= e3*/
double* unconstrained(double* x0,int n,double R);                           /*Finding solution for penalty function*/
double min(double a,double b);                                              /*function to find minimum between two numbers*/

int main()
{
	FILE *finput,*foutput;						/*Input and Output File*/
	finput= fopen("input_problem_3.txt","r");
	foutput= fopen("output_problem_3.txt","w");
    int n;
	fscanf(finput,"%d",&n);         /*Taking number of variables*/
    double* x0;
    int i;
	x0=malloc(sizeof(double*)*n);
	for(i=0;i<n;i++)							/*Step-1 : Entering x0 value*/
	{
		fscanf(finput,"%lf",&x0[i]);
	}
	double R=0.1;       /*Penalty paramater*/
    double e1=1.0e-3;   /*Termination parameters*/
    int k=1;            /*Parameter for Counter*/
    double* x_old;  x_old=malloc(sizeof(double*)*n);
    double* x_new;  x_new=malloc(sizeof(double*)*n);
    double value;
    double c=10;        /*Paramater to update R*/
    printf("************************************Problem_3************************************\n\n");
    fprintf(foutput,"************************************Problem_3************************************\n\n");
    printf("Sequence\tR\t\tx1\t\tx2\t\tf(x)\n");
    fprintf(foutput,"Sequence\tR\t\tx1\t\tx2\t\tf(x)\n");
    x_old=unconstrained(x0,n,R);        /*Step-2 : Finding solution x for penalty function using conjugate gradient method*/
    /*Step-3 : Comparing penalised function values. As this is first Iteration we do not have previous solution*/
    printf("%d\t\t%0.1lf\t\t",k,R);
    fprintf(foutput,"%d\t\t%0.1lf\t\t",k,R);
    for(i=0;i<n;i++)
    {
        printf("%lf\t",x_old[i]);
        fprintf(foutput,"%lf\t",x_old[i]);
    }
    printf("%lf\n",original_function(x_old,n));
    fprintf(foutput,"%lf\n",original_function(x_old,n));
    R=c*R;                                      /*Step-4 : Update R*/
    do
    {
        k++;                                        /*Increse Counter by 1*/
        x_new=unconstrained(x_old,n,R);             /*Step-2 : Finding solution x for penalty function using conjugate gradient method*/
        value = objective_function(x_new,n,R) - objective_function(x_old,n,(R/c));  /*Step-3 : Comparing penalised function values*/
        x_old=x_new;
        printf("%d\t\t%0.1lf\t\t",k,R);
        fprintf(foutput,"%d\t\t%0.1lf\t\t",k,R);
        for(i=0;i<n;i++)
        {
            printf("%lf\t",x_old[i]);
            fprintf(foutput,"%lf\t",x_old[i]);
        }
        printf("%lf\n",original_function(x_old,n));
        fprintf(foutput,"%lf\n",original_function(x_old,n));
        R=c*R;                                      /*Step-4 : Update R*/
    }while(fabs(value)>e1);
    printf("\n\nThe solution from the last Sequence is given by\n\n");
    fprintf(foutput,"\n\nThe solution from the last Sequence is given by\n\n");
	printf("Sequence\tR\t\tx1\t\tx2\t\tf(x)\n");
	fprintf(foutput,"Sequence\tR\t\tx1\t\tx2\t\tf(x)\n");
    printf("%d\t\t%0.1lf\t\t",k,R/c);
    fprintf(foutput,"%d\t\t%0.1lf\t\t",k,R/c);
    for(i=0;i<n;i++)
    {
        printf("%lf\t",x_old[i]);
        fprintf(foutput,"%lf\t",x_old[i]);
    }
    printf("%lf\n",original_function(x_old,n));
    fprintf(foutput,"%lf\n",original_function(x_old,n));
	fclose(finput);     /*Close input file*/
	fclose(foutput);    /*Close output file*/
	return 0;
}


double* unconstrained(double* x0,int n,double R)
{
    double e1=1.0e-10,e2=1.0e-10,e3=1.0e-10;		/*Termination parameters*/
	int k=0;									/*Iteration number*/
	double* x1;     x1=malloc(sizeof(double*)*n);
	double* x2;     x2=malloc(sizeof(double*)*n);
	double* x3;     x3=malloc(sizeof(double*)*n);
	double* s0;     s0=malloc(sizeof(double*)*n);
	double* s1;     s1=malloc(sizeof(double*)*n);
	double* s2;     s2=malloc(sizeof(double*)*n);
	int i;

double error=0;	
	do
	{
	    
	    s0=main_derivative(x0,n,R);												/*Step-2 : Compute s0*/
		x1=unidirectional_search(x0,s0,n,R);										/*Step-3 : Compute x1*/
	    for(i=0;i<=1000;i++){error=0;
		s1=fr_equation(x0,x1,s0,n,R);												/*Step-4 : Compute s1*/
		if(condition_1(s0,s1,n)<=e1 && condition_1(s0,s1,n)>=-e1){break;	}		/*linear independency check between two search directions*/
		x2=unidirectional_search(x1,s1,n,R);										/*Step-5 : Compute x2*/	
																			/*Increase counter k by 1*/								
	
		if(condition_2(x1,x2,n)<=e2 || condition_3(x2,n,R)<=e3){break; 	}			/*Step-6 : Checking conditions for termination*/

		s0=s1;

		x0=x1;

		x1=x2;
	
		error=1;
	    }
	    k++;
	    															/*Printing x2*/
	    if(error==0){break;}
		
	}while(k<100);
	
    return x2;
    
}





double condition_1(double* s0,double* s1,int n)		/*linear independency check between two search directions*/
{
	double sum=0;
	int i;
	double unit1=norm(s0,n);
	double unit2=norm(s1,n);
	for(i=0;i<n;i++)
	{
		sum = sum + (s0[i]*s1[i]/(unit1*unit2));
	}
	return (acos(sum));
}


double condition_2(double* x1,double* x2,int n)	/*To check ||x2-x1||/||x1|| <= e2*/			
{
	double* d;
	d=malloc(sizeof(double*)*n);
	int i;
	for(i=0;i<n;i++)
	{
		d[i] = x2[i] - x1[i];
	}
	double normm=norm(d,n)/norm(x1,n);
	return(normm);
}


double condition_3(double* x,int n,double R)		/*To check ||grad(f(x))|| <= e3*/
{
	double* main_derive=main_derivative(x,n,R);
	double normm=norm(main_derive,n);
	return(normm);
}


double condition_4(double* x1,double* x2,int n,double R)
{
    double ans;
    ans = objective_function(x2,n,R) - objective_function(x1,n,R);
    return (fabs(ans));
}

double norm(double* vector,int n)	/*To calculate norm of any vector*/
{
	int i;
	double sum=0;
	for(i=0;i<n;i++)
	{
		sum = sum + pow(vector[i],2);
	}
	return (sqrt(sum));
}




double* fr_equation(double* x0,double* x1,double* s0,int n,double R)			/*function for Fletcher & Reeves equation*/
{
    int i;
    double norm_s0=norm(s0,n);
	for(i=0;i<n;i++)
	{
	    s0[i] = s0[i]/norm_s0;
	}
	double* s1;
	s1=malloc(sizeof(double*)*n);
	double* d0; d0=malloc(sizeof(double*)*n);
	double* d1; d1=malloc(sizeof(double*)*n);
	
	d0=main_derivative(x0,n,R);
	d1=main_derivative(x1,n,R);
	
	for(i=0;i<n;i++)
	{
		s1[i] = d1[i] + (pow(norm(d1,n)/norm(d0,n),2)*s0[i]);
	}
	double norm_s1=norm(s1,n);
	for(i=0;i<n;i++)
	{
	    s1[i] = s1[i]/norm_s1;
	}
	
	return s1;
}


double* unidirectional_search(double* xx,double* s,int n,double R)			/*function for unidirectional search*/
{
    int i;
    double norm_s=norm(s,n);
	for(i=0;i<n;i++)
	{
	    s[i] = s[i]/norm_s;
	}
	double x0=1.0e-2,d=1.0e-5;		/*initial guess 'x0' and increment'd' */
	double x,y;					/*range (x,y) after bounding phase method*/
	bounding_phase(x0,d,&x,&y,xx,s,n,R);
	double alpha;
	
	alpha=bisection_method(x,y,xx,s,n,R);
	
	double* x_new;
	x_new=malloc(sizeof(double*)*n);
	
	for(i=0;i<n;i++)
	{
		x_new[i] = xx[i] + ((alpha)*(s[i]));
	}
	return x_new;
}


double bisection_method(double x1,double x2,double* xx,double* s,int n,double R)		/*function for bisection method*/
{
	double z,e=1.0e-3;	/*accuray epsilon 'e' */
	int i=1;
	if(derivative(x1,xx,s,n,R)>0)
	{
		double temp;
		temp=x1;
		x1=x2;
		x2=temp;	
	}
	z=(x1+x2)/2.00;
	while(fabs(derivative(z,xx,s,n,R))>e || i<1000)
	{i++;
		if(derivative(z,xx,s,n,R)>0)
		{
			x2=z;	
		}
		else
		{
			x1=z;
		}
		double k=z;
		z=(x1+x2)/2.00;
		if(fabs(k-z)<1.0e-14){break;}
	}
	//if(fabs(z)<1.0e-13){return 1.0e-14;}
	return z;
}


void bounding_phase(double x0,double d,double* xfinal,double* yfinal,double* xx,double* s,int n,double R)	/*function for bounding phase method*/
{
	double x[1000];
	x[0]=x0;
	if(function(x[0]-d,xx,s,n,R)>=function(x[0],xx,s,n,R)  &&  function(x[0],xx,s,n,R)>=function(x[0]+d,xx,s,n,R))
	{
		d=+d;
	}
	else
	{
		d=-d;
	}
	int k=0;
	x[k+1]=x[k]+(pow(2,k)*d);
	while(function(x[k+1],xx,s,n,R)<function(x[k],xx,s,n,R) )
	{
	    if(k>990){break;}
		k++;
		x[k+1]=x[k]+(pow(2,k)*d);
	}
	if(x[k+1]>x[k-1])
	{
		*xfinal=x[k-1];
		*yfinal=x[k+1];
	}
	else
	{
		*xfinal=x[k+1];
		*yfinal=x[k-1];
	}
}


double* main_derivative(double* xx,int n,double R)				/*To calculate -grad(f(x))*/
{
    double delta = 1.0e-6;
    double* derive_f;
    derive_f=malloc(sizeof(double*)*n);
    int i;
    for(i=0;i<n;i++)
    {
    	xx[i]=xx[i]-delta;
    	double y1=objective_function(xx,n,R);
    	xx[i]=xx[i]+(2*delta);
    	double y2=objective_function(xx,n,R);
    	derive_f[i]=-(y2-y1)/(2*delta);
    	xx[i]=xx[i]-(delta);
	}
    return derive_f;
}


double derivative(double x,double* xx,double* s,int n,double R)			/*Derivative of single variable function for unidirectional search*/
{
    double delta = 1.0e-6;
    return ((function(x+delta,xx,s,n,R)-function(x-delta,xx,s,n,R))/(2*delta));
}


double function(double alpha,double* xx,double* s,int n,double R)		/*Single variable function for unidirectional search*/
{
	int i;
	double pp[n];
	for(i=0;i<n;i++)
	{
		pp[i]=xx[i] + (alpha*s[i]);
	}
	return(objective_function(pp,n,R));
}


double original_function(double* x,int n)   /*main objective function*/
{
	double function;
    function = x[0] + x[1] + x[2] ;
	return function;
}

double objective_function(double* x, int n,double R)    /*penalty function*/
{
    double penalty_function;
    
    double g[7];
    g[1] = min(0,1 - (0.0025*(x[3]+x[5]))) ;
    g[2] = min(0,1 - (0.0025*(-x[3]+x[4]+x[6])) );
    g[3] = min(0,1 - (0.01*(-x[5]+x[7])));
    g[4] = min(0,-(100*x[0]) + (x[0]*x[5]) - (833.33252*x[3]) + 83333.333);
    g[5] = min(0,-(x[1]*x[3]) + (x[1]*x[6]) + (1250*x[3]) - (1250*x[4]));
    g[6] = min(0,-(x[2]*x[4]) + (x[2]*x[7]) + (2500*x[4]) - 1250000); 
    int i,k=11;
    for(i=3;i<=7;i++)
    {
        g[k] = min(0,x[i]-10);
        k++;
        g[k] = min(0,1000-x[i]);
        k++;
    }
    double sum=0;
    for(i=1;i<=6;i++)
    {
        sum += (R*pow(g[i],2));
    }
    penalty_function = original_function(x,n) +  sum ;
    return penalty_function;
}

double min(double a,double b)   /*function to find minimum between two numbers*/
{
    if(a>b)
    {
        return b;
    }
    else
    {
        return a;
    }
}


