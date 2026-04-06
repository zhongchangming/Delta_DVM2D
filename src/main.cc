
/* ****************upwind style**********************/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <omp.h>
#include <time.h>
//#define NUM_THREADS 4
#define L 1.0
#define M 21
#define M1 (M+1)
#define N 21
#define N1 (N+1)
#define Qx 81
#define Qy 81
#define Kn (10.0)
#define u_wall 0.16
#define sgn(x) (x>=0.0?1.0:-1.0)

double ex[Qx],ey[Qy],tpx[Qx], tpy[Qy];        //discrete velocity and weight
double Co_X[Qy],Co_W[Qy];       
double rho[M1][N1],u[M1][N1],v[M1][N1],T[M1][N1],T1[M1][N1],u1[M1][N1],v1[M1][N1]; //density x_velocity y_velocity
double f[M1][N1],g[M1][N1]; 
double qflux_x[M1][N1],qflux_y[M1][N1],stress[M1][N1];  
double fh[M1][N1],h[M1][N1];  

double xf[M1][N1],xfh[M1][N1]; 

double yf[M1][N1],yfh[M1][N1]; 
 
double rhoold[M1][N1];

double df[M1][N1],dfh[M1][N1]; 
double mu[M1][N1],mv[M1][N1],rho1[M1][N1],E1[M1][N1];  //mass momentum energy
double g1[Qx][Qy],g2[Qx][Qy];
double h1[Qx][Qy],h2[Qx][Qy];
double slg[M1][N1],slh[M1][N1],slgg[M1][N1],slhh[M1][N1]; 
void initial();
void boundary();
void slope();
void update();
void outputdata();
void Error();
void Cotes(double Vmin, double Vmax, int np);
double feq(int kx, int ky, double Ux, double Uy, double RHO,double TT);
double R,PI,dx,dy,tau,omega,gamma1,error,error1;
int i,j,kx,ky,inner,num,sumrun;
double err,err1,mu_ref,deltaf,deltaf1,df1rhs;
double tt[1000],rr[1000],rhs[1000],bushu[1000];
void Cotes(double Vmin, double Vmax, int np)
{
  int k,n,n1;
  double h, b, a;

  a=Vmin; b=Vmax;
  n=(np-1)/4;
  h=(b-a)/n;
  n1=np-1;

  for(k=0;k<np;k++) Co_X[k]=(k-0.5*n1)/n1*(b-a);

  for(k=0;k<n;k++)
  {
    Co_W[4*k]=14.0;
    Co_W[4*k+1]=32.0;
    Co_W[4*k+2]=12.0;
    Co_W[4*k+3]=32.0;
  }
  Co_W[0]=7.0;
  Co_W[np-1]=7.0;

  for(k=0;k<np;k++) Co_W[k]*=h/90;
}
void gauss_hermite();
void gauss_hermite()
{


  int k;
  //half-range Gauss-Hermit quadrature
  double GH_v[14]=  {
	 2.396724757843284e-002, 1.242911568927177e-001,  2.974477691473801e-001,  5.334998891822554e-001,   8.220981923307155e-001,
    1.154334829676328e+000, 1.523573280902749e+000,  1.925656879418240e+000,  2.358931078010347e+000,   2.824429077761555e+000,
    3.326604606373028e+000, 3.875436187413726e+000,  4.492761858309454e+000,  5.238744144080820e+000 };

  double GH_w[14]={
	 6.123724869886580e-002, 1.361115667521515e-001,  1.889014205320117e-001,  1.985861257727940e-001,  1.585835386099182e-001,
    9.276908911579565e-002, 3.789962589474001e-002,  1.024419454744042e-002,  1.720120272586066e-003,  1.656430867082354e-004,
    8.177599975023724e-006, 1.734291294554730e-007,  1.139604105901161e-009,  1.037777563554779e-012 };

  //Vmax=GH_v[13];

  for(k=0;k<14;k++)
  {
    ex[k]=-GH_v[13-k]; ex[14+k]=GH_v[k];
    tpx[k]=GH_w[13-k]; tpx[14+k]=GH_w[k];
  }

  for(k=0;k<Qy;k++) {ey[k]=ex[k]; tpy[k]=tpx[k];}
  for(k=0;k<Qy;k++) {tpy[k]=tpx[k]=tpx[k]*exp(ex[k]*ex[k]);} // Requirement: 2*R*T_ref=1.0


}
void Error()
{
	
	error=0;
	for(i=1;i<M;i++)
	{
		for(j=1;j<N;j++)
		{
 			  error=error+abs((T1[j][i]-T[j][i]) *(T1[j][i]-T[j][i]))+abs((u1[j][i]-u[j][i]) *(u1[j][i]-u[j][i]))+
                    abs((v1[j][i]-v[j][i]) *(v1[j][i]-v[j][i]))+abs((rhoold[j][i]-rho[j][i]) *(rhoold[j][i]-rho[j][i])); 
/* 			
 			if(error<fabs(T1[j][i]-T[j][i]))
		    {
		        error=fabs(T1[j][i]-T[j][i]);

		    }  */

		}
	}
	   error = sqrt(error/((M-1)*(N-1))); 


}


void slope()
{
	double s1,s2,s3,s4,x1,x2,x3,x4,q1,q2,q3,q4,ess;
	ess=1.0e-10;

	x1=0.5*(1+sgn(ex[kx]));
	x2=0.5*(1-sgn(ex[kx]));
	x3=0.5*(1+sgn(ey[ky]));
	x4=0.5*(1-sgn(ey[ky]));

	for(j=1;j<N;j++)
	{
		i=1;
		slg[i][j]=(f[i+1][j]-f[i][j])/dx;
		slh[i][j]=(fh[i+1][j]-fh[i][j])/dx;

		i=M-1;
		slg[i][j]=(f[i][j]-f[i-1][j])/dx;
		slh[i][j]=(fh[i][j]-fh[i-1][j])/dx;


		for(i=2;i<(M-1);i++)
		{
			s1=(f[i][j]-f[i-1][j])/dx;
			s2=(f[i+1][j]-f[i][j])/dx;
			slg[i][j]=(sgn(s1)+sgn(s2))*fabs(s1)*fabs(s2)/(fabs(s1)+fabs(s2)+1.0e-10);
			s1=(fh[i][j]-fh[i-1][j])/dx;
			s2=(fh[i+1][j]-fh[i][j])/dx;
			slh[i][j]=(sgn(s1)+sgn(s2))*fabs(s1)*fabs(s2)/(fabs(s1)+fabs(s2)+1.0e-10);

		}
		slg[0][j]=0.0;
		slg[M][j]=0.0;
		slh[0][j]=0.0;
		slh[M][j]=0.0;


		
	}
	for(i=1;i<M;i++)
	{

		j=1;
		slgg[i][j]=(f[i][j+1]-f[i][j])/dy;
		slhh[i][j]=(fh[i][j+1]-fh[i][j])/dy;

		j=N-1;
		slgg[i][j]=(f[i][j]-f[i][j-1])/dy;
		slhh[i][j]=(fh[i][j]-fh[i][j-1])/dy;

		for(j=2;j<(N-1);j++)
		{
			s1=(f[i][j+1]-f[i][j])/dy;
			s2=(f[i][j]-f[i][j-1])/dy;
			slgg[i][j]=(sgn(s1)+sgn(s2))*fabs(s1)*fabs(s2)/(fabs(s1)+fabs(s2)+1.0e-10);
			s1=(fh[i][j+1]-fh[i][j])/dy;
			s2=(fh[i][j]-fh[i][j-1])/dy;
			slhh[i][j]=(sgn(s1)+sgn(s2))*fabs(s1)*fabs(s2)/(fabs(s1)+fabs(s2)+1.0e-10);
		}
		
		slgg[i][0]=0.0;
		slgg[i][N]=0.0;
		slhh[i][0]=0.0;
		slhh[i][N]=0.0;

		
	}


	for(j=1;j<N;j++)
	{
		for(i=0;i<M;i++)
		{

			xf[i][j]=x1*(f[i][j]+dx/2*slg[i][j])+x2*(f[i+1][j]-dx/2*slg[i+1][j]);
			xfh[i][j]=x1*(fh[i][j]+dx/2*slh[i][j])+x2*(fh[i+1][j]-dx/2*slh[i+1][j]);

		}


	}


	for(i=1;i<M;i++)
	{
		for(j=0;j<N;j++)
		{

			yf[i][j]=x3*(f[i][j]+dy/2*slgg[i][j])+x4*(f[i][j+1]-dy/2*slgg[i][j+1]);
			yfh[i][j]=x3*(fh[i][j]+dy/2*slhh[i][j])+x4*(fh[i][j+1]-dy/2*slhh[i][j+1]);


		}

	}

}


double feq(int kx, int ky, double Ux, double Uy, double RHO,double TT)
{
	double x,eu;
	eu=(ex[kx]-Ux)*(ex[kx]-Ux)+(ey[ky]-Uy)*(ey[ky]-Uy);
	x=tpx[kx]*tpy[ky]*RHO/(2*PI*R*TT)*exp(-eu/(2.0*R*TT));
	return x;

}


int main()
{
  clock_t t1,t2;
  double factor,u_old,rho_old;
   double omega_ref; //reference parametere in HS model
  double alpha_ref; //reference coefficient in HS model
  int m,readdata,mmax,ss;
  //omp_set_num_threads(NUM_THREADS);
  R=0.5;
  num=1;
  df1rhs=1.0e-8;
  gamma1=5.0/3.0;
  PI=4.0*atan(1.0);

//  gauss_hermite();
  
  Cotes(-4.0,4.0,Qy);
  factor=0.0;
  for(m=0;m<Qx;m++)
  {
    ex[m]=ey[m]=Co_X[m];
    //tpx[m]=tpy[m]=Co_W[m]*exp(-ex[m]*ex[m]/(2*R))/sqrt(2*PI*R); 
    tpx[m]=tpy[m]=Co_W[m];
    //factor+=Co_W[m];
   

  } 


  


  dx=L/(M-1);
  dy=L/(N-1);

  
  initial();
  mu_ref=Kn*L/sqrt(0.5*PI*R)*R;

  omega=0.81;
  omega_ref=0.5;
  alpha_ref=1.0;
  mu_ref=5*(alpha_ref+1)*(alpha_ref+2)*sqrt(PI)/(4*alpha_ref*(5-2*omega_ref)*(7-2*omega_ref))*Kn;
  
  u_old=u[M/2][N/2];
  rho_old=rho[M/2][N/2];
  m=0;ss=0;
  mmax=99999;
AA:
  err=1.0; 
  err1=1.0;
  error=1.0;
  
  //printf("input mmax:\n");
  //scanf("%d",&mmax);
  mmax+=m;
  printf(" mmax=%d   \n", mmax);
  t1=clock();
  while((m<mmax) && error>3.0e-8)
	 {
	   m++;
	   update();
	   boundary();
	   
	   Error();
	   printf("deltaf  error =%e\t%e \n",deltaf,error);
	   if(m%5==0)
	   {
/* 	   	// Error();
	     err=fabs(u[M/2][N/2]-u_old)/u_wall;
	     err1=fabs(rho[M/2][N/2]-rho_old);
	     u_old=u[M/2][N/2];
	     rho_old=rho[M/2][N/2];
		 t2=clock();
		 rr[ss]=err;
		 rhs[ss]=error;
         tt[ss]=(double)(t2-t1)/CLOCKS_PER_SEC/60;
	     printf("err=%e  u=%e  m=%d  t=%f mins\n", err, u_old, m, (double)(t2-t1)/CLOCKS_PER_SEC/60);
	     printf("err1 error error1 =%e\t%e\t%e\n",err1,error,error1);
	     outputdata();
	     printf("rho u v T=%e\t%e\t%e\t%e\n",rho[M/2][N/2],u[M/2][N/2],v[M/2][N/2],T[M/2][N/2]);
	     sumrun=m;
         bushu[ss]=sumrun;
	     ss++; */

	   }
	 }
/* 	Error(); */
	t2=clock();
	/* rr[ss]=err;
	sumrun=m;
    bushu[ss]=sumrun;
    rhs[ss]=error;
    tt[ss]=(double)(t2-t1)/CLOCKS_PER_SEC/60;
    printf("err1 error error1 =%e\t%e\t%e\n",err1,error,error1);
    printf("err=%e  u=%e  m=%d  t=%f mins\n", err, u_old, m, (double)(t2-t1)/CLOCKS_PER_SEC/60); */
   
    outputdata();
/* 
   printf("rho u v T=%e\t%e\t%e\t%e\n",rho[M/2][N/2],u[M/2][N/2],v[M/2][N/2],T[M/2][N/2]);
 
  printf("Continue? (yes=1 no=0)\n");
  scanf("%d",&readdata);
  if(readdata) goto AA; */
  return 0;
}


void initial()
{
	for(i=0;i<=M;i++)
	{
		for(j=0;j<=N;j++)
		{
			u[i][j]=v[i][j]=0.0;
			rho[i][j]=1.0;
			T[i][j]=1.0;
			df[i][j]=0;
			dfh[i][j]=0;
		}
		u[i][0]=u_wall;
	}

	for(kx=0;kx<Qx;kx++)
	{
		for(ky=0;ky<Qy;ky++)
		{
			g1[kx][ky]=feq(kx,ky,0,0,1.0,1.0);
			g2[kx][ky]=feq(kx,ky,u_wall,0,1.0,1.0);
			h2[kx][ky]=R*g2[kx][ky];
			h1[kx][ky]=R*g1[kx][ky];   //h=R*T*g;
		}
	}
	
}
void update()
{
	double aa,bb,cc,dd,s1,s2;

	deltaf=0;
	
	for(i=0;i<=M;i++)
	{
		for(j=0;j<=N;j++)
		{
			mu[i][j]=0;
			rho1[i][j]=0;
			mv[i][j]=0;
			E1[i][j]=0;
			qflux_y[i][j]=0;
			qflux_x[i][j]=0;
			stress[i][j]=0;
		}
	}

	for(i=1;i<M;i++)
	{
		for(j=1;j<N;j++)
		{
			T1[i][j]=T[i][j];
			u1[i][j]=u[i][j];
			v1[i][j]=v[i][j];
      rhoold[i][j]=rho[i][j];
		}
	}
	//#pragma omp parallel for 
	for(kx=0;kx<Qx;kx++)
	{
		for(ky=0;ky<Qy;ky++)
		{
			aa=ex[kx]/dx;
			bb=ey[ky]/dy;
			
			
		
			for(i=0;i<=M;i++)
			{
				for(j=0;j<=N;j++)
				{
					g[i][j]=feq(kx,ky,u[i][j],v[i][j],rho[i][j],T[i][j]);
					h[i][j]=R*T[i][j]*g[i][j];
					f[i][j]=g[i][j];
					fh[i][j]=h[i][j];
				}
			}
			slope();
			
			if(ey[ky]>=0)//top to bottom
			{
				
				
				if(ex[kx]>=0)    //left to right
				{
					//dd=aa+bb+cc;
					for(inner=0;inner<num;inner++)
					{
						for(i=1;i<M;i++)
						{
							for(j=1;j<N;j++)
							{
								tau=mu_ref*exp(omega*log(T[i][j]))/(rho[i][j]*R*T[i][j]);
								cc=1.0/tau;
							//	df[i][j]=((g[i][j]-f[i][j])*cc-aa*(f[i][j]-f[i-1][j]+dx/2.0*(slg[i][j]-slg[i-1][j]))-bb*(f[i][j]-f[i][j-1]+dy/2.0*(slgg[i][j]-slgg[i][j-1]))+df[i-1][j]*aa+bb*df[i][j-1])/(aa+bb+cc);
							//	dfh[i][j]=((h[i][j]-fh[i][j])*cc-aa*(fh[i][j]-fh[i-1][j]+dx/2.0*(slh[i][j]-slh[i-1][j]))-bb*(fh[i][j]-fh[i][j-1]+dy/2.0*(slhh[i][j]-slhh[i][j-1]))+dfh[i-1][j]*aa+bb*dfh[i][j-1])/(aa+bb+cc);
								

								df[i][j]=((g[i][j]-f[i][j])*cc-aa*(xf[i][j]-xf[i-1][j])-bb*(yf[i][j]-yf[i][j-1])+df[i-1][j]*aa+bb*df[i][j-1])/(aa+bb+cc);
								dfh[i][j]=((h[i][j]-fh[i][j])*cc-aa*(xfh[i][j]-xfh[i-1][j])-bb*(yfh[i][j]-yfh[i][j-1])+dfh[i-1][j]*aa+bb*dfh[i][j-1])/(aa+bb+cc);
							

							}
							
						}
						for(i=1;i<M;i++)
						{
							for(j=1;j<N;j++)
							{
								f[i][j]+=df[i][j];
								fh[i][j]+=dfh[i][j];
							}
						}
						slope();

						

						deltaf1=0;

						for(i=1;i<M;i++)
						{
							for(j=1;j<N;j++)
							{
								if(deltaf1<fabs(df[i][j]))
								{
									deltaf1=fabs(df[i][j]);
								}

								if(deltaf1<fabs(dfh[i][j]))
								{
									deltaf1=fabs(dfh[i][j]);
								}

							}
						}
/* 						if(deltaf1<df1rhs)
						{
							num=1;
						}
						else
						{
							num=inner+2;
						} */
					


					}
					for(i=1;i<M;i++)
					{
						//f[i][N]=2*f[i][N-1]-f[i][N-2];
					//	f[i][N]=f[i][N-1]+slgg[i][N-1]*dy/2;
					//	fh[i][N]=fh[i][N-1]+slhh[i][N-1]*dy/2;

						f[i][N]=yf[i][N-1];
						fh[i][N]=yfh[i][N-1];

						j=N;
						mv[i][N]+=ey[ky]*f[i][N];
						//E1[i][j]+=0.5*((ey[ky]*ey[ky]+ex[kx]*ex[kx])*f[i][j]+fh[i][j]);
					}
					for(j=1;j<N;j++)
					{
						
						//f[M][j]=2*f[M-1][j]-f[M-2][j];
					//	f[M][j]=f[M-1][j]+slg[M-1][j]*dx/2;
					//	fh[M][j]=fh[M-1][j]+slh[M-1][j]*dx/2;

						f[M][j]=xf[M-1][j];
						fh[M][j]=xfh[M-1][j];


						mu[M][j]+=ex[kx]*f[M][j];
						//i=M;
						//E1[i][j]+=0.5*((ey[ky]*ey[ky]+ex[kx]*ex[kx])*f[i][j]+fh[i][j]);
					}

				}

				else   //right to left
				{
					//dd=-aa+bb+cc;
					for(inner=0;inner<num;inner++)
					{
						for(i=(M-1);i>=1;i--)
						{
							for(j=1;j<N;j++)
							{
								tau=mu_ref*exp(omega*log(T[i][j]))/(rho[i][j]*R*T[i][j]);
								cc=1.0/tau;
							//	df[i][j]=((g[i][j]-f[i][j])*cc-aa*(f[i+1][j]-f[i][j]-dx/2.0*(slg[i+1][j]-slg[i][j]))-bb*(f[i][j]-f[i][j-1]+dy/2.0*(slgg[i][j]-slgg[i][j-1]))-df[i+1][j]*aa+bb*df[i][j-1])/(-aa+bb+cc);
							//	dfh[i][j]=((h[i][j]-fh[i][j])*cc-aa*(fh[i+1][j]-fh[i][j]-dx/2.0*(slh[i+1][j]-slh[i][j]))-bb*(fh[i][j]-fh[i][j-1]+dy/2.0*(slhh[i][j]-slhh[i][j-1]))-dfh[i+1][j]*aa+bb*dfh[i][j-1])/(-aa+bb+cc);
						

								df[i][j]=((g[i][j]-f[i][j])*cc-aa*(xf[i][j]-xf[i-1][j])-bb*(yf[i][j]-yf[i][j-1])-df[i+1][j]*aa+bb*df[i][j-1])/(-aa+bb+cc);
								dfh[i][j]=((h[i][j]-fh[i][j])*cc-aa*(xfh[i][j]-xfh[i-1][j])-bb*(yfh[i][j]-yfh[i][j-1])-dfh[i+1][j]*aa+bb*dfh[i][j-1])/(-aa+bb+cc);
						
							}
							
							
						}
						for(i=1;i<M;i++)
						{
							for(j=1;j<N;j++)
							{
								f[i][j]+=df[i][j];
								fh[i][j]+=dfh[i][j];
							}
						}
						slope();

						

						deltaf1=0;

						for(i=1;i<M;i++)
						{
							for(j=1;j<N;j++)
							{
								if(deltaf1<fabs(df[i][j]))
								{
									deltaf1=fabs(df[i][j]);
								}

								if(deltaf1<fabs(dfh[i][j]))
								{
									deltaf1=fabs(dfh[i][j]);
								}

							}
						}
/* 						if(deltaf1<df1rhs)
						{
							num=1;
						}
						else
						{
							num=inner+2;
						} */

						

					}
					for(i=(M-1);i>=1;i--)
					{
						//f[i][N]=2*f[i][N-1]-f[i][N-2];
					//	f[i][N]=f[i][N-1]+slgg[i][N-1]*dy/2;
					//	fh[i][N]=fh[i][N-1]+slhh[i][N-1]*dy/2;

						f[i][N]=yf[i][N-1];
						fh[i][N]=yfh[i][N-1];


						mv[i][N]+=ey[ky]*f[i][N];
						//j=N;
						//E1[i][j]+=0.5*((ey[ky]*ey[ky]+ex[kx]*ex[kx])*f[i][j]+fh[i][j]);
					}

					
					for(j=1;j<N;j++)
					{
						
						//f[0][j]=2*f[1][j]-f[2][j];
					//	f[0][j]=f[1][j]-slg[1][j]*dx/2;
					//	fh[0][j]=fh[1][j]-slh[1][j]*dx/2;

						f[0][j]=xf[0][j];
						fh[0][j]=xfh[0][j];


						mu[0][j]+=ex[kx]*f[0][j];
						//i=0;
						//E1[i][j]+=0.5*((ey[ky]*ey[ky]+ex[kx]*ex[kx])*f[i][j]+fh[i][j]);

					}

				}
				
			}


			else      //bottom to top
			{

				if(ex[kx]>=0)         //left to right
				{
					//dd=aa-bb+cc;
					for(inner=0;inner<num;inner++)
					{

					
						
						for(i=1;i<M;i++)
						{
							for(j=(N-1);j>=1;j--)
							{
								tau=mu_ref*exp(omega*log(T[i][j]))/(rho[i][j]*R*T[i][j]);
								cc=1.0/tau;
							//	df[i][j]=((g[i][j]-f[i][j])*cc-aa*(f[i][j]-f[i-1][j]+dx/2.0*(slg[i][j]-slg[i-1][j]))-bb*(f[i][j+1]-f[i][j]-dy/2.0*(slgg[i][j+1]-slgg[i][j]))+df[i-1][j]*aa-bb*df[i][j+1])/(aa-bb+cc);
							//	dfh[i][j]=((h[i][j]-fh[i][j])*cc-aa*(fh[i][j]-fh[i-1][j]+dx/2.0*(slh[i][j]-slh[i-1][j]))-bb*(fh[i][j+1]-fh[i][j]-dy/2.0*(slhh[i][j+1]-slhh[i][j]))+dfh[i-1][j]*aa-bb*dfh[i][j+1])/(aa-bb+cc);
						

								df[i][j]=((g[i][j]-f[i][j])*cc-aa*(xf[i][j]-xf[i-1][j])-bb*(yf[i][j]-yf[i][j-1])+df[i-1][j]*aa-bb*df[i][j+1])/(aa-bb+cc);
								dfh[i][j]=((h[i][j]-fh[i][j])*cc-aa*(xfh[i][j]-xfh[i-1][j])-bb*(yfh[i][j]-yfh[i][j-1])+dfh[i-1][j]*aa-bb*dfh[i][j+1])/(aa-bb+cc);
						
							}
						}


						for(i=1;i<M;i++)
						{
							for(j=1;j<N;j++)
							{
								f[i][j]+=df[i][j];
								fh[i][j]+=dfh[i][j];
							}
						}


						slope();

						
						deltaf1=0;

						for(i=1;i<M;i++)
						{
							for(j=1;j<N;j++)
							{
								if(deltaf1<fabs(df[i][j]))
								{
									deltaf1=fabs(df[i][j]);
								}

								if(deltaf1<fabs(dfh[i][j]))
								{
									deltaf1=fabs(dfh[i][j]);
								}

							}
						}

						
/* 						if(deltaf1<df1rhs)
						{
							num=1;
						}
						else
						{
							num=inner+2;
						} */

						


					}


					for(i=1;i<M;i++)
					{
						//ftb[0][i][kx][ky]=2*f[i][1]-f[i][2];
						//f[i][0]=2*f[i][1]-f[i][2];
						f[i][0]=yf[i][0];
						fh[i][0]=yfh[i][0];
						mv[i][0]+=ey[ky]*f[i][0];
						//j=0;
						//E1[i][j]+=0.5*((ey[ky]*ey[ky]+ex[kx]*ex[kx])*f[i][j]+fh[i][j]);
					}
					for(j=1;j<N;j++)
					{
						
						//flr[1][j][kx][ky]=2*f[M-1][j]-f[M-2][j];
						//f[M][j]=2*f[M-1][j]-f[M-2][j];
						f[M][j]=xf[M-1][j];
						fh[M][j]=xfh[M-1][j];
						mu[M][j]+=ex[kx]*f[M][j];
						//i=M;
						//E1[i][j]+=0.5*((ey[ky]*ey[ky]+ex[kx]*ex[kx])*f[i][j]+fh[i][j]);
					}
				}

				else            //right to left
				{

					//dd=-aa-bb+cc;
					//dd=-aa-bb+cc/2;
					for(inner=0;inner<num;inner++)
					{

					
						for(i=(M-1);i>=1;i--)
						{
							for(j=(N-1);j>=1;j--)
							{
								tau=mu_ref*exp(omega*log(T[i][j]))/(rho[i][j]*R*T[i][j]);
								cc=1.0/tau;
							//	df[i][j]=((g[i][j]-f[i][j])*cc-aa*(f[i+1][j]-f[i][j]-dx/2.0*(slg[i+1][j]-slg[i][j]))-bb*(f[i][j+1]-f[i][j]-dy/2.0*(slgg[i][j+1]-slgg[i][j]))-df[i+1][j]*aa-bb*df[i][j+1])/(-aa-bb+cc);
							//	dfh[i][j]=((h[i][j]-fh[i][j])*cc-aa*(fh[i+1][j]-fh[i][j]-dx/2.0*(slh[i+1][j]-slh[i][j]))-bb*(fh[i][j+1]-fh[i][j]-dy/2.0*(slhh[i][j+1]-slhh[i][j]))-dfh[i+1][j]*aa-bb*dfh[i][j+1])/(-aa-bb+cc);
								

								df[i][j]=((g[i][j]-f[i][j])*cc-aa*(xf[i][j]-xf[i-1][j])-bb*(yf[i][j]-yf[i][j-1])-df[i+1][j]*aa-bb*df[i][j+1])/(-aa-bb+cc);
								dfh[i][j]=((h[i][j]-fh[i][j])*cc-aa*(xfh[i][j]-xfh[i-1][j])-bb*(yfh[i][j]-yfh[i][j-1])-dfh[i+1][j]*aa-bb*dfh[i][j+1])/(-aa-bb+cc);
						
							}
						}
						for(i=1;i<M;i++)
						{
							for(j=1;j<N;j++)
							{
								f[i][j]+=df[i][j];
								fh[i][j]+=dfh[i][j];
							}
						}
						slope();

						
						deltaf1=0;

						for(i=1;i<M;i++)
						{
							for(j=1;j<N;j++)
							{
								if(deltaf1<fabs(df[i][j]))
								{
									deltaf1=fabs(df[i][j]);
								}

								if(deltaf1<fabs(dfh[i][j]))
								{
									deltaf1=fabs(dfh[i][j]);
								}

							}
						}
/* 						if(deltaf1<df1rhs)
						{
							num=1;
						}
						else
						{
							num=inner+2;
						} */

						


					}
						
					for(i=1;i<M;i++)
					{
						//ftb[0][i][kx][ky]=2*f[i][1]-f[i][2];
						//f[i][0]=2*f[i][1]-f[i][2];
						f[i][0]=yf[i][0];
						fh[i][0]=yfh[i][0];
						mv[i][0]+=ey[ky]*f[i][0];
						//j=0;
						//E1[i][j]+=0.5*((ey[ky]*ey[ky]+ex[kx]*ex[kx])*f[i][j]+fh[i][j]);
					}
					
					for(j=1;j<N;j++)
					{	
						//flr[0][j][kx][ky]=2*f[1][j]-f[2][j];
						//f[0][j]=2*f[1][j]-f[2][j];
						f[0][j]=xf[0][j];
						fh[0][j]=xfh[0][j];
						mu[0][j]+=ex[kx]*f[0][j];
						//i=0;
						//E1[i][j]+=0.5*((ey[ky]*ey[ky]+ex[kx]*ex[kx])*f[i][j]+fh[i][j]);
					}
				}

			}


			for(i=1;i<M;i++)
			{
				for(j=1;j<N;j++)
				{
					if(deltaf<fabs(df[i][j]))
					{
						deltaf=fabs(df[i][j]);
					}

					if(deltaf<fabs(dfh[i][j]))
					{
						deltaf=fabs(dfh[i][j]);
					}

				}
			}

			for(i=1;i<M;i++)
			{
				for(j=1;j<N;j++)
				{
					rho1[i][j]+=f[i][j];
					mu[i][j]+=ex[kx]*f[i][j];
					mv[i][j]+=ey[ky]*f[i][j];
					E1[i][j]+=0.5*((ey[ky]*ey[ky]+ex[kx]*ex[kx])*f[i][j]+fh[i][j]);
					qflux_x[i][j]+=0.5*(ex[kx]-u[i][j])*(((ex[kx]-u[i][j])*(ex[kx]-u[i][j])+(ey[ky]-v[i][j])*(ey[ky]-v[i][j]))*f[i][j]+fh[i][j]);
					qflux_y[i][j]+=0.5*(ey[ky]-v[i][j])*(((ex[kx]-u[i][j])*(ex[kx]-u[i][j])+(ey[ky]-v[i][j])*(ey[ky]-v[i][j]))*f[i][j]+fh[i][j]);
				}
			}
		}
	}
	for(i=1;i<M;i++)
	{
		for(j=1;j<N;j++)
		{
			rho[i][j]=rho1[i][j];
			u[i][j]=mu[i][j]/rho1[i][j];
			v[i][j]=mv[i][j]/rho1[i][j];
			T[i][j]=(E1[i][j]-0.5*rho[i][j]*(u[i][j]*u[i][j]+v[i][j]*v[i][j]))/(1.5*R*rho1[i][j]);
		}
	}

}
void boundary()         //maxwellian reflection on the cell center
{
	double s1,s2,density;
	
	for(i=1;i<M;i++)
	{
		
		s2=0.0;
		s1=mv[i][N];
		for(ky=0;ky<Qy;ky++)  //bottom
		{
			if(ey[ky]>=0)
			{
				for(kx=0;kx<Qx;kx++)
				{
					
					//s1+=ey[ky]*ftb[1][i][kx][ky];
					s2+=ey[Qy-1-ky]*g1[kx][Qy-1-ky];//feq(kx,Qy-1-ky,0,0,1.0);
				}
			}
		}
		
		density=-s1/s2;
		rho[i][N]=density;
		
	

		s2=0.0;
		s1=mv[i][0];
		for(ky=0;ky<Qy;ky++)  //top
		{
			if(ey[ky]<0)
			{
				for(kx=0;kx<Qx;kx++)
				{
					
					//s1+=ey[ky]*ftb[0][i][kx][ky];
					s2+=ey[Qy-1-ky]*g2[kx][Qy-1-ky];//feq(kx,Qy-1-ky,u_wall,0,1.0);
				}
			}
			

		}
	
		density=-s1/s2;
		rho[i][0]=density;
		
	}
	for(j=1;j<N;j++)
	{
		
		s2=0.0;
		s1=mu[0][j];
		for(kx=0;kx<Qx;kx++)
		{
			if(ex[kx]<0)      //left
			{
				for(ky=0;ky<Qy;ky++)
				{
					//s1+=ex[kx]*flr[0][j][kx][ky];
					s2+=ex[Qx-1-kx]*g1[Qx-1-kx][ky];//feq(Qx-1-kx,ky,0,0,1.0);
				}
			}
		}
	
		density=-s1/s2;
		rho[0][j]=density;
		
		s2=0.0;
		s1=mu[M][j];
		for(kx=0;kx<Qx;kx++)
		{
			if(ex[kx]>=0)      //right
			{
				for(ky=0;ky<Qy;ky++)
				{
					//s1+=ex[kx]*flr[1][j][kx][ky];
					s2+=ex[Qx-1-kx]*g1[Qx-1-kx][ky];//feq(Qx-1-kx,ky,0,0,1.0);
				}
			}

		}
		density=-s1/s2;
		rho[M][j]=density;
		
	}
}


void outputdata()
{
	  FILE *fp;
	 fp=fopen("rhow.dat","w");
	 for(j=0;j<=N;j++)
	 {
	 	for(i=0;i<=M;i++)
	 	{
	 		fprintf(fp,"%e ",rho[i][j]);
	 	}
	 	fprintf(fp,"\n");
	 }
	  fclose(fp);
	 fp=fopen("uw.dat","w");
	 for(j=0;j<=N;j++)
	 {
	 	for(i=0;i<=M;i++)
	 	{
	 		fprintf(fp,"%24e ",u[i][j]);
	 	}
	 	fprintf(fp,"\n");
	 }
	  fclose(fp);
	 fp=fopen("vw.dat","w");
	 for(j=0;j<=N;j++)
	 {
	 	for(i=0;i<=M;i++)
	 	{
	 		fprintf(fp,"%24e ",v[i][j]);
	 	}
	 	fprintf(fp,"\n");
	 }
	  fclose(fp);
	  fp=fopen("Tw.dat","w");
	 for(j=0;j<=N;j++)
	 {
	 	for(i=0;i<=M;i++)
	 	{
	 		fprintf(fp,"%e ",T[i][j]);
	 	}
	 	fprintf(fp,"\n");
	 }
	  fclose(fp);
	  fp=fopen("ttw.dat","w");
	 for(j=0;j<1000;j++)
	 {
	 	fprintf(fp,"%e ",tt[j]);
	 	
	 }
	  fclose(fp);
	   fp=fopen("rrw.dat","w");
	 for(j=0;j<1000;j++)
	 {
	 	fprintf(fp,"%e ",rr[j]);
	 	
	 }
	  fclose(fp);

	  fp=fopen("rhsw.dat","w");
	 for(j=0;j<1000;j++)
	 {
	 	fprintf(fp,"%e ",rhs[j]);
	 	
	 }
	  fclose(fp);

	  fp=fopen("qyw.dat","w");
	 for(j=0;j<=N;j++)
	 {
	 	for(i=0;i<=M;i++)
	 	{
	 		fprintf(fp,"%e ",qflux_y[i][j]);
	 	}
	 	fprintf(fp,"\n");
	 }
	  fclose(fp);
	  fp=fopen("qxw.dat","w");
	 for(j=0;j<=N;j++)
	 {
	 	for(i=0;i<=M;i++)
	 	{
	 		fprintf(fp,"%e ",qflux_x[i][j]);
	 	}
	 	fprintf(fp,"\n");
	 }
	  fclose(fp);
	   fp=fopen("bushu.dat","w");
   for(j=0;j<1000;j++)
   {
    fprintf(fp,"%e ",bushu[j]);
    
   }
    fclose(fp);
}

