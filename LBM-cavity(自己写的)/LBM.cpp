#include<iostream>
# include<cmath>
# include<cstdlib>
# include<iomanip>
# include<fstream>
# include<sstream>
# include<string>
using namespace std;

const int Q = 9;
const int NX = 256;
const int NY = 256;
const double U = 0.1;

int e[Q][2] = {{0,0},{1,0},{0,1},{-1,0},{0,-1},{1,1},{-1,1},{-1,-1},{1,-1}};
double w[Q] = {4.0/9,1.0/9,1.0/9,1.0/9,1.0/9,1.0/36,1.0/36,1.0/36,1.0/36};
double rho[NX+1][NY+1],u[NX+1][NY+1][2],u0[NX+1][NY+1][2],f[NX+1][NY+1][Q],F[NX+1][NY+1][Q];
int i,j,k,ip,jp,n;
double c,Re,dx,dy,Lx,Ly,dt,rho0,P0,tau_f,niu,error;

double feq(int k,double rho,double u[2]);

double feq(int k,double rho,double u[2])
{
	double eu,uv,feq;
	
	eu = (e[k][0]*u[0]+e[k][1]*u[1]);
	uv = (u[0]*u[0]+u[1]*u[1]);
	feq = w[k]*rho*(1.0+3.0*eu+4.5*eu*eu-1.5*uv);
	
	
	return feq;
}







int main()
{
	
	dx = 1.0;
	dy = 1.0;
	Lx = dx*double(NX);
	Ly = dy*double(NY);
	dt = dx;
	c = dx/dt;
	rho0 = 1.0;
	Re = 1000;
	niu = (U*Lx)/Re;
	tau_f = 3.0*niu+0.5;
	
	
	for( i=0;i<=NX;i++)
	{
		for( j=0;j<=NY;j++)
		{
			u[i][j][0] = 0;
			u[i][j][1] = 0;
			rho[i][j] = rho0;
			u[i][NY][0] = U;
			for(int k=0;k<Q;k++)
			{
				f[i][j][k] = feq(k,rho[i][j],u[i][j]); 
				
			}
			
		}
		
	}
	
	
for(int t=0;t<10000;t++)
{
	
	
	for( i=1;i<NX;i++)
	{
		
		for( j=1;j<NY;j++)
		{
			
			
			for( k=0;k<Q;k++)
			{
				
				ip = i-e[k][0];
				jp = j-e[k][1];
				
				F[i][j][k] = f[ip][jp][k]+(feq(k,rho[ip][jp],u[ip][jp])-f[ip][jp][k])/tau_f;
				
				
			}
	
		}
		

	}
	
	
	for( i=1;i<NX;i++)
	{
		
		for( j=1;j<NY;j++)
		{
			u0[i][j][0] = u[i][j][0];
			u0[i][j][1] = u[i][j][1];
			rho[i][j] = 0;
			u[i][j][0] = 0; 
			u[i][j][1]  =0;
			
			for( k=0;k<Q;k++)
			{
				
				f[i][j][k] = F[i][j][k];
				rho[i][j] +=f[i][j][k];
				u[i][j][0] +=e[k][0]*f[i][j][k];
				u[i][j][1] +=e[k][1]*f[i][j][k];
			}
			
		 u[i][j][0] /=rho[i][j];
		 u[i][j][1] /=rho[i][j];
	
		}
		
	
	}
	
	for( j=1;j<NY;j++)
	{
		for( k=0;k<Q;k++)
		{
			rho[NX][j] = rho[NX-1][j];
			f[NX][j][k] = feq(k,rho[NX][j],u[NX][j])+f[NX-1][j][k]-feq(k,rho[NX-1][j],u[NX-1][j]);
			
			rho[0][j] = rho[1][j];
			f[0][j][k] = feq(k,rho[0][j],u[0][j])+f[1][j][k]-feq(k,rho[1][j],u[1][j]);
			
		
		}
		
		
	}
	
	
	for( i=0;i<=NX;i++)
	{
		for( k=0;k<Q;k++)
		{
			
			rho[i][0] = rho[i][1];
			f[i][0][k] = feq(k,rho[i][0],u[i][0])+f[i][1][k]-feq(k,rho[i][1],u[i][1]);
			
			rho[i][NY] = rho[i][NY-1];
			f[i][NY][k] = feq(k,rho[i][NY],u[i][NY])+f[i][NY-1][k]-feq(k,rho[i][NY-1],u[i][NY-1]);
		 	u[i][NY][0] = U;	
		 
	
		}

		
	}
}
	ostringstream name;
	  name<<"cavity_"<<3<<".dat";
	  ofstream out(name.str().c_str());
	  out<< "Title= \"LBM Lid Driven Flow\"\n" << "VARIABLES=\"X\",\"Y\",\"U\",\"V\"\n" << "ZONE T=\"BOX\",I=" << NX+1 << ",J=" << NY+1 << ",F=POINT" << endl;
	  for(int j=0;j<=NY;j++)
	     for(int i=0;i<=NX;i++)
	     {
	     	out<<double(i)/Lx<<" "<<double(j)/Ly<<" "<<u[i][j][0]<<" "<<u[i][j][1] <<endl;
		 }
	
	
	
	
	
	
	
	
	

	return 0;
 } 
