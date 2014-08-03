#include "config.h"

#ifdef HAVE_MALLOC_H
# include<malloc.h>
#endif
#include<math.h>
#include<stdio.h>
#include<stdlib.h>

#define FREE_ARG char*
#define NR_END 1
#define PI 3.141592653589793

float *vector(nl,nh)
long nh,nl;
/* allocate a float vector with subscript range v[nl..nh] */
{
        float *v;

        v=(float *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(float)));
        return v-nl+NR_END;
}

int *ivector(nl,nh)
long nh,nl;
/* allocate an int vector with subscript range v[nl..nh] */
{
        int *v;

        v=(int *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(int)));
        return v-nl+NR_END;
}

float **matrix(nrl,nrh,ncl,nch)
int nrl,nrh,ncl,nch;
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
    int i;
    float **m;

        /*allocate pointers to rows */
        m=(float **) malloc((unsigned) (nrh-nrl+1)*sizeof(float*));
    m -= nrl;

   /*allocate rows and set pointers to them */
      for(i=nrl;i<=nrh;i++) {
                      m[i]=(float *) malloc((unsigned) (nch-ncl+1)*sizeof(float)
);
            m[i] -= ncl;
    }
      /* return pointer to array of pointers to rows */
      return m;
}

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

void fourn(float data[], int nn[], int ndim, int isign)
{
	int idim;
	unsigned long i1,i2,i3,i2rev,i3rev,ip1,ip2,ip3,ifp1,ifp2;
	unsigned long ibit,k1,k2,n,nprev,nrem,ntot;
	float tempi,tempr;
	double theta,wi,wpi,wpr,wr,wtemp;

	for (ntot=1,idim=1;idim<=ndim;idim++)
		ntot *= nn[idim];
	nprev=1;
	for (idim=ndim;idim>=1;idim--) {
		n=nn[idim];
		nrem=ntot/(n*nprev);
		ip1=nprev << 1;
		ip2=ip1*n;
		ip3=ip2*nrem;
		i2rev=1;
		for (i2=1;i2<=ip2;i2+=ip1) {
			if (i2 < i2rev) {
				for (i1=i2;i1<=i2+ip1-2;i1+=2) {
					for (i3=i1;i3<=ip3;i3+=ip2) {
						i3rev=i2rev+i3-i2;
						SWAP(data[i3],data[i3rev]);
						SWAP(data[i3+1],data[i3rev+1]);
					}
				}
			}
			ibit=ip2 >> 1;
			while (ibit >= ip1 && i2rev > ibit) {
				i2rev -= ibit;
				ibit >>= 1;
			}
			i2rev += ibit;
		}
		ifp1=ip1;
		while (ifp1 < ip2) {
			ifp2=ifp1 << 1;
			theta=isign*6.28318530717959/(ifp2/ip1);
			wtemp=sin(0.5*theta);
			wpr = -2.0*wtemp*wtemp;
			wpi=sin(theta);
			wr=1.0;
			wi=0.0;
			for (i3=1;i3<=ifp1;i3+=ip1) {
				for (i1=i3;i1<=i3+ip1-2;i1+=2) {
					for (i2=i1;i2<=ip3;i2+=ifp2) {
						k1=i2;
						k2=k1+ifp1;
						tempr=(float)wr*data[k2]-(float)wi*data[k2+1];
						tempi=(float)wr*data[k2+1]+(float)wi*data[k2];
						data[k2]=data[k1]-tempr;
						data[k2+1]=data[k1+1]-tempi;
						data[k1] += tempr;
						data[k1+1] += tempi;
					}
				}
				wr=(wtemp=wr)*wpr-wi*wpi+wr;
				wi=wi*wpr+wtemp*wpi+wi;
			}
			ifp1=ifp2;
		}
		nprev *= n;
	}
}
#undef SWAP

int
main()
{    float dum,delta,alpha,fact,*w,delrho;
     int lattice_size_x,lattice_size_y,i,j;
     int *nn;
     FILE *fp1,*fp2;

     fp1=fopen("load2dandes","r");
     fp2=fopen("load2dandesdeflect50","w");
     lattice_size_x=2048;
     lattice_size_y=4096;
     delta=1.1;   /* km */
     alpha=50.0;  /* (4D/((rho_m-rho_c)*g))^0.25, in units of multiples of delta */
     delrho=0.27; /* (rho_m-rho_c)/rho_c */
     nn=ivector(1,2);
     nn[1]=lattice_size_x;
     nn[2]=lattice_size_y;
     w=vector(1,2*lattice_size_x*lattice_size_y);
     for (j=1;j<=lattice_size_y;j++)
      for (i=1;i<=lattice_size_x;i++)
       {fscanf(fp1,"%f",&dum);
        if (dum<1000) w[2*(i-1)*lattice_size_y+2*j-1]=0; 
         else w[2*(i-1)*lattice_size_y+2*j-1]=dum;
        w[2*(i-1)*lattice_size_y+2*j]=0.0;}
     fourn(w,nn,2,1);
     w[1]*=1/delrho;
     w[2]*=1/delrho;
     for (j=1;j<=lattice_size_y/2;j++)
      {fact=1/(delrho+4*delrho*pow(alpha*j*PI/lattice_size_y,4.0));
       w[2*j+1]*=fact;
       w[2*j+2]*=fact;
       w[2*lattice_size_y-2*j+1]*=fact;
       w[2*lattice_size_y-2*j+2]*=fact;}
     for (i=1;i<=lattice_size_x/2;i++)
      {fact=1/(delrho+4*delrho*pow(alpha*i*PI/lattice_size_x,4.0));
       w[2*i*lattice_size_y+1]*=fact;
       w[2*i*lattice_size_y+2]*=fact;
       w[2*lattice_size_x*lattice_size_y-2*i*lattice_size_y+1]*=fact;
       w[2*lattice_size_x*lattice_size_y-2*i*lattice_size_y+2]*=fact;}
     for (i=1;i<=lattice_size_x/2;i++)
      for (j=1;j<=lattice_size_y/2;j++)
       {fact=1/(delrho+4*delrho*pow(alpha*sqrt((PI*PI*i*i)/(lattice_size_x*lattice_size_x)+
         (PI*PI*j*j)/(lattice_size_y*lattice_size_y)),4.0));
        w[2*i*lattice_size_y+2*j+1]*=fact;
        w[2*i*lattice_size_y+2*j+2]*=fact;
        w[2*i*lattice_size_y+2*lattice_size_y-2*j+1]*=fact;
        w[2*i*lattice_size_y+2*lattice_size_y-2*j+2]*=fact;
        w[2*lattice_size_x*lattice_size_y-2*(i-1)*lattice_size_y-
          2*(lattice_size_y-j)+1]*=fact;
        w[2*lattice_size_x*lattice_size_y-2*(i-1)*lattice_size_y-
          2*(lattice_size_y-j)+2]*=fact;
        w[2*lattice_size_x*lattice_size_y-2*lattice_size_y-
          2*(i-1)*lattice_size_y+2*(lattice_size_y-j)+1]*=fact;         
		w[2*lattice_size_x*lattice_size_y-2*lattice_size_y-
          2*(i-1)*lattice_size_y+2*(lattice_size_y-j)+2]*=fact;}
     fourn(w,nn,2,-1);
     for (j=1;j<=lattice_size_y;j++)
      for (i=1;i<=lattice_size_x;i++)
       fprintf(fp2,"%f\n",w[2*(i-1)*lattice_size_y+2*j-1]
        /(lattice_size_x*lattice_size_y));
     fclose(fp1);
     fclose(fp2);
}
