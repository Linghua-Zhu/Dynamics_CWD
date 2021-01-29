//
//  main.cpp
//  FirstMethod
//
//  Created by Zhu Linghua on 8/1/16.
//  Copyright (c) 2016 Zhu Linghua. All rights reserved.
//

#include <stdio.h>
#include <math.h>

#define NN 512
#define Pi 3.14159265
#define PP 2001

void initial_kx(double kx[], int N)
{
    for(int i = 0; i < N; i++)
    {
        kx[i]=2*Pi/N*(i+1)-Pi;
    }
}
void initial_ky(double ky[], int N)
{
    for(int i = 0; i < N; i++)
    {
        ky[i]=2*Pi/N*(i+1)-Pi;
    }
}
void Set_Elk(double Elk[2][NN][NN], double kx[], double ky[], int N, double t, double lm, double u0)
{
    for(int i = 0; i < N; i++)
    {
        kx[i]=2*Pi/N*(i+1)-Pi;
    }
    for(int i = 0; i < N; i++)
    {
        ky[i]=2*Pi/N*(i+1)-Pi;
    }
    for(int q=0; q<2; q++)
    {
        for(int i=0; i<N; i++)
        {
            for(int j=0; j<N; j++)
            {
                Elk[q][i][j]=pow(-1,q)*sqrt(4*t*t*(cos(kx[i])+cos(ky[j]))*(cos(kx[i])+cos(ky[j]))+lm*lm*u0*u0);
            }
        }
    }
}
void Set_Elkp(double Elkp[2][NN][NN], double kx[], double ky[], int N, double t, double lm, double u0)
{
    for(int i = 0; i < N; i++)
    {
        kx[i]=2*Pi/N*(i+1)-Pi;
    }
    for(int i = 0; i < N; i++)
    {
        ky[i]=2*Pi/N*(i+1)-Pi;
    }
    for(int q=0; q<2; q++)
    {
        for(int i=0; i<N; i++)
        {
            for(int j=0; j<N; j++)
            {
                Elkp[q][i][j]=pow(-1,q)*sqrt(4*t*t*(cos(kx[i])+cos(ky[j]))*(cos(kx[i])+cos(ky[j]))+lm*lm*u0*u0);
            }
        }
    }
}
void give_print(double En[], double g[], double Dos[], int Ne, char *Str)
{
    FILE *fp=fopen(Str,"w");
    for(int i=0;i<Ne;i++)
    {
        if(Dos[i]!=0.0)
        {
            fprintf(fp,"%13.12f    %13.12f\n",En[i],g[i]);
        }
    }
    fclose(fp);
}
int main()
{
    int i,j,k,q,p,Np,N=NN,Nk=1500,Ne=3600,P=PP;
    int j1,j2,j3,j4,jmin,jmax;
    int Flag[Ne];
    double m=0.0033212,t,lm,K;
    double gm=0.009,gm0,gm3,gm4,eps;
    double kb=0.000086173324;
    double Emax,Emin,IDos[Ne],Dos[Ne],delE;
    double dt,Tt,Sum;
    double v0,v,w,w0;
    double kx[N],ky[N];
    double Elk[2][NN][NN],En[Ne+1];
    double fe[Ne],fk[2][N][N],Elkp[2][NN][NN];
    double feold[Ne],fkold[2][N][N],f_r[Ne];
    double u0,u,u_eq;
    double d0,DE,W,dfmax,gee_const;
    double df[2][N][N];
    double T,Te;
    double WD,Ep[20],Dpb[20],Dp[20],GM,Alf,Sp,gep[Ne],g[Ne],bw[Np],x,y,z;
    double gee[Ne],Kee,Kep;
    double dE2,dE1,dE0,Etot;
    double a,b;
    double dn,eve[NN][NN];
    double up,Eel0,Nlkupper,Nlklower;
    char *Str;
    char Name[12];
    double Eel;
    FILE *fp,*fp2,*fp3,*fp1;
    
    initial_kx(kx,NN);
    initial_ky(ky,NN);
    
    dt=0.035;
    
    t=0.5;
    lm=1.2*0.78;
    K=0.85;
    
    DE=2.0;
    W=0.02;
    dfmax=0.4*0.0;
    
    Kee=651.0*3.0;
    Kep=0.93*0.0;
    
    T=260.0;
    Te=kb*T;
    
    u0=0.000001;
    v0=0.0;
    
    eps=0.0000000001;
    
    fp=fopen("0Parameters.txt","w");
    fprintf(fp,"N=%d\nNe=%d\nm=%f\ngm=%f\nt=%f\nlm=%f\nK=%f\ngm=%f\nTini=%f\ndfmax=%f\nKee=%f\nW=%f\nu0=%f\nv0=%f\nKep=%f\n",N,Ne,m,gm,t,lm,K,gm,T,dfmax,Kee,W,u0,v0,Kep);
    fclose(fp);

    for(p=0;p<P;p++)
    {
        Tt=p*dt;
        Set_Elk(Elk,kx,ky,NN,t,lm,u0);
        
        if(p==0)
        {
            for(q=0;q<2;q++)
            {
                for(i=0;i<N;i++)
                {
                    for(j=0;j<N;j++)
                    {
                        fk[q][i][j]=1.0/(exp(Elk[q][i][j]/Te)+1.0);
                    }
                }
            }
        }
        if(p==0)
        {
            for(q=0;q<2;q++)
            {
                for(i=0;i<N;i++)
                {
                    for(j=0;j<N;j++)
                    {
                        df[q][i][j]=pow(-1,q)*dfmax*exp(-(Elk[q][i][j]-pow(-1,q)*DE/2)*(Elk[q][i][j]-pow(-1,q)*DE/2)/(2*W*W));
                    }
                }
            }
            for(i=0;i<N;i++)
            {
                for(j=0;j<N;j++)
                {
                    fk[0][i][j]+=df[0][i][j];
                    fk[1][i][j]+=df[1][i][j];
                }
            }
          
            fp=fopen("Start_Distribution_k.txt","w");
            for(q=0;q<2;q++)
	        {
                for(i=0;i<N;i++)
                {
                     for(j=0;j<N;j++)
                     {
                         fprintf(fp,"%16.15f      %32.31f\n",Elk[q][i][j],fk[q][i][j]);
                     }
                 }
	        }
	        fclose(fp);
         }

        dE0=0.0;
        for(i=0;i<N;i++)
        {
            for(j=0;j<N;j++)
            {
                dE0+=lm*lm*u0/Elk[0][i][j]*fk[0][i][j]-lm*lm*u0/Elk[0][i][j]*fk[1][i][j];
            }
        }
        dE1=2*K*u0+dE0/N/N/2;
        
        dE0=0.0;
        for(i=0;i<N;i++)
        {
            for(j=0;j<N;j++)
            {
                dE0+=(lm*lm*Elk[0][i][j]-lm*lm*lm*lm*u0*u0/Elk[0][i][j])/(4*t*t*(cos(kx[i])+cos(ky[j]))*(cos(kx[i])+cos(ky[j]))+lm*lm*u0*u0)*(fk[0][i][j]-fk[1][i][j]);
            }
        }
        dE2=2*K+dE0/N/N/2;
        
        a=dE2;
        b=dE1-a*u0;
        
        printf("%d      u0=%16.15f\n",p,u0);
        if(a==0.0)
        {
            printf("a is zero!");
        }
        
        u_eq=-b/a;
        
        d0=gm*gm/(4*m*m)-2.0*a/m;
        
        if(d0<0.0)
        {
            w0=sqrt(a/2/m);
            gm0=gm/4/m;
            w=sqrt(w0*w0-gm0*gm0);
            u=sin(w*dt)/w*exp(-gm0*dt)*v0+(gm0/w*sin(w*dt)+cos(w*dt))*exp(-gm0*dt)*(u0-u_eq)+u_eq;
            v=(-gm0/w*sin(w*dt)+cos(w*dt))*exp(-gm0*dt)*v0-(gm0*gm0/w+w)*exp(-gm0*dt)*sin(w*dt)*(u0-u_eq);
        }
        if(d0==0.0)
        {
            w0=sqrt(a/2/m);
            gm0=gm/4/m;
            u=exp(-gm0*dt)*(1+gm0*dt)*(u0-u_eq)+exp(-gm0*dt)*v0*dt+u_eq;
            v=-gm0*gm0*dt*exp(-gm0*dt)*(u0-u_eq)+(1-gm0*dt)*exp(-gm0*dt)*v0;
        }
        if(d0>0.0)
        {
            gm3=-gm/4/m+sqrt(gm*gm-4*a*2*m)/4/m;
            gm4=-gm/4/m-sqrt(gm*gm-4*a*2*m)/4/m;
            u=(-gm4/(gm3-gm4)*(u0-u_eq)+1/(gm3-gm4)*v0)*exp(gm3*dt)+(gm3/(gm3-gm4)*(u0-u_eq)-1/(gm3-gm4)*v0)*exp(gm4*dt)+u_eq;
            v=(-gm3*gm4/(gm3-gm4)*(u0-u_eq)+gm3/(gm3-gm4)*v0)*exp(gm3*dt)+(gm3*gm4/(gm3-gm4)*(u0-u_eq)-gm4/(gm3-gm4)*v0)*exp(gm4*dt);
        }
        
        fp2=fopen("2u.txt","a");
        fprintf(fp2,"%f    %16.15f\n",Tt,u0);
        fclose(fp2);
        
        Emax=sqrt(16*t*t+lm*lm*u0*u0);                       /*Density of States*/
        Emin=-Emax;
        delE=(Emax-Emin)/(Ne-1);
        
        for(k=0;k<Ne;k++)
        {
            En[k]=Emin+k*delE;
        }
        for(i=0;i<Ne;i++)
        {
            IDos[i]=0.0;
        }
        for(q=0;q<2;q++)
        {
            for(i=0;i<N;i++)
            {
                for(j=0;j<N;j++)
                {
                    k=int((Elk[q][i][j]-Emin+0.5*delE)/delE);
                    IDos[k]+=1/delE;
                }
            }
        }
        for(k=0;k<Ne;k++)
        {
            Dos[k]=IDos[k]/(2*N*N);
        }
      
        /*Sum=0.0;
        for(i=0;i<Ne+1;i++)
        {
            Sum+=Dos[i]*delE;
        }
        printf("%f\n",Sum);*/
        for(i=0;i<Ne;i++)
        {
            fe[i]=0.0;
        }
        for(q=0;q<2;q++)                  /*Distribution of Energy*/
        {
            for(i=0;i<N;i++)
            {
                for(j=0;j<N;j++)
                {
                    k=int((Elk[q][i][j]-Emin+0.5*delE)/delE);
                    fe[k]+=fk[q][i][j]/(Dos[k]*N*N*delE)/2;
                }
            }
        }
        /*Sum=0.0;
        for(i=0;i<Ne+1;i++)
        {
            Sum+=fe[i]*Dos[i]*delE;
        }
        printf("%f\n",Sum);*/
        
        for(k=0;k<Ne;k++)
        {
            if(fe[k]>eps && fe[k+1]<eps)
            {
                jmax=k;
            }
        }
        
        jmin=Ne-jmax-1;

        if(p==0)
        {
            Str="Start_Distribution_E.txt";
            give_print(En,fe,Dos,Ne,Str);
        }
        for(i=0;i<Ne;i++)                                                /*Boltzmann Scattering*/
        {
            if(Dos[i]<0.000000000000001)
            {
                Flag[i]=0;
            }
            if(Dos[i]>0.000000000000001)
            {
                Flag[i]=1;
            }
        }
        for(i=0;i<Ne;i++)
        {
            g[i]=0.0;
            gee[i]=0.0;
            gep[i]=0.0;
        }
        gee_const = delE * delE * Kee * 0.5;
        for(j1=jmin;j1<jmax+1;j1++)
        {
            for(j2=jmin;j2<jmax+1;j2++)
            {
                for(j3=jmin;j3<jmax+1;j3++)
                {
                    j4 = j1 + j3 - j2;
                    if (j4>=jmin&&j4<jmax+1)
                    {
                        if (j1==j2||j1==j4) continue;
                        if (Flag[j1]==1 && Flag[j2]==1 && Flag[j3]==1 && Flag[j4]==1)
                        {
                            gee[j1]+=gee_const*Dos[j2]*Dos[j3]*Dos[j4]*(fe[j1]*fe[j3]*(fe[j2]+fe[j4]-1.0)-fe[j2]*fe[j4]*(fe[j1]+fe[j3]-1.0));
                        }
                    }
                }
            }
        }
        
        for(i=0;i<Ne+1;i++)
        {
            g[i]=gee[i]+gep[i];
        }
        for(i=0;i<Ne+1;i++)
        {
            fe[i]+=g[i]*dt;
        }
        
        for(q=0;q<2;q++)                            /*ReGetting lk from energy*/
        {
            for(i=0;i<N;i++)
            {
                for(j=0;j<N;j++)
                {
                    k=int((Elk[q][i][j]-Emin+0.5*delE)/delE);
                    fk[q][i][j]=fe[k];
                }
            }
        }
        
        u0=u;
        v0=v;
    }
    
    Str="New_Distribution_E.txt";
    give_print(En,fe,Dos,Ne,Str);

    
    fp=fopen("Distribution_k_ForSave.txt","w");
    for(q=0;q<2;q++)
    {
        for(i=0;i<N;i++)
        {
            for(j=0;j<N;j++)
            {
                fprintf(fp,"%32.31f\n",fk[q][i][j]);
            }
        }
    }
    fclose(fp);
    
    fp=fopen("New_Distribution_k.txt","w");
    for(q=0;q<2;q++)
    {
        for(i=0;i<N;i++)
        {
            for(j=0;j<N;j++)
            {
                fprintf(fp,"%16.15f      %32.31f\n",Elk[q][i][j],fk[q][i][j]);
            }
        }
    }
    fclose(fp);

    
    return(0);
}
