//Rescaling time TM and TM/TD = 1; Double exponential

    
#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <iomanip>  // for setting precision, see http://stackoverflow.com/questions/5907031/printing-the-correct-number-of-decimal-points-with-cout
#include <math.h>       /* round, floor, ceil, trunc */
#include <stdlib.h>     /* atof */
#include <stdio.h> /* for printf */
#include <cmath>       /* sqrt */
#include <algorithm>  /* std::max*/
#include <time.h>

#include <random>

// Split string at delimiter
// Source: http://code.runnable.com/VHb0hWMZp-ws1gAr/splitting-a-string-into-a-vector-for-c%2B%2B
#include <vector>
typedef std::vector<double> vec;
const int NUM_SECONDS =  1800;//half hour;

bool myfn(double i, double j) { return i<j; }
vec Evolution(vec V, vec Rand, int dim, vec Tau, vec Gamma);


int main(int argc, char** argv){
    
    //Print time
    double time_counter = 0;
    clock_t this_time = clock();
    clock_t last_time = this_time;
    
    
    int i, j, z;
    double q;
    int N;
    const int NExp=2;
    int Ntau;
    
    N = 2*(NExp + 2);
    Ntau = NExp;


    
    //Random number
    double r0;
    double r1;
    
    const int U0 = 3.;
    const float e = 2.71828183;
    
    
    time_t start, end;
    int elapsed_time, h, min, sec;
    start = time(NULL);
    
    vec v(N,0.);
    vec v0(N,0.);
    vec dv(N,0.);
    vec rand(2*NExp,0.);
    vec gamma(NExp,0.);
    vec tau(Ntau,0.);   //tau[0] = tauD/tauM tau[1] = tauD/tauG1  tau[2] = tauD/tauG2
    vec tauG2(8, 0.);

    
    tau[0] = 1/(0.1);
    tau[1] = 1/(10.);
    //tau[2] = 1/(0.1);
    
    
    //open file
    FILE *fv;
    fv = fopen("V3.txt", "w");
    FILE *fx;
    fx = fopen("X3.txt", "w");
    FILE *fy;
    fy = fopen("Y3.txt", "w");
    FILE *fxy;
    fxy = fopen("XY3.txt", "w");
    
    //FILE *fV;
    //fV = fopen("velocity1.txt", "w");
    int alpha =-100;
    gamma[0] = 1-alpha;  //=gamma_i = G_i/G    G_i= g_i/Nexp    G = sum_i G_i
    gamma[1] = alpha;  // gamma_1/gamma_2 = g_1/g_2
    
    printf("\n %f \n",gamma[0]);
    printf("\n %f \n",gamma[1]);
    
   ////////////////////////Time
    
    
    
    const double T1 = (1000.)*(pow(e,U0))*((pow(U0,-1)/(1 + (10)*U0*pow(tau[1],-1))) + pow(U0,-1)*pow(tau[0],-1) + 2*sqrt(pow(U0,-1)*pow(tau[0],-1)) + e*(pow(tau[1],-1)*pow(tau[1],-1)));
    
    const double T2 = (1000.)*(pow(e,U0))*((pow(U0,-1)/(1 + (10)*U0*pow(tau[2],-1))) + pow(U0,-1)*pow(tau[0],-1) + 2*sqrt(pow(U0,-1)*pow(tau[0],-1)) + e*(pow(tau[2],-1)*pow(tau[2],-1)));
    
    //printf("\nT1: %f\n",T1);
   // printf("\nT2: %f\n",T2);
    
    double t_final = 100000; //fmax(T1,T2);
    printf("\nt_final: %f\n",t_final);
    
    
    
    double t_step;
    double t;
    //printf("\nt_final: %f\n",T);
    
    //////////////////////////////////////////////////////////////////////////
    //Simulations
    //////////////////////////////////////////////////////////////////////////
    
        //fprintf(fA, "\t");
        //fprintf(fB, "\t");
        t = 0;
        t_step = 0.1;
    
        
      //printf("\nTau: %f\n",1./tau[1]);
    
    
        //check on t_step (tauD/tauM)*(1/t_step)<<1 & (tauD/tauG)*(1/t_step)<<1
        for (i = 0; i<(Ntau); i++){
            printf("\nTau: %f\n",1./tau[i]);
            if((tau[i])*(t_step)>0.01){
                t_step = 0.01*(1/tau[i]);
            }
        }
        
        printf("\ndt: %f\n",t_step);
    
        
        r0 = std::sqrt(2.*gamma[0]*(1./t_step));
        r1 = std::sqrt(2.*abs(gamma[1])*(1./t_step));
    
        std::random_device rd;
        std::mt19937 gen(rd());
    
        std::normal_distribution<> Gauss0(0, r0);
        std::normal_distribution<> Gauss1(0, r1);
        
     printf("\n tstep: %f\n",t_step);
     printf("\n tfinal: %f\n",t_final);
    
    double step = t_final/t_step;
    printf("\n step: %f\n",step);


        
        //initial condition
        for (i = 0; i<N; i++){
            v0[i] = 0.;
        }
        
        
        //Runge Kutta
        for (q = 0; q < step; q++){
            
            t = t + t_step;
            
            //Print t
            this_time = clock();
            time_counter += (double)(this_time - last_time);
            last_time = this_time;
            
            if(time_counter > (double)(NUM_SECONDS * CLOCKS_PER_SEC))
            {
                time_counter -= (double)(NUM_SECONDS * CLOCKS_PER_SEC);
                printf("\nTime:%f\t", t);
                printf("Step:%f\n", q);
            }
            
            
                rand[0] = Gauss0(gen);
                rand[1] = Gauss0(gen);
                rand[2] = Gauss1(gen);
                rand[3] = Gauss1(gen);
    
            
            //rand = sqrt(2./t_step);
            
            vec K1(N,0.);
            vec K2(N,0.);
            vec K3(N,0.);
            vec K4(N,0.);
            vec v1(N,0.);
            //std::cout<<"step\n"<<j;
            
            //K1
            K1 = Evolution(v0, rand, N, tau, gamma);
            
            //K2
            for (i = 0; i<N; i++){
                
                v1[i] = v0[i] + (t_step*K1[i])/2.;
                // std::cout<<"\nv1\n"<<v1[i]<<std::endl;
            }
            K2 = Evolution(v1, rand, N, tau, gamma);
            
            //K3
            for (i = 0; i<N; i++){
                
                v1[i] = v0[i] + (t_step*K2[i])/2.;
            }
            K3 = Evolution(v1, rand, N, tau, gamma);
            
            //K4
            for (i = 0; i<N; i++){
                
                v1[i] = v0[i] + t_step*K3[i];
            }
            K4 = Evolution(v1, rand, N, tau, gamma);
            
            
            for (i = 0; i<N; i++){
                
                v[i] = v0[i] + t_step*(K1[i] + 2*K2[i] + 2*K3[i] + K4[i])/6.;
            }
            
    
            fprintf(fv, "%f\n", v[1]);
            fprintf(fx, "%f\n", v[0]);
            fprintf(fy, "%f\n", v[3]);
            fprintf(fxy, "%f\n", v[3]*v[3] +v[0]*v[0]);
            //fprintf(fV, "%f\n", v[1]);
            for (i = 0; i<N; i++){
                
                v0[i] = v[i];
            }

            //printf("x0:%f\n",v0[0]);
            //printf("v0:%f\n",v0[1]);

        }
    
         printf("\nFinish Simulation:\n");

    
    //close file
    
    fclose(fv);
    fclose(fx);
    fclose(fy);
    fclose(fxy);
    //fclose(fV);
    
    
    /* Elapsed time */
    end = time(NULL);
    elapsed_time = (int)(end-start);
    h=(int)(elapsed_time/3600);
    min=(int)((elapsed_time-h*3600)/60);
    sec=elapsed_time-h*3600-min*60;
    printf("\nElapsed time  %d h  %d min  %d sec \n \n",h, min, sec);
    
    return 0;
}


//******Function****//
vec Evolution(vec V, vec Rand, int dim, vec Tau, vec Gamma){
    vec DV(dim);
    int i;
    //printf("r1:%f\n",Rand[0]);
    //printf("R1:%f\n",Rand[1]);
    //printf("r2:%f\n",Rand[2]);
    //printf("R2:%f\n",Rand[3]);
    //printf("tau1:%f\n",Tau[0]);
    //printf("tau2:%f\n",Tau[1]);
    
    //Evolution
    DV[2] = - Tau[0]*(V[2] + V[1]*(Gamma[0]) - Rand[0]);
    DV[5] = - Tau[0]*(V[5] + V[4]*(Gamma[0]) - Rand[1]);

    DV[6] = - Tau[1]*(V[6] + V[1]*(Gamma[1]) - Rand[2]);
    DV[7] = - Tau[1]*(V[7] + V[4]*(Gamma[1]) - Rand[3]);
    //printf("\nsum:%f\n",Sum);
    DV[0] = V[1];
    DV[1] = V[2]+V[6];
    
    DV[3] = V[4] ;
    DV[4] = V[5]+V[7];
    
    //printf("x:%f\n",V[0]);
    //printf("v:%f\n",V[1]);
    //printf("r1:%f\n",V[2]);
    //printf("y:%f\n",V[3]);
    //printf("w:%f\n",V[4]);
    //printf("R1:%f\n",V[5]);
    //printf("r2:%f\n",V[6]);
    //printf("R2:%f\n",V[7]);
    return DV;
}



