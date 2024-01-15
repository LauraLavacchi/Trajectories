//Vic code con gamma diverse

    
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
const int NUM_SECONDS = 1800;// half hour;

bool myfn(double i, double j) { return i<j; }
vec Evolution(vec V, double Rand, int dim, vec Tau, double B, double U);


int main(int argc, char** argv){
    
    //Print time
    double time_counter = 0;
    clock_t this_time = clock();
    clock_t last_time = this_time;
    
    int i, j, z;
    double q;
    
    int N;
    const int NExp=1;
    int Ntau;
    
    N = NExp + 2;
    Ntau = NExp + 1;
    
    //Random number
    double r;
    double rand;
    
    const float e = 2.71828183;
    const int U0 = 3.;
    double b = 0.1;
    double Max = -(1./8)*(3.*b - sqrt(9.*b*b-10.*b+17.) +1.);
    double Min = -(1./8)*(3.*b +sqrt(9.*b*b-10.*b+17.) +1.);
    
    printf("\n b:%f\n", b);
    printf("\n Max:%f\n", Max);
    printf("\n Min:%f\n", Min);
    
    
    time_t start, end;
    int elapsed_time, h, min, sec;
    start = time(NULL);
    
    vec v(N,0.);
    vec v0(N,0.);
    vec dv(N,0.);
    vec tau(Ntau,0.);   //tau[0] = tauD/tauM tau[1] = tauD/tauG1  tau[2] = tauD/tauG2

    
    tau[0] = 1./(0.1);
    tau[1] = 1./(1.);
    //open file
    //
    //FILE *fdiff1;
    //fdiff1 = fopen("Diff1G1b01.txt", "w");
    //FILE *fdiff2;
    //fdiff2 = fopen("Diff2G1b01.txt", "w");
    //printf("\n DiffG1b01.txt \n");
    //FILE *fX;
    //fX = fopen("Xb1.txt", "w");
    FILE *fX;
    fX = fopen("TXb01.txt", "w");

   ////////////////////////Time
   
    double t_final = 500000; //(100000)*(pow(e,U0))*((pow(U0,-1)/(1 + (10)*U0*pow(tau[1],-1))) + pow(U0,-1)*pow(tau[0],-1) + 2*sqrt(pow(U0,-1)*pow(tau[0],-1)) + e*(pow(tau[1],-1)*pow(tau[1],-1)));;
    
    double t_step;
    double t;
  
        t = 0;
        t_step = 0.1;
    
   
        for (i = 0; i<(Ntau); i++){
            printf("\n Tau: %f\n",1./tau[i]);
            if((tau[i])*(t_step)>0.01){
                t_step = 0.01*(1/tau[i]);
            }
        }
    
        
        r = std::sqrt(2.*(1/t_step));
    
        std::random_device rd;
        std::mt19937 gen(rd());
    
        std::normal_distribution<> Gauss(0, r);
    
     printf("\n tstep: %f\n",t_step);
     printf("\n tfinal: %f\n",t_final);
    
    double step = t_final/t_step;
        printf("\n step: %f\n",step);

    
        //int step = 100;
        vec A(0, 0.);
        vec B(0, 0.);
        vec C(0, 0.);
    
        
        //initial condition
        for (i = 0; i<N; i++){
            v0[i] = 0.;
        }
        
        
        //Runge Kutta
        for (q = 0; q < step; q++){
            
            t = t + t_step;
            
            this_time = clock();
            time_counter += (double)(this_time - last_time);
            last_time = this_time;
            
            if(time_counter > (double)(NUM_SECONDS * CLOCKS_PER_SEC))
            {
                time_counter -= (double)(NUM_SECONDS * CLOCKS_PER_SEC);
                printf("Time:%f\t", t);
                printf("Step:%f\n", q);
                printf("Time/2:%f\n", (step)/(q));
            }
            
                rand = Gauss(gen);
            
            //rand = sqrt(2./t_step);
            
            vec K1(N,0.);
            vec K2(N,0.);
            vec K3(N,0.);
            vec K4(N,0.);
            vec v1(N,0.);
            //std::cout<<"step\n"<<j;
            
            //K1
            K1 = Evolution(v0, rand, N, tau, b, U0);
            
            //K2
            for (i = 0; i<N; i++){
                
                v1[i] = v0[i] + (t_step*K1[i])/2.;
                // std::cout<<"\nv1\n"<<v1[i]<<std::endl;
            }
            K2 = Evolution(v1, rand, N, tau, b, U0);
            
            //K3
            for (i = 0; i<N; i++){
                
                v1[i] = v0[i] + (t_step*K2[i])/2.;
            }
            K3 = Evolution(v1, rand, N, tau, b, U0);
            
            //K4
            for (i = 0; i<N; i++){
                
                v1[i] = v0[i] + t_step*K3[i];
            }
            K4 = Evolution(v1, rand, N, tau, b, U0);
            
            
            for (i = 0; i<N; i++){
                
                v[i] = v0[i] + t_step*(K1[i] + 2*K2[i] + 2*K3[i] + K4[i])/6.;
            }
            
            
            if ((v0[0] > 1 && v[0] < 1) || (v0[0] < 1 && v[0] > 1)){
                A.push_back(t);
                //std::cout << "Size A: " << A.size() << '\n';
                //fprintf(fA, "%f\n", t);
            }
            
            else if((v0[0] > Max && v[0] < Max) || (v0[0] < Max && v[0] > Max)){
                
                B.push_back(t);
                //std::cout << "Size B: " << B.size() << '\n';
                //fprintf(fB, "%f\n", t);
                //printf("\ntime %f\n",t);
                //printf("\nFalse \n");
            }
            
            else if((v0[0] > Min && v[0] < Min) || (v0[0] < Min && v[0] > Min)){
                
                C.push_back(t);
                //fprintf(fC, "%f\n", t);
            }

    
            for (i = 0; i<N; i++){
                
                v0[i] = v[i];
            }

            fprintf(fX,"%f\t",t);
            fprintf(fX,"%f\n",v0[0]);
            //fprintf(fX,"%f\n",v0[1]);

        }
    
    
         printf("\nFinish Simulation:\n");
    

    vec Buffer(0.,0.);
    
    vec Diff1(0.,0.);
    vec Diff2(0.,0.);

    int ka=0;
    int kb=0;
    int kc=0;
    
    int MA;
    int MB;
    int MC;
    MA = A.size();
    MB = B.size();
    MC = C.size();
    
    printf("\n sizaA:%d\n", MA);
    printf("\n sizaB:%d\n", MB);
    printf("\n sizaC:%d\n", MC);

   while(ka<MA && kb<MB && kc<MC){
  
       if(A[ka]==B[kb] or B[kb]==C[kc]){
            printf("Error");
        }
       else if(A[ka]<B[kb]){
             //printf("Error1");
           while(A[ka]<B[kb] && ka<MA){
               //printf("A:%f\n",A[ka]);
               Buffer.push_back(A[ka]);
               A.pop_back();
               //printf("\n sizaA:%d\n", A.size());
               //printf("\n sizaB:%d\n", B.size());
               ka = ka+1;
               //printf("ka:%d\n",ka);
                }
           
            if(A[ka]>B[kb]){
               for(j=0; j<Buffer.size(); j++){
                   Diff1.push_back(B[kb]-Buffer[j]);
                   //printf("\n sizeDiff1:%d\n", Diff1.size());
                   //printf("DiffA:%f \n",B[kb]-Buffer[j]);
                }
                Buffer.clear();
            }
           //printf("kb:%d\n",kb);
           kb = kb +1;
        }
       
        else if(C[kc]<B[kb]) {
            while(C[kc]<B[kb] && kc<MC){
                //printf("C:%f\n",C[kc]);
                Buffer.push_back(C[kc]);
                C.pop_back();
                //printf("\n sizaC:%d\n", C.size());
                //printf("\n sizaB:%d\n", B.size());
                kc = kc+1;
                //printf("kc:%d\n",kc);
            }
            
            if(C[kc]>B[kb]){
                for(j=0; j<Buffer.size(); j++){
                    Diff2.push_back(B[kb]-Buffer[j]);
                    //printf("\n sizeDiff2:%d\n", Diff2.size());
                    //printf("DiffC:%f \n",B[kb]-Buffer[j]);
                }
                Buffer.clear();
            }
        //printf("kb:%d\n",kb);
        kb = kb +1;
        }
        else{
            kb = kb+1;
            
        }
    }


    
    for (i = 0; i<Diff1.size(); i++){
        //printf("Diff:%f \n",Diff[i]);
        //fprintf(fdiff1, "%f\n", Diff1[i]);
    }
    
    for (i = 0; i<Diff2.size(); i++){
        //printf("Diff:%f \n",Diff[i]);
        //fprintf(fdiff2, "%f\n", Diff2[i]);
    }
    
    printf("size Diff1:%lu\n",Diff1.size());
    printf("size Diff2:%lu\n",Diff2.size());
    
    
    //close file
   fclose(fX);
   //fclose(fdiff1);
   //fclose(fdiff2);


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
vec Evolution(vec V, double Rand, int dim, vec Tau, double B, double U){
    vec DV(dim);
    int i;
    double Sum=0;
    double F = -U*((4/B)*(V[0]*V[0]*V[0])+(3-3/B)*(V[0]*V[0]) -(2+2/B)*V[0]-1+(1/B)); //Derivation of the potential
    
    //Evolution
    DV[2] = - Tau[1]*(V[2] + V[1] - Rand);
    
    DV[0] = V[1];
    DV[1] = (Tau[0])*(V[2]+ F);
    return DV;
}



