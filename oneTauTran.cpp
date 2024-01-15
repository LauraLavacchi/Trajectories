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
vec Evolution(vec V, double Rand, int dim, vec Tau);


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
    double r0;
    //double r1;
    
    double rand;
    
    const int U0 = 3.;
    const float e = 2.71828183;
    
    
    time_t start, end;
    int elapsed_time, h, min, sec;
    start = time(NULL);
    
    vec v(N,0.);
    vec v0(N,0.);
    vec dv(N,0.);
  
    vec gamma(NExp,0.);
    vec tau(Ntau,0.);   //tau[0] = tauD/tauM tau[1] = tauD/tauG1  tau[2] = tauD/tauG2

    
    tau[0] = 1./(0.001);
    tau[1] = 1./(1);

    //open file
    //FILE *fdiff;
    //fdiff = fopen("DifTm001G001.txt", "w");
    //printf("\n DiffTm0G001.txt \n");
    FILE *fX;
    fX = fopen("XM0001T1L.txt", "w");
    
    
   ////////////////////////Time
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
    
        
        r0 = std::sqrt((2./t_step));
    
        std::random_device rd;
        std::mt19937 gen(rd());
    
        std::normal_distribution<> Gauss0(0, r0);
        
     printf("\n tstep: %f\n",t_step);
    
    double step = (5)*pow(10,7);
        printf("\n step: %f\n",step);

        //int step = 100;
        vec A(0, 0.);
        vec B(0, 0.);
        

        
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
            
                rand = Gauss0(gen);
            
            vec K1(N,0.);
            vec K2(N,0.);
            vec K3(N,0.);
            vec K4(N,0.);
            vec v1(N,0.);
            //std::cout<<"step\n"<<j;
            
            //K1
            K1 = Evolution(v0, rand, N, tau);
            
            //K2
            for (i = 0; i<N; i++){
                
                v1[i] = v0[i] + (t_step*K1[i])/2.;
                // std::cout<<"\nv1\n"<<v1[i]<<std::endl;
            }
            K2 = Evolution(v1, rand, N, tau);
            
            //K3
            for (i = 0; i<N; i++){
                
                v1[i] = v0[i] + (t_step*K2[i])/2.;
            }
            K3 = Evolution(v1, rand, N, tau);
            
            //K4
            for (i = 0; i<N; i++){
                
                v1[i] = v0[i] + t_step*K3[i];
            }
            K4 = Evolution(v1, rand, N, tau);
            
            
            for (i = 0; i<N; i++){
                
                v[i] = v0[i] + t_step*(K1[i] + 2*K2[i] + 2*K3[i] + K4[i])/6.;
            }
            
            
            //printf("x:%f\n",v[0]);
            //printf("v:%f\n",v[1]);
            
            if(q>pow(10,7)){
                if ((v0[0] > 0 && v[0] <0) || (v0[0] < 0 && v[0] > 0)){
                    A.push_back(t);
                //std::cout << "Size A: " << A.size() << '\n';
                
                //fprintf(fA, "%f\n", t);
                }
            
                else if((v0[0] > -1 && v[0] < -1) || (v0[0] < -1 && v[0] > -1)){
                
                    B.push_back(t);
                //std::cout << "Size B: " << B.size() << '\n';
             
                //fprintf(fB, "%f\n", t);
                //printf("\ntime %f\n",t);
                //printf("\nFalse \n");
                }
            }

    
            for (i = 0; i<N; i++){
                
                v0[i] = v[i];
            }

            if(q>pow(10,7)){
                fprintf(fX,"%f\t",t);
                fprintf(fX,"%f\n",v0[0]);
            //printf("x:%f\n",v0[0]);
            }

        }
    
    
         printf("\nFinish Simulation:\n");
    

    vec Buffer(0.,0.);
    vec Diff(0.,0.);
    double Summ=0;
    int Norm=0;
    int ka=0;
    int kb=0;
    
    int MA;
    int MB;
    MA = A.size();
    MB = B.size();
    
    printf("\n sizaA:%d\n", MA);
    printf("\n sizaB:%d\n", MB);
    
    //printf("size A:%lu\n",A.size());
    //printf("size B:%lu\n",B.size());
   while(ka<MA && kb<MB){
  
       if(A[ka]==B[kb]){
            printf("Error");
        }
    
       else if(A[ka]<B[kb]){
             //printf("Error1");
           while(A[ka]<B[kb] && ka<MA){
                //printf("A:%f\n",A[ka]);
               //printf("\n sizaA:%d\n", A.size());
               //printf("\n sizaB:%d\n", B.size());
               
               ka = ka+1;
                //printf("ka:%d\n",ka);

                }
          
            if(A[ka]>B[kb]){
                    
                   Diff.push_back(B[kb]-A[ka-1]);
                   //printf("\n sizaDiff:%d\n", Diff.size());
                    //printf("DiffA:%f \n",B[kb]-Buffer[j]);
              
            }
          
        }
       
        else{
           
            while(B[kb]<A[ka] && kb<MB){
                //printf("B:%f\n",B[kb]);
                //printf("\n sizaA:%d\n", A.size());
                //printf("\n sizaB:%d\n", B.size());
                
                kb = kb+1;
                //printf("kb:%d\n",kb);
            }
            
            if(B[kb]>A[ka]){
                for(j=0; j<Buffer.size(); j++){
                    //:Diff.push_back(A[ka]-Buffer[j]);
                   //printf("\n sizaDiff:%d\n", Diff.size());
                      //printf("DiffB:%f \n",A[ka]-Buffer[j]);
                }
                Buffer.clear();
            }
        }
    }


    for (i = 0; i<Diff.size(); i++){
        //printf("Diff:%f \n",Diff[i]);
        //fprintf(fdiff, "%f\n", Diff[i]);
        
    }
    
    printf("size Diff:%lu\n",Diff.size());
    
    
    //close file
   //fclose(fdiff);
   fclose(fX);


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
vec Evolution(vec V, double Rand, int dim, vec Tau){
    vec DV(dim);
    double F = -12*(V[0]*V[0]*V[0] - V[0]); //Derivation of the potential
    
    //Evolution

    DV[0] = V[1];
    DV[1] = (Tau[0])*(V[2]+ F);
    DV[2] = - Tau[1]*(V[2] + V[1] - Rand);
    return DV;
}



