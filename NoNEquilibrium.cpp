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
vec Evolution(vec V, double Rand, int dim, vec Tau, double U);


int main(int argc, char** argv){
    
    //Print time
    double time_counter = 0;
    clock_t this_time = clock();
    clock_t last_time = this_time;
    
    int i, j, z;
    double q;
    
    int N=4;
    
    //Random number
    double r;
    
    const float e = 2.71828183;
    
    time_t start, end;
    int elapsed_time, h, min, sec;
    start = time(NULL);
    
    vec v(4,0.);
    vec v0(4,0.);
    vec dv(4,0.);
    double rand;
    vec tau(3,0.);   //tau[0] = tauD/tauM tau[1] = tauD/tauV  tau[2] = tauD/tauR

    tau[0] = 1./(0.1);
    tau[1] = 1./(1.);
    tau[2] = 1./(1.);
    double U0 = 3.;
    double u=4*U0;
    
    //open file
    //FILE *fdiff;
    //fdiff = fopen("DiffTm01T108TG1U6.txt", "w");
    //printf("\n DiffTm01T108TG1U6\n");
    
    //FILE *fX;
    //fX = fopen("XTm01T103TG01.txt", "w");
    //printf("\n XTm01T103TG01\n");
    
    FILE *fV;
    fV = fopen("VTm01T11TG1L3.txt", "w");
    printf("\n VTm01T11TG1L3\n");
    

   ////////////////////////Time
    
    const double T1 = (5000.)*(pow(e,U0))*((pow(U0,-1)/(1 + (10)*U0*pow(tau[1],-1))) + pow(U0,-1)*pow(tau[0],-1) + 2*sqrt(pow(U0,-1)*pow(tau[0],-1)) + e*(pow(tau[1],-1)*pow(tau[1],-1)));
    
    const double T2 = (5000.)*(pow(e,U0))*((pow(U0,-1)/(1 + (10)*U0*pow(tau[2],-1))) + pow(U0,-1)*pow(tau[0],-1) + 2*sqrt(pow(U0,-1)*pow(tau[0],-1)) + e*(pow(tau[2],-1)*pow(tau[2],-1)));
  
    
    double t_final = 150;//fmax(T1,T2);
    
    printf("\nt_final: %f\n",t_final);
    
    double t_step;
    double t;
    
    //////////////////////////////////////////////////////////////////////////
    //Simulations
    //////////////////////////////////////////////////////////////////////////
  
        t = 0;
        t_step = 0.1;

        //check on t_step (tauD/tauM)*(1/t_step)<<1 & (tauD/tauG)*(1/t_step)<<1
        for (i = 0; i<(3); i++){
            printf("\n Tau: %f\n",1./tau[i]);
            if((tau[i])*(t_step)>0.01){
                t_step = 0.01*(1/tau[i]);
            }
        }
     printf("\n U_0=%f \n",u/4);
    
        r = std::sqrt(2./t_step);
    
        std::random_device rd;
        std::mt19937 gen(rd());
    
        std::normal_distribution<> Gauss(0, r);
        std::normal_distribution<> Gauss0(0, 1);
        
     printf("\n tstep: %f\n",t_step);
    
    double step = t_final/t_step;
        printf("\n step: %f\n",step);

        //int step = 100;
        vec A(0, 0.);
        vec B(0, 0.);
    
        //initial condition
        for (i = 0; i<N; i++){
            v0[i] = Gauss0(gen);
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
                
            vec K1(N,0.);
            vec K2(N,0.);
            vec K3(N,0.);
            vec K4(N,0.);
            vec v1(N,0.);
            //std::cout<<"step\n"<<j;
            
            //K1
            K1 = Evolution(v0, rand, N, tau,u);
            
            //K2
            for (i = 0; i<N; i++){
                
                v1[i] = v0[i] + (t_step*K1[i])/2.;
                // std::cout<<"\nv1\n"<<v1[i]<<std::endl;
            }
            K2 = Evolution(v1, rand, N, tau,u);
            
            //K3
            for (i = 0; i<N; i++){
                
                v1[i] = v0[i] + (t_step*K2[i])/2.;
            }
            K3 = Evolution(v1, rand, N, tau,u);
            
            //K4
            for (i = 0; i<N; i++){
                
                v1[i] = v0[i] + t_step*K3[i];
            }
            K4 = Evolution(v1, rand, N, tau,u);

            
            for (i = 0; i<N; i++){
                
                v[i] = v0[i] + t_step*(K1[i] + 2*K2[i] + 2*K3[i] + K4[i])/6.;
            }
            
            //printf("v:%f\n",v[1]);
            if ((v0[0] > 1 && v[0] < 1) || (v0[0] < 1 && v[0] > 1)){
                A.push_back(t);
            }
            
            else if((v0[0] > -1 && v[0] < -1) || (v0[0] < -1 && v[0] > -1)){
             B.push_back(t);
            }

    
            for (i = 0; i<N; i++){
                
                v0[i] = v[i];
            }

            //if (q>100 && q<6000){
                //fprintf(fX,"%f\n",v[1]);
            //}
            //fprintf(fX, "%f\t", t);
            //fprintf(fX, "%f\n", v[0]);
            fprintf(fV, "%f\n", v[1]);
    
        }
    
         printf("\nFinish Simulation:\n");

    vec Buffer(0.,0.);
    vec Diff(0.,0.);
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
               Buffer.push_back(A[ka]);

               A.pop_back();
               //printf("\n sizaA:%d\n", A.size());
               //printf("\n sizaB:%d\n", B.size());
               
               ka = ka+1;
                //printf("ka:%d\n",ka);

                }
          
            if(A[ka]>B[kb]){
               for(j=0; j<Buffer.size(); j++){
                    
                   Diff.push_back(B[kb]-Buffer[j]);
                   //printf("\n sizaDiff:%d\n", Diff.size());
                    //printf("DiffA:%f \n",B[kb]-Buffer[j]);
              
                }
                Buffer.clear();
            }
          
        }
       
        else{
           
            while(B[kb]<A[ka] && kb<MB){
                //printf("B:%f\n",B[kb]);
                Buffer.push_back(B[kb]);
                
                B.pop_back();
                //printf("\n sizaA:%d\n", A.size());
                //printf("\n sizaB:%d\n", B.size());
                
                kb = kb+1;
                //printf("kb:%d\n",kb);
            }
            
            if(B[kb]>A[ka]){
                for(j=0; j<Buffer.size(); j++){
                    Diff.push_back(A[ka]-Buffer[j]);
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
   //fclose(fX);
   fclose(fV);
   //fclose(fdiff);
  


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
vec Evolution(vec V, double Rand, int dim, vec Tau,double U){
    vec DV(dim);
    int i;
    double F = -U*(V[0]*V[0]*V[0] - V[0]); //Derivation of the potential
    
    //Evolution
    DV[2] = -Tau[1]*(V[2] + V[1]);
    DV[3] = -Tau[2]*(V[3] - Rand);
    
    DV[1] = (Tau[0])*(V[2] + V[3] + F);
    DV[0] = V[1];
    
    return DV;
}



