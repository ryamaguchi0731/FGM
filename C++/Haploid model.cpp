//
//  Multiple runs_haploid.cpp
//  Fisher's geometrical model
//
//  Created by Ryo Yamaguchi on 2018/07/07.
//  Copyright © 2018年 Ryo Yamaguchi. All rights reserved.
//

#include <iostream>
#include <stdio.h>
#include <fstream>
#include <time.h>
#include <vector>
#include <numeric>
#include <algorithm>
#include <boost/random.hpp>
#include <boost/progress.hpp>
#include <cstdio>

using namespace std;
using namespace boost;

double UniformRandom() //uniform random variable [0,1]
{
    static mt19937 gen( static_cast<unsigned int>(time(0)) );
    static uniform_real<> dst(0,1);
    static variate_generator< mt19937&, uniform_real<> > rand(gen, dst);
    return rand();
}

int UniformRandom_int(int x, int y) //uniform random integer [x,y]
{
    static mt19937 gen( static_cast<unsigned int>(time(0)) );
    static uniform_smallint<> dst( x, y );
    static variate_generator< mt19937&, uniform_smallint<> > rand( gen, dst );
    return rand();
}

double StandardNormalRandom() //standard normal random variables (mean=0,SD=1)
{
    static lagged_fibonacci1279 gen( static_cast<unsigned int>(time(0)) );
    static normal_distribution<> dst( 0, 1 );
    static variate_generator< lagged_fibonacci1279&, normal_distribution<> > rand( gen, dst );
    return rand();
}

double NormalRandom(double x, double y) //standard normal random variables (mean=x,SD=y)
{
    static lagged_fibonacci1279 gen( static_cast<unsigned int>(time(0)) );
    static normal_distribution<> dst( x, y );
    static variate_generator< lagged_fibonacci1279&, normal_distribution<> > rand( gen, dst );
    return rand();
}

double ExponentialRandom(double x) //random variable from exponential distribution (mean=1/x,SD=1/x)
{
    static mt19937 gen( static_cast<unsigned int>(time(0)) );
    static exponential_distribution<> dst(x);
    static variate_generator< mt19937&, exponential_distribution<> > rand(gen, dst);
    return rand();
}

double GammaRandom(double x, double y) //x:shape parameter、y:scale parameter (mean=x*y, variance=x*y^2)
{
    static lagged_fibonacci1279 gen( static_cast<unsigned int>(time(0)) );
    gamma_distribution<> dst( x, y );
    variate_generator< lagged_fibonacci1279&, gamma_distribution<> > rand( gen, dst );
    return rand();
}

double Getfitness(double x, double y, double z) //fitness function
{
    double fitness;
    fitness=exp(-x*pow(z,y));
    return(fitness);
}

int main(int argc, const char * argv[]) {
    
    progress_timer t;
    char fname1[1000];
    char fname2[1000];
    
    for(double mutsize=0.04; mutsize<=0.04; mutsize=mutsize+0.04){
    
    for(int run=0; run<50; run++){
    
        cout << mutsize << " " << run << endl;
    
    sprintf(fname1, "Scenario1.1_Haploid_r=%.2f_pop1_%d.dat", mutsize, run );
    sprintf(fname2, "Scenario1.1_Haploid_r=%.2f_pop2_%d.dat", mutsize, run );
    ofstream file1( fname1 ) ;
    ofstream file2( fname2 ) ;
    

    int numt=10; //number of trait
    int numgen=500; //number of generations
    double totdif=0.0,totdiftemp=0.0; //total sqared difference from the optimum
    double distopt=0.0,distopttemp=0.0; //distance from the optimum
    double totdif2=0.0,totdiftemp2=0.0; //total sqared difference from the optimum in pop2
    double distopt2=0.0,distopttemp2=0.0; //distance from the optimum in pop2
    double q=1.0; // parameter of fitness function
    double k=2.0; // parameter of fitness function (degree of epistasis)
    double s=0.0; //selection coefficient (initialized as 0)
    double popsize=2000.0; // population size
    double r; // mutation size
    double sum;
        
        if(mutsize==0.04){
            numgen=500;
        }
        else if (mutsize==0.08){
            numgen=350;
        }
        else{
            numgen=250;
        }
        
    
    //parameter for mutation size
    double alpha=1.0; //shape parameter of gamma
    double theta=mutsize; //scale parameter of gamma
    
    //condition of moving optimum
    //if you do not want to move the optimum, set tau1>numgen
    int numtmov=1; //number of moving traits
    int tau1=10000; //generation when optima start moving
    int tau2=0; //generation when optima stop moving
        
        //Scenario2.1 Slow
        //tau1=0;
        //tau2=200;
        
        //Scenario2.2 Medium
        //tau1=0;
        //tau2=100;
        
        //Scenario2.3 Fast
        //tau1=0;
        //tau2=40;
        
    double nextopt=2.0; // optima value after movement (endpoint)
    double speed=nextopt/(double)(tau2-tau1); //speed of environnmental change
    if (speed>0.0) {
        cout << "spped is " << speed <<endl;
    }
    
    vector<double> trait(numt,0.0);  // initial value of trait
    vector<double> traittemp(numt,0.0);  // initial value of trait
    vector<double> trait2(numt,0.0);  // initial value of trait in population 2
    vector<double> traittemp2(numt,0.0);  // initial value of trait in population 2
    vector<double> traitopt1(numt,0.0);  // optimum 1
    vector<double> traitopt2(numt,0.0);  // optimum 2
    vector<double> mutvalue(numt,0.0);  // mutation value for each trait (initialized as 0)
    vector<double> mutvalue2(numt,0.0);  // mutation value for each trait (initialized as 0)
    vector<double> nextopttemp(numt,0.0); // temporal optima value (initialized as 0)
    
    /* fixed optimum setting */
    //    Scenario 1.1
        traitopt1[0]=2.0;
        traitopt2[0]=2.0;
    
    //    Scenario 3.1 (pi/6)
    //    traitopt1[0]=-sqrt(2.0-sqrt(3.0));
    //    traitopt2[0]=sqrt(2.0-sqrt(3.0));
    //    traitopt1[1]=-sqrt(2.0+sqrt(3.0));
    //    traitopt2[1]=-sqrt(2.0+sqrt(3.0));
    
    //    Scenario 3.2 (pi/4)
    //    traitopt1[0]=-sqrt(2.0-sqrt(2.0));
    //    traitopt2[0]=sqrt(2.0-sqrt(2.0));
    //    traitopt1[1]=-sqrt(2.0+sqrt(2.0));
    //    traitopt2[1]=-sqrt(2.0+sqrt(2.0));
    
    //    Scenario 3.3 (pi/2)
    //    traitopt1[0]=-sqrt(2.0);
    //    traitopt2[0]=sqrt(2.0);
    //    traitopt1[1]=-sqrt(2.0);
    //    traitopt2[1]=-sqrt(2.0);
    
    //    Scenario 3.4 (pi)
    //traitopt1[0]=-2.0;
    //traitopt2[0]=2.0;
    
    file1 << "generation " << "s " << "r " << "distopt ";
    for (int i=0; i<numt; i++) {
        file1 << "trait[" << i << "] ";
    }
    file1 << endl;
    
    file2 << "generation " << "s " << "r " << "distopt ";
    for (int i=0; i<numt; i++) {
        file2 << "trait[" << i << "] ";
    }
    file2 << endl;
        
        
    //Main loop
    for (int loop=0; loop<numgen; loop++) {
        int ct=0;
        
        while (ct<1) {
            
            if (UniformRandom()<0.5) {
                
                /*Population 1*/
                
                /* mutation phase */
                sum=0.0;
                for (int i=0; i<numt; i++) {
                    mutvalue[i]=StandardNormalRandom();
                    sum=sum+pow(mutvalue[i],2);
                }
                
                r=GammaRandom(alpha, theta);
                
                for (int i=0; i<numt; i++) {
                    traittemp[i]=trait[i]+r*(mutvalue[i]/sqrt(sum));
                }
                
                /* selection phase */
                if (loop<=tau1) {  //without changing optima
                    totdif=0.0;
                    totdiftemp=0.0;
                    for (int i=0; i<numt; i++) {
                        totdif=totdif+pow((trait[i]-traitopt1[i]), 2); //squared distance before mutation
                        totdiftemp=totdiftemp+pow((traittemp[i]-traitopt1[i]), 2); //squared distance after mutation
                    }
                }
                else if ((tau1<loop)&&(loop<=tau2)){  //optima start changing
                    totdif=0.0;
                    totdiftemp=0.0;
                    
                    for (int i=0; i<numt; i++) {
                        if (i<numtmov){
                            totdif=totdif+pow((trait[i]-(traitopt1[i]+speed*(double)(loop-tau1))),2); //squared distance before mutation
                            totdiftemp=totdiftemp+pow((traittemp[i]-(traitopt1[i]+speed*(double)(loop-tau1))),2); //squared distance after mutation
                        }
                        else{
                            totdif=totdif+pow((trait[i]-traitopt1[i]), 2); //squared distance before mutation
                            totdiftemp=totdiftemp+pow((traittemp[i]-traitopt1[i]), 2); //squared distance after mutation
                        }
                    }
                }
                else{
                    totdif=0.0;
                    totdiftemp=0.0;
                    for (int i=0; i<numt; i++) {
                        if (i<numtmov){
                            totdif=totdif+pow((trait[i]-(traitopt1[i]+nextopt)),2); //squared distance before mutation
                            totdiftemp=totdiftemp+pow((traittemp[i]-(traitopt1[i]+nextopt)),2); //squared distance after mutation
                        }
                        else{
                            totdif=totdif+pow((trait[i]-traitopt1[i]), 2); //squared distance before mutation
                            totdiftemp=totdiftemp+pow((traittemp[i]-traitopt1[i]), 2); //squared distance after mutation
                        }
                    }
                }
                
                distopt=sqrt(totdif);
                distopttemp=sqrt(totdiftemp);
                
                s=Getfitness(q, k, distopttemp)/Getfitness(q, k, distopt)-1.0;
                
                if (UniformRandom()<(1.0 - exp(-2.0*s))/(1.0 - exp(-2.0*popsize*s))) {
                    
                    file1 << loop+1 << " " << s << " " << r << " " << distopt;
                    for (int i=0; i<numt; i++) {
                        file1 << " " << traittemp[i]-trait[i] << " ";
                    }
                    file1 << endl;
                    
                    for (int i=0; i<numt; i++) {
                        trait[i]=traittemp[i];
                    }
                    ct++;
                    
                }
                
            }
            
            
            else{
                /*Population 2*/
                
                /* mutation phase */
                sum=0.0;
                for (int i=0; i<numt; i++) {
                    mutvalue2[i]=StandardNormalRandom();
                    sum=sum+pow(mutvalue2[i],2);
                }
                
                r=GammaRandom(alpha, theta);
                
                for (int i=0; i<numt; i++) {
                    traittemp2[i]=trait2[i]+r*(mutvalue2[i]/sqrt(sum));
                }
                
                /* selection phase */
                if (loop<=tau1) {  //without changing optima
                    totdif2=0.0;
                    totdiftemp2=0.0;
                    for (int i=0; i<numt; i++) {
                        totdif2=totdif2+pow((trait2[i]-traitopt2[i]), 2); //squared distance before mutation
                        totdiftemp2=totdiftemp2+pow((traittemp2[i]-traitopt2[i]), 2); //squared distance after mutation
                    }
                }
                else if ((tau1<loop)&&(loop<=tau2)){  //optima start changing
                    totdif2=0.0;
                    totdiftemp2=0.0;

                    for (int i=0; i<numt; i++) {
                        if (i<numtmov){
                            totdif2=totdif2+pow((trait2[i]-(traitopt2[i]+speed*(double)(loop-tau1))),2); //squared distance before mutation
                            totdiftemp2=totdiftemp2+pow((traittemp2[i]-(traitopt2[i]+speed*(double)(loop-tau1))),2); //squared distance after mutation
                        }
                        else{
                            totdif2=totdif2+pow((trait2[i]-traitopt2[i]), 2); //squared distance before mutation
                            totdiftemp2=totdiftemp2+pow((traittemp2[i]-traitopt2[i]), 2); //squared distance after mutation
                        }
                    }
                }
                else{
                    totdif2=0.0;
                    totdiftemp2=0.0;
                    for (int i=0; i<numt; i++) {
                        if (i<numtmov){
                            totdif2=totdif2+pow((trait2[i]-(traitopt2[i]+nextopt)),2); //squared distance before mutation
                            totdiftemp2=totdiftemp2+pow((traittemp2[i]-(traitopt2[i]+nextopt)),2); //squared distance after mutation
                        }
                        else{
                            totdif2=totdif2+pow((trait2[i]-traitopt2[i]), 2); //squared distance before mutation
                            totdiftemp2=totdiftemp2+pow((traittemp2[i]-traitopt2[i]), 2); //squared distance after mutation
                        }
                    }
                }
                
                distopt2=sqrt(totdif2);
                distopttemp2=sqrt(totdiftemp2);
                
                s=Getfitness(q, k, distopttemp2)/Getfitness(q, k, distopt2)-1.0;
                
                if (UniformRandom()<(1.0 - exp(-2.0*s))/(1.0 - exp(-2.0*popsize*s))) {
                    
                    file2 << loop+1 << " " << s << " " << r << " " << distopt2;
                    for (int i=0; i<numt; i++) {
                        file2 << " " << traittemp2[i]-trait2[i] << " ";
                    }
                    file2 << endl;
                    
                    for (int i=0; i<numt; i++) {
                        trait2[i]=traittemp2[i];
                    }
                    
                    ct++;
                }
                
            }
            
        }
        
    }
     
        file1.close();
        file2.close();
        
    }
        
    }
    
    return 0;
}
