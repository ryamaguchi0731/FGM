//
//  Multiple runs_haploid modurality.cpp
//  Fisher's geometrical model
//
//  Created by Ryo Yamaguchi on 2018/07/09.
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
    uniform_smallint<> dst( x, y );
    variate_generator< mt19937&, uniform_smallint<> > rand( gen, dst );
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
    normal_distribution<> dst( x, y );
    variate_generator< lagged_fibonacci1279&, normal_distribution<> > rand( gen, dst );
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
    
    for(int run=0; run<50; run++){
        
        cout << run << endl;
        
        sprintf(fname1, "Modurality_Haploid_sphere=1_trait=10_alpha=1.0_theta=0.04_pop1_%d.dat", run );
        sprintf(fname2, "Modurality_Haploid_sphere=1_trait=10_alpha=1.0_theta=0.04_pop2_%d.dat", run );
        ofstream file1( fname1 ) ;
        ofstream file2( fname2 ) ;
    
    int nums=1; //number of spheres
    int spheremut; // a sphere chosen for a mutation process
    int numt=10; //number of trait
    int numgen=500; //number of generations
    double sig=1.0; // parameter of fitness function
    double k=2.0; // parameter of fitness function (degree of epistasis)
    double s=0.0; //selection coefficient (initialized as 0)
    double popsize=2000.0; // population size
    double r; // mutation size
    double sum;
    double fit,fittemp;
    
    //parameter for mutation size
    double alpha=2.0; //shape parameter of gamma
    double theta=0.04; //scale parameter of gamma
    
    //parameter for exponential distribution
    
    vector<vector<double>> trait(nums,vector<double>(numt,0.0));  // initial value of trait
    vector<vector<double>> traittemp(nums,vector<double>(numt,0.0));  // initial value of trait
    vector<vector<double>> traitopt(nums,vector<double>(numt,0.0));  // optimum value of trait is 0
    vector<vector<double>> trait2(nums,vector<double>(numt,0.0));  // initial value of trait
    vector<vector<double>> traittemp2(nums,vector<double>(numt,0.0));  // initial value of trait
    vector<vector<double>> traitopt2(nums,vector<double>(numt,0.0));  // optimum value of trait is 0
    
    vector<double> totdif(nums,0.0);  //total sqared difference from the optimum
    vector<double> totdiftemp(nums,0.0);
    vector<double> distopt(nums,0.0);//distance from the optimum
    vector<double> distopttemp(nums,0.0);
    vector<double> fitness(nums,0.0);
    vector<double> fitnesstemp(nums,0.0);
    vector<double> mutvalue(numt,0.0);  // mutation value for each trait (initialized as 0)
    vector<double> totdif2(nums,0.0);  //total sqared difference from the optimum
    vector<double> totdiftemp2(nums,0.0);
    vector<double> distopt2(nums,0.0);//distance from the optimum
    vector<double> distopttemp2(nums,0.0);
    vector<double> fitness2(nums,0.0);
    vector<double> fitnesstemp2(nums,0.0);
    vector<double> mutvalue2(numt,0.0);  // mutation value for each trait (initialized as 0)
    
    /* fixed optimum setting */
    traitopt[0][0]=2.0;
    traitopt2[0][0]=2.0;
    
    file1 << "generation " << "s " << "r ";
    for (int j=0; j<nums; j++) {
        for (int i=0; i<numt; i++) {
            file1 << "trait[" << j << "]" << "[" << i << "] ";
        }
    }
    file1 << endl;
    
    file2 << "generation " << "s " << "r ";
    for (int j=0; j<nums; j++) {
        for (int i=0; i<numt; i++) {
            file2 << "trait[" << j << "]" << "[" << i << "] ";
        }
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
            spheremut=UniformRandom_int(0, nums-1);
            
            for (int j=0; j<nums; j++) {
                for (int i=0; i<numt; i++) {
                    traittemp[j][i]=0.0;
                    if (j==spheremut) {
                        traittemp[j][i]=trait[j][i]+r*(mutvalue[i]/sqrt(sum));
                    }
                    else{
                        traittemp[j][i]=trait[j][i];
                    }
                }
            }
            
            /* selection phase */
            for (int j=0; j<nums; j++) {
                totdif[j]=0.0;
                totdiftemp[j]=0.0;
            }
            
            for (int j=0; j<nums; j++) {
                for (int i=0; i<numt; i++) {
                    totdif[j]=totdif[j]+pow((trait[j][i]-traitopt[j][i]), 2); //squared distance before mutation
                    totdiftemp[j]=totdiftemp[j]+pow((traittemp[j][i]-traitopt[j][i]), 2); //squared distance after mutation
                }
            }
            
            for (int j=0; j<nums; j++) {
                distopt[j]=sqrt(totdif[j]);
                distopttemp[j]=sqrt(totdiftemp[j]);
                fitness[j]=Getfitness(sig, k, distopt[j]);
                fitnesstemp[j]=Getfitness(sig, k, distopttemp[j]);
            }
            
            fit=1.0;
            fittemp=1.0;
            for (int j=0; j<nums; j++) {
                fit=fit*fitness[j];
                fittemp=fittemp*fitnesstemp[j];
            }
            
            s=fittemp/fit-1.0;
            
            if (UniformRandom()<(1.0 - exp(-2.0*s))/(1.0 - exp(-2.0*popsize*s))) {
                
                file1 << loop+1 << " " << s << " " << r;
                for (int j=0; j<nums; j++) {
                    for (int i=0; i<numt; i++) {
                        file1 << " " << traittemp[j][i]-trait[j][i] << " ";
                        trait[j][i]=traittemp[j][i];
                    }
                }
                file1 << endl;
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
            spheremut=UniformRandom_int(0, nums-1);
            
            for (int j=0; j<nums; j++) {
                for (int i=0; i<numt; i++) {
                    traittemp2[j][i]=0.0;
                    if (j==spheremut) {
                        traittemp2[j][i]=trait2[j][i]+r*(mutvalue2[i]/sqrt(sum));
                    }
                    else{
                        traittemp2[j][i]=trait2[j][i];
                    }
                }
            }
            
            /* selection phase */
            for (int j=0; j<nums; j++) {
                totdif2[j]=0.0;
                totdiftemp2[j]=0.0;
            }
            
            for (int j=0; j<nums; j++) {
                for (int i=0; i<numt; i++) {
                    totdif2[j]=totdif2[j]+pow((trait2[j][i]-traitopt2[j][i]), 2); //squared distance before mutation
                    totdiftemp2[j]=totdiftemp2[j]+pow((traittemp2[j][i]-traitopt2[j][i]), 2); //squared distance after mutation
                }
            }
            
            for (int j=0; j<nums; j++) {
                distopt2[j]=sqrt(totdif2[j]);
                distopttemp2[j]=sqrt(totdiftemp2[j]);
                fitness2[j]=Getfitness(sig, k, distopt2[j]);
                fitnesstemp2[j]=Getfitness(sig, k, distopttemp2[j]);
            }
            
            fit=1.0;
            fittemp=1.0;
            for (int j=0; j<nums; j++) {
                fit=fit*fitness2[j];
                fittemp=fittemp*fitnesstemp2[j];
            }
            
            s=fittemp/fit-1.0;
            
            if (UniformRandom()<(1.0 - exp(-2.0*s))/(1.0 - exp(-2.0*popsize*s))) {
                
                file2 << loop+1 << " " << s << " " << r;
                for (int j=0; j<nums; j++) {
                    for (int i=0; i<numt; i++) {
                        file2 << " " << traittemp2[j][i]-trait2[j][i] << " ";
                        trait2[j][i]=traittemp2[j][i];
                        }
                }
                file2 << endl;
                ct++;
            }
            
        }
        
    }
    
    
    }

        file1.close();
        file2.close();
        
    }
    
    return 0;
}
