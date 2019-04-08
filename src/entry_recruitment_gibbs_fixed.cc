////////////////////////////////////////////////////////////////////
// entry_recruitment_gibbs_fixed.cc
//
// entry_recruitment_gibbs_fixed.cc samples from the posterior distribution of a
// Poisson linear regression model
//
////////////////////////////////////////////////////////////////////
//
// Original code by Ghislain Vieilledent, June 2011
// CIRAD UR B&SEF
// ghislain.vieilledent@cirad.fr / ghislainv@gmail.com
//
////////////////////////////////////////////////////////////////////
//
// This software is distributed under the terms of the GNU GENERAL
// PUBLIC LICENSE Version 2, June 1991.  See the package LICENSE
// file for more information.
//
// Copyright (C) 2011 Ghislain Vieilledent
// 
////////////////////////////////////////////////////////////////////
//
// Revisions: 
// - G. Vieilledent, on 9 Sept 2011
//
////////////////////////////////////////////////////////////////////


// Scythe libraries
#include "matrix.h"
#include "distributions.h"
#include "stat.h"
#include "la.h"
#include "ide.h"
#include "smath.h"
#include "rng.h"
#include "mersenne.h"
// R libraries
#include <R.h> // needed to use Rprintf()
#include <R_ext/Utils.h> // needed to allow user interrupts

using namespace scythe;
using namespace std;

extern "C"{

    /* Gibbs sampler function */
    void entry_recruitment_gibbs_fixed (
	
	// Constants and data
	const int *ngibbs, const int *nthin, const int *nburn, // Number of iterations, burning and samples
	const int *nobs, // Constants
	const int *np, // Number of fixed and random covariates
	const double *Y_vect, // Observed response variable
	const double *X_vect, // Covariate for fixed effects
	const double *Int_vect, // Time interval (in year)
	const double *SCell_vect, // Cell surface (in m2)
        // Parameters to save
	double *beta_vect, // Fixed effects
	double *V, // Variance of residuals
	// Defining priors
	const double *mubeta_vect, const double *Vbeta_vect,
	const double *s1_V, const double *s2_V,
	// Diagnostic
	double *Deviance,
	double *log_lambda_pred, // log(Annual recruitment rate)
	// Seeds
	const int *seed,
	// Verbose
	const int *verbose,
        // Overdispersion
	const int *FixOD
	
	) {
	
	////////////////////////////////////////////////////////////////////////////////
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// Defining and initializing objects

	///////////////////////////
	// Redefining constants //
	const int NGIBBS=ngibbs[0];
	const int NTHIN=nthin[0];
	const int NBURN=nburn[0];
	const int NSAMP=(NGIBBS-NBURN)/NTHIN;
	const int NOBS=nobs[0];
	const int NP=np[0];

	///////////////
	// Constants //

	// Small fixed matrices indexed on k for data access
	Matrix<double> *Xi_arr = new Matrix<double>[NOBS];
	for(int i=0; i<NOBS; i++) {
	    Xi_arr[i] = Matrix<double>(1,NP);
	    for (int p=0; p<NP; p++) {
		Xi_arr[i](p)=X_vect[p*NOBS+i];
	    }
	} 

	////////////////////
	// Priors objects //
	Matrix<double> mubeta(NP,1,mubeta_vect);
	Matrix<double> Vbeta(NP,NP,Vbeta_vect);

	/////////////////////////////////////
	// Initializing running parameters //
	Matrix<double> beta_run(NP,1,false); // Unicolumn matrix of fixed effects
	for (int p=0; p<NP; p++) {
	    beta_run(p)=beta_vect[p*NSAMP];
	}
	double V_run=V[0];
	Matrix<double> *log_lambdai_run = new Matrix<double>[NOBS];
	for (int i=0;i<NOBS;i++) {
	    log_lambdai_run[i] = Matrix<double>(1,1);
	    log_lambdai_run[i](0)=log_lambda_pred[i];
	    log_lambda_pred[i]=0; // We reinitialize lambda_pred to zero to compute mean of MCMC
	}
	double Deviance_run=Deviance[0];

        ////////////////////////////////////////////////////////////
	// Proposal variance and acceptance for adaptive sampling //
	double *sigmap = new double[NOBS];
	double *nA = new double[NOBS];
	for (int n=0; n<NOBS; n++) {
	    nA[n]=0;
	    sigmap[n]=1;
	}
	double *Ar = new double[NOBS]; // Acceptance rate 
	for (int n=0; n<NOBS; n++) {
	    Ar[n]=0;
	}

	////////////
	// Message//
	Rprintf("\nRunning the Gibbs sampler. It may be long, please keep cool :)\n\n");
	R_FlushConsole();
	//R_ProcessEvents(); for windows	

	/////////////////////
	// Set random seed //
	mersenne myrng;
	myrng.initialize(*seed);
	
	///////////////////////////////////////////////////////////////////////////////////////
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// Gibbs sampler (see Chib et Carlin, 1999 p.4 "Blocking for Gaussian mixed models")

	for (int g=0; g<NGIBBS; g++) {

            ////////////////////////////////////////////////
            // vector log_lambda : Metropolis algorithm // 

	    for (int i=0; i<NOBS; i++) {

		// Mean of the prior
		Matrix<double> log_lambda_hat=Xi_arr[i]*beta_run;

		// Proposal
		double log_lambda_prop=myrng.rnorm(log_lambdai_run[i](0),sigmap[i]);
		
		// lambda_prim
		double lambda_prim_prop=exp(log_lambda_prop)*Int_vect[i]*SCell_vect[i];
		double lambda_prim_run=exp(log_lambdai_run[i](0))*Int_vect[i]*SCell_vect[i];
		
		// Ratio of probabilities
		double p_prop=log(dpois(Y_vect[i],lambda_prim_prop))+
		    log(dnorm(log_lambda_prop,log_lambda_hat(0),sqrt(V_run)));
		
		double p_now=log(dpois(Y_vect[i],lambda_prim_run))+
		    log(dnorm(log_lambdai_run[i](0),log_lambda_hat(0),sqrt(V_run)));
		
		double r=exp(p_prop-p_now); // ratio
		double z=myrng.runif();
		
		// Actualization
		if (z < r) {
		    log_lambdai_run[i](0)=log_lambda_prop;
		    nA[i]++;
		}
	    }


	    //////////////////////////////////
	    // vector beta: Gibbs algorithm //

	    // invVi, sum_V and sum_v
	    Matrix<double> sum_V(NP,NP);
	    Matrix<double> sum_v(NP,1);
	    for (int i=0; i<NOBS; i++) {
	    	sum_V += crossprod(Xi_arr[i]);
	    	sum_v += t(Xi_arr[i])*log_lambdai_run[i];
	    }

	    // big_V
	    Matrix<double> big_V=invpd(invpd(Vbeta)+sum_V/V_run);
	    
	    // small_v
	    Matrix<double> small_v=invpd(Vbeta)*mubeta+sum_v/V_run;

	    // Draw in the posterior distribution
	    beta_run=myrng.rmvnorm(big_V*small_v,big_V);


	    ////////////////
	    // variance V //
	    
	    // e
	    Matrix<double> e(1,1,true,0.0);
	    for (int i=0; i<NOBS; i++) {
	    	e+=crossprod(log_lambdai_run[i]-Xi_arr[i]*beta_run);
	    }

	    // Parameters
	    double S1=*s1_V+(NOBS/2); //shape
	    double S2=*s2_V+0.5*e(0); //rate

	    // Draw in the posterior distribution
	    if (FixOD[0]==1) {
	    	V_run=V[0];
	    }
	    else {
	    	V_run=1/myrng.rgamma(S1,S2);
	    }


	    //////////////////////////////////////////////////
	    //// Deviance

	    // logLikelihood
	    double logLk=0;
	    for (int i=0; i<NOBS; i++) {
		// lambda_prim_run
		double lambda_prim_run=exp(log_lambdai_run[i](0))*Int_vect[i]*SCell_vect[i];
		// L
		logLk+=log(dpois(Y_vect[i],lambda_prim_run));
	    }

	    // Deviance
	    Deviance_run=-2*logLk;


	    //////////////////////////////////////////////////
	    // Output
	    if(((g+1)>NBURN) && (((g+1)%(NTHIN))==0)) {
		int isamp=((g+1)-NBURN)/(NTHIN);
		for (int p=0; p<NP; p++) {
		    beta_vect[p*NSAMP+(isamp-1)]=beta_run(p);
		}
		V[isamp-1]=V_run;
		Deviance[isamp-1]=Deviance_run;
	    	for (int i=0;i<NOBS;i++){
		    Matrix<double> log_lambda_hat=Xi_arr[i]*beta_run;
		    log_lambda_pred[i]+=log_lambda_hat(0)/NSAMP; // We compute the mean of NSAMP values
	    	}
	    }
	    

	    ///////////////////////////////////////////////////////
	    // Adaptive sampling (on the burnin period)
	    int DIV=0;
	    if (NGIBBS >=1000) DIV=100;
	    else DIV=NGIBBS/10;
	    if((g+1)%DIV==0 && (g+1)<=NBURN){
	    	for (int n=0; n<NOBS; n++) {
		    const double ropt=0.44;
	    	    Ar[n]=nA[n]/(DIV*1.0);
	    	    if(Ar[n]>=ropt) sigmap[n]=sigmap[n]*(2-(1-Ar[n])/(1-ropt));
	    	    else sigmap[n]=sigmap[n]/(2-Ar[n]/ropt);
	    	    nA[n]=0.0; // We reinitialize the number of acceptance to zero
	    	}
	    }
	    if((g+1)%DIV==0 && (g+1)>NBURN){
	    	for (int n=0; n<NOBS; n++) {
	    	    Ar[n]=nA[n]/(DIV*1.0);
	    	    nA[n]=0.0; // We reinitialize the number of acceptance to zero
	    	}
	    }

    
	    //////////////////////////////////////////////////
	    // Progress bar
	    double Perc=100*(g+1)/(NGIBBS);
	    if(((g+1)%(NGIBBS/100))==0 && (*verbose==1)){  
	    	Rprintf("*");
	    	R_FlushConsole();
	    	//R_ProcessEvents(); for windows
	    	if(((g+1)%(NGIBBS/10))==0){
		    double mAr=0; // Mean acceptance rate
		    for (int n=0; n<NOBS; n++) {
			mAr+=Ar[n]/(NOBS);
		    }
		    Rprintf(":%.1f%%, mean accept. rate=%.3f\n",Perc,mAr);
      	    	    R_FlushConsole();
	    	    //R_ProcessEvents(); for windows
	    	}
	    } 

            //////////////////////////////////////////////////
	    // User interrupt
	    R_CheckUserInterrupt(); // allow user interrupt 
	    
	} // Gibbs sampler


	///////////////
	// Delete memory allocation (see new)
	delete[] sigmap;
	delete[] nA;
	delete[] Xi_arr;
	delete[] log_lambdai_run;
	delete[] Ar;
	
    } //  end entry_recruitment_gibbs_fixed function

} // end extern "C"

////////////////////////////////////////////////////////////////////
// END
////////////////////////////////////////////////////////////////////
