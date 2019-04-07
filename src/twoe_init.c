#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .C calls */
extern void entry_growth_gibbs_fixed(const int *ngibbs, const int *nthin, const int *nburn, const int *nobs, const int *np, const double *Y_vect, const double *X_vect, double *beta_vect, double *V, const double *mubeta_vect, const double *Vbeta_vect, const double *s1_V, const double *s2_V, double *Deviance, double *Y_pred, const int *seed, const int *verbose);
extern void entry_growth_gibbs_mixed(const int *ngibbs, const int *nthin, const int *nburn, const int *nobs, const int *ngroup,	const int *np, const int *nq, const double *Y_vect, const double *W_vect, double *beta_vect, double *b_vect, double *Vb_vect, double *V, const double *mubeta_vect, const double *Vbeta_vect, const double *r, const double *R_vect, const double *s1_V, const double *s2_V, double *Deviance, double *Y_pred, const int *seed, const int *verbose);
extern void entry_mortality_gibbs_fixed(const int *ngibbs, const int *nthin, const int *nburn, const int *nobs, const int *np, const double *Y_vect, const double *X_vect, const double *Int_vect, double *beta_vect, double *V, const double *mubeta_vect, const double *Vbeta_vect, const double *s1_V, const double *s2_V, double *Deviance, double *logit_theta_pred, const int *seed, const int *verbose, const int *FixOD);	
extern void entry_mortality_gibbs_mixed(const int *ngibbs, const int *nthin, const int *nburn, const int *nobs, const int *ngroup, const int *np, const int *nq, const int *IdentGroup, const double *Y_vect, const double *X_vect, const double *W_vect, const double *Int_vect, double *beta_vect, double *b_vect, double *Vb_vect, double *V, const double *mubeta_vect, const double *Vbeta_vect, const double *r, const double *R_vect, const double *s1_V, const double *s2_V, double *Deviance, double *logit_theta_pred, const int *seed, const int *verbose, const int *FixOD);
extern void entry_ngrowth_gibbs_fixed(const int *ngibbs, const int *nthin, const int *nburn, const int *nobs, const int *np, const double *Y_vect, const double *X_vect, const double *D_vect, const double *a_sd, const double *b_sd, double *beta_vect, double *V, double *Y_true, const double *mubeta_vect, const double *Vbeta_vect, const double *s1_V, const double *s2_V, double *Deviance, double *log_Y_pred, const int *seed, const int *verbose);
extern void entry_ngrowth_gibbs_mixed(const int *ngibbs, const int *nthin, const int *nburn, const int *nobs, const int *ngroup, const int *np, const int *nq, const int *IdentGroup, const double *Y_vect, const double *X_vect, const double *W_vect, const double *D_vect, const double *a_sd, const double *b_sd, double *beta_vect, double *b_vect, double *Vb_vect, double *V, double *Y_true, const double *mubeta_vect, const double *Vbeta_vect, const double *r, const double *R_vect, const double *s1_V, const double *s2_V, double *Deviance, double *log_Y_pred, const int *seed, const int *verbose);
extern void entry_recruitment_gibbs_fixed(const int *ngibbs, const int *nthin, const int *nburn, const int *nobs, const int *np, const double *Y_vect, const double *X_vect, const double *Int_vect,const double *SCell_vect,double *beta_vect,double *V,const double *mubeta_vect, const double *Vbeta_vect, const double *s1_V, const double *s2_V, double *Deviance, double *lambda_pred, const int *seed, const int *verbose, const int *FixOD);
extern void entry_recruitment_gibbs_mixed(const int *ngibbs, const int *nthin, const int *nburn, const int *nobs, const int *ngroup, const int *np, const int *nq, const int *IdentGroup, const double *Y_vect, const double *X_vect, const double *W_vect, const double *Int_vect, const double *SCell_vect, double *beta_vect, double *b_vect, double *Vb_vect, double *V, const double *mubeta_vect, const double *Vbeta_vect, const double *r, const double *R_vect, const double *s1_V, const double *s2_V, double *Deviance, double *lambda_pred, const int *seed, const int *verbose, const int *FixOD);

static const R_CMethodDef CEntries[] = {
    {"entry_growth_gibbs_fixed",      (DL_FUNC) &entry_growth_gibbs_fixed,      17},
    {"entry_growth_gibbs_mixed",      (DL_FUNC) &entry_growth_gibbs_mixed,      25},
    {"entry_mortality_gibbs_fixed",   (DL_FUNC) &entry_mortality_gibbs_fixed,   19},
    {"entry_mortality_gibbs_mixed",   (DL_FUNC) &entry_mortality_gibbs_mixed,   27},
    {"entry_ngrowth_gibbs_fixed",     (DL_FUNC) &entry_ngrowth_gibbs_fixed,     21},
    {"entry_ngrowth_gibbs_mixed",     (DL_FUNC) &entry_ngrowth_gibbs_mixed,     29},
    {"entry_recruitment_gibbs_fixed", (DL_FUNC) &entry_recruitment_gibbs_fixed, 20},
    {"entry_recruitment_gibbs_mixed", (DL_FUNC) &entry_recruitment_gibbs_mixed, 28},
    {NULL, NULL, 0}
};

void R_init_twoe(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}