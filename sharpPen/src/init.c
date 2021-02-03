/*
  *  Part of R package sharpPen
*  Copyright (C) D. Wang (2020)
*
  *  Unlimited use and distribution (see LICENCE).
*/
  
  #include <stddef.h>
  #include <R.h>
  #include <Rinternals.h>
  #include <R_ext/Rdynload.h>
  
  void F77_SUB(blkest46)(double *x, double *y, int *n, int *q, int *qq,
                       int *nval, double *xj, double *yj, double *coef,
                       double *xmat, double *wk, double *qraux,
                       double *sigsqe, double *th44e, double *th46e);

void F77_SUB(cpsd)(double *x, double *y, int *n,
                  int *nmax, double *rss, double *xj,
                 double *yj, double *coef, double *xmat, double *wk,
                 double *qraux, double *cpvals);

void F77_SUB(linbinsd)(double *x, int *n, double *a,
                     double *b, int *m, double *gcounts);


void F77_SUB(rlbinsd)(double *x, double *y, int *n,
                    double *a, double *b, int *m, double *xcounts,
                    double *ycounts);

void F77_SUB(sdiagsd)(double *xcounts, double *delta,
                    double *hdisc, int *lvec, int *indic, int *midpts,
                    int *m,  double *fkap, int *ipp, int *ippp,
                    double *ss, double *smat, double *work, double *et,
                    int *ipvt, double *sd);

void F77_SUB(sstdgsd)(double *xcounts, double *delta,
                    double *hdisc, int *lvec, int *indic, int *midpts,
                    int *m,  double *fkap, int *ipp, int *ippp,
                    double *ss, double *uu, double *smat, double *umat,
                    double *work, double *det, int *ipvt, double *sstd);

static const R_FortranMethodDef FortEntries[] = {
  {"blkest46", (DL_FUNC) &F77_SUB(blkest46), 15},
  {"cpsd",     (DL_FUNC) &F77_SUB(cpsd),     12},
  {"linbinsd", (DL_FUNC) &F77_SUB(linbinsd),  6},
  {"rlbinsd",  (DL_FUNC) &F77_SUB(rlbinsd),   8},
  {"sdiagsd",  (DL_FUNC) &F77_SUB(sdiagsd),  16},
  {"sstdgsd",  (DL_FUNC) &F77_SUB(sstdgsd),  18},
  {NULL, NULL, 0}
};


void R_init_sharpPen(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, NULL, FortEntries, NULL);
  R_useDynamicSymbols(dll, FALSE);
}

