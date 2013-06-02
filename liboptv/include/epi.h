#ifndef EPI_H
#define EPI_H

#include <math.h>
#include <stdio.h>
#include <string.h>


#include "calibration.h"
#include "tracking_frame_buf.h"
#include "parameters.h"
#include "ray_tracing.h"
#include "typedefs.h"
#include "lsqadj.h"

typedef struct {
  int  	pnr;
  double  tol, corr;
} candidate;

double epi_line(double xl, double yl, Exterior Ex1, Interior I1, Glass G1,
    Exterior Ex2, Interior I2, Glass G2);
    
int  epi_mm(double xl, double yl, Exterior Ex1, Interior I1, Glass G1,
    Exterior Ex2, Interior I2, Glass G2, mm_np mmp, volume_par *vpar,
    control_par *cp, double *xmin, double *ymin, double *xmax, double *ymax);
    
int  epi_mm_2D(double xl, double yl, Exterior Ex1, Interior I1, Glass G1,
    mm_np mmp, volume_par *vpar, double *xout, double *yout, double *zout);
    
void find_candidate_plus_msg (coord_2d crd[], target pix[], int num, double xa,\
double ya,double xb,double yb, int n, int nx, int ny, int sumg,\
candidate cand[], int *count, int i12, volume_par *vpar, control_par *cp, \
Interior I[], ap_52 ap[]);
    
void find_candidate_plus (coord_2d crd[], target pix[], int num, double xa, \
double ya, double xb, double yb, int n, int nx, int ny, int sumg,\
candidate cand[], int *count, int nr, volume_par *vpar, control_par *cp, \
Interior *I, ap_52 ap[]);
    


void correct_brown_affin (double x, double y, ap_52 ap, double *x1, double *y1);

void img_xy_mm_geo (double X,double Y,double Z, Exterior Ex, Interior I, 
Glass G, mm_np mm, double *x, double *y, int n_img);


void  multimed_nlay_v2 (Exterior ex, Exterior ex_o, mm_np mm, double X,\
double Y, double Z,double *Xq, double *Yq, int n_img, Exterior *Ex);

void back_trans_Point(double X_t,double Y_t,double Z_t,mm_np mm,Glass G,\
double cross_p[3], double cross_c[3], double *X, double *Y, double *Z);


void trans_Cam_Point(Exterior ex, mm_np mm,Glass gl,\
double X,double Y, double Z,Exterior *ex_t,double *X_t,\
double *Y_t, double *Z_t,double *cross_p[3], double *cross_c[3]);

#endif
