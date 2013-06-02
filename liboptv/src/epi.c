#include <math.h>
#include <stdio.h>
#include <string.h>

#include "epi.h"
#include "lsqadj.h"

int dumbbell_pyptv;



double epi_line (double xl, double yl, Exterior Ex1, Interior I1, Glass G1, \
Exterior Ex2, Interior I2, Glass G2){

  int i,j;
  double m2;
  double vect1[3], vect2[3], vect3[3], nk[3], n2[3],
    p1l[3], K2[3], k2[3], D2t[3][3];

  /* base O1 -> O2 */
  vect1[0] = Ex2.x0 - Ex1.x0;
  vect1[1] = Ex2.y0 - Ex1.y0;
  vect1[2] = Ex2.z0 - Ex1.z0;

  /* coordinates of arbitrary point P1 in image1 */
  p1l[0] = xl;  p1l[1] = yl;	p1l[2] = - I1.cc;

  /* beam O1 -> P in space */
  matmul (vect2, (double *) Ex1.dm, p1l, 3,3,1);

  /* normale to epipolar plane */
  crossprod (vect1,vect2,nk);

  /* normale to image2 */
  vect3[0] = 0;	vect3[1] = 0;	vect3[2] = - I2.cc;

  /* normale to image 2, in space */
  matmul (n2, (double *) Ex2.dm, vect3, 3,3,1);

  /* epipolar line in image2, in space */
  crossprod (nk,n2,K2);

  /* epipolar line in image2 */
  for (i=0; i<3; i++)  for (j=0; j<3; j++)  D2t[i][j] = Ex2.dm[j][i];
  matmul (k2, (double *) D2t, K2, 3,3,1);
  m2 = k2[1] / k2[0];
  return (m2);
}



int epi_mm (double x1, double y1, Exterior Ex1, Interior I1, Glass G1, \
Exterior Ex2, Interior I2, Glass G2, mm_np mmp, volume_par *vpar, \
control_par *cp, double *xmin,\
double *ymin, double *xmax, double *ymax){
  /*  ray tracing gives the point of exit and the direction
      cosines at the waterside of the glass;
      min. and max. depth give window in object space,
      which can be transformed into _2 image
      (use img_xy_mm because of comparison with img_geo)  */

  double a, b, c, xa,ya,xb,yb;
  double X1,Y1,Z1, X, Y, Z;
  double Zmin, Zmax;

  //ray_tracing    (x1,y1, Ex1, I1,     mmp, &X1, &Y1, &Z1, &a, &b, &c);
  ray_tracing_v2 (x1,y1, Ex1, I1, G1, mmp, &X1, &Y1, &Z1, &a, &b, &c);

  /* calculate min and max depth for position (valid only for one setup) */
  Zmin = vpar->Zmin_lay[0]
    + (X1 - vpar->X_lay[0]) * (vpar->Zmin_lay[1] - vpar->Zmin_lay[0]) / 
    (vpar->X_lay[1] - vpar->X_lay[0]);
  Zmax = vpar->Zmax_lay[0]
    + (X1 - vpar->X_lay[0]) * (vpar->Zmax_lay[1] - vpar->Zmax_lay[0]) / 
    (vpar->X_lay[1] - vpar->X_lay[0]);

  Z = Zmin;   X = X1 + (Z-Z1) * a/c;   Y = Y1 + (Z-Z1) * b/c;
  //img_xy_mm_geo_old (X,Y,Z, Ex2, I2,     mmp, &xa, &ya);
  img_xy_mm_geo     (X,Y,Z, Ex2, I2, G2, mmp, &xa, &ya, cp->num_cams);

  Z = Zmax;   X = X1 + (Z-Z1) * a/c;   Y = Y1 + (Z-Z1) * b/c;
  //img_xy_mm_geo_old (X,Y,Z, Ex2, I2,     mmp, &xb, &yb);
  img_xy_mm_geo(X,Y,Z, Ex2, I2, G2, mmp, &xb, &yb, cp->num_cams);

  /*  ==> window given by xa,ya,xb,yb  */

  *xmin = xa;  *ymin = ya;  *xmax = xb;  *ymax = yb;

  return (0);
}

int epi_mm_2D (double x1, double y1, Exterior Ex1, Interior I1, Glass G1, \
mm_np mmp, volume_par *vpar, double *xp, double *yp, double *zp){
  /*  ray tracing gives the point of exit and the direction
      cosines at the waterside of the glass;
      min. and max. depth give window in object space,
      which can be transformed into _2 image
      (use img_xy_mm because of comparison with img_geo)  */

  double a, b, c;
  double X1,Y1,Z1,X,Y,Z;
  
  double Zmin, Zmax;

  ray_tracing_v2 ( x1,y1, Ex1, I1, G1, mmp, &X1, &Y1, &Z1, &a, &b, &c);

  /* calculate min and max depth for position (valid only for one setup) */
  Zmin = vpar->Zmin_lay[0]
    + (X1 - vpar->X_lay[0]) * (vpar->Zmin_lay[1] - vpar->Zmin_lay[0]) / 
    (vpar->X_lay[1] - vpar->X_lay[0]);
  Zmax = vpar->Zmax_lay[0]
    + (X1 - vpar->X_lay[0]) * (vpar->Zmax_lay[1] - vpar->Zmax_lay[0]) /
    (vpar->X_lay[1] - vpar->X_lay[0]);

  Z = 0.5*(Zmin+Zmax);   
  X = X1 + (Z-Z1) * a/c;   
  Y = Y1 + (Z-Z1) * b/c;
  
  *xp=X; *yp=Y; *zp=Z;

  return (0);
}

void find_candidate_plus (coord_2d crd[], target pix[], int num, double xa, \
double ya, double xb, double yb, int n, int nx, int ny, int sumg,\
candidate cand[], int *count, int nr, volume_par *vpar, control_par *cp, \
Interior *I, ap_52 ap[]){
/*  binarized search in a x-sorted coord-set, exploits shape information  */
  register int	j;
  int dummy;
  int	       	j0, dj, p2;
  double      	m, b, d, temp, qn, qnx, qny, qsumg, corr;
  double       	xmin, xmax, ymin, ymax,particle_size;
  int dumbbell=0;
  double tol_band_width;
  
//Beat Mai 2010 for dumbbell

    if (dumbbell_pyptv==1){    dumbbell=1;
  }

  if (dumbbell==0){
	
	  /////here is new Beat version of April 2010
	  if (nx>ny) particle_size=nx;
	  else       particle_size=ny;
	  tol_band_width = vpar->eps0*0.5*(cp->pix_x + cp->pix_y)*particle_size;
  }
  else{
      tol_band_width = vpar->eps0;
  }
  if(tol_band_width<0.05){
       tol_band_width=0.05;
  }

  /* define sensor format for search interrupt */
  xmin = (-1) * cp->pix_x * cp->imx/2;	
  xmax = cp->pix_x * cp->imx/2;
  ymin = (-1) * cp->pix_y * cp->imy/2;	
  ymax = cp->pix_y * cp->imy/2;
  xmin -= I[nr].xh;	
  ymin -= I[nr].yh;
  xmax -= I[nr].xh;	
  ymax -= I[nr].yh;
  correct_brown_affin (xmin,ymin, ap[nr], &xmin,&ymin);
  correct_brown_affin (xmax,ymax, ap[nr], &xmax,&ymax);

  for (j=0; j<4; j++)
    {
      cand[j].pnr = -999;  cand[j].tol = -999;  cand[j].corr = -999;
    }


  /* line equation: y = m*x + b */
  if (xa == xb)	xa += 1e-10;
  m = (yb-ya)/(xb-xa);  b = ya - m*xa;

  if (xa > xb)
    {
      temp = xa;  xa = xb;  xb = temp;
    }
  if (ya > yb)
    {
      temp = ya;  ya = yb;  yb = temp;
    }

  if ( (xb>xmin) && (xa<xmax) && (yb>ymin) && (ya<ymax))  /* sensor area */
    {
      /* binarized search for start point of candidate search */
      for (j0=num/2, dj=num/4; dj>1; dj/=2)
	{
	  if (crd[j0].x < (xa - tol_band_width))  j0 += dj;
	  else  j0 -= dj;
	}
      j0 -= 12;  if (j0 < 0)  j0 = 0;		       	/* due to trunc */

      for (j=j0, *count=0; j<num; j++)			/* candidate search */
	{
	  if (crd[j].x > xb+tol_band_width)  return;		/* finish search */

	  if ((crd[j].y > ya-tol_band_width) && (crd[j].y < yb+tol_band_width))
	    {
	      if ((crd[j].x > xa-tol_band_width) && (crd[j].x < xb+tol_band_width))
		{
		  d = fabs ((crd[j].y - m*crd[j].x - b) / sqrt(m*m+1));
          
		  /* Beat: modified in April 2010 to allow for better treatment of 
		  //different sized traced particles, in particular colloids and tracers
		  if ( d < eps ){
          */
          /////here is new Beat version of April 2010
		   //if (nx>ny) particle_size=nx;
		   //else       particle_size=ny;
		   if ( d < tol_band_width ){
		   ///////end of new Beat version

		      p2 = crd[j].pnr;
		      if (n  < pix[p2].n)      	qn  = (double) n/pix[p2].n;
		      else		       	qn  = (double) pix[p2].n/n;
		      if (nx < pix[p2].nx)	qnx = (double) nx/pix[p2].nx;
		      else		       	qnx = (double) pix[p2].nx/nx;
		      if (ny < pix[p2].ny)	qny = (double) ny/pix[p2].ny;
		      else		       	qny = (double) pix[p2].ny/ny;
		      if (sumg < pix[p2].sumg)
			        qsumg = (double) sumg/pix[p2].sumg;
		      else	qsumg = (double) pix[p2].sumg/sumg;

		      // empirical correlation coefficient
			  // from shape and brightness parameters 
		      corr = (4*qsumg + 2*qn + qnx + qny);
		      // create a tendency to prefer those matches
			  // with brighter targets 
		      corr *= ((double) (sumg + pix[p2].sumg));

		      if (qn >= vpar->cn && qnx >= vpar->cnx && \
                 qny >= vpar->cny && qsumg > vpar->csumg) {

				 if ( *count < maxcand) {
			        cand[*count].pnr = j;
			        cand[*count].tol = d;
			        cand[*count].corr = corr;
			        (*count)++;
		         } else {
			        dummy=(int)maxcand;
			        printf("in find_candidate_plus: count > maxcand\n");}
			     }
		      }
		   }
           
           
	    }
	}
    }

  else  *count = -1;	 	/* out of sensor area */
}


void find_candidate_plus_msg (coord_2d crd[], target pix[], int num, double xa,\
double ya,double xb,double yb, int n, int nx, int ny, int sumg,\
candidate cand[], int *count, int i12, volume_par *vpar, control_par *cp, \
Interior I[], ap_52 ap[])
{

/*  binarized search in a x-sorted coord-set, exploits shape information  */
/*  gives messages (in examination)  */

  register int	j;
  int	       	j0, dj, p2;
  double        m, b, d, temp, qn, qnx, qny, qsumg, corr;
  double       	xmin, xmax, ymin, ymax;
  double tol_band_width,particle_size;

  /* define sensor format for search interrupt */
  xmin = (-1) * cp->pix_x * cp->imx/2;	
  xmax = cp->pix_x * cp->imx/2;
  ymin = (-1) * cp->pix_y * cp->imy/2;	
  ymax = cp->pix_y * cp->imy/2;
  
  xmin -= I[i12].xh;	ymin -= I[i12].yh;
  xmax -= I[i12].xh;	ymax -= I[i12].yh;
  correct_brown_affin (xmin,ymin, ap[i12], &xmin,&ymin);
  correct_brown_affin (xmax,ymax, ap[i12], &xmax,&ymax);

  if (nx>ny) particle_size=nx;
  else       particle_size=ny;
  tol_band_width = vpar->eps0*0.5*(cp->pix_x + cp->pix_y)*particle_size;

  for (j=0; j<4; j++)
    {
      cand[j].pnr = -999;  cand[j].tol = 999;
    }
  m = (yb-ya)/(xb-xa);  b = ya - m*xa;   /* line equation: y = m*x + b */

  if (xa > xb)
    {
      temp = xa;  xa = xb;  xb = temp;
    }
  if (ya > yb)
    {
      temp = ya;  ya = yb;  yb = temp;
    }

  if ( (xb>xmin) && (xa<xmax) && (yb>ymin) && (ya<ymax)) /* sensor area */
    {
      /* binarized search for start point of candidate search */
      for (j0=num/2, dj=num/4; dj>1; dj/=2)
	{
	  if (crd[j0].x < (xa - tol_band_width))  j0 += dj;
	  else  j0 -= dj;
	}
      j0 -= 12;  if (j0 < 0)  j0 = 0;  	/* due to trunc */

      for (j=j0, *count=0; j<num; j++) 	/* candidate search */
	{
	  if (crd[j].x > xb+tol_band_width)  return;      	/* finish search */

	  if ((crd[j].y > ya-tol_band_width) && (crd[j].y < yb+tol_band_width))
	    {
	      if ((crd[j].x > xa-tol_band_width) && (crd[j].x < xb+tol_band_width))
		{
		  d = fabs ((crd[j].y - m*crd[j].x - b) / sqrt(m*m+1));
          if ( d < tol_band_width ){
		      p2 = crd[j].pnr;
		      if (n  < pix[p2].n)      	qn  = (double) n/pix[p2].n;
		      else		       	qn  = (double) pix[p2].n/n;
		      if (nx < pix[p2].nx)	qnx = (double) nx/pix[p2].nx;
		      else		       	qnx = (double) pix[p2].nx/nx;
		      if (ny < pix[p2].ny)	qny = (double) ny/pix[p2].ny;
		      else		       	qny = (double) pix[p2].ny/ny;
		      if (sumg < pix[p2].sumg)
			qsumg = (double) sumg/pix[p2].sumg;
		      else	qsumg = (double) pix[p2].sumg/sumg;


		      /* empirical correlation coefficient
			 from shape and brightness parameters */
		      corr = (4*qsumg + 2*qn + qnx + qny);
		      /* create a tendency to prefer those matches
			 with brighter targets */
		      corr *= ((double) (sumg + pix[p2].sumg));

            if (qn >= vpar->cn && qnx >= vpar->cnx && qny >= vpar->cny && 
                qsumg > vpar->csumg)
			{
			  if (*count>=maxcand)
			    { printf("More candidates than (maxcand): %d\n",*count); return; }
			  cand[*count].pnr = p2;
			  cand[*count].tol = d;
 			  cand[*count].corr = corr;
			  (*count)++;
			  printf ("%d %3.0f/%3.1f \n", p2, corr, d*1000);
			}
		    }
		}
	    }
	}
      if (*count == 0)  puts ("- - -");
    }
  else  *count = -1;		       	       /* out of sensor area */

}


/* copied from imgcoord.c */
void img_xy_mm_geo (double X,double Y,double Z, Exterior Ex, Interior I, 
Glass G, mm_np mm, double *x, double *y, int n_img){
  double deno;
  Exterior Ex_t;
  double X_t,Y_t,Z_t;
  double *cross_p, *cross_c;
  double Xh,Yh,Zh;

  trans_Cam_Point(Ex,mm,G,X,Y,Z,&Ex_t,&X_t,&Y_t,&Z_t,&cross_p,&cross_c);
  multimed_nlay_v2 (Ex_t,Ex,mm,X_t,Y_t,Z_t,&X_t,&Y_t, n_img, Ex);
  back_trans_Point(X_t,Y_t,Z_t,mm, G,cross_p,cross_c,&X,&Y,&Z);

  deno = Ex.dm[0][2] * (X-Ex.x0)
    + Ex.dm[1][2] * (Y-Ex.y0)
    + Ex.dm[2][2] * (Z-Ex.z0);

  *x = - I.cc *  (Ex.dm[0][0] * (X-Ex.x0)
		  + Ex.dm[1][0] * (Y-Ex.y0)
		  + Ex.dm[2][0] * (Z-Ex.z0)) / deno;

  *y = - I.cc *  (Ex.dm[0][1] * (X-Ex.x0)
		  + Ex.dm[1][1] * (Y-Ex.y0)
		  + Ex.dm[2][1] * (Z-Ex.z0)) / deno;
}

/* copied from multimed.c */
void  multimed_nlay_v2 (Exterior ex,Exterior ex_o, mm_np mm, double X,\
double Y, double Z,double *Xq, double *Yq, int n_img,\
Exterior *Ex){
  
  //Beat Lüthi, Nov 2007 comment actually only Xq is affected since all Y and Yq are always zero
  int		i, it=0;
  double	 beta1, beta2[32], beta3, r, rbeta, rdiff, rq, mmf;
  
  // interpolation in mmLUT, if selected (requires some global variables) 
  if (mm.lut)
    {    
      // check, which is the correct image 
      for (i=0; i<n_img; i++)
	if (Ex[i]->x0 == ex_o->x0  &&  Ex[i]->y0 == ex_o->y0  &&  Ex[i]->z0 == ex_o->z0)
	  break;
      
      mmf = get_mmf_from_mmLUT (i, X,Y,Z);
      
      if (mmf > 0)
	{
	  *Xq = ex.x0 + (X-ex.x0) * mmf;
	  *Yq = ex.y0 + (Y-ex.y0) * mmf;
	  return;
	}
    }
  
  // iterative procedure (if mmLUT does not exist or has no entry) 
  r = sqrt ((X-ex.x0)*(X-ex.x0)+(Y-ex.y0)*(Y-ex.y0));
  rq = r;
  
  do
    {
      beta1 = atan (rq/(ex.z0-Z));
      for (i=0; i<mm.nlay; i++)	beta2[i] = asin (sin(beta1) * mm.n1/mm.n2[i]);
      beta3 = asin (sin(beta1) * mm.n1/mm.n3);
      
      rbeta = (ex.z0-mm.d[0]) * tan(beta1) - Z * tan(beta3);
      for (i=0; i<mm.nlay; i++)	rbeta += (mm.d[i] * tan(beta2[i]));
      rdiff = r - rbeta;
      rq += rdiff;
      it++;
    }
  while (((rdiff > 0.001) || (rdiff < -0.001))  &&  it < 40);
  
  if (it >= 40)
    {
      *Xq = X; *Yq = Y;
      puts ("Multimed_nlay stopped after 40. Iteration");	return;
    }
    
  if (r != 0)
    {
      *Xq = ex.x0 + (X-ex.x0) * rq/r;
      *Yq = ex.y0 + (Y-ex.y0) * rq/r;
    }
  else
    {
      *Xq = X;
      *Yq = Y;
    }
	
}




/* from multimed.c */

void back_trans_Point(double X_t,double Y_t,double Z_t,mm_np mm, Glass G,\
double cross_p[3], double cross_c[3], double *X, double *Y, double *Z){
    
    double nVe,nGl;
	nGl=sqrt(pow(G.vec_x,2.)+pow(G.vec_y,2.)+pow(G.vec_z,2.));
	nVe=sqrt( pow(cross_p[0]-(cross_c[0]-mm.d[0]*G.vec_x/nGl),2.)
		     +pow(cross_p[1]-(cross_c[1]-mm.d[0]*G.vec_y/nGl),2.)
			 +pow(cross_p[2]-(cross_c[2]-mm.d[0]*G.vec_z/nGl),2.));
	

	*X=cross_c[0]-mm.d[0]*G.vec_x/nGl+X_t*(cross_p[0]-(cross_c[0]-mm.d[0]*G.vec_x/nGl))/nVe+Z_t*G.vec_x/nGl;
	*Y=cross_c[1]-mm.d[0]*G.vec_y/nGl+X_t*(cross_p[1]-(cross_c[1]-mm.d[0]*G.vec_y/nGl))/nVe+Z_t*G.vec_y/nGl;
	*Z=cross_c[2]-mm.d[0]*G.vec_z/nGl+X_t*(cross_p[2]-(cross_c[2]-mm.d[0]*G.vec_z/nGl))/nVe+Z_t*G.vec_z/nGl;

}



void trans_Cam_Point(Exterior ex, mm_np mm, Glass gl,\
double X, double Y, double Z, Exterior *ex_t, double *X_t,\
double *Y_t, double *Z_t,double *cross_p, double *cross_c){
  //--Beat Lüthi June 07: I change the stuff to a system perpendicular to the interface
  double dist_cam_glas,dist_point_glas,dist_o_glas; //glas inside at water 
  
  dist_o_glas = sqrt(gl.vec_x*gl.vec_x+gl.vec_y*gl.vec_y+gl.vec_z*gl.vec_z);
  
  dist_cam_glas =   ex.x0*gl.vec_x/dist_o_glas +\
                    ex.y0*gl.vec_y/dist_o_glas +\
                    ex.z0*gl.vec_z/dist_o_glas -\
                    dist_o_glas-mm.d[0];
                    
  dist_point_glas = X*gl.vec_x/dist_o_glas+\
                    Y*gl.vec_y/dist_o_glas+\
                    Z*gl.vec_z/dist_o_glas-\
                    dist_o_glas; 

  cross_c[0]=ex.x0-dist_cam_glas*gl.vec_x/dist_o_glas;
  cross_c[1]=ex.y0-dist_cam_glas*gl.vec_y/dist_o_glas;
  cross_c[2]=ex.z0-dist_cam_glas*gl.vec_z/dist_o_glas;
  cross_p[0]=X    -dist_point_glas*gl.vec_x/dist_o_glas;
  cross_p[1]=Y    -dist_point_glas*gl.vec_y/dist_o_glas;
  cross_p[2]=Z    -dist_point_glas*gl.vec_z/dist_o_glas;

  ex_t->x0 = 0.0;
  ex_t->y0 = 0.0;
  ex_t->z0 = dist_cam_glas + mm.d[0];

  *X_t=sqrt( pow(cross_p[0]-(cross_c[0]-mm.d[0]*gl.vec_x/dist_o_glas),2.)
		    +pow(cross_p[1]-(cross_c[1]-mm.d[0]*gl.vec_y/dist_o_glas),2.)
			+pow(cross_p[2]-(cross_c[2]-mm.d[0]*gl.vec_z/dist_o_glas),2.));
  *Y_t=0;
  *Z_t=dist_point_glas;
      
}


/* copied from trafo.c */
/*  correct crd to geo with Brown + affine  */
/* Arguments:
*   double x,y
*   additional orientation parameters (skewness, etc.) ap_52
* Returns:
*   corrected x1,y1 (double)
*/  
void correct_brown_affin (double x, double y, ap_52 ap, double *x1, double *y1){
  double  r, xq, yq;
	

  r = sqrt (x*x + y*y);
  if (r != 0)
    {
      xq = (x + y*sin(ap.she))/ap.scx
	- x * (ap.k1*r*r + ap.k2*r*r*r*r + ap.k3*r*r*r*r*r*r)
	- ap.p1 * (r*r + 2*x*x) - 2*ap.p2*x*y;
      yq = y/cos(ap.she)
	- y * (ap.k1*r*r + ap.k2*r*r*r*r + ap.k3*r*r*r*r*r*r)
	- ap.p2 * (r*r + 2*y*y) - 2*ap.p1*x*y;
    }
  r = sqrt (xq*xq + yq*yq);		/* one iteration */
  if (r != 0)
    {
      *x1 = (x + yq*sin(ap.she))/ap.scx
	- xq * (ap.k1*r*r + ap.k2*r*r*r*r + ap.k3*r*r*r*r*r*r)
	- ap.p1 * (r*r + 2*xq*xq) - 2*ap.p2*xq*yq;
      *y1 = y/cos(ap.she)
	- yq * (ap.k1*r*r + ap.k2*r*r*r*r + ap.k3*r*r*r*r*r*r)
	- ap.p2 * (r*r + 2*yq*yq) - 2*ap.p1*xq*yq;
    }
}

