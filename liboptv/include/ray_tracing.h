
void modu(double a[3], double *m);
void norm_cross(double a[3], double b[3], double *n1, double *n2, double *n3);

void ray_tracing (double x,double y,Exterior Ex,Interior I,mm_np mm,\
double *Xb2,double *Yb2,double *Zb2, double *a3, double *b3,double *c3);

void point_line_line(Exterior Ex0, Interior I0, Glass G0, mm_np mm, double gX0, 
double gY0, double gZ0, double a0, double b0, double c0, Exterior Ex1, 
Interior I1, Glass G1, double gX1, double gY1, double gZ1,double a1, double b1,
double c1, double * x,double * y,double * z);

void norm_cross(double a[3],double b[3],double * n1,double * n2,double * n3);
void dot(double a[3], double b[3],double * d);
void modu(double a[3], double * m);
void ray_tracing_v2 (double x,double y,Exterior Ex,Interior I,Glass G,mm_np mm,\
double *Xb2,double *Yb2,double *Zb2, double *a3,double *b3,double *c3);
	                 
