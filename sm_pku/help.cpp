#include"help.h"
int phase(int n)
{
  return n%2?-1:1;
}

int delta(int a,int b)
{
  return (a==b)?1:0;
}

double cg(int jj1,int jj2,int jj3,int mm1,int mm2,int mm3)
{
  return gsl_sf_coupling_3j(jj1,jj2,jj3,mm1,mm2,-mm3)*sqrt(jj3+1.)*phase((jj1-jj2+mm3)/2);
}

int mapab(int a,int b,int dim)
{
  return dim*a-((a+3)*a)/2+b-1;
}


int mapab_e(int a,int b,int dim)
{
  return dim*a-(a*(a+1))/2+b;
}


int Cnr(int iN, int iR)
{
  if (iR < 0 || iR > iN)
    {
      return 0;
    }
  int iComb = 1;
  int i = 0;
  while (i < iR)
    {
      ++i;
      iComb *= iN - i + 1;
      iComb /= i;
    }
  return iComb;
}




complexd sqrt(complexd c)
{
  double r=std::abs(c);
  double phi=std::arg(c);
  phi/=2;
  return std::polar(sqrt(r),phi);
}

bool operator < (const complexd & c1,const complexd &c2)
{
  return c1.real()<c2.real();
}



double ho_k(int n,int l,double b,double k)
{
  double z=b*k;
  double temp=phase(n) * exp(-0.5*z*z) * pow(z,l) * sqrt(2.0) *\
    exp( 0.5*(gsl_sf_lngamma(n+1)-gsl_sf_lngamma(n+l+1.5)) )*\
    gsl_sf_laguerre_n(n,l+0.5,z*z)*pow(b,1.5);
  return temp;
}

double ho_r(int n,int l,double b,double r)
{
  double z=r/b;
  double temp=exp(-0.5*z*z) * pow(z,l) * sqrt(2.0) *\
    exp( 0.5*(gsl_sf_lngamma(n+1)-gsl_sf_lngamma(n+l+1.5)) )*\
    gsl_sf_laguerre_n(n,l+0.5,z*z)/pow(b,1.5);
  return temp;
}

double ho_f(double x,void *params)
{
  nlb*p=(nlb*)params;
  return ho_r(p->n1,p->l1,1,x)*ho_r(p->n2,p->l2,1,x)*x*x*pow(x,p->L);
}
double ho_L(int n1,int l1,int n2,int l2,int L)
{
  gsl_integration_workspace *w=gsl_integration_workspace_alloc(500);
  gsl_function F;
  nlb p;
  p.n1=n1;p.l1=l1;p.n2=n2;p.l2=l2;p.L=L;
  F.function=&ho_f;
  F.params=&p;
  double result=0,error=0;
  gsl_integration_qagiu(&F,0,1e-4,1e-4,500,w,&result,&error);
  gsl_integration_workspace_free(w);
  return result;
}



void gauleg(const double x1, const double x2, vector<double> &x, vector<double> &w)
{
  const double EPS=1.0e-14;
  double z1,z,xm,xl,pp,p3,p2,p1;
  int n=x.size();
  int m=(n+1)/2;
  xm=0.5*(x2+x1);
  xl=0.5*(x2-x1);
  for (int i=0;i<m;i++)
    {
      z=cos(Pi*(i+0.75)/(n+0.5));
      do {
	p1=1.0;
	p2=0.0;
	for (int j=0;j<n;j++)
	  {
	    p3=p2;
	    p2=p1;
	    p1=((2.0*j+1.0)*z*p2-j*p3)/(j+1);
	  }
	pp=n*(z*p1-p2)/(z*z-1.0);
	z1=z;
	z=z1-p1/pp;
      } while (abs(z-z1) > EPS);
      x[i]=xm-xl*z;
      x[n-1-i]=xm+xl*z;
      w[i]=2.0*xl/((1.0-z*z)*pp*pp);
      w[n-1-i]=w[i];
    }
}


double reducedYMatEle(int la,int jja,int lb,int jjb,int L)
{
  //parity check
  if( (la+lb+L)%2 ) return 0;
  return phase( (jja+1)/2 )*sqrt((jja+1)*(jjb+1)*(2*L+1)/(4.0*Pi))*gsl_sf_coupling_3j(jja,2*L,jjb,1,0,-1);
}




complexd csphbessel_jn(int n,complexd z)
{
  vcomplexd jn,yn,djn,dyn;
  csphbessel_jy(n+1,z,jn,yn,djn,dyn);
  return jn[n];
}
complexd csphbessel_yn(int n,complexd z)
{
  vcomplexd jn,yn,djn,dyn;
  csphbessel_jy(n,z,jn,yn,djn,dyn);
  return yn[n];
}
complexd csphbessel_djn(int n,complexd z)
{
  vcomplexd jn,yn,djn,dyn;
  csphbessel_jy(n,z,jn,yn,djn,dyn);
  return djn[n];
}
complexd csphbessel_dyn(int n,complexd z)
{
  vcomplexd jn,yn,djn,dyn;
  csphbessel_jy(n,z,jn,yn,djn,dyn);
  return dyn[n];
}
double envj(int n,double x)
{
  return 0.5*log10(6.28*n)-n*log10(1.36*x/n);
}
int msta1(double x,int mp)
{
  double a0=abs(x);
  int n0=int(1.1*a0)+1;
  double f0=envj(n0,a0)-mp;
  int n1=n0+5;
  double f1=envj(n1,a0)-mp;
  int nn;
  for(int it=1;it<=20;it++)
    {
      nn=n1-(n1-n0)/(1.-f0/f1);
      double f=envj(nn,a0)-mp;
      if(abs(nn-n1)<1) break;
      n0=n1;
      f0=f1;
      n1=nn;
      f1=f;
    }
  return nn;
}

int msta2(double x,int n,int mp)
{
  double a0=abs(x);
  double hmp=0.5*mp;
  double ejn=envj(n,a0);
  double obj;
  int n0;
  if(ejn <= hmp)
    {
      obj=mp;
      n0=int(1.1*a0)+1;
    }
  else
    {
      obj=hmp+ejn;
      n0=n;
    }
  double f0=envj(n0,a0)-obj;
  int n1=n0+5;
  double f1=envj(n1,a0)-obj;
  int nn;
  for(int it=1;it<=20;it++)
    {
      nn=n1-(n1-n0)/(1-f0/f1);
      double f=envj(nn,a0)-obj;
      if(abs(nn-n1) < 1) break;
      n0=n1;
      f0=f1;
      n1=nn;
      f1=f;
    }
  return nn+10;
}

void csphbessel_jy(const int n,complexd z,vcomplexd&jn,vcomplexd&yn,vcomplexd&djn,vcomplexd&dyn)
{
  jn.resize(n+1);
  yn.resize(n+1);
  djn.resize(n+1);
  dyn.resize(n+1);
  double a0=abs(z);
  int nmax=n;
  if(a0 < 1e-8)
    {
      for(int i=0;i<=n;i++)
	{
	  jn[i]=0.;
	  djn[i]=0.;
	  yn[i]=-1e300;
	  dyn[i]=1e300;
	  jn[0]=1.;
	  djn[1]=1./3;
	}
      return;
    }
  jn[0]=sin(z)/z;
  jn[1]=( jn[0]-cos(z) )/z;
  if(n >= 2)
    {
      complexd a=jn[0];
      complexd b=jn[1];
      int m=msta1(a0,200);
      if(m<n)
	nmax=m;
      else
	m=msta2(a0,n,15);
      complexd f0=0.;
      complexd f1=1e-100;
      complexd f;
      for(int k=m;k>=0;k--)
	{
	  f=(2.*k+3.)*f1/z-f0;
	  if(k<=nmax) jn[k]=f;
	  f0=f1;
	  f1=f;
	}
      complexd s;
      if(abs(a) > abs(b)) s=a/f;
      if(abs(a)<=abs(b)) s=b/f0;
      for(int k=0;k<=nmax;k++)
	jn[k]*=s;
    }
  djn[0]=(cos(z)-sin(z)/z)/z;
  for(int k=1;k<=nmax;k++)
    djn[k]=jn[k-1]-(k+1.0)*jn[k]/z;
  yn[0]=-cos(z)/z;
  yn[1]=( yn[0]-sin(z) )/z;
  dyn[0]=(sin(z)+cos(z)/z)/z;
  dyn[1]=(2.*dyn[0]-cos(z))/z;
  for(int k=2;k<=nmax;k++)
    {
      if(abs(jn[k-1]) > abs(jn[k-2]) )
	yn[k]=( jn[k]*yn[k-1]-1.0/(z*z) )/jn[k-1];
      else
	yn[k]=(jn[k]*yn[k-2]-(2.*k-1.)/pow(z,3))/jn[k-2];
    }
  for(int k=2;k<=nmax;k++)
    dyn[k]=yn[k-1]-(k+1.)*yn[k]/z;
}



complexd cmplx_laguerre(int n,double alpha,complexd z)
{
  double x=z.real();
  double y=z.imag();

  if(abs(y)<1e-8) return gsl_sf_laguerre_n(n,alpha,x);
  double u(0),v(0);
  //real part
  int summax=n/2;
  for(int i=0;i<=summax;i++)
    {
      u+=phase(i) * pow(y,2*i) * exp(-gsl_sf_lngamma(2*i+1) ) * gsl_sf_laguerre_n(n-2*i,2*i+alpha,x);
    }

  //imag part
  if(n>=1)
    {
      summax=(n-1)/2;
      for(int i=0;i<=summax;i++)
	{
	  v+=phase(i-1) * pow(y,2*i+1) * exp(-gsl_sf_lngamma(2*i+2) ) * gsl_sf_laguerre_n(n-1-2*i,2*i+alpha+1,x);
	}
    }

  return complexd(u,v);
}

complexd cmplx_ho_k(int n,int l,double b,complexd k)
{
  complexd z=b*k;
  complexd temp=double(phase(n)) * exp(-0.5*z*z) * pow(z,l) * sqrt(2.0) * \
    exp( 0.5*(gsl_sf_lngamma(n+1)-gsl_sf_lngamma(n+l+1.5)) )*\
    cmplx_laguerre(n,l+0.5,z*z)*pow(b,1.5);
  return temp;
}


