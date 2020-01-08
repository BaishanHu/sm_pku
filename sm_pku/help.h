#ifndef _HELP_H
#define _HELP_H
#include<cmath>
#include<vector>
#include<complex>
#include<bitset>
#include<limits>//numeric_limits

//gsl special functions, 3j,6j,9j coupling coeff. etc.
#include<gsl/gsl_sf.h>
#include<gsl/gsl_integration.h>


//eigen library
#include<eigen3/Eigen/Eigen>//matrix library

#include<iostream>
#include"constants.h"

using std::vector;
using std::abs;//important for abs taken double type as parameter.
using std::bitset;
using std::complex;

using Eigen::Matrix;
using Eigen::Dynamic;

typedef complex<double> complexd;

int phase(int n);///< return -1 if n is odd, 1 if n is even.

int delta(int a,int b);

///cg. coeff
double cg(int jj1,int jj2,int jj3,int mm1,int mm2,int mm3);


///map a<b matrix indice to one dimension indice, a must less than b
int mapab(int a,int b,int dim);
///map a<=b matrix indice to one dimension indice, a must not greater than b
int mapab_e(int a,int b,int dim);


///sqrt of complex
complexd sqrt(complexd c);
//compare between complex
bool operator < (const complexd & c1,const complexd &c2);



///HO. orbits in k space, radial part: R_nl(k), b is oscillator length. b=sqrt( hbar/(m omega) )
///note that the phase is (-i)^(2n+l) compare to ho_r, we throw the (-i)^l away, we should be careful and add the phase explicitly in our expression.
double ho_k(int n,int l,double b,double k);

///HO. orbits in r space, radial part: R_nl(r)
double ho_r(int n,int l,double b,double r);

/*
for numeric integration of <HO|r^\lambda|HO>,using gsl.
*/
struct nlb
{
  int n1,l1,n2,l2;
  int L;
};
double ho_f(double x,void *params);

double ho_L(int n1,int l1,int n2,int l2,int L);///<return numeric integration <HO1|r^L|HO2>,considering only the radial part. The length of HO. is 1.


///gauss-legendre quadrature
void gauleg(const double x1, const double x2, vector<double> &x, vector<double> &w);



///return reduced Matrix element<la ja|| Y^L ||lb jb>, jja is twice of ja, jja is odd, same with jb. 
double reducedYMatEle(int la,int jja,int lb,int jjb,int L);


using std::complex;
typedef complex<double> complexd;
typedef vector<complexd> vcomplexd;
///complex spherical bessel function jn(z),yn(z)
void csphbessel_jy(int n,complexd z,vcomplexd&jn,vcomplexd&yn,vcomplexd&djn,vcomplexd&dyn);
complexd csphbessel_jn(int n,complexd z);
complexd csphbessel_yn(int n,complexd z);
complexd csphbessel_djn(int n,complexd z);
complexd csphbessel_dyn(int n,complexd z);


///complex laguerre based on gsl_sf_laguerre
complexd cmplx_laguerre(int n,double alpha,complexd z);

///complex HO. orbital in k space
complexd cmplx_ho_k(int n,int l,double b,complexd k);


//****************************************************************
///C(n,r). combination number
int Cnr(int iN, int iR);

///some operator for bit operation
template<size_t N>
bitset<N> operator + (const bitset<N>&a, const bitset<N>&b)
{
  bitset<N> ans,atemp(a),btemp(b);
  while(btemp.any())
    {
      ans = atemp^btemp;
      btemp = (atemp&btemp)<<1;
      atemp = ans;
    }
  return ans;
}
template<size_t N>
bitset<N> operator - (const bitset<N>&a)
{
  return ~a + bitset<N>(1);
}
template<size_t N>
bitset<N> operator - (const bitset<N>&a,const bitset<N>&b)
{
  return a + (-b);
}
///number of zeros after first 1, number of trailing zeros
template<size_t N>
int ctz(bitset<N> a)
{
  if(a.none()) return 0;
  int count=0;
  while( !a.test(0) )
    {
      a>>=1;
      ++count;
    }
  return count;
}
template<size_t N>
bitset<N> next(const bitset<N>&a)
{
  bitset<N> one(1);
  bitset<N> t= a | (a-one);
  return (t + one) | (((~t & -~t) - one) >> (ctz(a) + 1));
}

/// find the indices of lowest num ones
template<size_t N>
void findones(const bitset<N> &a,vector<int>&pos,int num)
{
  pos.resize(num);
  int count=0;
  for(int i=0;i<N;i++)
    {
      if( a.test(i) )
	{
	  pos[count]=i;
	  ++count;
	}
      if(count==num) break;
    }
}

///return the number of ones between position region (n,m] of a bit, n>m, not including n,including m, position indice start from 0
template<size_t N>
int count(const bitset<N> & a,int n,int m)
{
  if(n<=m) return 0;
  bitset<N> one(1);
  return (( (one<<n) - (one<<m) )&a).count();
}
//****************************************************************



template <class ScalarType>
ScalarType pythag(const ScalarType a,const ScalarType b)
{
  double absa=abs(a),absb=abs(b);
  return (absa>absb? a * sqrt(1.+pow(b/a,2)):
  (absb< std::numeric_limits<double>::epsilon() ? 0.:b*sqrt(1.+ pow(a/b,2))));
  //  return sqrt(pow(a,2)+pow(b,2));//the vesion up if for overflow considration.
}

///QL algorithm taken from numerical recipe c++
template <class ScalarType>
void TriDiagQL(Matrix<ScalarType,Dynamic,1> & diag,Matrix<ScalarType,Dynamic,1> & subdiag,Matrix<ScalarType,Dynamic,Dynamic> & Q,int end=-1,bool yesvecs=true)  
{
  const int maxIter=200;
  if(end==-1) end=diag.rows();
  int m,l,iter,i,k;
  ScalarType s,r,p,g,f,c,b;
  double dd;
  //  const double EPS=std::numeric_limits<double>::epsilon();
  const double EPS=1e-8;
  subdiag[end-1]=0.;
  for(l=0;l<end;l++)
    {
      iter=0;
      do{
	for(m=l;m<end-1;m++)
	  {
	    dd=abs(diag[m]) + abs(diag[m+1]);
	    //	    if(abs(subdiag[m]) <= EPS*dd ) break;
	    if(abs(subdiag[m]) <= EPS ) break;
	  }
	if(m!=l)
	  {
	    if(iter++ == maxIter) throw("Too many iterations!");
	    //	    cout<<iter<<endl;
	    g=(diag[l+1]-diag[l])/(2.0*subdiag[l]);
	    r=pythag(g,ScalarType(1.0));
	    ScalarType t1=subdiag[l]/(g+r);
	    ScalarType t2=subdiag[l]/(g-r);
	    g=diag[m]-diag[l] +  ( (abs(t1)<abs(t2))?t1:t2 );
		    
	    s=c=1.0;
	    p=0.0;
	    for(i=m-1;i>=l;i--)
	      {
		f=s*subdiag[i];
		b=c*subdiag[i];
		subdiag[i+1]=(r=pythag(f,g));
		if(abs(r)<EPS)
		  {
		    diag[i+1]-=p;
		    subdiag[m]=0.0;
		    break;
		  }
		s=f/r;
		c=g/r;
		g=diag[i+1]-p;
		r=(diag[i]-g)*s+2.0*c*b;
		diag[i+1]=g+(p=s*r);
		g=c*r-b;

		if(yesvecs)
		  {
		    for(k=0;k<Q.rows();k++)
		      {
			f=Q(k,i+1);
			Q(k,i+1)=s*Q(k,i)+c*f;
			Q(k,i)=c*Q(k,i)-s*f;
		      }
		  }
	      }
	    if(abs(r)<EPS && i>=1) continue;
	    diag[l]-=p;
	    subdiag[l]=g;
	    subdiag[m]=0.0;
	  }
      } while(m!=l);
    }
}

//****************************************************************



#endif
