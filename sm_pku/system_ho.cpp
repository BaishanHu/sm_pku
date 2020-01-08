#include"system_ho.h"
#include<fstream>
#include<cmath>

using std::ifstream;
using std::abs;
using std::swap;

void System_Tz_HO::setupOrbitals()
{
  //read HO. orbitals from file
  ifstream OrbIn(OrbFile.c_str());
  if(!OrbIn) std::cerr<<"can't open Orbital Input file!\n";
  int totalOrbitals;
  OrbIn.ignore(128,'\n');
  OrbIn.ignore(128,'\n');
  OrbIn.ignore(128,':');
  OrbIn>>lengthb>>hbar_omega;
  OrbIn.ignore(128,'\n');
  OrbIn.ignore(128,'\n');
  OrbIn.ignore(128,'\n');
  OrbIn.ignore(128,'\n');
  OrbIn.ignore(128,'\n');
  OrbIn.ignore(128,'\n');
  OrbIn.ignore(45,'\n');
  OrbIn>>totalOrbitals;
  OrbIn.ignore(128,'\n');
  OrbIn.ignore(128,'\n');
  Orbitals.resize(totalOrbitals);
  int num,n,l,jj,tz,nnl;
  double e,eval;
  for(int i=0;i<totalOrbitals;i++)
    {
      OrbIn.ignore(128,':');
      OrbIn>>num>>n>>l>>jj>>tz>>nnl>>e>>eval;
      OrbIn.ignore(128,'\n');
      Orbitals[i]=OrbitalType(n,l,jj,tz,e);
    }
  setupGroups();
  int numGroups=Groups.size();
  // int nmax=12;
  // int newtotalOrbitals=numGroups*nmax;
  // vector<OrbitalType> TempOrbs=Orbitals;
  // for(int n=0;n<nmax;n++)
  //   {
  //     for(int i=0;i<numGroups;i++)
  // 	{
  // 	  int dim=Groups[i].size();
  // 	  int l=TempOrbs[ Groups[i][0] ].l;
  // 	  int jj=TempOrbs[ Groups[i][0] ].jj;
  // 	  int tz=TempOrbs[ Groups[i][0] ].tz;
  // 	  int j=0;
  // 	  for(j=0;j<dim;j++)
  // 	    {
  // 	      if(TempOrbs[Groups[i][j]].n==n) break;
  // 	    }
  // 	  if(j>=dim) Orbitals.push_back( OrbitalType(n,l,jj,tz,(2*n+l+1.5)*hbar_omega) );
  // 	}
  //   }
  // using std::cout;
  // using std::endl;
  // if(Orbitals.size()!=newtotalOrbitals) cout<<"something wrong"<<endl;
  // setupGroups();
}

double System_Tz_HO::get1B(int A,int bra,int ket) const
{
  double fac=(1-1./A);
  return fac*0.5*k2(bra,ket)*hbar_omega;
}

double System_Tz_HO::get1B_p2(int bra,int ket) const
{
  return 0.5*k2(bra,ket)*hbar_omega;
}

double System_Tz_HO::get2B(int A,int a,int b,int c,int d,int J) const
{
  TwoBodyMatEle MatEle(System_Tz_HO::MatEledim);
  get2BmatAt(a,b,c,d,J,MatEle);
  if(MatEle.empty())
    return 0;
   double V=(MatEle[0]-MatEle[3]*hbar_omega/A);
  //    double V=(MatEle[0]+0.5*relk2(a,b,c,d,J)*hbar_omega/A);
  //  for JISP (p1-p2)^2
  //   double V=MatEle[0]+MatEle[1]*hbar_omega/A*2;
  const double sqrt2=sqrt(2.);
  V *= ( (a==b)? sqrt2 : 1.0 );
  V *= ( (c==d)? sqrt2 : 1.0 );
  return V;
}

double System_Tz_HO::getV(int a,int b,int c,int d,int J) const
{
  TwoBodyMatEle MatEle(System_Tz_HO::MatEledim);
  get2BmatAt(a,b,c,d,J,MatEle);
  if(MatEle.empty())
    return 0;
   double V=MatEle[0];
  const double sqrt2=sqrt(2.);
  V *= ( (a==b)? sqrt2 : 1.0 );
  V *= ( (c==d)? sqrt2 : 1.0 );
  return V;
}

const int System_Tz_HO::MatEledim=4;//initialize of static variable
void System_Tz_HO::setupTwoBodyMat()
{
  setupIndiceIn2BStates();
  int totalChannels=TwoBodyChannels.size();
  TwoBodyMat.resize(totalChannels);
  for(int channel=0;channel<totalChannels;channel++)
    {
      int num2BInChannel=TwoBodyChannels[channel].size();
      TwoBodyMat[channel].clear();
      TwoBodyMat[channel].resize( (num2BInChannel*(num2BInChannel+1))/2 );
    }
  ifstream MatIn(MatFile.c_str());
  if(!MatIn) std::cerr<<"can't open Mat Input file!\n";
  MatIn.ignore(128,'\n');
  MatIn.ignore(128,'\n');
  MatIn.ignore(128,'\n');
  MatIn.ignore(128,'\n');
  MatIn.ignore(128,'\n');
  MatIn.ignore(128,':');

  int totalLines,ppLines,pnLines,nnLines;
  MatIn>>totalLines>>ppLines>>pnLines>>nnLines;
  MatIn.ignore(128,'\n');
  MatIn.ignore(128,'\n');
  MatIn.ignore(128,'\n');
  int Tz,Par,JJ,a,b,c,d;
  TwoBodyMatEle MatEle(MatEledim);
  /*MatEledim=4,the 4 doubles are:
    0  interaction
    1  (p1 \dot p2)/(m hbar omega) + m omega r1 \dot r2 / hbar
    2  (r1-r2)^2 m omega/(2 hbar)
    3  (p1 \dot p2)/(m hbar omega)
  */
  for(int i=0;i<ppLines;i++)
    {
      MatIn>>Tz>>Par>>JJ>>a>>b>>c>>d;
      for(int j=0;j<MatEledim;j++)
	MatIn>>MatEle[j];
      --a;--b;--c;--d;//the indice start from zero here, the input interaction from jensen start from 1.
      set2BmatAt(a,b,c,d,JJ/2,MatEle);
    }
  for(int i=0;i<pnLines;i++)
    {
      MatIn>>Tz>>Par>>JJ>>a>>b>>c>>d;
      for(int j=0;j<MatEledim;j++)
	MatIn>>MatEle[j];
      --a;--b;--c;--d;
      set2BmatAt(a,b,c,d,JJ/2,MatEle);
    }
  for(int i=0;i<nnLines;i++)
    {
      MatIn>>Tz>>Par>>JJ>>a>>b>>c>>d;
      for(int j=0;j<MatEledim;j++)
	MatIn>>MatEle[j];
      --a;--b;--c;--d;
      set2BmatAt(a,b,c,d,JJ/2,MatEle);
    }

}

double System_Tz_HO::k2(int bra,int ket) const
{
  double t=0;
  //have same quantum number
  if(Orbitals[bra].isInAGroupWith(Orbitals[ket]) )
    {
      int n_bra=Orbitals[bra].n;
      int n_ket=Orbitals[ket].n;
      int l_bra=Orbitals[bra].l;
      int l_ket=Orbitals[ket].l;
      if(n_bra==n_ket)
	t=(2*n_bra+l_bra+1.5);
      else if(n_bra-n_ket==1)
	t=0.5*sqrt(2.0*n_bra)*sqrt(2.0*n_bra+2.0*l_bra+1.0);
      else if(n_ket-n_bra==1)
	t=0.5*sqrt(2.0*n_ket)*sqrt(2.0*n_ket+2.0*l_ket+1.0);
    }
  return t;
}

double System_Tz_HO::r2(int bra,int ket) const
{
  double t=0;
  //have same quantum number
  if(Orbitals[bra].isInAGroupWith(Orbitals[ket]) )
    {
      int n_bra=Orbitals[bra].n;
      int n_ket=Orbitals[ket].n;
      int l_bra=Orbitals[bra].l;
      int l_ket=Orbitals[ket].l;
      if(n_bra==n_ket)
	t=(2*n_bra+l_bra+1.5);
      else if(n_bra-n_ket==1)
	t=-0.5*sqrt(2.0*n_bra)*sqrt(2.0*n_bra+2.0*l_bra+1.0);
      else if(n_ket-n_bra==1)
	t=-0.5*sqrt(2.0*n_ket)*sqrt(2.0*n_ket+2.0*l_ket+1.0);
    }
  return t;
}



double System_Tz_HO::rL(int bra,int ket,int L) const
{
  int nbra=Orbitals[bra].n;
  int lbra=Orbitals[bra].l;
  int nket=Orbitals[ket].n;
  int lket=Orbitals[ket].l;
  double result=0;

  int QuadPoints=x.size();
  for(int i=0;i<QuadPoints;i++)
    {
      double r=x[i];
      result+=ho_r(nbra,lbra,lengthb,r) * r*r * pow(r,L)  * ho_r(nket,lket,lengthb,r) * w[i];
    }

  // result=ho_L(nbra,lbra,nket,lket,L);
  // result*=pow(lengthb,L);
  return result;
}

double System_Tz_HO::kL(int bra,int ket,int L) const
{
  int nbra=Orbitals[bra].n;
  int lbra=Orbitals[bra].l;
  int nket=Orbitals[ket].n;
  int lket=Orbitals[ket].l;

  double result=0;
  int QuadPoints=x.size();
  for(int i=0;i<QuadPoints;i++)
    {
      double k=x[i];
      result+=ho_k(nbra,lbra,lengthb,k) * k*k * pow(k,L)  * ho_k(nket,lket,lengthb,k) * w[i];
    }

  // result=ho_L(nbra,lbra,nket,lket,L);
  // result/=pow(lengthb,L);
  // result*=phase(nbra+nket);//take care of the phase of HO. in k space
  return result;
}

double System_Tz_HO::k1k2(int a ,int b,int c,int d,int J) const
{
  int Jminab=abs(Orbitals[a].jj-Orbitals[b].jj)/2;
  int Jmincd=abs(Orbitals[c].jj-Orbitals[d].jj)/2;
  int Jmaxab=(Orbitals[a].jj+Orbitals[b].jj)/2;
  int Jmaxcd=(Orbitals[c].jj+Orbitals[d].jj)/2;

  if( J<max(Jminab,Jmincd) ) return 0;
  if( J>min(Jmaxab,Jmaxcd) ) return 0;

  int la=Orbitals[a].l;
  int jja=Orbitals[a].jj;
  int tza=Orbitals[a].tz;
	  
  int lb=Orbitals[b].l;
  int jjb=Orbitals[b].jj;
  int tzb=Orbitals[b].tz;
	  
  int lc=Orbitals[c].l;
  int jjc=Orbitals[c].jj;
  int tzc=Orbitals[c].tz;
	    
  int ld=Orbitals[d].l;
  int jjd=Orbitals[d].jj;
  int tzd=Orbitals[d].tz;

  if((la+lb+lc+ld)%2) return 0;//parity
  double result=0;
  //direct term
  if( (tza==tzc)&&(tzb==tzd) )
    {
      result+=phase((jjb+jjc)/2+J)*\
	gsl_sf_coupling_6j(jja,jjb,2*J,jjd,jjc,2)*
	4*Pi/3*kL(a,c,1)*kL(b,d,1)*\
	reducedYMatEle(la,jja,lc,jjc,1)*\
	reducedYMatEle(lb,jjb,ld,jjd,1);
    }
  //exchange term
  if( (tzb==tzc)&&(tza==tzd) )
    {
      result+=phase((jjb+jjc)/2)*\
	gsl_sf_coupling_6j(jja,jjb,2*J,jjc,jjd,2)*\
	4*Pi/3*kL(a,d,1)*kL(b,c,1)*\
	reducedYMatEle(la,jja,ld,jjd,1)*\
	reducedYMatEle(lb,jjb,lc,jjc,1);
    }
  if(a==b) result/=sqrt(2.0);
  if(c==d) result/=sqrt(2.0);

  //the final phase is the result of the phase (-i)^l, see "ho_k" in file "system.h"
  return result*lengthb*lengthb*phase( (la+lb-lc-ld)/2 );
}

double System_Tz_HO::relk2(int a,int b,int c,int d,int J) const
{
  int Jminab=abs(Orbitals[a].jj-Orbitals[b].jj)/2;
  int Jmincd=abs(Orbitals[c].jj-Orbitals[d].jj)/2;
  int Jmaxab=(Orbitals[a].jj+Orbitals[b].jj)/2;
  int Jmaxcd=(Orbitals[c].jj+Orbitals[d].jj)/2;

  if( J<max(Jminab,Jmincd) ) return 0;
  if( J>min(Jmaxab,Jmaxcd) ) return 0;

	  
  int jjc=Orbitals[c].jj;
  int jjd=Orbitals[d].jj;

  double result=0;
  //k1^2
  result+=k2(a,c)*delta(b,d)-phase((jjc+jjd)/2-J)*k2(a,d)*delta(b,c);
  //k2^2
  result+=k2(b,d)*delta(a,c)-phase((jjc+jjd)/2-J)*k2(b,c)*delta(a,d);
  if(a==b) result/=sqrt(2.0);
  if(c==d) result/=sqrt(2.0);
  //-2k1k2
  result-=2*k1k2(a,b,c,d,J);
  return result;
}

double System_Tz_HO::r1r2(int a,int b,int c,int d,int J) const
{
  int Jminab=abs(Orbitals[a].jj-Orbitals[b].jj)/2;
  int Jmincd=abs(Orbitals[c].jj-Orbitals[d].jj)/2;
  int Jmaxab=(Orbitals[a].jj+Orbitals[b].jj)/2;
  int Jmaxcd=(Orbitals[c].jj+Orbitals[d].jj)/2;

  if( J<max(Jminab,Jmincd) ) return 0;
  if( J>min(Jmaxab,Jmaxcd) ) return 0;

	  
  int la=Orbitals[a].l;
  int jja=Orbitals[a].jj;
  int tza=Orbitals[a].tz;
	  
  int lb=Orbitals[b].l;
  int jjb=Orbitals[b].jj;
  int tzb=Orbitals[b].tz;
	  
  int lc=Orbitals[c].l;
  int jjc=Orbitals[c].jj;
  int tzc=Orbitals[c].tz;
	    
  int ld=Orbitals[d].l;
  int jjd=Orbitals[d].jj;
  int tzd=Orbitals[d].tz;

  double t1=0,t2=0;
  double result=0;
  //direct term
  if( (tza==tzc)&&(tzb==tzd) )
    {
      result+=phase((jjb+jjc)/2+J)*\
	gsl_sf_coupling_6j(jja,jjb,2*J,jjd,jjc,2)*
	4*Pi/3*rL(a,c,1)*rL(b,d,1)*\
	reducedYMatEle(la,jja,lc,jjc,1)*\
	reducedYMatEle(lb,jjb,ld,jjd,1);
    }
  //exchange term
  if( (tzb==tzc)&&(tza==tzd) )
    {
      result+=phase((jjb+jjc)/2)*\
	gsl_sf_coupling_6j(jja,jjb,2*J,jjc,jjd,2)*\
	4*Pi/3*rL(a,d,1)*rL(b,c,1)*\
	reducedYMatEle(la,jja,ld,jjd,1)*\
	reducedYMatEle(lb,jjb,lc,jjc,1);
    }
  if(a==b) result/=sqrt(2.0);
  if(c==d) result/=sqrt(2.0);

  return result/(lengthb*lengthb);
}

double System_Tz_HO::relr2(int a,int b,int c,int d,int J) const
{
  int Jminab=abs(Orbitals[a].jj-Orbitals[b].jj)/2;
  int Jmincd=abs(Orbitals[c].jj-Orbitals[d].jj)/2;
  int Jmaxab=(Orbitals[a].jj+Orbitals[b].jj)/2;
  int Jmaxcd=(Orbitals[c].jj+Orbitals[d].jj)/2;

  if( J<max(Jminab,Jmincd) ) return 0;
  if( J>min(Jmaxab,Jmaxcd) ) return 0;
	  
  int jjc=Orbitals[c].jj;
  int jjd=Orbitals[d].jj;

  double result=0;
  //r1^2
  result+=r2(a,c)*delta(b,d)-phase((jjc+jjd)/2-J)*r2(a,d)*delta(b,c);
  //r2^2
  result+=r2(b,d)*delta(a,c)-phase((jjc+jjd)/2-J)*r2(b,c)*delta(a,d);
  if(a==b) result/=sqrt(2.0);
  if(c==d) result/=sqrt(2.0);
  //-2r1r2
  result-=2*r1r2(a,b,c,d,J);
  return result;
}



void System_Tz_HO::printOrbitals() const
{
  for(int i=0;i<Orbitals.size();i++)
    {
      std::cout<<Orbitals[i].l<<"\t"<<Orbitals[i].jj<<"\t"<<Orbitals[i].tz<<"\t"<<Orbitals[i].e<<std::endl;
    }
}




//****************************************************************************************
//****************************************************************************************


void System_T_HO::setupOrbitals()
{
  ifstream OrbIn(OrbFile.c_str());
  if(!OrbIn) std::cerr<<"can't open Orbital Input file!\n";
  int totalOrbitals;
  OrbIn.ignore(128,'\n');
  OrbIn>>totalOrbitals;
  Orbitals.resize(totalOrbitals);
  lengthb=sqrt(hbar_c*hbar_c/(hbar_omega*mc2));
  int number,n,l,jj,nnl;
  for(int i=0;i<totalOrbitals;i++)
    {
      OrbIn>>number>>n>>l>>jj>>nnl;
      double e=(2*n+l+1.5)*hbar_omega;
      Orbitals[i]=OrbitalType(n,l,jj,e);
    }
  setupGroups();  
}

const int System_T_HO::MatEledim=6;//initialize of static variable
void System_T_HO::setupTwoBodyMat()
{
  setupIndiceIn2BStates();
  int totalChannels=TwoBodyChannels.size();
  TwoBodyMat.resize(totalChannels);
  for(int channel=0;channel<totalChannels;channel++)
    {
      int num2BInChannel=TwoBodyChannels[channel].size();
      TwoBodyMat[channel].clear();
      TwoBodyMat[channel].resize( (num2BInChannel*(num2BInChannel+1))/2 );
    }
  ifstream MatIn(MatFile.c_str());
  if(!MatIn) std::cerr<<"can't open Mat Input file!\n";
  int totalLines;
  MatIn>>totalLines;
  MatIn.ignore(128,'\n');    
  int a,b,c,d,J,T,Par;
  char ac[4],bc[4],cc[4],dc[4];//to read JISP with fixed width integer given by fortran program,some number are not separated by space,but every integer has a width of 4
	  
  for(int i=0;i<totalLines;i++)
    {
      MatIn.get(ac,4);
      MatIn.get(bc,4);
      MatIn.get(cc,4);
      MatIn.get(dc,4);
      a=atoi(ac);b=atoi(bc);c=atoi(cc);d=atoi(dc);
	      
      MatIn>>J>>T;
      TwoBodyMatEle MatEle(MatEledim);
      /*MatEledim=6,the 6 double are:
	0  (p1-p2)^2/(4m hbar omega)
	1  (p1-p2)^2/(4m habr omega)+ (r1-r2)^2 m omega/(4 hbar), relative hamiltion
	2  coloumb interaction
	3  nuclear interaction for Tz=0
	4  nuclear interaction for Tz=1
	5  nuclear interaction for Tz=-1
      */
      for(int j=0;j<MatEledim;j++)
	MatIn>>MatEle[j];
      MatIn.ignore(128,'\n');
      MatEle[2]*=sqrt(mc2/938.093);
      --a;--b;--c;--d;//the indice start from zero here, the input interaction start from 1.
      set2BmatAt(a,b,c,d,J,T,MatEle);
    }
}


//****************************************************************************************
//****************************************************************************************


void TtoTz(const System_T_HO&HT,System_Tz_HO&HTz)
{
  //copy some parameter
  HTz.lengthb=HT.lengthb;
  HTz.hbar_omega=HT.hbar_omega;
  //double Orbitals
  int TtotalOrbitals=HT.Orbitals.size();
  HTz.Orbitals.resize(2*TtotalOrbitals);
  for(int i=0;i<TtotalOrbitals;i++)
    {
      int n=HT.Orbitals[i].n;
      int l=HT.Orbitals[i].l;
      int jj=HT.Orbitals[i].jj;
      double e=HT.Orbitals[i].e;
      //first proton, second neutron
      HTz.Orbitals[2*i]=System_Tz_HO::OrbitalType(n,l,jj,-1,e);
      HTz.Orbitals[2*i+1]=System_Tz_HO::OrbitalType(n,l,jj,1,e);
    }
  HTz.setupGroups();
  //setup twobody states
  HTz.setupTwoBodyStates();
  //setup IndiceIn2BStates
  HTz.setupIndiceIn2BStates();
  //setup TwoBodyMat
  int TztotalChannels=HTz.TwoBodyChannels.size();
  HTz.TwoBodyMat.resize(TztotalChannels);
  for(int channel=0;channel<TztotalChannels;channel++)
    {
      int num2BInChannel=HTz.TwoBodyChannels[channel].size();
      HTz.TwoBodyMat[channel].clear();
      HTz.TwoBodyMat[channel].resize( (num2BInChannel*(num2BInChannel+1))/2 );
    }
#pragma omp parallel for
  for(int channel=0;channel<TztotalChannels;channel++)
    {
      int numInAChan=HTz.TwoBodyChannels[channel].size();
      for(int bra=0;bra<numInAChan;bra++)
	for(int ket=bra;ket<numInAChan;ket++)
	  {
	    int BRA=HTz.TwoBodyChannels[channel][bra];
	    int KET=HTz.TwoBodyChannels[channel][ket];
	    const System_Tz_HO::TwoBodyStateType & Bra2B=HTz.TwoBodyStates[BRA];
	    const System_Tz_HO::TwoBodyStateType & Ket2B=HTz.TwoBodyStates[KET];
	    int a=Bra2B.a;
	    int b=Bra2B.b;
	    int c=Ket2B.a;
	    int d=Ket2B.b;
	    int J=Bra2B.J;
	    int Tz=Bra2B.Tz;
	    int pha=1;
	    //change to pn order
	    if(HTz.Orbitals[a].tz>HTz.Orbitals[b].tz)
	      {
		swap(a,b);
		//phase for changing to pn order
		pha*=phase((HTz.Orbitals[a].jj+HTz.Orbitals[b].jj)/2-J+1);
	      }
	    if(HTz.Orbitals[c].tz>HTz.Orbitals[d].tz)
	      {
		swap(c,d);
		//phase for changing to pn order
		pha*=phase((HTz.Orbitals[c].jj+HTz.Orbitals[d].jj)/2-J+1);
	      }
	    System_T_HO::TwoBodyMatEle TMatEle;
	    /*dim of TMatEle=6,the 6 double are:
	      0  (p1-p2)^2/(4m hbar omega)
	      1  (p1-p2)^2/(4m habr omega)+ (r1-r2)^2 m omega/(4 hbar), relative hamiltion
	      2  coloumb interaction
	      3  nuclear interaction for Tz=0
	      4  nuclear interaction for Tz=1
	      5  nuclear interaction for Tz=-1
	    */
	    HT.get2BmatAt(a/2,b/2,c/2,d/2,J,1,TMatEle);
	    System_Tz_HO::TwoBodyMatEle MatEle(System_Tz_HO::MatEledim);
	    /*dim of MatEle=4,the 4 double are:
	      dim of MatEle=4,the 4 double are:
	      0  interaction
	      1  (p1 \dot p2)/(m hbar omega) + m omega r1 \dot r2 / hbar
	      2  (r1-r2)^2 m omega/(2 hbar)
	      3  (p1 \dot p2)/(m hbar omega)
	    */
	    //convert
	    MatEle.setZero();
	    if(Tz==1)
	      {
		if(TMatEle.empty()) continue;
		MatEle[0]=TMatEle[4];
		MatEle[1]=TMatEle[0];
	      }
	    else if(Tz==-1)
	      {
		if(TMatEle.empty()) continue;
		MatEle[0]=TMatEle[5]+TMatEle[2];
		MatEle[1]=TMatEle[0];
	      }
	    else if(Tz==0)
	      {
		if(!TMatEle.empty())
		  {
		    MatEle[0]=TMatEle[3];
		    MatEle[1]=TMatEle[0];
		  }
		HT.get2BmatAt(a/2,b/2,c/2,d/2,J,0,TMatEle);
		if(!TMatEle.empty())
		  {
		    MatEle[0]+=TMatEle[3];
		    MatEle[1]+=TMatEle[0];
		  }
		MatEle*=0.5;
		if(a/2==b/2) MatEle*=sqrt(2.);
		if(c/2==d/2) MatEle*=sqrt(2.);
	      }
	    //some terms can calculate directly, not transform from System_T_HO
	    //	    double p1p2=HTz.k1k2(a,b,c,d,J);
	    //	    double r1r2=HTz.r1r2(a,b,c,d,J);
	    //	    double relr2=HTz.relr2(a,b,c,d,J)*0.5;
	    //	    MatEle[1]=p1p2+r1r2;
	    //	    MatEle[2]=relr2;
	    //	    MatEle[3]=p1p2;
		    
	    MatEle*=pha;
	    HTz.set2Bmat(channel,bra,ket,MatEle);
	  }
    }
}

//TODO
void TztoT(const System_Tz_HO&,System_T_HO&)
{
	  
}
