#include"sm_system.h"
#include"help.h"
#include<fstream>
#include<iostream>
#include<cmath>//abs,sqrt
#include<algorithm>//min,max

using std::ifstream;
using std::abs;
using std::min;
using std::max;

namespace SM
{
  void System_Tz_SM::setupOrbitals()
  {
    //read sp orbitals from file
    ifstream OrbIn(OrbFile.c_str());
    if(!OrbIn) std::cerr<<"can't open Orbital Input file!\n";
    int totalOrbitals;
    OrbIn>>totalOrbitals;
    OrbIn.ignore(128,'\n');
    Orbitals.resize(totalOrbitals);
    int num,type,n,l,jj,tz,min,max;;
    DataType e;
    double ereal;
    double eimag;
    for(int i=0;i<totalOrbitals;i++)
      {
	OrbIn>>num>>type>>min>>max>>n>>l>>jj>>tz>>ereal>>eimag;
	//Orbitals[num-1]=OrbitalType(l,jj,tz,ereal);	
	Orbitals[num-1]=OrbitalType(l,jj,tz,complexd(ereal,eimag),min,max,type);
      }
    setupGroups();
    int channels=Groups.size();
    OneBodyMat.resize(channels);
    for(int i=0;i<channels;i++)
      {
        int dim=Groups[i].size();
        OneBodyMat[i].resize(dim*(dim+1)/2);
      }
    int a,b;
    while(OrbIn>>a>>b>>ereal>>eimag)
      {
        a--;b--;
        int channel=GroupIndice[a];
        int channelb=GroupIndice[b];
        if(channel == channelb)
        {
          int aa=OrbIndice[a],bb=OrbIndice[b];
          int dim=Groups[channel].size();
          OneBodyMat[channel][mapab_e(aa,bb,dim)]=complexd(ereal,eimag);
        }
      }
  }

  const int System_Tz_SM::MatEledim=1;//initialize of static variable
  void System_Tz_SM::setupTwoBodyMat()
  {
    int totalChannels=TwoBodyChannels.size();
    TwoBodyMat.resize(totalChannels);
    for(int channel=0;channel<totalChannels;channel++)
      {
	TwoBodyMat[channel].clear();
	int num2BInChannel=TwoBodyChannels[channel].size();
	TwoBodyMat[channel].resize( (num2BInChannel*(num2BInChannel+1))/2 );
      }
    ifstream MatIn(MatFile.c_str());
    if(!MatIn) std::cerr<<"can't open Mat Input file!\n";
    int totalLines;
    MatIn>>totalLines;

    int J,a,b,c,d;
    double vreal,vimag;    
    TwoBodyMatEle MatEle(MatEledim);
    /*MatEledim=1,the value is interaction*/
    for(int i=0;i<totalLines;i++)
      {
	MatIn>>a>>b>>c>>d>>J;
	for(int j=0;j<MatEledim;j++)
	  {
	    MatIn>>vreal>>vimag;
	    //	    MatEle[j]=vreal;	    
	    MatEle[j]=complexd(vreal,vimag);
	  }
	--a;--b;--c;--d;//the indice start from zero here, the input interaction orbital order start from 1.
	set2BmatAt(a,b,c,d,J,MatEle);
      }
  }

  void System_Tz_SM::setupOper1B(const string &file)
  {
    if(file=="*")
      {
	lambda1B=-1;
	return;
      }
	
    ifstream fin(file.c_str());
    char ch;
    do
      {
	fin>>ch;
	if(ch=='!')
	  fin.ignore(128,'\n');
      }while(ch=='!');
      fin.putback(ch);
      fin>>lambda1B;

      Oper1B.clear();
      int a,b;
      double oper;
      int dim=Orbitals.size();
      while(fin>>a>>b>>oper)
	{
	  a--;b--;
	  Oper1B[a*dim+b]=oper;
	}
      //std::cout<<Oper1B<<endl;
  }
  void System_Tz_SM::setupOper2B(const string &file)
  {
    if(file=="*")
      {
	lambda2B=-1;
	return;
      }    
    ifstream fin(file.c_str());
    char ch;
    do
      {
	fin>>ch;
	if(ch=='!')
	  fin.ignore(128,'\n');
      }while(ch=='!');
    fin.putback(ch);
    fin>>lambda2B;

    Oper2B.clear();
    int a,b,c,d,Jab,Jcd;
    double oper;
    int dim=TwoBodyStates.size();
    while(fin>>a>>b>>c>>d>>Jab>>Jcd>>oper)
      {
	double phase=1.;
	a--;b--;c--;d--;
	if(a>b)
	  {
	    swap(a,b);
	    phase*=( (Orbitals[a].jj+Orbitals[b].jj)/2 -Jab + 1 )%2?-1:1;
	  }
	if(c>d)
	  {
	    swap(c,d);
	    phase*=( (Orbitals[c].jj+Orbitals[d].jj)/2 -Jcd + 1 )%2?-1:1;
	  }
	int Jminab=abs(Orbitals[a].jj-Orbitals[b].jj)/2;
	int Jmincd=abs(Orbitals[c].jj-Orbitals[d].jj)/2;
	int ab=IndiceIn2BStates[a][b][Jab-Jminab];
	int cd=IndiceIn2BStates[c][d][Jcd-Jmincd];
	int dim=TwoBodyStates.size();
	if(ab>cd)
	  {
	    swap(ab,cd);
	    phase*=(Jab-Jcd)%2?-1:1;
	  }
	Oper2B[ab*dim+cd]=oper*phase;
      }
  }
  
  DataType System_Tz_SM::get1BOper(int a,int b) const
  {
    int dim=Orbitals.size();
    int ab=a*dim+b;
    if( Oper1B.count(ab) )
      return Oper1B.at(ab);      
    else
      return 0;
    
    // if(a==b)
    //   return Orbitals[a].e*sqrt(Orbitals[a].jj+1.);
    // else
    //   return 0;
  }

  DataType System_Tz_SM::get2BOper(int a,int b,int c,int d,int Jab,int Jcd) const
  {
    double phase=1.;
    if(a>b)
      {
	swap(a,b);
	phase*=( (Orbitals[a].jj+Orbitals[b].jj)/2 -Jab + 1 )%2?-1:1;
      }
    if(c>d)
      {
	swap(c,d);
	phase*=( (Orbitals[c].jj+Orbitals[d].jj)/2 -Jcd + 1 )%2?-1:1;
      }
    int Jminab=abs(Orbitals[a].jj-Orbitals[b].jj)/2;
    int Jmincd=abs(Orbitals[c].jj-Orbitals[d].jj)/2;
    int ab=IndiceIn2BStates[a][b][Jab-Jminab];
    int cd=IndiceIn2BStates[c][d][Jcd-Jmincd];
    int dim=TwoBodyStates.size();
    if(ab>cd)
      {
	swap(ab,cd);
	phase*=(Jab-Jcd)%2?-1:1;
      }    
    int abcd=ab*dim+cd;
    if( Oper2B.count(abcd) )
      return Oper2B.at(abcd)*phase;
    else
      return 0;
    
    // if(Jab==Jcd)
    //   {
    // 	TwoBodyMatEle MatEle(MatEledim);
    // 	get2BmatAt(a,b,c,d,Jab,MatEle);
    // 	return MatEle[0]*sqrt(2*Jab+1.);
    //   }
    // else
    //   return 0;
  }
  
  void System_Tz_SM::printOrbitals() const
  {
    for(int i=0;i<Orbitals.size();i++)
      {
	std::cout<<Orbitals[i].l<<"\t"<<Orbitals[i].jj<<"\t"<<Orbitals[i].tz<<"\t"<<Orbitals[i].e<<std::endl;
      }
  }


  void System_Tz_SM::printTwoBodyStates() const
  {
    for(int i=0;i<TwoBodyStates.size();i++)
      {
	std::cout<<TwoBodyStates[i].a<<"\t"<<TwoBodyStates[i].b<<"\t"<<TwoBodyStates[i].Tz<<"\t"<<TwoBodyStates[i].Par<<"\t"<<TwoBodyStates[i].J<<std::endl;
      }
  }  
  //****************************************************************




  
  void System_Tz_SM_M::setupOrbitals()
  {
    vector<OrbitalType> OrbitalsTemp;
    vector<int> OrdersTemp;
    //split different m orbitals from System_Tz_SM
    for(int i=0;i<pSystem->Orbitals.size();i++)
      {
	int l=pSystem->Orbitals[i].l;
	int jj=pSystem->Orbitals[i].jj;
	int tz=pSystem->Orbitals[i].tz;
	DataType e=pSystem->Orbitals[i].e;
	for(int mm=-jj;mm<=jj;mm+=2)
	  {
	    OrbitalsTemp.push_back(OrbitalType(l,jj,mm,tz,e));
	    OrdersTemp.push_back(i);
	  }
      }
    //store the Orbitals so that protons and neutrons are separated, proton first.
    //notice that we use this property when generating configurations in Configs.
    //the orbitals with same j but adjacent m is put together, this is important to simplify the cal. of diagonal element of J^2 between configs, see getJ2 in Configs.
    Orbitals.clear();
    Orders.clear();
    Orbitals.reserve( OrbitalsTemp.size() );
    Orders.reserve( OrdersTemp.size() );
    for(int i=0;i<OrbitalsTemp.size();i++)
      {
	if(OrbitalsTemp[i].tz==-1)
	  {
	    Orbitals.push_back(OrbitalsTemp[i]);
	    Orders.push_back(OrdersTemp[i]);
	  }
      }
    maxProtons=Orbitals.size();
    for(int i=0;i<OrbitalsTemp.size();i++)
      {
	if(OrbitalsTemp[i].tz==1)
	  {
	    Orbitals.push_back(OrbitalsTemp[i]);
	    Orders.push_back(OrdersTemp[i]);
	  }
      }
    maxNeutrons=Orbitals.size()-maxProtons;
  }
  
  void System_Tz_SM_M::setupTwoBodyStates()
  {
    TwoBodyStates.clear();
    for(int i=0;i<Orbitals.size();i++)
      for(int j=i+1;j<Orbitals.size();j++)
	{
	  int Tz=(Orbitals[i].tz + Orbitals[j].tz)/2;
	  int Par=(Orbitals[i].par + Orbitals[j].par)%2;
	  int M=(Orbitals[i].mm + Orbitals[j].mm)/2;;
	  TwoBodyStates.push_back(TwoBodyStateType(i,j,Tz,Par,M));
	}
    setupTwoBodyChannels();//noting that the twobody states are sorted after this calling.
  }

  const int  System_Tz_SM_M::MatEledim=System_Tz_SM_M::VSystem::MatEledim;//initialize of static variable
  void System_Tz_SM_M::setupTwoBodyMat()
  {
    setupIndiceIn2BStates();
    int totalChannels=TwoBodyChannels.size();
    TwoBodyMat.resize(totalChannels);
    #pragma omp parallel for
    for(int channel=0;channel<totalChannels;channel++)
      {
	int num2BInChannel=TwoBodyChannels[channel].size();
	TwoBodyMat[channel].clear();
	TwoBodyMat[channel].resize( (num2BInChannel*(num2BInChannel+1))/2 );
	for(int bra=0;bra<num2BInChannel;bra++)
	  for(int ket=bra;ket<num2BInChannel;ket++)
	    {
	      int ab=TwoBodyChannels[channel][bra];
	      int cd=TwoBodyChannels[channel][ket];
	      int a=TwoBodyStates[ab].a;
	      int b=TwoBodyStates[ab].b;
	      int c=TwoBodyStates[cd].a;
	      int d=TwoBodyStates[cd].b;
	      TwoBodyMatEle MatEle(MatEledim);
	      get2BmatOnNeed(a,b,c,d,MatEle);
	      set2Bmat(channel,bra,ket,MatEle);
	    }
      }
  }

  void System_Tz_SM_M::setupIndiceIn2BStates()
  {
    int totalOrbitals=Orbitals.size();
    IndiceIn2BStates.resize((totalOrbitals*(totalOrbitals-1))/2);
    for(int i=0;i<TwoBodyStates.size();i++)
      {
	int a=TwoBodyStates[i].a;
	int b=TwoBodyStates[i].b;
	//a<b
	IndiceIn2BStates[mapab(a,b)]=i;
      }
  }

  void System_Tz_SM_M::get2BmatAt(int a,int b,int c,int d,TwoBodyMatEle & MatEle) const
  {
    MatEle.setZero();
    if(a==b || c==d) return;//pauli principle
    int phase=1;
    //make sure a<b,c<d,(ab)<=(cd)
    if(a>b)
      {
	swap(a,b);
	phase*=-1;
      }
    if(c>d)
      {
	swap(c,d);
	phase*=-1;
      }
    int Mab=(Orbitals[a].mm + Orbitals[b].mm)/2;
    int Mcd=(Orbitals[c].mm + Orbitals[d].mm)/2;

    if( Mab!=Mcd ) return;
    int ab=IndiceIn2BStates[mapab(a,b)];
    int cd=IndiceIn2BStates[mapab(c,d)];
    if(ab==-1||cd==-1) return;
    int channelab=ChannelIndice[ab];
    int channelcd=ChannelIndice[cd];
    //channelab should = channelcd
    if(channelab!=channelcd) return;
    int bra=TwoBodyStateIndice[ab];
    int ket=TwoBodyStateIndice[cd];
    if(bra>ket)
      swap(bra,ket);
    get2Bmat(channelab,bra,ket,MatEle);
    MatEle*=phase;
  }

  void System_Tz_SM_M::get2BmatOnNeed(int a,int b,int c,int d,TwoBodyMatEle & MatEle) const
  {
    MatEle.resize(MatEledim);
    MatEle.setZero();

    if(a==b||c==d) return;//pauli
    int jja=Orbitals[a].jj;
    int mma=Orbitals[a].mm;
    int jjb=Orbitals[b].jj;
    int mmb=Orbitals[b].mm;
    int Tzab=(Orbitals[a].tz+Orbitals[b].tz)/2;;
    int MMab=mma+mmb;
    
    int jjc=Orbitals[c].jj;
    int mmc=Orbitals[c].mm;    
    int jjd=Orbitals[d].jj;
    int mmd=Orbitals[d].mm;    
    int Tzcd=(Orbitals[c].tz+Orbitals[d].tz)/2;
    int MMcd=mmc+mmd;

    int Jmin=max(abs(jja-jjb)/2,abs(jjc-jjd)/2);
    int Jmax=min((jja+jjb)/2,(jjc+jjd)/2);

    if(Tzab!=Tzcd) return;
    if(MMab!=MMcd) return;
    if((Orbitals[a].l+Orbitals[b].l+Orbitals[c].l+Orbitals[d].l)%2) return;//parity check
    if(Jmin>Jmax) return;

    int MM=MMab;
    int ordera=Orders[a];
    int orderb=Orders[b];
    int orderc=Orders[c];
    int orderd=Orders[d];
    VSystem::TwoBodyMatEle Temp;
    for(int J=Jmin;J<=Jmax;J++)
      {
	pSystem->get2BmatAt(ordera,orderb,orderc,orderd,J,Temp);
	if(Temp.empty()) continue;
	if(ordera==orderb) Temp*=sqrt(2.);
	if(orderc==orderd) Temp*=sqrt(2.);
	MatEle+=Temp*cg(jja,jjb,J*2,mma,mmb,MM)*cg(jjc,jjd,J*2,mmc,mmd,MM);
      }
  }

  DataType System_Tz_SM_M::get1BOper(int a,int b) const
  {
    if(Orbitals[a].mm==Orbitals[b].mm)
      {
	int lambda=pSystem->lambda1B;
	int jja=Orbitals[a].jj;
	int mma=Orbitals[a].mm;
	int jjb=Orbitals[b].jj;
	int mmb=Orbitals[b].mm;
	if(abs(jja-jjb)>2*lambda||(jja+jjb)<2*lambda) return 0;
	double fac=cg(jjb,2*lambda,jja,mmb,0,mma)/sqrt(jja+1.);
//bhu
       if(lambda==0) fac=1.0;
//bhu
	return pSystem->get1BOper(Orders[a],Orders[b]) * fac;
      }
    else
      return 0;    
  }

  DataType System_Tz_SM_M::get2BOper(int a,int b,int c,int d) const
  {    
    if(a==b||c==d) return 0;//pauli
    int jja=Orbitals[a].jj;
    int mma=Orbitals[a].mm;
    int jjb=Orbitals[b].jj;
    int mmb=Orbitals[b].mm;
    int MMab=mma+mmb;
    
    int jjc=Orbitals[c].jj;
    int mmc=Orbitals[c].mm;    
    int jjd=Orbitals[d].jj;
    int mmd=Orbitals[d].mm;    
    int MMcd=mmc+mmd;

    if(MMab!=MMcd) return 0;    
    int Jabmin=abs(jja-jjb)/2;
    int Jcdmin=abs(jjc-jjd)/2;
    int Jabmax=(jja+jjb)/2;
    int Jcdmax=(jjc+jjd)/2;


    int lambda=pSystem->lambda2B;
    int MM=MMab;
    int ordera=Orders[a];
    int orderb=Orders[b];
    int orderc=Orders[c];
    int orderd=Orders[d];
    DataType Temp,Oper(0);
    for(int Jab=Jabmin;Jab<=Jabmax;Jab++)
      {
	for(int Jcd=Jcdmin;Jcd<=Jcdmax;Jcd++)
	  {
	    if(abs(Jab-Jcd)>lambda || (Jab+Jcd)<lambda) continue;
	    Temp=pSystem->get2BOper(ordera,orderb,orderc,orderd,Jab,Jcd);
	    double fac=cg(jja,jjb,Jab*2,mma,mmb,MM)*cg(jjc,jjd,Jcd*2,mmc,mmd,MM);
//	    fac*=cg(Jcd*2,lambda*2,Jab*2,MM,0,MM)/sqrt(Jab*2+1.);
//bhu
            if(lambda!=0) fac*=cg(Jcd*2,lambda*2,Jab*2,MM,0,MM)/sqrt(Jab*2+1.);
//bhu

	    if(ordera==orderb) Temp*=sqrt(2.);
            if(orderc==orderd) Temp*=sqrt(2.);
	    Oper+=Temp*fac;
	  }
      }
    return Oper;
  }

  
  double System_Tz_SM_M::jplus(int a,int b) const
  {
    if(Orders[a]!=Orders[b])
      return 0;
    if( (Orbitals[a].mm-Orbitals[b].mm)/2 ==1 )
      {
	int jj=Orbitals[b].jj;
	int mm=Orbitals[b].mm;
	double j=0.5*jj;
	double m=0.5*mm;
	return sqrt( j*(j+1)-m*(m+1) );
      }
    else
      return 0;
  }

  double System_Tz_SM_M::jminus(int a,int b) const
  {
    if(Orders[a]!=Orders[b])
      return 0;
    if( (Orbitals[b].mm-Orbitals[a].mm)/2 == 1 )
      {
	int jj=Orbitals[b].jj;
	int mm=Orbitals[b].mm;
	double j=0.5*jj;
	double m=0.5*mm;
	return sqrt( j*(j+1)-m*(m-1) );
      }
    else
      return 0;
  }
  
  void System_Tz_SM_M::printOrbitals() const
  {
    for(int i=0;i<Orbitals.size();i++)
      {
	std::cout<<i<<"\t"<<Orbitals[i].l<<"\t"<<Orbitals[i].jj<<"\t"<<Orbitals[i].mm<<"\t"<<Orbitals[i].tz<<"\t"<<Orbitals[i].e<<std::endl;
      }
  }


  void System_Tz_SM_M::printTwoBodyStates() const
  {
    for(int i=0;i<TwoBodyStates.size();i++)
      {
	std::cout<<TwoBodyStates[i].a<<"\t"<<TwoBodyStates[i].b<<"\t"<<TwoBodyStates[i].Tz<<"\t"<<TwoBodyStates[i].Par<<"\t"<<TwoBodyStates[i].M<<std::endl;
      }
  }
  //****************************************************************

  void Configs::addConfigs(int vN,int vZ,int MM,int Par,const vector<int> &restriction,double Emax)
  {
    setvNZ(vN,vZ);
    //generate partitions first;
    genPartitions(Par,restriction,Emax);
    addConfigs(MM);//from partitions, generate configs.
    dim=configs.size();
    // for(int i=0;i<dim;i++)
    //   std::cout<<i<<"\t"<<configs[i]<<std::endl;
  }

  void Configs::addConfigs(int MM)
  {
    vector<int> blocks;//every block has same quantum number except m
    vector<int> shift;//bits should be shifted of every block
    int bits=0;
    for(int i=0;i<Orbitals.size();)
      {
	shift.push_back(bits);
	int jj1=Orbitals[i].jj+1;
	blocks.push_back(jj1);
	bits+=jj1;
	i+=jj1;
      }
    vector< vector<Config> > BlockConfigs(blocks.size());//configs of every block
    for(int i=0;i<partitions.size();i++)
      {
	#pragma omp parallel for
	for(int block=0;block<blocks.size();block++)
	  {
	    BlockConfigs[block].clear();
	    //"partitions[i][block]" particles in  "blocks[block]" orbitals
	    Config start(0),last(0);
	    for(int j=0;j<partitions[i][block];j++)
	      {
		start[j]=1;
		last[blocks[block]-1-j]=1;
	      }
	    Config config(start);
	    while(true)
	      {
		BlockConfigs[block].push_back(config);
		if(config==last)
		  break;
		config=next(config);
	      }
	  }
	//combine configs between blocks
	vector<Config> configs_partition;//configs in a partition
	configs_partition.push_back(Config(0));
	for(int block=0;block<blocks.size();block++)
	  {
	    vector<Config> temp=configs_partition;
	    configs_partition.clear();
	    for(int k=0;k<BlockConfigs[block].size();k++)
	      for(int j=0;j<temp.size();j++)
		{
		  configs_partition.push_back( (BlockConfigs[block][k]<<shift[block])^temp[j] );
		}
	  }
	for(int j=0;j<configs_partition.size();j++)
	  {
	    if(getMM(configs_partition[j])==MM)
	      configs.push_back(configs_partition[j]);
	  }
      }
  }
  

  void Configs::genPartitions(int Par,const vector<int> &restriction,double Emax)
  {
    const vector<System_Tz_SM_M::VSystem::OrbitalType> & JOrbitals=pSystem->pSystem->Orbitals;
    int totalOrbitals=JOrbitals.size();
    int ProOrbitals=0,NeuOrbitals=0;//num of proton and neutron J scheme orbitals;
    for(int i=0;i<totalOrbitals;i++)
      {
	if(JOrbitals[i].tz==-1) ++ProOrbitals;
	else if(JOrbitals[i].tz==1) ++NeuOrbitals;
      }
    //generate partitions
    vector<int> n_low_limits(NeuOrbitals),n_up_limits(NeuOrbitals);
    vector<int> p_low_limits(ProOrbitals),p_up_limits(ProOrbitals);
    vector<int> n_type(NeuOrbitals),p_type(ProOrbitals);
    vector<int> n_l(NeuOrbitals),p_l(ProOrbitals);
    vector<DataType> n_e(NeuOrbitals),p_e(ProOrbitals);    
    //setup limits
    int in=0,ip=0;
    for(int i=0;i<totalOrbitals;i++)
      {
	if(JOrbitals[i].tz==1)
	  {
	    n_low_limits[in]=JOrbitals[i].low_limit;
	    n_up_limits[in]=JOrbitals[i].up_limit;
	    n_type[in]=JOrbitals[i].type;
	    n_l[in]=JOrbitals[i].l;
	    n_e[in]=JOrbitals[i].e;
	    //set to default value if input is absurd
	    if(n_low_limits[in]<0)  n_low_limits[in]=0;
	    if(n_up_limits[in]<0)  n_up_limits[in]=JOrbitals[i].jj+1;
	    ++in;	    
	  }
	else if(JOrbitals[i].tz==-1)
	  {
	    p_low_limits[ip]=JOrbitals[i].low_limit;
	    p_up_limits[ip]=JOrbitals[i].up_limit;
	    p_type[ip]=JOrbitals[i].type;
	    p_l[ip]=JOrbitals[i].l;
	    p_e[ip]=JOrbitals[i].e;	    
	    //set to default value if input is absurd
	    if(p_low_limits[ip]<0)  p_low_limits[ip]=0;
	    if(p_up_limits[ip]<0)  p_up_limits[ip]=JOrbitals[i].jj+1;
	    ++ip;	    
	  }
      }
    vector<Partition> Npartitions,Ppartitions;
    part(vN,NeuOrbitals,Npartitions,n_low_limits.data(),n_up_limits.data());
    part(vZ,ProOrbitals,Ppartitions,p_low_limits.data(),p_up_limits.data());
    //combine neutron proton partitions togather
    const size_t MaxNum=1000000;
//hubsh    const size_t MaxNum=100000;
    partitions.clear();
    partitions.reserve(min(Npartitions.size()*Ppartitions.size(),MaxNum));
    for(int i=0;i<Npartitions.size();i++)
      {
	int sum_n=0;//sum of particles on scattering states
	for(int k=0;k<NeuOrbitals;k++)
	  {
	    if(n_type[k]==1)
	      sum_n+=Npartitions[i][k];
	    if(sum_n>restriction[0])
	      break;
	  }
	if(sum_n>restriction[0]) continue;
        int par_temp=Par;
	DataType E_n=0;
	for(int k=0;k<NeuOrbitals;k++)
	  {
	    par_temp^=( (Npartitions[i][k]*n_l[k])%2 );
	    E_n+=Npartitions[i][k]*1.*n_e[k];
	  }
	if(ProOrbitals==0)
	  {
	    if(par_temp) continue;
	    if(Emax<E_n) continue;	    
	    partitions.push_back(Npartitions[i]);
	  }
        int par_n=par_temp; //hu
	for(int j=0;j<Ppartitions.size();j++)
	  {
	    int sum_p=0;//sum of particles on scattering states
	    for(int k=0;k<ProOrbitals;k++)
	      {
		if(p_type[k]==1)
		  sum_p+=Ppartitions[j][k];
		if(sum_p>restriction[1])
		  break;
	      }
	    if(sum_p>restriction[1]) continue;
            par_temp=par_n; //hu
	    DataType E_p=0;
	    for(int k=0;k<ProOrbitals;k++)
	      {
		par_temp^=( (Ppartitions[j][k]*p_l[k])%2 );
		E_p+=Ppartitions[j][k]*1.*p_e[k];		
	      }
	    if(par_temp) continue;
	    if(Emax<(E_n+E_p)) continue;

	    Partition temp(totalOrbitals);
	    for(int k=0;k<ProOrbitals;k++)
	      temp[k]=Ppartitions[j][k];
	    for(int k=0;k<NeuOrbitals;k++)
	      temp[k+ProOrbitals]=Npartitions[i][k];
	    partitions.push_back(temp);
	  }
      }
  }
  
  void Configs::setup_simple(int MM,int Par)
  {
    configs.clear();
    //generate proton and neutron configs separately
    bitset<numBits> start_p,start_n;
    bitset<numBits> last_p,last_n;
    for(int i=0;i<vZ;i++)
      {
	start_p[i]=1;
	last_p[maxProtons-1-i]=1;
      }
    for(int i=0;i<vN;i++)
      {
	start_n[i]=1;
	last_n[maxNeutrons-1-i]=1;
      }
    vector< bitset<numBits> > configs_p;
    vector< bitset<numBits> > configs_n;
    bitset<numBits> config(start_p);
    while(true)
      {
	configs_p.push_back(config);
	if(config==last_p)
	  break;
	config=next(config);
      }
    config=start_n;
    while(true)
      {
	configs_n.push_back(config);
	if(config==last_n)
	  break;
	config=next(config);
	
      }
    //combine proton and neutron configs togather and choose those the total angular momentum projection is MM/2
    for(int i=0;i<configs_n.size();i++)
      for(int j=0;j<configs_p.size();j++)
	{
	  bitset<numBits> temp(configs_p[j]);
	  for(int k=0;k<maxNeutrons;k++)
	    temp[maxProtons+k]=configs_n[i][k];
	  if(getMM(temp)==MM && getPar(temp)==Par)
	    configs.push_back(temp);
	}
    dim=configs.size();
    // for(int i=0;i<configs.size();i++)
    //   std::cout<<i<<"\t"<<configs[i]<<std::endl;
  }

  int Configs::getMM(const bitset<numBits> & config) const
  {
    int MM=0;
    for(int i=0;i<maxProtons+maxNeutrons;i++)
      {
	if(config.test(i))
	  MM+=Orbitals[i].mm;
      }
    return MM;
  }

  int Configs::getPar(const bitset<numBits> & config) const
  {
    int Par=0;
    for(int i=0;i<maxProtons+maxNeutrons;i++)
      {
	if(config.test(i))
	  Par+=Orbitals[i].par;
      }
    return Par%2;
  }
  DataType Configs::getME(int i,int j) const
  {
    DataType temp(0);
    bitset<numBits> diff=configs[i]^configs[j];
    bitset<numBits> diffket=diff&configs[j];
    bitset<numBits> diffbra=diff&configs[i];
    int diffCount=diff.count();
    vector<int> pos;
    if(diffCount > 4) return temp;//we only have twobody interaction
    else if(diffCount==4)
      {
	findones(diffket,pos,2);
	int c=pos[1];
	int d=pos[0];
	findones(diffbra,pos,2);
	int a=pos[1];
	int b=pos[0];
	temp=get2B(a,b,c,d);
	int sign=count(configs[i],a,b) + count(configs[j],c,d);
	temp=sign%2? -temp:temp;
      }
    else if(diffCount==2)
      {
    	findones(diffket,pos,1);
    	int b=pos[0];
    	findones(diffbra,pos,1);
    	int a=pos[0];
    	int sign=count(configs[i],a,0) + count(configs[j],b,0);
    	temp=get1B(a,b);
    	for(int k=0;k<maxProtons+maxNeutrons;k++)
    	  {
    	    if(configs[i].test(k))
    	      temp+=get2B(a,k,b,k);
    	  }
    	temp=sign%2?-temp:temp;
      }
    else if(diffCount==0)
      {
	for(int k=0;k<maxProtons+maxNeutrons;k++)
	  {
	    if(configs[i].test(k))
	      {
		temp+=get1B(k,k);
		for(int l=0;l<maxProtons+maxNeutrons;l++)
		  {
		    if(configs[i].test(l))
		      {
			temp+=0.5*get2B(k,l,k,l);
		      }
		  }
	      }
	  }
      }
    return temp;
  }


  double Configs::getJ2(int i,int j) const
  {
    double temp(0);
    bitset<numBits> diff=configs[i]^configs[j];
    bitset<numBits> diffket=diff&configs[j];
    bitset<numBits> diffbra=diff&configs[i];
    int diffCount=diff.count();
    vector<int> pos;
    //J^2= J-J+ + Jz^2 +Jz
    if(diffCount > 4) return temp;//we only have twobody interaction
    else if(diffCount==4)
      {
	findones(diffket,pos,2);
	int c=pos[1];
	int d=pos[0];
	findones(diffbra,pos,2);
	int a=pos[1];
	int b=pos[0];
	temp=pSystem->jminus(b,d)*pSystem->jplus(a,c)+pSystem->jminus(a,c)*pSystem->jplus(b,d);
	temp-=pSystem->jminus(a,d)*pSystem->jplus(b,c)+pSystem->jminus(b,c)*pSystem->jplus(a,d);
	int sign=count(configs[i],a,b) + count(configs[j],c,d);
	temp=sign%2? -temp:temp;
      }
    // else if(diffCount==2)//it's not necessary if different m orbitals of every shell are put together
    //   {
    // 	findones(diffket,pos,1);
    // 	int b=pos[0];
    // 	findones(diffbra,pos,1);
    // 	int a=pos[0];
    // 	int sign=count(configs[i],a,0) + count(configs[j],b,0);
    // 	for(int k=0;k<maxProtons+maxNeutrons;k++)
    // 	  {
    // 	    if(configs[i].test(k))
    // 	      temp-=pSystem->jminus(k,b)*pSystem->jplus(a,k);
    // 	    else
    // 	      temp+=pSystem->jminus(a,k)*pSystem->jplus(k,b);
	    
    // 	  }
    // 	temp=sign%2?-temp:temp;
    //   }    
    else if(diffCount==0)
      {
	//J-J+ matrix element between configs = sum_{hp}(J-)_hp (J+)_ph , p(h) is orbital (not) occupied by this config. because the orbital with same j but with adjacent  m is adjacent, we simplify to use one loop only
	for(int ii=0;ii<maxNeutrons+maxProtons;ii++)
	  {
	    if(configs[i].test(ii))
	      {
		int jj=ii+1;
		if(jj>=maxNeutrons+maxProtons) continue;
		// for(int jj=0;jj<maxNeutrons+maxProtons;jj++)
		//   {
		if(configs[i].test(jj)) continue;
		double j=0.5*Orbitals[ii].jj;
		double m=0.5*Orbitals[ii].mm;
		temp+=j*(j+1)-m*(m+1);
		//		  }
	      }
	  }
	double M=getMM(configs[i])*0.5;
	//Jz^2+Jz
	temp+=M*(M+1);
      }
    return temp;
  }
  
  DataType Configs::getOper(int i,int j) const
  {
    DataType temp(0);
    bitset<numBits> diff=configs[i]^configs[j];
    bitset<numBits> diffket=diff&configs[j];
    bitset<numBits> diffbra=diff&configs[i];
    int diffCount=diff.count();
    vector<int> pos;
    if(diffCount==4)
      {
	findones(diffket,pos,2);
	int c=pos[1];
	int d=pos[0];
	findones(diffbra,pos,2);
	int a=pos[1];
	int b=pos[0];
	temp=pSystem->get2BOper(a,b,c,d);
	int sign=count(configs[i],a,b) + count(configs[j],c,d);
	temp=sign%2? -temp:temp;
      }
    else if(diffCount==2)
      {
    	findones(diffket,pos,1);
    	int b=pos[0];
    	findones(diffbra,pos,1);
    	int a=pos[0];
    	int sign=count(configs[i],a,0) + count(configs[j],b,0);
	temp=pSystem->get1BOper(a,b);
    	for(int k=0;k<maxProtons+maxNeutrons;k++)
    	  {
    	    if(configs[i].test(k))
    	      temp+=pSystem->get2BOper(a,k,b,k);
    	  }
    	temp=sign%2?-temp:temp;
      }
    else if(diffCount==0)
      {
	for(int k=0;k<maxProtons+maxNeutrons;k++)
	  {
	    if(configs[i].test(k))
	      {
		temp+=pSystem->get1BOper(k,k);		    		
		for(int l=0;l<maxProtons+maxNeutrons;l++)
		  {
		    if(configs[i].test(l))
		      {
			temp+=0.5*pSystem->get2BOper(k,l,k,l);
		      }
		  }
	      }
	  }
      }
    return temp;
  }

  
  //recursion algorithm to part 'num' particles into numorbs parts.
  void Configs::part(int num,int numorbs,vector<Partition> & partitions,int *low_limits,int *up_limits)
  {
    partitions.clear();
    //divide and conquer
    int numorbs1=numorbs/2;
    int numorbs2=numorbs-numorbs1;
    if(numorbs==0) return;
    if(numorbs==1)
      {
	if(num>=low_limits[0] && num<=up_limits[0])
	  partitions.push_back(Partition(1,num));
	return;
      }
    for(int num2=0;num2<=num;num2++)
      {
	int num1=num-num2;
	vector<Partition> temp1,temp2;
	part(num1,numorbs1,temp1,low_limits,up_limits);
	part(num2,numorbs2,temp2,low_limits+numorbs1,up_limits+numorbs1);
	for(auto part1:temp1)
	  for(auto part2:temp2)
	    {
	      Partition combined_part(numorbs);
	      for(int i=0;i<numorbs1;i++)
		combined_part[i]=part1[i];
	      for(int i=0;i<numorbs2;i++)
		combined_part[i+numorbs1]=part2[i];
	      partitions.push_back(combined_part);
	    }
      }
  }
  
  //****************************************************************


  

}

