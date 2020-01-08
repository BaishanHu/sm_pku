#include<iostream>
#include<fstream>
#include<boost/mpi.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
using namespace boost::mpi;
#include"sm_system.h"
#include"serialization.h"
#include"sm_solver.h"

using std::cout;
using std::cin;
using std::endl;
using std::ofstream;

using namespace SM;

int main()
{
  environment env;
  communicator world;
  System_Tz_SM effV;
  System_Tz_SM_M effVm;
  string outputfile;
  timer clock; 
  if(world.rank()==0)
    {
      string spfile,vfile;
      cout<<"input <sp file> and <int file>:\n";
      cin>>spfile>>vfile>>outputfile;
      cout<<spfile<<endl;
      cout<<vfile<<endl;
      effV.setup(spfile.c_str(),vfile.c_str());
      effVm.setup(&effV);
    }

  int vN,vZ,MM,Par,num;
  double Emax;
  bool printOccuNums,calOper;
  vector<int> restriction(2);
  if(world.rank()==0)
    {
      cout<<"input <vN> \t <vZ> \t Emax \t <restriction> \t <MM> \t <Par> \t <num>\n";
      cin>>vN>>vZ>>Emax>>restriction[0]>>restriction[1]>>MM>>Par>>num;
      cin>>printOccuNums;
      cin>>calOper;
      cout<<"\t"<<vN<<"\t"<<vZ<<"\t"<<Emax<<"\t[ "<<restriction[0]<<"\t"<<restriction[1]<<" ]\t"<<MM<<"\t"<<Par<<"\t"<<num<<endl;
      cout<<outputfile<<endl;
      if(calOper)
	{
	  string Operf1B,Operf2B;
	  cin>>Operf1B>>Operf2B;
	  effV.setupOper1B(Operf1B);
	  effV.setupOper2B(Operf2B);
	}
      
    }
  broadcast(world,effVm,0);
  broadcast(world,vN,0);
  broadcast(world,vZ,0);
  broadcast(world,Emax,0);  
  broadcast(world,restriction,0);
  broadcast(world,MM,0);
  broadcast(world,Par,0);
  broadcast(world,num,0);
  broadcast(world,printOccuNums,0);
  broadcast(world,calOper,0);  

/* hubsh-test
  Configs configs(&effVm);
  configs.addConfigs(vN,vZ,MM,Par,restriction,Emax);
  for(auto p:configs.partitions)
  {
    for(auto i:p)
    cout<<i<<",";
    cout<<endl;
  } 
  for(auto c:configs.configs)
  {
    cout<<c<<endl;
  }
  return 0;
hubsh-test */

  SMSolver smsolver(world,&effVm);//,vN,vZ
  smsolver.addConfigs(vN,vZ,MM,Par,restriction,Emax);//MM,Par,restriction

  bool savevecs=printOccuNums;
//  smsolver.diag(calOper);  
  smsolver.diag_lanczos(num,savevecs,calOper);//num,savevecs,step=1

  ofstream fout(outputfile.c_str());
  smsolver.printStates(num,fout);
  if(calOper)
    smsolver.printOperInfo(num,fout);
  if(printOccuNums)
    smsolver.printOccuNums(num);
  if(world.rank()==0)
    cout<<"total time: "<<clock.elapsed()/60.<<" min"<<endl;
  return 0;
}
