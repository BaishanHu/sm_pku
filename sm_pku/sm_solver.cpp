#include<iostream>
#include<iomanip>
#include"serialization.h"
#include"sm_solver.h"

namespace SM
{
  void SMSolver::addConfigs(int vN,int vZ,int MM,int Par,const vector<int> &restriction,double Emax)
  {
    int rank=world.rank();

    if(rank==0)
      {
	ConfigMM=MM;
	configs.addConfigs(vN,vZ,MM,Par,restriction,Emax);
      }
    broadcast(world,configs,0);
  
    if(rank==0)
      std::cout<<"total dim:\t"<<configs.dim<<std::endl;
  }
  
  void SMSolver::diag(bool calOper)
  {
    if(world.size()!=1)
      {
	std::cerr<<"!!!Direct diagonalization requires only 1 process!\n";
	environment::abort(0);
      }
    if(world.rank()==0)
      {
	int dim=configs.dim;
	numStates=dim;
	if(dim==0)
	  {
	    eigVals.resize(0);
	    eigVecs.resize(0,0);
	    J2s.resize(0);
	    Oper.resize(0,0);
	  }
	else
	  {
	    Mat mat(dim,dim);
	    for(int i=0;i<dim;i++)
	      for(int j=i;j<dim;j++)
		{
		  mat(i,j)=configs.getME(i,j);
		  mat(j,i)=mat(i,j);
		}
	    ComplexEigenSolver<Mat> es(mat);
	    //SelfAdjointEigenSolver<Mat> es(mat);
	    eigVals=es.eigenvalues();
	    eigVecs=es.eigenvectors();

	    //norm
	    for(int i=0;i<dim;i++)
	      {
		DataType norm=eigVecs.col(i).transpose()*eigVecs.col(i);
		eigVecs.col(i)/=sqrt(norm);
	      }

	    //    sort
	    quickSort(eigVals,eigVecs,0,dim-1);

	    for(int i=0;i<dim;i++)
	      for(int j=i;j<dim;j++)
		{
		  mat(i,j)=configs.getJ2(i,j);
		  mat(j,i)=mat(i,j);
		}
	    J2s=(eigVecs.transpose()*mat*eigVecs).diagonal();


	    if(calOper)
	      {
		Mat OperMat(dim,dim);
		for(int i=0;i<dim;i++)
		  for(int j=i;j<dim;j++)
		    {
		      OperMat(i,j)=configs.getOper(i,j);
		      OperMat(j,i)=OperMat(i,j);
		    }
		Oper=eigVecs.transpose()*OperMat*eigVecs;
		//std::cout<<Oper.diagonal()<<endl;
	      }
	  }
      }
    
    //needed for calOccuNums
    numOfparts=1;
    partdim=eigVals.size();
    diagcom=world;
  }

  void SMSolver::diag_lanczos(int num,bool savevecs,bool calOper,int step)
  {
    int rank=world.rank();
    totaldim=configs.dim;
    if(totaldim==0)
      {
	eigVecs.resize(0,0);
	eigVals.resize(0);
	J2s.resize(0);
	numStates=0;
	numOfparts=1;	
	return;
      }

    if(num>totaldim) num=totaldim;
    numStates=num;
    
    LanczosIter<Configs> lanczositer(configs,world);
    
    int lanczosdim=num;
    double sigma;
    Vec eigVals_last;
    lanczositer.compute(eigVals,lanczosdim,false);    
    do
      {
    	eigVals_last=eigVals;
	lanczosdim+=step;	
    	lanczositer.compute(eigVals,lanczosdim,false);

    	if(rank==0)
    	  {
    	    sigma=(eigVals.head(num)-eigVals_last.head(num)).norm()/num;
    	    //show convergance
	    std::cout<<lanczosdim<<"\t"<<eigVals.head(num).transpose()<<"\t"<<sigma<<std::endl;	    
    	  }
    	broadcast(world,sigma,0);
      }while(sigma>1e-5);
  

    lanczositer.compute(eigVals,lanczosdim,true);
    lanczositer.calJ2(J2s,numStates);

    if(calOper)
      {
	lanczositer.setupOper();
	Mat OperMat;
	Mat Q;
	lanczositer.calOperMat(OperMat);	
	lanczositer.getVecs(Q,numStates);
	if(rank==0)
	  {
	    Oper=Q.transpose()*OperMat*Q;
	    //std::cout<<Oper.diagonal()<<endl;
	  }
      }

    //save eigenvectors, needed by calOccuNums
    if(savevecs)
      {
	lanczositer.getEigVecs(eigVecs,numStates);
	numOfparts=lanczositer.numOfparts;
	partdim=lanczositer.partdim;
	diagcom=lanczositer.diagcom;
      }
  }
  
  void SMSolver::calOccuNums(int state)
  {
    int rank=world.rank();
    if(rank<numOfparts)
      {
	int totalOrbitals=configs.Orbitals.size();
	OccuNums.resize(totalOrbitals);
	for(int i=0;i<totalOrbitals;i++)
	  {
	    OccuNums[i]=0;
	  }

	for(int j=0;j<totalOrbitals;j++)
	  {
	    DataType temp=0;
	    for(int i=0;i<eigVecs.rows();i++)
	      {
		int index=rank*partdim;
		if( configs.configs[index+i].test(j) )
		  temp+=pow(eigVecs(i,state),2);
	      }
	    DataType sum;
	    reduce(diagcom,temp,sum,std::plus<DataType>(),0);
	    if(rank==0)
	      OccuNums[j]=sum;
	  }
      }
  }
  
  void SMSolver::printOccuNums(int num)
  {
    if(num==-1||num > numStates) num=numStates;//default num=-1, print all
    for(int state=0;state<num;state++)
      {
	calOccuNums(state);
	if(world.rank()==0)
	  {
	    using std::cout;
	    using std::setw;
	    using std::endl;

	    cout<<"for "<<state<<" th state\n";

	    // cout<<setw(4)<<"l"<<"\t"<<setw(4)<<"jj"<<"\t"<<setw(4)<<"mm"<<"\t"<<setw(4)<<"tz"<<"\t"<<setw(10)<<"OccuNums"<<endl;
	    // for(int i=0;i<OccuNums.size();i++)
	    //   {
	    //     cout<<setw(4)<<configs.Orbitals[i].l<<"\t"<<setw(4)<<configs.Orbitals[i].jj<<"\t"<<setw(4)<<configs.Orbitals[i].mm<<"\t"<<setw(4)<<configs.Orbitals[i].tz<<"\t"<<setw(10)<<OccuNums[i]<<endl;
	    //   }
	    
	    
	    cout<<setw(4)<<"l"<<"\t"<<setw(4)<<"jj"<<"\t"<<setw(4)<<"\t"<<setw(4)<<"tz"<<"\t"<<setw(10)<<"OccuNums"<<endl;
	    int i=0;
	    while(i<OccuNums.size())
	      {
		int jj=configs.Orbitals[i].jj;
		DataType OccuNum=0;
		for(int k=0;k<jj+1;k++)
		  OccuNum+=OccuNums[i+k];
		cout<<setw(4)<<configs.Orbitals[i].l<<"\t"<<setw(4)<<jj<<"\t"<<setw(4)<<configs.Orbitals[i].tz<<"\t"<<setw(10)<<OccuNum<<endl;
		i+=jj+1;
	      }
	  }
      }
  }



  void SMSolver::printStates(int num, std::ostream&fout)
  {
    if(world.rank()==0)
      {
	using std::cout;
	using std::setw;
	using std::endl;
	cout.precision(6);
	cout.setf( std::ios::fixed, std:: ios::floatfield );
	
	if(num==-1||num > numStates) num=numStates;//default num=-1, print all
      
	//cout<<"Num"<<"\t"<<setw(25)<<"Energy"<<"\t"<<setw(25)<<"J(J+1)"<<endl;
        cout<<"Num"<<"\t"<<setw(17)<<"J(J+1)"<<"\t"<<setw(12)<<"J"<<"\t"<<setw(22)<<"Energy"<<"\t"<<setw(28)<<"Ex-E0"<<endl;
        fout<<"Num"<<"\t"<<setw(17)<<"J(J+1)"<<"\t"<<setw(12)<<"J"<<"\t"<<setw(20)<<"Energy"<<"\t"<<setw(28)<<"Ex-E0"<<endl;
	
	for(int i=0;i<num;i++)
	  {
//	    cout<<i<<"\t"<<setw(25)<<eigVals(i)<<"\t"<<setw(25)<<(i==0?0:eigVals(i)-eigVals(0))<<"\t"<<setw(6)<<int(sqrt(4*J2s(i).real()+1)-1 + 0.5)/2.<<"\t"<<setw(25)<<J2s(i)<<endl;	    
//	    fout<<i<<"\t"<<setw(25)<<eigVals(i)<<"\t"<<setw(25)<<(i==0?0:eigVals(i)-eigVals(0))<<"\t"<<setw(6)<<int(sqrt(4*J2s(i).real()+1)-1 + 0.5)/2.<<"\t"<<setw(25)<<J2s(i)<<endl;	    
            cout<<i<<"\t"<<setw(12)<<J2s(i).real()<<"\t"<<setw(12)<<J2s(i).imag()<<" | \t"<<setw(3)<<int(sqrt(4*J2s(i).real()+1)-1 + 0.5)/2.<<" | \t"<<setw(12)<<eigVals(i).real()<<"\t"<<setw(12)<<eigVals(i).imag()<<" | \t"<<setw(12)<<eigVals(i).real()-eigVals(0).real()<<"\t"<<setw(12)<<eigVals(i).imag()-eigVals(0).imag()<<endl;
            fout<<i<<"\t"<<setw(12)<<J2s(i).real()<<"\t"<<setw(12)<<J2s(i).imag()<<"\t"<<setw(3)<<int(sqrt(4*J2s(i).real()+1)-1 + 0.5)/2.<<"     "<<setw(12)<<eigVals(i).real()<<"\t"<<setw(12)<<eigVals(i).imag()<<"     "<<setw(12)<<eigVals(i).real()-eigVals(0).real()<<"\t"<<setw(12)<<eigVals(i).imag()-eigVals(0).imag()<<endl;
	    //		cout<<i<<"\t"<<setw(10)<<eigVals(i)<<"\t"<<setw(8)<<int(sqrt(4*J2s(i)+1)-1 + 0.5)/2.<<endl;
	    //cout<<i<<"\t"<<setw(25)<<eigVals(i)<<"\t"<<setw(25)<<J2s(i)<<endl;
	    //	cout<<i<<"\t"<<setw(10)<<eigVals(i)<<"\t"<<setw(8)<<int(J2s(i).real() + 0.5)<<endl;
	  }
      }
  }

  void SMSolver::printOperInfo(int num, std::ostream&fout)
  {
    if(world.rank()==0)
      {
	using std::cout;
	using std::setw;
	using std::endl;
	using std::setprecision;
	if(num==-1||num > numStates) num=numStates;//default num=-1, print all
	int lambda1B=configs.pSystem->pSystem->lambda1B;
	int lambda2B=configs.pSystem->pSystem->lambda2B;
	int lambda;
	if( (lambda1B!=-1&&lambda2B!=-1) && (lambda1B!=lambda2B) )
	  {
	    cout<<"lambda in 1B and 2B doesn't match\n";
	    return;
	  }
	if(lambda1B==-1&&lambda2B==-1)
	  {
	    cout<<"please input operator file if you intend to cal. oper!\n";
	    return;
	  }
	else
	  {
	    lambda=(lambda1B==-1)?lambda2B:lambda1B;
	  }
	cout<<"Below are transition table\n";
        cout<<"state1  state2"<<"\t"<<setw(16)<<"dE"<<"\t"<<setw(32)<<"B(2>1) or u or Q"<<"\t"<<setw(22)<<"B(1>2) or u or Q"<<endl;
	fout<<"Below are transition table\n";
        fout<<"state1  state2"<<"\t"<<setw(16)<<"dE"<<"\t"<<setw(32)<<"B(2>1) or u or Q"<<"\t"<<setw(22)<<"B(1>2) or u or Q"<<endl;
	cout.precision(6);
	cout.setf( std::ios::fixed, std:: ios::floatfield );
	for(int i=0;i<num;i++)
	  {
	    int jji=int(sqrt(4*J2s(i).real()+1)-1 + 0.5);
	    for(int j=i;j<num;j++)
	      {
		int jjj=int(sqrt(4*J2s(j).real()+1)-1 + 0.5);
		if( (abs(jji-jjj)>2*lambda1B || (jji+jjj)<2*lambda1B) && (abs(jji-jjj)>2*lambda2B || (jji+jjj)<2*lambda2B) ) continue;
		double fac1=cg(jjj,2*lambda,jji,ConfigMM,0,ConfigMM)/sqrt(jji+1.);
		double fac2=cg(jji,2*lambda,jjj,ConfigMM,0,ConfigMM)/sqrt(jjj+1.);
		DataType qf,qb;
		if(abs(fac1)<1e-8||abs(fac2)<1e-8)
		  {
		    cout<<"Warning:   for "<<i<<"<=>"<<j<<" current cal. is invalid, change MM and try again\n";
		    qf=nan("");
		    qb=nan("");
		  }
		else
		  {
		    qf=Oper(i,j)/fac1;
		    qb=Oper(j,i)/fac2;
//hu
                    if(jji == jjj and lambda ==1)
                    {
		      qf=sqrt(4.0*Pi/3.0)*sqrt(jji/2./(jji/2.+1.0)/(jji+1.0))*qf;
		      qb=sqrt(4.0*Pi/3.0)*sqrt(jji/2./(jji/2.+1.0)/(jji+1.0))*qb;
                    }
                    else if(jji == jjj and lambda ==2)
                    {
		      qf=sqrt(16.0*Pi/5.0)*sqrt(jji/2.0*(jji-1.0)/(jji/2.0+1.0)/(jji+1.0)/(jji+3.0))*qf;
		      qb=sqrt(16.0*Pi/5.0)*sqrt(jji/2.0*(jji-1.0)/(jji/2.0+1.0)/(jji+1.0)/(jji+3.0))*qb;
                    }
                    else if(jji == jjj and lambda ==0)
                    {
		      qf=qf*fac1;
		      qb=qb*fac2;
                    }
                    else
                    {
		      qf=qf*qf/(jjj+1.);
		      qb=qb*qb/(jji+1.);
                    }
//hu
		  }
		DataType de=eigVals(j)-eigVals(i);
		cout<<i<<"  "
		    <<"\t"<<j<<"  "
		    <<"\t"<<setw(12)<<de.real()<<"  "
		    <<"\t"<<setw( 6)<<de.imag()<<"  "
		    <<"\t"<<setw(12)<<qf.real()<<"  "
		    <<"\t"<<setw( 6)<<qf.imag()<<"  "
		    <<"\t"<<setw(12)<<qb.real()<<"  "
		    <<"\t"<<setw( 6)<<qb.imag()
		    <<endl;
		fout<<i<<"  "
                    <<"\t"<<j<<"  "
                    <<"\t"<<setw(12)<<de.real()<<"  "
                    <<"\t"<<setw( 6)<<de.imag()<<"  "
                    <<"\t"<<setw(12)<<qf.real()<<"  "
                    <<"\t"<<setw( 6)<<qf.imag()<<"  "
                    <<"\t"<<setw(12)<<qb.real()<<"  "
                    <<"\t"<<setw( 6)<<qb.imag()
		    <<endl;
	      }
	  }
      }
  }


}
