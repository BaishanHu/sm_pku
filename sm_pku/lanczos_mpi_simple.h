#ifndef _LANCZOS_H_
#define _LANCZOS_H_
#include<boost/mpi.hpp>//boost mpi library
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include<eigen3/Eigen/Eigen>//matrix library
#include<eigen3/Eigen/Sparse>//matrix library
#include<cmath>
#include<limits>
#include"serialization.h"
using std::abs;
using std::numeric_limits;
using namespace boost::mpi;


using Eigen::Matrix;
using Eigen::Dynamic;
using Eigen::SparseMatrix;
using Eigen::Triplet;


namespace SM
{

  //****************************************************************
  //sort algorithm, sort according to eigvals and reorder the eigvecs respectively
  template<typename Scalar>
    void quickSort(Matrix<Scalar,Dynamic,1> & eigvals,Matrix<Scalar,Dynamic,Dynamic> &eigvecs,int startPos, int endPos) 
    { 
      int i,j; 
      Scalar key; 
      key=eigvals[startPos];
      Matrix<Scalar,Dynamic,1> Vec=eigvecs.col(startPos);
      i=startPos; 
      j=endPos;
      while(i<j) 
	{ 
	  while(!(eigvals[j]<key)&& i<j) --j; 
	  eigvals[i]=eigvals[j];
	  eigvecs.col(i)=eigvecs.col(j);
	  while(eigvals[i]<key && i<j) ++i; 
	  eigvals[j]=eigvals[i];
	  eigvecs.col(j)=eigvecs.col(i);
	} 
      eigvals[i]=key;
      eigvecs.col(i)=Vec; 
      if(i-1>startPos) quickSort(eigvals,eigvecs,startPos,i-1); 
      if(endPos>i+1) quickSort(eigvals,eigvecs,i+1,endPos); 
    } 

  template<typename Scalar>
    void quickSort(Matrix<Scalar,Dynamic,1> & eigvals,int startPos, int endPos) 
    { 
      int i,j; 
      Scalar key; 
      key=eigvals[startPos];
      i=startPos; 
      j=endPos;
      while(i<j) 
	{ 
	  while(!(eigvals[j]<key)&& i<j) --j; 
	  eigvals[i]=eigvals[j];
	  while(eigvals[i]<key && i<j) ++i; 
	  eigvals[j]=eigvals[i];
	} 
      eigvals[i]=key;
      if(i-1>startPos) quickSort(eigvals,startPos,i-1); 
      if(endPos>i+1) quickSort(eigvals,i+1,endPos); 
    } 


  //****************************************************************

  
  //Lanczos algorithm, MatGen has variable dim and method getME(i,j),getJ2(i,j),getOper(i,j)
  template<typename MatGen>
    class LanczosIter
    {
    public:
      typedef DataType Scalar;
      typedef Matrix<Scalar,Dynamic, Dynamic> MatrixType;
      typedef Matrix<Scalar,Dynamic, 1> VectorType;
      //typedef MatrixType MatrixAType;
      typedef SparseMatrix<Scalar> MatrixAType;//use sparse Matrix
      
      LanczosIter(const MatGen & _matgen,const communicator& _world);//:matgen(_matgen),dim(0),totaldim(_matgen.dim),world(_world)
	  
      void clear();
      void setV0(const VectorType & _V0)
      {
	if(world.rank()==0)
	  {
	    V0=_V0;
	    V0/=sqrt(V0.transpose()*V0);
	  }
	broadcast(world,V0,0);
	clear();      
      }
      void AmultiplyV(const MatrixAType& _partA,const VectorType & _partV,VectorType & _partAV);//finish AV=A*V using distributed partA and partV, and distribute AV into partAV in diagonal processors.
      void udotv(const VectorType & _partu,const VectorType & _partv,Scalar & _uv);//cal. u dot v using distributed partu,partv on diagonal processors, the result scalar is on process with rank 0
      bool next();//get next lanczos vector
      int getAlphaBeta(VectorType & diag, VectorType & subdiag,int _dim);//return the tridiagnal matrix in form of diag and subdiag
      vector<Scalar> alpha;//meaningful only for rank 0
      vector<Scalar> beta;//meaningful only for rank 0

      void compute(VectorType & m_eivalues,int _max_dim=-1,bool yesvecs=true);
      void calJ2(VectorType & J2s,int numStates);

      void setupOper();
      void calOperMat(MatrixType & OperMat);
      void getVecs(MatrixType& eigVecs,int numStates);//eigvecs under lanczos basis
      
      void getEigVecs(MatrixType & eigVecs,int numStates);//eigvecs under config. basis

      
      vector<VectorType> partV;//meaningful for rank 0~numOfparts-1
      VectorType partu;//the rest,meaningful for rank 0~numOfparts-1
      int dim;//the actual dim of alpha,beta, once next() is called successfully, dim increases 1

      MatrixAType partA;//part of the whole matrix, meaningful for all rank
      MatrixAType partJ2;//for J2 cal.
      MatrixAType partOper;//for other Opers cal.
      int partdim;//dim of partA
      int numOfparts;

      void reportMemory();
      unsigned long int numAnonZeros;//total number of nonzero matrix elements of A, meaningful in rank 0.
      unsigned long int numJ2nonZeros;

      const int totaldim;
      const MatGen & matgen;
      VectorType V0;


      MatrixType MEJ2;//J2 matrix element under lanczos vectors, only meaningful on rank 0 processor

      MatrixType Q;//only meaningful on rank 0 processor, wavefun under lanczos vectors
    
      communicator rowcom;
      communicator colcom;
      communicator diagcom;
     
      const communicator & world;
    };



  template<typename MatGen>
    LanczosIter<MatGen>::LanczosIter(const MatGen & _matgen,const communicator& _world):matgen(_matgen),dim(0),totaldim(_matgen.dim),world(_world)
  {
    //setup start vector V0
    if(world.rank()==0)
      {
	srand((unsigned int) time(0));
	V0=VectorType::Random(totaldim);
	//V0=MatrixType::Ones(totaldim,1);//another choice
	V0/=sqrt(V0.transpose()*V0);
      }
    broadcast(world,V0,0);
      
    //divide task
    int totalProcessors=world.size();
    numOfparts=(int(sqrt(1.+8*totalProcessors))-1)/2;//totalProcessors=numOfparts*(numOfparts+1)/2
    if(world.rank()==0 && (totalProcessors!=numOfparts*(numOfparts+1)/2) )
      {
	std::cerr<<"==============================================================================\n";
	std::cerr<<"!!!!Please set total number of processors as n(n+1)/2 with n a integer!!!!\n";
	std::cerr<<"look at the table below:\n";
	std::cerr<<"i"<<"\t"<<"i*(i+1)/2\n";
	for(int i=1;i<40;i++)
	  std::cerr<<i<<"\t"<<i*(i+1)/2<<std::endl;
	std::cerr<<"==============================================================================\n";	
	environment::abort(0);
      }
    partdim=totaldim/numOfparts;

    //setup communicators for row,col and diagonal, partA maxtrix is also setup for every processor
    for(int col=0;col<numOfparts;col++)
      for(int row=col;row<numOfparts;row++)
	{
	  int rank=col+((2*numOfparts+1-row+col)*(row-col))/2;
	  if(world.rank()==rank)
	    {
	      int startrow=row*partdim;
	      int endrow=(row+1)*partdim;
	      if(row==numOfparts-1)
		endrow=totaldim;
	      int startcol=col*partdim;
	      int endcol=(col+1)*partdim;
	      if(col==numOfparts-1)
		endcol=totaldim;

	      /*
	      //setup matrix
	      partA.resize(endrow-startrow,endcol-startcol);
	      partJ2.resize(endrow-startrow,endcol-startcol);
	      if(world.rank()<numOfparts)
	        {
	          for(int i=startrow;i<endrow;i++)
	            for(int j=i;j<endcol;j++)
		      {
			partA(i-startrow,j-startcol)=matgen.getME(i,j);
			partA(j-startcol,i-startrow)=partA(i-startrow,j-startcol);
			  
			partJ2(i-startrow,j-startcol)=matgen.getJ2(i,j);
			partJ2(j-startcol,i-startrow)=partJ2(i-startrow,j-startcol);
		      }
	        }
	      else
	        {
	          for(int i=startrow;i<endrow;i++)
	            for(int j=startcol;j<endcol;j++)
		      {
			partA(i-startrow,j-startcol)=matgen.getME(i,j);
			  
			partJ2(i-startrow,j-startcol)=matgen.getJ2(i,j);
		      }
	        }
	      */
	      
	      //setup matrix,using sparse matrix
	      partA.resize(endrow-startrow,endcol-startcol);
	      partJ2.resize(endrow-startrow,endcol-startcol);
	      if(world.rank()<numOfparts)
		{
		  vector< Triplet<Scalar> > nonzeros;
		  for(int i=startrow;i<endrow;i++)
		    for(int j=i;j<endcol;j++)
		      {
			Scalar val=matgen.getME(i,j);
			if(abs(val)<1e-10) continue;	  			
			int ii=i-startrow;
			int jj=j-startcol;
			nonzeros.push_back(Triplet<Scalar>(ii,jj,val));
			if(ii!=jj)
			  nonzeros.push_back(Triplet<Scalar>(jj,ii,val));			  
		      }
		  partA.setFromTriplets(nonzeros.begin(),nonzeros.end());
		      
		  nonzeros.clear();
		  for(int i=startrow;i<endrow;i++)
		    for(int j=i;j<endcol;j++)
		      {
			Scalar val=matgen.getJ2(i,j);
			if(abs(val)<1e-10) continue;	  			
			int ii=i-startrow;
			int jj=j-startcol;
			nonzeros.push_back(Triplet<Scalar>(ii,jj,val));
			if(ii!=jj)
			  nonzeros.push_back(Triplet<Scalar>(jj,ii,val));			  
		      }
		  partJ2.setFromTriplets(nonzeros.begin(),nonzeros.end());
		}
	      else
		{
		  vector< Triplet<Scalar> > nonzeros;		      
		  for(int i=startrow;i<endrow;i++)
		    for(int j=startcol;j<endcol;j++)
		      {
			Scalar val=matgen.getME(i,j);
			if(abs(val)<1e-10) continue;	  			
			int ii=i-startrow;
			int jj=j-startcol;
			nonzeros.push_back(Triplet<Scalar>(ii,jj,val));
		      }
		  partA.setFromTriplets(nonzeros.begin(),nonzeros.end());
		      
		  nonzeros.clear();
		  for(int i=startrow;i<endrow;i++)
		    for(int j=startcol;j<endcol;j++)
		      {
			Scalar val=matgen.getJ2(i,j);
			if(abs(val)<1e-10) continue;	  			
			int ii=i-startrow;
			int jj=j-startcol;
			nonzeros.push_back(Triplet<Scalar>(ii,jj,val));
		      }
		  partJ2.setFromTriplets(nonzeros.begin(),nonzeros.end());
		}
	      
		
		
	      colcom=world.split(col);
	      rowcom=world.split(row);
	    }
	}
    
    if(world.rank()<numOfparts)
      diagcom=world.split(numOfparts);//0-numOfparts-1 belong to diagcom.
    else
      diagcom=world.split(MPI_UNDEFINED);
  }


  template<typename MatGen>
    void LanczosIter<MatGen>::setupOper()
    {
      for(int col=0;col<numOfparts;col++)
	for(int row=col;row<numOfparts;row++)
	  {
	    int rank=col+((2*numOfparts+1-row+col)*(row-col))/2;
	    if(world.rank()==rank)
	      {
		int startrow=row*partdim;
		int endrow=(row+1)*partdim;
		if(row==numOfparts-1)
		  endrow=totaldim;
		int startcol=col*partdim;
		int endcol=(col+1)*partdim;
		if(col==numOfparts-1)
		  endcol=totaldim;

		/*
		//setup matrix		
		partOper.resize(endrow-startrow,endcol-startcol);
		if(world.rank()<numOfparts)
		{
		for(int i=startrow;i<endrow;i++)
		for(int j=i;j<endcol;j++)
		{
		partOper(i-startrow,j-startcol)=matgen.getOper(i,j);
		partOper(j-startcol,i-startrow)=partOper(i-startrow,j-startcol);
			  
		}
		}
		else
		{
		for(int i=startrow;i<endrow;i++)
		for(int j=startcol;j<endcol;j++)
		{
		partOper(i-startrow,j-startcol)=matgen.getOper(i,j);
		}
		}
		*/
		
		//setup matrix,using sparse matrix
		partOper.resize(endrow-startrow,endcol-startcol);
		if(world.rank()<numOfparts)
		  {
		    vector< Triplet<Scalar> > nonzeros;
		    for(int i=startrow;i<endrow;i++)
		      for(int j=i;j<endcol;j++)
			{
			  Scalar val=matgen.getOper(i,j);
			  if(abs(val)<1e-10) continue;	  			
			  int ii=i-startrow;
			  int jj=j-startcol;
			  nonzeros.push_back(Triplet<Scalar>(ii,jj,val));
			  if(ii!=jj)
			    nonzeros.push_back(Triplet<Scalar>(jj,ii,val));			  
			}
		    partOper.setFromTriplets(nonzeros.begin(),nonzeros.end());
		      
		  }
		else
		  {
		    vector< Triplet<Scalar> > nonzeros;		      
		    for(int i=startrow;i<endrow;i++)
		      for(int j=startcol;j<endcol;j++)
			{
			  Scalar val=matgen.getOper(i,j);
			  if(abs(val)<1e-10) continue;	  			
			  int ii=i-startrow;
			  int jj=j-startcol;
			  nonzeros.push_back(Triplet<Scalar>(ii,jj,val));
			}
		    partOper.setFromTriplets(nonzeros.begin(),nonzeros.end());
		  }
	      }
	  }
    }
  
  
  template<typename MatGen>
    void LanczosIter<MatGen>::reportMemory()
    {
      //cal total number of no-zero matrix elements
      if(world.rank()==0)
	{
	  numAnonZeros=0;
	  numJ2nonZeros=0;
	  reduce(world,static_cast<unsigned long int>( partA.nonZeros() ),numAnonZeros,std::plus<unsigned long int>(),0);
	  reduce(world,static_cast<unsigned long int>( partJ2.nonZeros() ),numJ2nonZeros,std::plus<unsigned long int>(),0);	  
	  std::cout<<"total number of non-zero A matrix elements: "<<numAnonZeros<<std::endl;
	  std::cout<<"total number of non-zero J2 matrix elements: "<<numJ2nonZeros<<std::endl;
	  std::cout<<"estimated memory: "<<(numAnonZeros+numJ2nonZeros)*sizeof(Scalar)/double(1<<30)<<"G\n";
	}
      else
	{
	  reduce(world,partA.nonZeros(),std::plus<unsigned long int>(),0);
	  reduce(world,partJ2.nonZeros(),std::plus<unsigned long int>(),0);
	}      
    }
  
  
  template<typename MatGen>
    void LanczosIter<MatGen>::clear()
    {
      int rank=world.rank();
      if(rank==0)
	{
	  alpha.clear();
	  beta.clear();
	}
      if(rank<numOfparts)
	{
	  partV.clear();
	}
      dim=0;
    }

  template<typename MatGen>
    void LanczosIter<MatGen>::AmultiplyV(const MatrixAType& _partA,const VectorType & _partV,VectorType & _partAV)
    {
      VectorType AV,Vtemp;
      int rank=world.rank();
      if(rank<numOfparts)
	{
	  Vtemp=_partV;
	}
      broadcast(colcom,Vtemp,0);
      if(rank<numOfparts)
	{
	  VectorType Temp=_partA*Vtemp;
	  reduce(rowcom,Temp,AV,std::plus<VectorType>(),0);
	  _partAV=AV;
	}
      else
	{
	  VectorType Temp=_partA*Vtemp;
	  reduce(rowcom,Temp,std::plus<VectorType>(),0);
	}

      broadcast(rowcom,Vtemp,0);
      if(rank<numOfparts)
	{
	  int partrows=partdim;
	  if(rank==numOfparts-1)
	    partrows=totaldim-rank*partdim;
	  VectorType Temp=VectorType::Zero(partrows);
	  reduce(colcom,Temp,AV,std::plus<VectorType>(),0);
	  _partAV+=AV;
	}
      else
	{
	  VectorType Temp=_partA.transpose()*Vtemp;
	  reduce(colcom,Temp,std::plus<VectorType>(),0);
	}      
    }


  
  template<typename MatGen>
    void LanczosIter<MatGen>::udotv(const VectorType & _partu,const VectorType & _partv,Scalar & _uv)
    {
      int rank=world.rank();
      if(rank < numOfparts)
	{
	  Scalar partuv=_partu.transpose()*_partv;
	  reduce(diagcom,partuv,_uv,std::plus<Scalar>(),0);
	}
    }

  template<typename MatGen>
    bool LanczosIter<MatGen>::next()
    {
      int rank=world.rank();
      VectorType partVtemp;
      if(dim==0)
	{
	  if(rank<numOfparts)
	    {
	      int start=rank*partdim;
	      int partrows=partdim;
	      if(rank==numOfparts-1)
		partrows=totaldim-rank*partdim;
	      partV.push_back(V0.segment(start,partrows));
	    }
	  AmultiplyV(partA,partV[partV.size()-1],partu);
	  if(rank<numOfparts)
	    {
	      Scalar a;
	      udotv(partu,partV[partV.size()-1],a);
	      broadcast(diagcom,a,0);
	      if(rank==0)
		{
		  alpha.push_back(a);
		}
	      partu-=a*partV[partV.size()-1];
	    }
	  dim++;
	}
      else
	{
	  bool succeed=true;
	  Scalar norm;	  
	  if(rank<numOfparts)
	    {
	      udotv(partu,partu,norm);
	      if(rank==0)
		{
		  norm=sqrt(norm);
		  if(abs(norm)<1e-8)
		    {
		      succeed=false;
		    }
		  else
		    {
		      beta.push_back( norm );
		    }
		}
	    }
	  broadcast(world,succeed,0);
	  if(succeed)
	    {
	      if(rank<numOfparts)
		{
		  broadcast(diagcom,norm,0);
		  partV.push_back( partu/norm );
		}
	    }
	  else
	    return false;
	  
	  AmultiplyV(partA,partV[partV.size()-1],partu);

	  if(rank<numOfparts)
	    {
	      Scalar a;
	      udotv(partu,partV[partV.size()-1],a);
	      //broadcast(diagcom,a,0);
	      //Scalar b;
	      if(rank==0)
		{
		  //b=beta[beta.size()-1];
		  alpha.push_back(a);
		}
	      //broadcast(diagcom,b,0);
	      //partu-=a*partV[partV.size()-1]+b*partV[partV.size()-2];

	      //reorthogonal in case of numeric error accumulating
	      for(int i=partV.size()-1;i>=0;i--)
		{
		  Scalar t;
		  udotv(partu,partV[i],t);
		  broadcast(diagcom,t,0);
		  partu-=t*partV[i];
		}
	    }

	  dim++;
	}
      return true;
    }



  template<typename MatGen>
    int LanczosIter<MatGen>::getAlphaBeta(VectorType & diag, VectorType & subdiag,int _dim)
    {
      if(dim<_dim)
	{
	  while(next())
	    {
	      if(dim>=_dim)
		break;
	    }
	}
      int rank=world.rank();
      if(rank==0)
	{
	  if(dim<_dim)
	    {
	      _dim=dim;
	      std::cout<<"warning: The lanczos iteration terminated, you may need to change the starting vec.\n";
	    }
	  diag.resize(_dim);
	  subdiag.resize(_dim);
	  for(int i=0;i<_dim;i++)
	    {
	      diag[i]=alpha[i];
	    }
	  for(int i=0;i<_dim-1;i++)
	    {
	      subdiag[i]=beta[i];
	    }
	  subdiag[_dim-1]=0;
	}
      return _dim;
    }

  template<typename MatGen>
    void LanczosIter<MatGen>::compute(VectorType & m_eivalues,int max_dim,bool yesvecs)
    {
      if(max_dim==0) return;
      if(max_dim>totaldim||max_dim==-1) max_dim=totaldim;

      int rank=world.rank();
      
      VectorType &diag=m_eivalues;
      VectorType subdiag;
      int n=getAlphaBeta(diag,subdiag,max_dim);
      
      if(rank==0)
	{
	  double scale=abs(diag[0]);
	  for(int i=0;i<n;i++)
	    {
	      double temp=abs(diag[i]);
	      if(scale<temp)
		scale=temp;
	      if(i<n-1)
		{
		  temp=abs(subdiag[i]);
		  if(scale<temp)
		    scale=temp;
		}
	    }
	  if(scale==double(0)) scale = double(1);
	  diag/=scale;
	  subdiag/=scale;

	  
	  Q=MatrixType::Identity(n,n); 
	  TriDiagQL(diag,subdiag,Q,n,yesvecs);
	  
	  if(yesvecs)
	    quickSort(m_eivalues,Q,0,n-1);
	  else
	    quickSort(m_eivalues,0,n-1);
	  // scale back the eigen values
	  m_eivalues *= scale;
	}
    }
  
  template<typename MatGen>
    void LanczosIter<MatGen>::getEigVecs(MatrixType & eigVecs,int numStates)
    {
      int rank=world.rank();
      int n=Q.rows();
      if(rank<numOfparts)
	{
	  broadcast(diagcom,n,0);
	  int partrows=partdim;
	  if(rank==numOfparts-1)
	    partrows=totaldim-rank*partdim;
	  eigVecs.resize(partrows,n);
	  for(int i=0;i<n;i++)
	    {
	      eigVecs.col(i)=partV[i];
	    }
	  broadcast(diagcom,Q,0);
      	  eigVecs=(eigVecs*Q).leftCols(numStates);
	}
    }
  
  template<typename MatGen>
    void LanczosIter<MatGen>::calJ2(VectorType & J2s,int numStates)
    {
      int rank=world.rank();
      int n=Q.rows();
      broadcast(world,n,0);
      if(rank==0)
	{
	  J2s.resize(n);
	  MEJ2.resize(n,n);
	}
      for(int i=0;i<n;i++)
	{
	  for(int j=i;j<n;j++)
	    {
	      VectorType Vtemp;
	      Scalar t;
	      AmultiplyV(partJ2,partV[i],Vtemp);
	      udotv(partV[j],Vtemp,t);
	      if(rank==0)
		{
		  MEJ2(i,j)=t;
		  MEJ2(j,i)=t;		  
		}
	    }
	}
      if(rank==0)
      	{
      	  J2s=(Q.transpose()*MEJ2*Q).diagonal().head(numStates);
      	}
    }

  template<typename MatGen>
    void LanczosIter<MatGen>::calOperMat(MatrixType & OperMat)
    {
      int rank=world.rank();
      int n=Q.rows();
      broadcast(world,n,0);
      if(rank==0)
	{
	  OperMat.resize(n,n);
	}
      for(int i=0;i<n;i++)
	{
	  for(int j=i;j<n;j++)
	    {
	      VectorType Vtemp;
	      Scalar t;
	      AmultiplyV(partOper,partV[i],Vtemp);
	      udotv(partV[j],Vtemp,t);
	      if(rank==0)
		{
		  OperMat(i,j)=t;
		  OperMat(j,i)=t;	  
		}
	    }
	}
    }

  template<typename MatGen>
    void LanczosIter<MatGen>::getVecs(MatrixType & eigVecs,int numStates)
    {
      if(world.rank()==0)
	eigVecs=Q.leftCols(numStates);
    }
}
#endif
