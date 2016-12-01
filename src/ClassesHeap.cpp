/*
 *  ClassesHeap.cpp
 *  PseudoMatrix
 *
 *  Created by Michel Koskas on 06/11/15.
 *  Copyright 2015 INRA, INA. All rights reserved.
 *
 */

#include "ClassesHeap.h"
#include <iostream>
#include <fstream>
#include <Rcpp.h>
using namespace Rcpp;

ClassesHeap::ClassesHeap()
{
  HeapSize = 0;
  MaxSize = 0;
  MyMatrix = NULL;
  Which = NULL;
  NumFusionnedClasses = 0;
  OutputStream = NULL;
}

ClassesHeap::ClassesHeap(PseudoMatrix *M, int NbWC, std::string OutputPath)
{
  Initialize(M, NbWC, OutputPath);
}

void ClassesHeap::Initialize(PseudoMatrix *M, int NbWC, std::string OutputPath)
{

  assert(NbWC <= M->p);
  const clock_t Begin = clock();

  if (OutputPath != "")
  {
    std::ofstream *f = new std::ofstream;
    f->open(OutputPath.c_str());
    OutputStream = f;
  }
  else
    //    OutputStream = &(std::cout);
    OutputStream = &(Rcpp::Rcout);

  Output = new double[3 * (M->p - NbWC)];

  NumFusionnedClasses = 0;
  HeapSize = 0;
  MyMatrix = M;
  MaxSize = MyMatrix->p - 1;
  assert(NbWC > 0);
  NbWantedClasses = NbWC;
  CurrentNbClasses = 1;
  Heap = new int[MaxSize];
  Which = new int[MaxSize + 1];

  const clock_t EndInit = clock();

  //std::cout << "Heap Init time = " << ((double) (EndInit - Begin)) / CLOCKS_PER_SEC << std::endl;
  Rcpp::Rcout << "Heap Init time = " << ((double) (EndInit - Begin)) / CLOCKS_PER_SEC << std::endl;
  for (int i = 0; i < MaxSize; i++)
  {
    AddNode();
    // assert(CheckMe());
    // std::cout << HeapSize << " ";
  }
  Which[MaxSize] = MaxSize;


  while ((CurrentNbClasses > NbWantedClasses) && (HeapSize > 0))
  {
    // std::cout << "------------------------ >>>>>>   Before fusion:" << std::endl;
    // std::cout << *this;
    MakeAFusion();
    // std::cout << "------------------------ >>>>>>   After fusion:" << std::endl;
    // std::cout << *this;
		// std::cout << "NbClasses = " << CurrentNbClasses << std::endl;

    // assert(CheckMe());
    // std::cout << HeapSize << " ";
  }
/*  if ((HeapSize == 1) && (CurrentNbClasses > NbWantedClasses))
  {
    std::cerr << "Pretty strange" << std::endl;
    std::cerr << "HeapSize = 1" << std::endl;
    std::cerr << "NbWantedClasses = " << NbWantedClasses << std::endl;
    std::cerr << "CurrentNbClasses = " << CurrentNbClasses << std::endl;
  }*/

  const clock_t End = clock();
  // std::cout << "Time computation for the heap = " << ((double) (End - EndInit)) / CLOCKS_PER_SEC << std::endl;
  Rcpp::Rcout << "Time computation for the heap = " << ((double) (End - EndInit)) / CLOCKS_PER_SEC << std::endl;
}

ClassesHeap::~ClassesHeap()
{
  delete[] Heap;
  delete[] Which;
  delete MyMatrix;
}

bool ClassesHeap::CheckMe()
{
  for (int i = 0; i < HeapSize; i++)
  {
    FusionnedClasses *C = MyMatrix->MyClasses + Heap[i];
    if (C->MyIndex != C->MyAvailableIndex)
      return false;
    if (Which[Heap[i]] != i)
      return false;
    if (Heap[Which[i]] != i)
      return false;
    if (C->NextAvailableIndex == -1)
      return false;
    if (2 * i + 1 < HeapSize)
      if (MyMatrix->MyClasses[Heap[2 * i + 1]] < MyMatrix->MyClasses[Heap[i]])
        return false;
    if (2 * i + 2 < HeapSize)
      if (MyMatrix->MyClasses[Heap[2 * i + 2]] < MyMatrix->MyClasses[Heap[i]])
        return false;
  }
  return true;
}


ClassesHeap::ClassesHeap(double *V, int P, int H, int NbWC, std::string OutputPath)
{
  PseudoMatrix *M = new PseudoMatrix;
  M->Initialize(V, P, H);
  Initialize(M, NbWC, OutputPath);
}

void ClassesHeap::Swap(int i, int j)
{
  MySwap<int>(Heap[i], Heap[j]);
  MySwap<int>(Which[Heap[i]], Which[Heap[j]]);

}

bool ClassesHeap::RebalanceToDown(int Index)
{
  int CurrentIndex = Index;
  bool Res = false;
  while (2 * CurrentIndex + 2 < HeapSize)
  {
    FusionnedClasses *Left = MyMatrix->MyClasses + Heap[2 * CurrentIndex + 1], *Right = MyMatrix->MyClasses + Heap[2 * CurrentIndex + 2];
    int IndexMinSon = 2 * CurrentIndex + 1;
    if (*Right < *Left)
      IndexMinSon++;
    if (MyMatrix->MyClasses[Heap[IndexMinSon]] < MyMatrix->MyClasses[Heap[CurrentIndex]])
    {
      Swap(CurrentIndex, IndexMinSon);
      CurrentIndex = IndexMinSon;
      Res = true;
    }
    else
      return Res;
  }
  if (2 * CurrentIndex + 1 < HeapSize)
    if (MyMatrix->MyClasses[Heap[2 * CurrentIndex + 1]] < MyMatrix->MyClasses[Heap[CurrentIndex]])
    {
      Swap(CurrentIndex, 2 * CurrentIndex + 1);
      Res = true;
    }
  return Res;
}

bool ClassesHeap::RebalanceToUp(int IndexInHeap)
{
  bool Res = false;

  int CurrentIndex = IndexInHeap;
  if (CurrentIndex == -1)
    CurrentIndex = HeapSize - 1;
  while ((CurrentIndex > 0) && (MyMatrix->MyClasses[Heap[CurrentIndex]] < MyMatrix->MyClasses[Heap[(CurrentIndex - 1) / 2]]))
  {
    Res = true;
    Swap(CurrentIndex, (CurrentIndex - 1) / 2);
    CurrentIndex = (CurrentIndex - 1) / 2;
  }
  return Res;
}


void ClassesHeap::AddNode()
{
  Heap[HeapSize] = HeapSize;
  Which[HeapSize] = HeapSize;
  HeapSize++;
  RebalanceToUp();
  CurrentNbClasses++;
}

void ClassesHeap::MakeAFusion()
{
  if (HeapSize <= 0)
    return;
  int IndexNextClass = MyMatrix->MyClasses[Heap[0]].NextAvailableIndex;
	assert(IndexNextClass > 0);
  int NextToNextClass = MyMatrix->MyClasses[IndexNextClass].NextAvailableIndex;
  int IndexCurrentClass = Heap[0];
  int IndexPrecClass;
  if (Heap[0] > 0)
    IndexPrecClass = MyMatrix->MyClasses[Heap[0] - 1].MyAvailableIndex;
  else
    IndexPrecClass = -1;
  MyMatrix->Fusion(IndexCurrentClass, NumFusionnedClasses, this);
	MyMatrix->MyClasses[IndexCurrentClass].ComputeMyFusionCost();
	if (IndexPrecClass > -1)
	{
		MyMatrix->MyClasses[IndexPrecClass].ComputeMyFusionCost();
		FullRebalance(Which[IndexPrecClass]);
	}
  if (NextToNextClass == -1)
  {
    Swap(0, HeapSize - 1);
    HeapSize--;
    RebalanceToDown();
  }
  else
  {
    int ModifiedIndex = Which[IndexNextClass];
    Swap(ModifiedIndex, HeapSize - 1);
    HeapSize--;
    FullRebalance(ModifiedIndex);
  }
	CurrentNbClasses--;
}

void ClassesHeap::FullRebalance(int IndexInHeap)
{
  if (IndexInHeap >= HeapSize - 1)
    return;
  if (IndexInHeap == 0)
  {
    RebalanceToDown();
    return;
  }
  if (!RebalanceToUp(IndexInHeap))
    RebalanceToDown(IndexInHeap);
}

std::ostream &operator<<(std::ostream &s, const ClassesHeap &H)
{
  s << "HeapSize = " << H.HeapSize << std::endl;
  s << "MaxSize = " << H.MaxSize << std::endl;
  s << "CurrentNbClasses = " << H.CurrentNbClasses << std::endl;
  s << "NbWantedClasses = " << H.NbWantedClasses << std::endl;
  s << "Heap = (";
  for (int i = 0; i < H.HeapSize - 1; i++)
    s << H.Heap[i] << ", ";
  if (H.HeapSize  > 0)
    s << H.Heap[H.HeapSize - 1]<< ")" << std::endl;
  s << "Which = (";
  for (int i = 0; i < H.MyMatrix->p - 1; i++)
    s << H.Which[i] << ", ";
  s << H.Which[H.MyMatrix->p - 1]<< ")" << std::endl;
  s << "My Pseudo Matrice is: " << *(H.MyMatrix);

  for (int i = 0; i < H.MyMatrix->p - H.NbWantedClasses; i++)
    s << H.Output[3 * i] << '\t' << H.Output[3 * i + 1] << '\t' << H.Output[3 * i + 2] <<std::endl;
  return s;
}


