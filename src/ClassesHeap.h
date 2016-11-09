/*
 *  ClassesHeap.h
 *  PseudoMatrix
 *
 *  Created by Michel Koskas on 06/11/15.
 *  Copyright 2015 INRA, INA. All rights reserved.
 *
 */

#ifndef ClassesHeap_h_
#define ClassesHeap_h_

#include "PseudoMatrix.h"
#include "FusionnedClasses.h"
#include "GeneralFunctions.h"
#include <iostream>

class ClassesHeap
{
public:
  int HeapSize;
  int MaxSize;
  int CurrentNbClasses;
  int NbWantedClasses;
  int NumFusionnedClasses;
  double *Output;
  std::ostream *OutputStream;
  int *Heap;
    // Which is an array which associates an absolute index to the index int he heap. 
  int *Which;
  PseudoMatrix *MyMatrix;
  ClassesHeap();
  ClassesHeap(PseudoMatrix *M, int NbWC, std::string OutputPath = "");
  void Initialize(PseudoMatrix *M, int NbWC, std::string OutputPath);
  ClassesHeap(double *V, int p, int h, int NbWC, std::string OutputPath="");
  ~ClassesHeap();
  void AddNode();
  bool CheckMe();
  void MakeAFusion();
private:
  void Swap(int i, int j);
  bool RebalanceToDown(int Index = 0);
  bool RebalanceToUp(int Index = -1);
  void FullRebalance(int);
  bool RebalanceUp(int Index = -1);
  bool RebalanceDown(int Index = 0);
};

std::ostream &operator<<(std::ostream &s, const ClassesHeap &H);










#endif



