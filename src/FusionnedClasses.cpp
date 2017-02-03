/*
 *  FunsionnedClasses.cpp
 *  PseudoMatrix
 *
 *  Created by Michel Koskas on 02/11/15.
 *  Copyright 2015 INRA, INA. All rights reserved.
 *
 */

#include "FusionnedClasses.h"
#include "ClassesHeap.h"
#include <string>
#include <sstream>

#include <iostream>


/* (default is ok (no allocation is ever made).
FusionnedClasses::~FusionnedClasses() TrucTruc
{
  
}
*/
 
FusionnedClasses::FusionnedClasses()
{
  NbFusions = 0;
  FusionCost = -1;
  MyIndex = -1;
  MyAvailableIndex = -1;
  NextAvailableIndex = -1;
  PrevAvailableIndex = -1;
  MyMatrix = NULL;
  WhoIAm = 0;
}

FusionnedClasses::FusionnedClasses(PseudoMatrix *M, int ClassIndex)
{
  Initialize(M, ClassIndex);
}

void FusionnedClasses::Initialize(PseudoMatrix *M, int ClassIndex)
{
  MyMatrix = M;
  MyIndex = ClassIndex;
  MyAvailableIndex = ClassIndex;
  PrevAvailableIndex = ClassIndex - 1;
  NextAvailableIndex = ClassIndex + 1;
  if (ClassIndex == MyMatrix->p - 1)
    NextAvailableIndex = -1;
  WhoIAm = -(ClassIndex + 1);
}

void FusionnedClasses::InitializeFusionCost()
{
  if (MyIndex == MyMatrix->p - 1)
    NextAvailableIndex = -1;
  else
    NextAvailableIndex = MyIndex + 1;
  if (NextAvailableIndex < 0)
    FusionCost = -1;
  else
    {
      double Coefficient = 0.5;
      FusionCost = Coefficient * (MyMatrix->Value(MyIndex, MyIndex) +  MyMatrix->Value(MyIndex + 1, MyIndex + 1) - 2 * MyMatrix->Value(MyIndex, MyIndex + 1));
    }
}

double FusionnedClasses::MyValue() const
{
  assert(MyIndex == MyAvailableIndex);
  return MyMatrix->Value(MyAvailableIndex, MyAvailableIndex);
}

void FusionnedClasses::Swallow(int &NumFusionnedClass)
{
  assert(MyAvailableIndex == MyIndex);
  assert(NextAvailableIndex > -1);
  FusionnedClasses *Next = (MyMatrix->MyClasses) + NextAvailableIndex;
  NbFusions += 1 + Next->NbFusions;
  Next->MyAvailableIndex = MyIndex;
  NextAvailableIndex = Next->NextAvailableIndex;
  if (NextAvailableIndex > -1)
  {
    Next = (MyMatrix->MyClasses) + NextAvailableIndex;
    Next->PrevAvailableIndex = MyIndex;
  }
  WhoIAm = NumFusionnedClass + 1;
  NumFusionnedClass++;
}


void FusionnedClasses::ComputeMyFusionCost()
{
  assert(MyIndex == MyAvailableIndex);
  if (NextAvailableIndex == -1)
    {
      FusionCost = -1;
      return;
    }
  double a = MyCardinal();
  double b = MyMatrix->MyClasses[NextAvailableIndex].MyCardinal();
  double Coefficient = (a * b) / (a + b);
  
  FusionCost = Coefficient * (MyValue() / (a * a)  - 2 * MyMatrix->Value(MyIndex, NextAvailableIndex) / (a * b) + MyMatrix->Value(NextAvailableIndex, NextAvailableIndex) / (b * b));
}

bool FusionnedClasses::Exist()
{
  return MyIndex == MyAvailableIndex;
}

int FusionnedClasses::MyCardinal() const
{
  return NbFusions + 1;
  /*assert(MyIndex == MyAvailableIndex);
  if (NextAvailableIndex == -1)
    return MyMatrix->p - MyAvailableIndex;
  return NextAvailableIndex - MyAvailableIndex;
   */
}

std::ostream &operator<<(std::ostream &s, const FusionnedClasses &C)
{
    // s << "(" << &C << "-> ";
  s << "WhoIAm = " << C.WhoIAm;
  // s << C.WhoIAm;
  s << "MyCardinal = " << C.MyCardinal() << ", ";
  s << "MyIndex = " << C.MyIndex << ", ";
  s << "FusionCost = " << C.FusionCost << ", ";
  s << "MyAvailableIndex = " << C.MyAvailableIndex << ", ";
  s << "NextAvailableIndex = " << C.NextAvailableIndex << ", ";
  s << "MyValue = " << C.MyValue() << ", ";
  s << "NbFusions = " << C.NbFusions << ", ";
  s << "MyMatrix = " << C.MyMatrix << ".";
  s << ")  ";
  return s; 
}


bool operator<(const FusionnedClasses &Left, const FusionnedClasses &Right)
{
  if (Left.FusionCost == Right.FusionCost)
    return (Left.MyAvailableIndex < Right.MyAvailableIndex);
  return (Left.FusionCost < Right.FusionCost);
}

