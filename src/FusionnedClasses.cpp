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
FusionnedClasses::~FusionnedClasses()
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
  if (MyIndex == MyMatrix->p - 1)
    NextAvailableIndex = -1;
  else
    NextAvailableIndex = ClassIndex + 1;
  WhoIAm = -(ClassIndex + 1);
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


void FusionnedClasses::Swallow(int &NumFusionnedClass, ClassesHeap *H)
{
  assert(MyAvailableIndex == MyIndex);
  assert(NextAvailableIndex > -1);
    // int NextClassIndex = NextAvailableIndex;
  FusionnedClasses *Next = (MyMatrix->MyClasses) + NextAvailableIndex;
  NbFusions += 1 + Next->NbFusions;
  for (int i = MyIndex; i < Next->MyIndex; i++)
    MyMatrix->MyClasses[i].NextAvailableIndex = Next->NextAvailableIndex;
  int LastIndex = MyMatrix->p;
  if (Next->NextAvailableIndex > -1)
    LastIndex = Next->NextAvailableIndex;
  for (int i = Next->MyIndex; i < LastIndex; i++)
    (*MyMatrix).MyClasses[i].MyAvailableIndex = MyAvailableIndex;

  H->Output[3 * NumFusionnedClass] = WhoIAm;
  H->Output[3 * NumFusionnedClass + 1] = Next->WhoIAm;
  H->Output[3 * NumFusionnedClass + 2] = FusionCost;
    // std::cout << "( " << WhoIAm << '\t' << Next->WhoIAm << ")" << std::endl;


  WhoIAm = NumFusionnedClass + 1;
  NumFusionnedClass++;
  if (NextAvailableIndex > -1)
    {
      ComputeMyFusionCost();
    }
  else
    FusionCost = -1;

  if (MyAvailableIndex > 0)
  {
    FusionnedClasses *Prec = MyMatrix->MyClasses + MyMatrix->MyClasses[MyAvailableIndex - 1].MyAvailableIndex;
    Prec->ComputeMyFusionCost();
  }
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

  // BEWARE: It's really MyIndex + 1 because coefficients have been stored even on unavailable columns!

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













