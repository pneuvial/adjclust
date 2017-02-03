//
//  ColsArray.cpp
//  MPseudoMatrix
//
//  Created by Michel Koskas on 10/03/2016.
//  Copyright © 2016 INRA. All rights reserved.
//

#include "ColsArray.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include "GeneralFunctions.h"
#include "PseudoMatrix.h"
#include <Rcpp.h>
using namespace Rcpp;

/*


 class ColsArray
 {
 public:
 double *Data;
 int *ColIndex;
 int ArraySize;
 PseudoMatrix *MyMatrix;
 ColsArray();
 ColsArray(PseudoMatrix *M);
 void Set(int Col, double Value);
 double Value(int Col);
 void Restructure();
 private:
 };

*/


void ColsArray::Restructure()
{
  int WriteIndex = 0;
  for (int ReadIndex = 0; ReadIndex < ArraySize; ReadIndex++)
  {
    if (MyMatrix->MyClasses[ColIndex[ReadIndex]].Exist())
    {
      ColIndex[WriteIndex] = ColIndex[ReadIndex];
      Data[WriteIndex] = Data[ReadIndex];
      WriteIndex++;
    }
  }
  ArraySize = WriteIndex;
  if (ArraySize > MyMatrix->h + 1)
    Rcpp::Rcerr << "Strange restructuration..." << std::endl;
}

bool ColsArray::CheckMe()
{
  if (ArraySize == 0)
    return true;
  for (int i = 0; i < ArraySize - 1; i++)
  {

    int First = ColIndex[i];
    int Last = ColIndex[i + 1];
    for (int i = First + 1; i < Last; i++)
      if (MyMatrix->MyClasses[i].Exist())
        return false;
    First = Last;
  }
  return true;
}



ColsArray::ColsArray()
{
  Data = NULL;
  ColIndex = NULL;
  MyMatrix = NULL;
  ArraySize = 0;
}

ColsArray::ColsArray(PseudoMatrix *M)
{
  Initialize(M);
}

void ColsArray::Initialize(PseudoMatrix *M)
{
  Data = new double[M->h + 1];
  ColIndex = new int[M->h + 1];
  MyMatrix = M;
}

void ColsArray::Set(int Col, double Value, int ForceIndex, bool Assumption)
{
  if ((Assumption) && (ColIndex[ForceIndex] == Col))
    Assumption = false;
  if (!Assumption)
  {
    ColIndex[ForceIndex] = Col;
    Data[ForceIndex] = Value;
    if (ArraySize < ForceIndex)
      ArraySize = ForceIndex;
    return;
  }
  Set(Col, Value);
  return;
}



void ColsArray::Set(int Col, double Value)
{
  if (Col >= MyMatrix->p)
    return;
  if (Col < 0)
    return;
  if (ArraySize >= MyMatrix->h)
    Restructure();
  assert(MyMatrix->MyClasses[Col].Exist());
  assert(ArraySize <= MyMatrix-> h + 1);
  int Where;
  bool IsPresent = Find<int>(ColIndex, ArraySize, Col, Where);
  if (IsPresent)
  {
    Data[Where] = Value;
    return;
  }
  if (Where == ArraySize)
  {
    ColIndex[ArraySize] = Col;
    Data[ArraySize] = Value;
    ArraySize++;
    return;
  }
  // Odd case (hopefully unlike)
  for (int i = ArraySize; i > Where; i--)
  {
    ColIndex[i] = ColIndex[i - 1];
    Data[i] = ColIndex[i - 1];
  }
  ColIndex[Where] = Col;
  Data[Where] = Value;
  ArraySize++;
  /*Rcpp::cerr << "<<<<<-------------------Value not found ------------->>>>>" << std::endl;
  Rcpp::cerr << *this << std::endl;;
  Rcpp::cerr << "Wanted Column: " << Col << std::endl;
  exit(17);*/
}

double ColsArray::Value(int Col, int SupposedPlace)
{
  if (ColIndex[SupposedPlace] == Col)
    return Data[SupposedPlace];
  if (ColIndex[SupposedPlace] < Col)
    for (int i = SupposedPlace + 1; (i < ArraySize); i++)
    {
      if (ColIndex[i] == Col)
        return Data[i];
      if (ColIndex[i] > Col)
        return 0;
    }
  return Value(Col);
}

double ColsArray::Value(int Col)
{
  int Where;
  bool IsPresent = Find<int>(ColIndex, ArraySize, Col, Where);
  if (IsPresent)
    return Data[Where];
  else
    return 0;
}

std::ostream & operator<<(std::ostream &s, const ColsArray &C)
{
  s << "ArraySize = " << C.ArraySize << std::endl;
  for (int i = 0; i < C.ArraySize; i++)
  {
    bool Present = C.MyMatrix->MyClasses[C.ColIndex[i]].Exist();
    if (!Present)
      s << "(";
    s << C.ColIndex[i];
    if (!Present)
      s << ")";
    s << '\t';
  }
  s << std::endl;
  for (int i = 0; i < C.ArraySize; i++)
  {
    bool Present = C.MyMatrix->MyClasses[C.ColIndex[i]].Exist();
    if (!Present)
      s << "(";
    s << C.Data[i];
    if (!Present)
      s << ")";
    s << '\t';
  }
  s << std::endl;
  return s;
}





