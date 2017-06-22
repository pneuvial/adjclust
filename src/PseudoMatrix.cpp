/*
 *  PseudoMatrix.cpp
 *  PseudoMatrix
 *
 *  Created by Michel Koskas on 02/11/15.
 *  Copyright 2015 INRA, INA. All rights reserved.
 *
 */

#include "PseudoMatrix.h"
#include <iostream>
#include <iomanip>
#include <utility>
#include "GeneralFunctions.h"
#include <math.h>
#include <Rcpp.h>
using namespace Rcpp;

PseudoMatrix::PseudoMatrix()
{
  h = 0;
  p = 0;
  MyClasses = NULL;
  NbClasses = 0;
  MyData = NULL;
}

PseudoMatrix::PseudoMatrix(double *V, int P, int H)
{
  Initialize(V, P, H);
}

void PseudoMatrix::Initialize(double *V, int P, int H)
{
  clock_t Begin = clock();
  p = P;
  h = H;
  MyClasses = new FusionnedClasses[p];
  for (int i = 0; i < p; i++)
    MyClasses[i].Initialize(this, i);
  MyData = new ColsArray[p];
  for (int Line = 0; Line < p; Line++)
  {
    MyData[Line].Initialize(this);
    for (int Col = 0; Col < h + 1; Col++)
      MyData[Line].Set(Line + Col, V[Line * (h + 1) + Col], Col);
  }
  for (int i = 0; i < p; i++)
    MyClasses[i].InitializeFusionCost();

  clock_t End = clock();
//  Rcpp::Rcout << "Time needed to Initialize the pseudo matrix = " << ((double) End - Begin) / CLOCKS_PER_SEC << std::endl;
}

double PseudoMatrix::Value(int LineInA, int ColumnInA, int SupposedPlace) const
{
  return MyData[LineInA].Value(ColumnInA, SupposedPlace);
}

double PseudoMatrix::Value(int LineInA, int ColumnInA) const
{
  return MyData[LineInA].Value(ColumnInA);
}

void PseudoMatrix::Set(int LineInA, int ColumnInA, double v)
{
  MyData[LineInA].Set(ColumnInA, v);
}

void PseudoMatrix::Set(int LineInA, int ColumnInA, double v, int ForceIndex, bool Assumption)
{
  MyData[LineInA].Set(ColumnInA, v, ForceIndex, Assumption);
}


PseudoMatrix::~PseudoMatrix()
{
  delete[] MyClasses;
  delete[] MyData;
}


void PseudoMatrix::Fusion(int LineIndex, int &NumFusionnedClass)
{
  assert(MyClasses[LineIndex].MyAvailableIndex == LineIndex);
  int NextAvailableIndex = MyClasses[LineIndex].NextAvailableIndex;
  assert(NextAvailableIndex > -1);

  MyData[LineIndex].Restructure();

  // Taking care of the diagonal value
  double v = Value(LineIndex, LineIndex, 0) + Value(NextAvailableIndex, NextAvailableIndex, 0) + 2 * Value(LineIndex, NextAvailableIndex, 1);
  Set(LineIndex, LineIndex, v, 0);

  // DisplayMatrixA(Rcpp::Rcout);

  // Taking care of the line:
  int CurrentIndex = MyClasses[NextAvailableIndex].NextAvailableIndex;
  for (int d = 2; ((d<= h + 1) && (CurrentIndex > -1)); d++)
  {
		v = Value(LineIndex, CurrentIndex, d) + Value(NextAvailableIndex, CurrentIndex, d - 1);
		Set(LineIndex, CurrentIndex, v, d - 1);
		CurrentIndex = MyClasses[CurrentIndex].NextAvailableIndex;
    MyData[LineIndex].ArraySize = d;
		// DisplayMatrixA(Rcpp::Rcout);
  }

  // Now taking care of the column:
  CurrentIndex = MyClasses[LineIndex].PrevAvailableIndex;
  for (int i = 1; ((i <= h) && (CurrentIndex > -1)); i++)
  {
    v = Value(CurrentIndex, LineIndex, i) + Value(CurrentIndex, NextAvailableIndex, i + 1);
    Set(CurrentIndex, LineIndex, v, i, true);
    CurrentIndex = MyClasses[CurrentIndex].PrevAvailableIndex;
  }
  MyClasses[LineIndex].Swallow(NumFusionnedClass);
  delete[] MyData[NextAvailableIndex].ColIndex;
  delete[] MyData[NextAvailableIndex].Data;
}

void PseudoMatrix::DisplayMatrixA(std::ostream &s) const
{
  long int CurrentPrecision = s.precision();

  for (int i = 0; i < p; i++)
  {
		if (MyClasses[i].MyAvailableIndex == MyClasses[i].MyIndex)
		{
			for (int j = 0; j < p; j++)
				if (MyClasses[j].MyAvailableIndex == MyClasses[j].MyIndex)
				{
					double V;
					if (j >= i)
						V = Value(i, j);
					else
						V =Value(j, i);
					MyPrint(V, s);
				}
      s << std::endl << std::endl;
		}
  }
  s << std::endl << std::endl;
  s.precision(CurrentPrecision);
}


std::ostream &operator<<(std::ostream &s, const PseudoMatrix &M)
{
  s << "Pseudo Matrix " << &M << "-> [";
  s << "h = " << M.h << ", ";
  s << "p = " << M.p << ", ";
	s << "Displaying WhichColumn: " << std::endl;
    // s << "MyValues = (" << std::endl;
  int p = M.p;
    // int h = M.h;
  long int CurrentPrecision = s.precision();
  s.precision(4);
  /*int CurrentIndex = 0;
  while (CurrentIndex > -1)
  {

    for (int j = 0; j < h; j++)
      s << std::setw(5) << M.MyValues[CurrentIndex * h + j] << " ";
    s << std::endl;
    CurrentIndex = M.MyClasses[CurrentIndex].NextAvailableIndex;
  }
  s << ")" << std::endl;*/
  s << "MyClasses : " << std::endl;
  for (int i = 0; i < p; i++)
    if (M.MyClasses[i].Exist())
      s << M.MyClasses[i] << std::endl;
  s << std::endl;
  s << "Displaying ColsArrays. "  << std::endl;
  for (int i = 0; i < p; i++)
    if (M.MyClasses[i].Exist())
      s << M.MyData[i] << std::endl;
  s << "Now displaying Matrix A:" << std::endl;
  M.DisplayMatrixA(s);



  s.precision(CurrentPrecision);
  return s;
}


