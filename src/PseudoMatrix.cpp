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
#include <utility>

PseudoMatrix::PseudoMatrix()
{
  h = 0;
  p = 0;
  MyClasses = NULL;
  NbClasses = 0;
  NbFusions = 0;
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
  for (int Line = 0; Line < p; Line++)
    for (int Col = 0; Col <= h; Col++)
      Set(Line, Line + Col, V[Line * (h + 1) + Col]);
  MyClasses = new FusionnedClasses[p];
  for (int i = 0; i < p; i++)
    MyClasses[i].Initialize(this, i);
  clock_t End = clock();
  std::cout << "Time needed to Initialize the pseudo matrix = " << ((double) End - Begin) / CLOCKS_PER_SEC << std::endl; 
}

double PseudoMatrix::Value(int LineInA, int ColumnInA) const
{
  if (LineInA >= p)
    return 0;
  if (ColumnInA >= p)
    return 0;
  if (MyClasses[LineInA].MyIndex != MyClasses[LineInA].MyAvailableIndex)
    return 0;
  std::map<std::pair<int, int>, double>::const_iterator El = MyValues.find(std::make_pair(LineInA, ColumnInA));
  if (El == MyValues.end())
    return 0;
  else
    return El->second;
}

void PseudoMatrix::Set(int LineInA, int ColumnInA, double v)
{
  if (LineInA >= p)
    return;
  if (LineInA < 0)
    return;
  if (ColumnInA >= p)
    return;
  if (ColumnInA < 0)
    return;
  MyValues[std::make_pair(LineInA, ColumnInA)] = v;
}

void PseudoMatrix::Erase(int LineIndex, int ColumnIndex)
{
  std::map<std::pair<int, int>, double>::iterator El = MyValues.find(std::make_pair(LineIndex, ColumnIndex));
  if (El != MyValues.end())
  {
    MyValues.erase(El);
    return;
  }
  // std::cerr << "Warning: you're trying to erase an inexisting value" << std::endl;
}


PseudoMatrix::~PseudoMatrix()
{
  delete[] MyClasses;
}


void PseudoMatrix::Fusion(int LineIndex, int &NumFusionnedClass, ClassesHeap *H)
{
  assert(MyClasses[LineIndex].MyAvailableIndex == LineIndex);
  int NextAvailableIndex = MyClasses[LineIndex].NextAvailableIndex;
  assert(NextAvailableIndex > -1);

  // Taking care of the diagonal value
  double v = Value(LineIndex, LineIndex) + Value(NextAvailableIndex, NextAvailableIndex) + 2 * Value(LineIndex, NextAvailableIndex);
  Set(LineIndex, LineIndex, v);
	Erase(LineIndex, NextAvailableIndex);
	Erase(NextAvailableIndex, NextAvailableIndex);
	
  // DisplayMatrixA(std::cout);
  
  // Taking care of the line:
  int CurrentIndex = MyClasses[NextAvailableIndex].NextAvailableIndex;
  for (int d = 2; ((d<= h + 1) && (CurrentIndex > -1)); d++)
  {
		v = Value(LineIndex, CurrentIndex) + Value(NextAvailableIndex, CurrentIndex);
		Erase(NextAvailableIndex, CurrentIndex);
		Set(LineIndex, CurrentIndex, v);
		CurrentIndex = MyClasses[CurrentIndex].NextAvailableIndex;
		// DisplayMatrixA(std::cout);
  }
	
  // Now taking care of the column:
  if (LineIndex > 0)
    CurrentIndex = MyClasses[LineIndex - 1].MyAvailableIndex;
  else
    CurrentIndex = -1;
  for (int i = 1; ((i <= h) && (CurrentIndex > -1)); i++)
  {
    v = Value(CurrentIndex, LineIndex) + Value(CurrentIndex, NextAvailableIndex);
    Set(CurrentIndex, LineIndex, v);
    if (i < h)
      Erase(CurrentIndex, NextAvailableIndex);
    if (CurrentIndex > 0)
      CurrentIndex = MyClasses[CurrentIndex - 1].MyAvailableIndex;
    else
      CurrentIndex = -1;
  }
  NbFusions += 1 + MyClasses[NextAvailableIndex].NbFusions;
  MyClasses[LineIndex].Swallow(NumFusionnedClass, H);
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
		}
		s << std::endl << std::endl;
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
  for (std::map<std::pair<int, int>, double>::const_iterator I = M.MyValues.begin(); I != M.MyValues.end(); I++)
  {
    s << "V[" << I->first.first << ", " << I->first.second << "] = " << I->second << std::endl;
  }
  
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
  s << "Now displaying Matrix A:" << std::endl;
  M.DisplayMatrixA(s);
  
  
  
  s.precision(CurrentPrecision);
  return s;
}


