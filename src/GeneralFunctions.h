/*
 *  GeneralFunctions.h
 *  PseudoMatrix
 *
 *  Created by Michel Koskas on 06/11/15.
 *  Copyright 2015 INRA, INA. All rights reserved.
 *
 */


#ifndef GeneralFunctions_h_
#define GeneralFunctions_h_


#include <iostream>



void MyPrint(double &f, std::ostream &s);

template<typename T>
void MySwap(T &A, T &B)
{
  T C = A;
  A = B;
  B = C;
}

template <typename T>
bool Find(T *Array, int ArraySize, T x, int &Res)
{
  int First =0, Last = ArraySize;
  if (ArraySize == 0)
  {
    Res = 0;
    return false;
  }
  if (x < Array[0])
  {
    Res = 01;
    return false;
  }
  if (x > Array[ArraySize - 1])
  {
    Res = ArraySize;
    return false;
  }
  while (Last - First > 1)
  {
    int m = (First + Last) / 2;
    if (x < Array[m])
      Last = m;
    else
      First = m;
  }
  Res = First;
  if (x == Array[First])
    {Res = First; return true;}
  else
    {Res = First + 1; return false;}
}

















#endif
