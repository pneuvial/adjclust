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




















#endif
