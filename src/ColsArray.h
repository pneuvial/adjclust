//
//  ColsArray.hpp
//  MPseudoMatrix
//
//  Created by Michel Koskas on 10/03/2016.
//  Copyright © 2016 INRA. All rights reserved.
//

#ifndef ColsArray_hpp
#define ColsArray_hpp

#include <stdio.h>
#include <iostream>
#include <iomanip>


class PseudoMatrix;

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
  void Set(int Col, double Value, int ForceIndex, bool Assumption = false);
  double Value(int Col);
  double Value(int Col, int SupposedPlace);
  void Restructure();
  void Initialize(PseudoMatrix *M);
private:
  bool CheckMe();
};

std::ostream & operator<<(std::ostream &s, const ColsArray &C);


#endif /* ColsArray_hpp */
