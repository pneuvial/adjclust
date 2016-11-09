//
//  GeneralFunctions.cpp
//  MPseudoMatrix
//
//  Created by Michel Koskas on 19/01/16.
//  Copyright © 2016 INRA. All rights reserved.
//

#include <stdio.h>

#include "GeneralFunctions.h"


void MyPrint(double &f, std::ostream &s)
{
  s << f << '\t';
  int n = f;
  if (n == f)
    s << '\t' << '\t';
  int m = 100 * f;
  if ((m > 0) && (n != f))
    s << '\t';
}
