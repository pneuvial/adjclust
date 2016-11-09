//
//  Permutation.cpp
//  MPseudoMatrix
//
//  Created by Michel Koskas on 22/02/2016.
//  Copyright © 2016 INRA. All rights reserved.
//

#include "Permutation.hpp"
#include "GeneralFunctions.h"
#include <assert.h>

/*

 class Permutation
 {
 public:
 int *P;
 int *IP;
 int *Free;
 bool *Used;
 int NbFree;
 Permutation();
 Permutation(int);
 void Initialize(int h);
 void Kill(int i);
 void Add(int i);
 void Displace(int i);
 };


*/

Permutation::Permutation()
{
  P = NULL;
  IP = NULL;
  Free = NULL;
  Used =NULL;
  NbFree = 0;
}

Permutation::Permutation(int hMax)
{
  Initialize(hMax);
}

void Permutation::Initialize(int hMax)
{
  P = new int[hMax + 1];
  IP = new int[hMax + 1];
  Free = new int[hMax + 1];
  for (int i = 0; i < hMax + 1; i++)
    P[i] = i;
  for (int i = 0; i < hMax + 1; i++)
    IP[i] = i;
  NbFree = 0;
  Used = new bool[hMax + 1];
  for (int i = 0; i < hMax + 1; i++)
    Used[i] = true;
}

void Permutation::Kill(int i)
{
  assert(Used[i]);
  if (!Used[i])
  {
    std::cerr << "I Can NOT kill a corpse..."  << std::endl << "Please reconsider! " << std::endl << "I'm history." << std::endl;
    exit(257);
  }
  Free[NbFree] = IP[i];
  NbFree++;
  Used[i] = false;
}

void Permutation::Displace(int i)
{
  if (!Used[i])
  {
    std::cerr << "I Can NOT displace a corpse..."  << std::endl << "Please reconsider! " << std::endl << "I'm history." << std::endl;
    exit(258);
  }
  if (NbFree == 0)
    return;
  if (Free[NbFree - 1] > IP[i])
    return;
  MySwap<int>(i, P[Free[NbFree - 1]]);
  MySwap<int>(IP[i], Free[NbFree - 1]);
}

void Permutation::Add(int i)
{
  if ((Used[i]) || (NbFree <= 0))
  {
    std::cerr << "I Can NOT give birth to a living body..."  << std::endl << "Please reconsider! " << std::endl << "I'm history." << std::endl;
    exit(259);
  }
  Used[i] = true;
  P[Free[NbFree - 1]] = i;
  IP[i] = Free[NbFree - 1];
  NbFree--;
}

void Permutation::Decremente(int i)
{
  assert(!Used[i - 1]);
  MySwap<int>(P[i], P[i - 1]);
  MySwap<int>(IP[P[i]], IP[P[i - 1]]);
  MySwap<bool>(Used[i], Used[i - 1]);
}

Permutation::~Permutation()
{
  if (P != NULL)
  {
    delete[] P;
    delete [] IP;
    delete[] Free;
    delete[] Used;
  }
}




