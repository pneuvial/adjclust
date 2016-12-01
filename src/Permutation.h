//
//  Permutation.hpp
//  MPseudoMatrix
//
//  Created by Michel Koskas on 22/02/2016.
//  Copyright © 2016 INRA. All rights reserved.
//

#ifndef Permutation_hpp
#define Permutation_hpp

#include <stdio.h>
#include <iostream>


// Free contains the list of free storage places.
// P says the order in which the neibourghs are stored
// IP is the inverse permutation of P, that is to say that IP[i] says where the i-th neibourgh is stored
// USed says if the i-th neibourgh is alive.
// At each moment P is a permutation of [0, h]. 
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
  void Decremente(int i);
  ~Permutation();
};







#endif /* Permutation_hpp */






