/*
 *  FunsionnedClasses.h
 *  PseudoMatrix
 *
 *  Created by Michel Koskas on 02/11/15.
 *  Copyright 2015 INRA, INA. All rights reserved.
 *
 */

#ifndef FusionnedClasses_h_
#define FusionnedClasses_h_

#include <string>
#include <assert.h>
#include "PseudoMatrix.h"
#include <iostream>


class PseudoMatrix;
class ClassesHeap;

class FusionnedClasses
{
public: 
  int NbFusions;
  double MyValue() const; // somme des elements diagonaux de la classe
  double FusionCost;
  int MyIndex;  // indice (absolu) du premier element de la classe: ne change jamais
  PseudoMatrix *MyMatrix;
  
    // MyAvailableIndex is the smallest index of the class I belong to
  int MyAvailableIndex; // indice actualise de la classe
    // NextAvailableIndex is the smallest index of the next class (-1 if I bezlon to the last class) 
  int NextAvailableIndex;
  int PrevAvailableIndex;
  int WhoIAm;
  FusionnedClasses(PseudoMatrix *M, int ClassIndex);
  FusionnedClasses();
  void Initialize(PseudoMatrix *M, int ClassIndex);
    // A FusionnedClass should only swallow next one
  void Swallow(int &NumFusionnedClass); // avale la voisine de droite
  void Disappear();
  bool Exist();
  void DisplayMatrixA();
    // ~FusionnedClasses();
  int MyCardinal() const;
  void ComputeMyFusionCost();
  void InitializeFusionCost();
  
private:
};

std::ostream &operator<<(std::ostream &s, const FusionnedClasses &C);
bool operator<(const FusionnedClasses &Left, const FusionnedClasses &Right);








#endif
