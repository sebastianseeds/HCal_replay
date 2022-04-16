#ifndef BLIND_H
#define BLIND_H
#include "TRandom.h"
#include "TMath.h"

double blind( double yield )
{
  double RNG = gRandom->Uniform( 0.9, 1.1 );
  double blindedYield = yield*RNG;
  return blindedYield;
}

#endif // BLIND_H
