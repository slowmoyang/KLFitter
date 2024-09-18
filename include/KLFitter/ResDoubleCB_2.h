#ifndef KLFITTER_RESDOUBLECB_2_H_
#define KLFITTER_RESDOUBLECB_2_H_

#include <cmath>
#include <iostream>
#include <vector>

#include "KLFitter/ResDoubleCBBase.h"

// ---------------------------------------------------------

namespace KLFitter {
class ResDoubleCB_2 : public ResDoubleCBBase {
 public:
  explicit ResDoubleCB_2(const char * filename);
  virtual ~ResDoubleCB_2();

  double GetMean(double x);
  double GetSigma(double x);
  double GetAlphaL(double x);
  double GetAlphaR(double x);
  double GetNL(double x);
  double GetNR(double x);

  double evaluate(const double x, const size_t start);

}; // ResDoubleCB_2
}  // namespace KLFitter

#endif  // KLFITTER_RESDOUBLECB_2_H_
