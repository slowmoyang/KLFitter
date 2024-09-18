#ifndef KLFITTER_RESDOUBLECB_3_H_
#define KLFITTER_RESDOUBLECB_3_H_

#include <cmath>
#include <iostream>
#include <vector>

#include "KLFitter/ResDoubleCBBase.h"

// ---------------------------------------------------------

namespace KLFitter {
class ResDoubleCB_3 : public ResDoubleCBBase {
 public:
  explicit ResDoubleCB_3(const char * filename);
  virtual ~ResDoubleCB_3();

  double GetMean(double x);
  double GetSigma(double x);
  double GetAlphaL(double x);
  double GetAlphaR(double x);
  double GetNL(double x);
  double GetNR(double x);

  double evaluate(const double x, const size_t start);

}; // ResDoubleCB_3
}  // namespace KLFitter

#endif  // KLFITTER_RESDOUBLECB_3_H_
