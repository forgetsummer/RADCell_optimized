#ifndef StandardNormalDistribution_h
#define StandardNormalDistribution_h
#include <cmath>
#include <vector>


class StandardNormalDistribution
{
 public:
  StandardNormalDistribution();
  ~StandardNormalDistribution();

  // Distribution functions
  virtual double pdf(const double& x) const;
  virtual double cdf(const double& x) const;

  // Inverse cumulative distribution function (aka the probit function)
  virtual double inv_cdf(const double& quantile) const;
  
  
  double RationalApproximation(double t);
  double NormalCDFInverse(double p);

};

#endif