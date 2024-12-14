#ifndef NA_NEOHOOKEAN_H
#define NA_NEOHOOKEAN_H

#include "na/macros.h"
#include "Eigen/Core"

// Neo-Hookean [Smith et al. 2018] materials (without the logarithmic term)
class NeoHookean
{
public:
    NeoHookean() = default;
    NeoHookean(const double mu, const double lambda);

    double psi(const Eigen::Matrix3d &F) const;

    Eigen::Matrix3d gradF(const Eigen::Matrix3d &F) const;

    Eigen::Matrix9d hessF(const Eigen::Matrix3d &F) const;

private:
    double m_mu;
    double m_lambda;
    double m_ratio;
    double m_alpha;
};

#endif // !NA_NEOHOOKEAN_H
