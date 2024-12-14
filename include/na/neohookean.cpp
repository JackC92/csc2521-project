#include "na/neohookean.h"
#include "Eigen/Core"
#include "Eigen/Geometry"

NeoHookean::NeoHookean(const double mu, const double lambda)
{
    m_mu = mu;
    m_lambda = lambda;
    m_ratio = mu / lambda;
    m_alpha = 1.0 + m_ratio;
}

double NeoHookean::psi(const Eigen::Matrix3d &F) const
{
    // First Right Cauchy-Green invariant
    const double Ic = F.squaredNorm();
    const double Jminus1 = F.determinant() - m_alpha;
    return 0.5 * (m_mu * (Ic - 3.0) + m_lambda * Jminus1 * Jminus1);
}

static Eigen::Matrix3d partialJ_partialF(const Eigen::Matrix3d &F)
{
    Eigen::Matrix3d pJpF;

    pJpF.col(0) = F.col(1).cross(F.col(2));
    pJpF.col(1) = F.col(2).cross(F.col(0));
    pJpF.col(2) = F.col(0).cross(F.col(1));

    return pJpF;
}

Eigen::Matrix3d NeoHookean::pk1(const Eigen::Matrix3d &F) const
{
    const Eigen::Matrix3d pJpF = partialJ_partialF(F);
    const double Jminus1 = F.determinant() - m_alpha;
    return m_mu * F + m_lambda * Jminus1 * pJpF;
}
