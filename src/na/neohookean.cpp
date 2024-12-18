#include "na/neohookean.h"
#include "Eigen/Core"
#include "Eigen/Geometry"

NeoHookean::NeoHookean(const double mu, const double lambda)
{
    m_mu = mu;
    m_lambda = lambda;
    m_ratio = m_mu / m_lambda;
    m_alpha = 1.0 + m_ratio;
}

double NeoHookean::psi(const Eigen::Matrix3d &F) const
{
    const double Ic = F.squaredNorm();
    const double Jminus1 = F.determinant() - m_alpha;
    return 0.5 * (m_mu * (Ic - 3.0) + m_lambda * Jminus1 * Jminus1);
}

Eigen::Matrix3d partialJ_partialF(const Eigen::Matrix3d &F)
{
    Eigen::Matrix3d pJpF;

    pJpF.col(0) = F.col(1).cross(F.col(2));
    pJpF.col(1) = F.col(2).cross(F.col(0));
    pJpF.col(2) = F.col(0).cross(F.col(1));

    return pJpF;
}

Eigen::Vector9d NeoHookean::gradF(const Eigen::Matrix3d &F) const
{
    const double Jminus1 = F.determinant() - m_alpha;

    Eigen::Vector9d grad;
    grad(0) = +(F(1, 1) * F(2, 2) - F(2, 1) * F(1, 2)) * m_lambda * Jminus1 + m_mu * F(0, 0);
    grad(1) = -(F(0, 1) * F(2, 2) - F(2, 1) * F(0, 2)) * m_lambda * Jminus1 + m_mu * F(1, 0);
    grad(2) = +(F(0, 1) * F(1, 2) - F(1, 1) * F(0, 2)) * m_lambda * Jminus1 + m_mu * F(2, 0);
    grad(3) = -(F(1, 0) * F(2, 2) - F(2, 0) * F(1, 2)) * m_lambda * Jminus1 + m_mu * F(0, 1);
    grad(4) = +(F(0, 0) * F(2, 2) - F(2, 0) * F(0, 2)) * m_lambda * Jminus1 + m_mu * F(1, 1);
    grad(5) = -(F(0, 0) * F(1, 2) - F(1, 0) * F(0, 2)) * m_lambda * Jminus1 + m_mu * F(2, 1);
    grad(6) = +(F(1, 0) * F(2, 1) - F(2, 0) * F(1, 1)) * m_lambda * Jminus1 + m_mu * F(0, 2);
    grad(7) = -(F(0, 0) * F(2, 1) - F(2, 0) * F(0, 1)) * m_lambda * Jminus1 + m_mu * F(1, 2);
    grad(8) = +(F(0, 0) * F(1, 1) - F(1, 0) * F(0, 1)) * m_lambda * Jminus1 + m_mu * F(2, 2);
    return grad;
}

Eigen::Matrix9d NeoHookean::hessF(const Eigen::Matrix3d &F) const
{
    const double Jminus1 = F.determinant() - m_alpha;

    Eigen::Vector9d C;
    C(0) = +(F(1, 1) * F(2, 2) - F(2, 1) * F(1, 2));
    C(1) = -(F(0, 1) * F(2, 2) - F(2, 1) * F(0, 2));
    C(2) = +(F(0, 1) * F(1, 2) - F(1, 1) * F(0, 2));
    C(3) = -(F(1, 0) * F(2, 2) - F(2, 0) * F(1, 2));
    C(4) = +(F(0, 0) * F(2, 2) - F(2, 0) * F(0, 2));
    C(5) = -(F(0, 0) * F(1, 2) - F(1, 0) * F(0, 2));
    C(6) = +(F(1, 0) * F(2, 1) - F(2, 0) * F(1, 1));
    C(7) = -(F(0, 0) * F(2, 1) - F(2, 0) * F(0, 1));
    C(8) = +(F(0, 0) * F(1, 1) - F(1, 0) * F(0, 1));

    Eigen::Matrix9d hess = C * C.transpose();
    hess(0, 4) += +F(2, 2) * Jminus1;
    hess(0, 8) += +F(1, 1) * Jminus1;
    hess(0, 5) += -F(1, 2) * Jminus1;
    hess(0, 7) += -F(2, 1) * Jminus1;

    hess(1, 3) += -F(2, 2) * Jminus1;
    hess(1, 8) += -F(0, 1) * Jminus1;
    hess(1, 5) += +F(0, 2) * Jminus1;
    hess(1, 6) += +F(2, 1) * Jminus1;

    hess(2, 3) += +F(1, 2) * Jminus1;
    hess(2, 7) += +F(0, 1) * Jminus1;
    hess(2, 4) += -F(0, 2) * Jminus1;
    hess(2, 6) += -F(1, 1) * Jminus1;

    hess(3, 1) += -F(2, 2) * Jminus1;
    hess(3, 8) += -F(1, 0) * Jminus1;
    hess(3, 2) += +F(1, 2) * Jminus1;
    hess(3, 7) += +F(2, 0) * Jminus1;

    hess(4, 0) += +F(2, 2) * Jminus1;
    hess(4, 8) += +F(0, 0) * Jminus1;
    hess(4, 2) += -F(0, 2) * Jminus1;
    hess(4, 6) += -F(2, 0) * Jminus1;

    hess(5, 0) += -F(1, 2) * Jminus1;
    hess(5, 7) += -F(0, 0) * Jminus1;
    hess(5, 1) += +F(0, 2) * Jminus1;
    hess(5, 6) += +F(1, 0) * Jminus1;

    hess(6, 1) += +F(2, 1) * Jminus1;
    hess(6, 5) += +F(1, 0) * Jminus1;
    hess(6, 2) += -F(1, 1) * Jminus1;
    hess(6, 4) += -F(2, 0) * Jminus1;

    hess(7, 0) += -F(2, 1) * Jminus1;
    hess(7, 5) += -F(0, 0) * Jminus1;
    hess(7, 2) += +F(0, 1) * Jminus1;
    hess(7, 3) += +F(2, 0) * Jminus1;

    hess(8, 0) += +F(1, 1) * Jminus1;
    hess(8, 4) += +F(0, 0) * Jminus1;
    hess(8, 1) += -F(0, 1) * Jminus1;
    hess(8, 3) += -F(1, 0) * Jminus1;

    hess *= m_lambda;
    hess.diagonal().array() += m_mu;
    return hess;
}
