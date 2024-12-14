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

static Eigen::Matrix3d partialJ_partialF(const Eigen::Matrix3d &F)
{
    Eigen::Matrix3d pJpF;

    pJpF.col(0) = F.col(1).cross(F.col(2));
    pJpF.col(1) = F.col(2).cross(F.col(0));
    pJpF.col(2) = F.col(0).cross(F.col(1));

    return pJpF;
}

Eigen::Matrix3d NeoHookean::gradF(const Eigen::Matrix3d &F) const
{
    const Eigen::Matrix3d pJpF = partialJ_partialF(F);
    const double Jminus1 = F.determinant() - m_alpha;
    return m_mu * F + m_lambda * Jminus1 * pJpF;
}

static Eigen::Matrix3d hat_matrix(const Eigen::Vector3d& v)
{
    Eigen::Matrix3d hat;
    hat << 0.0, -v(2), v(1), v(2), 0.0, -v(0), -v(1), v(0), 0.0;
    return hat;
}

Eigen::Matrix9d NeoHookean::hessF(const Eigen::Matrix3d &F) const
{
    const double scale = m_lambda * (F.determinant() - m_alpha);

    Eigen::Matrix3d ahat = hat_matrix(scale * F.col(0));
    Eigen::Matrix3d bhat = hat_matrix(scale * F.col(1));
    Eigen::Matrix3d chat = hat_matrix(scale * F.col(2));

    Eigen::Matrix9d FJ;
    FJ.block<3, 3>(0, 0).setZero();
    FJ.block<3, 3>(0, 3) = -chat;
    FJ.block<3, 3>(0, 6) = bhat;

    FJ.block<3, 3>(3, 0) = chat;
    FJ.block<3, 3>(3, 3).setZero();
    FJ.block<3, 3>(3, 6) = -ahat;

    FJ.block<3, 3>(6, 0) = -bhat;
    FJ.block<3, 3>(6, 3) = ahat;
    FJ.block<3, 3>(6, 6).setZero();

    const Eigen::Vector9d pJpF = partialJ_partialF(F).reshaped();
    return m_mu * Eigen::Matrix9d::Identity() + m_lambda * pJpF * pJpF.transpose() + FJ;
}
