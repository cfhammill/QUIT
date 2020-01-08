// #define QI_DEBUG_BUILD 1
#include "Macro.h"
#include "mupa_model_1c.h"
#include "mupa_ss.hpp"

using AugMat = Eigen::Matrix<double, 4, 4>;
using AugVec = Eigen::Vector<double, 4>;

AugMat RF_1c(double const &B1x, double const &B1y) {
    AugMat rf;
    rf << 0, 0, -B1y, 0, //
        0, 0, B1x, 0,    //
        B1y, -B1x, 0, 0, //
        0, 0, 0, 0;
    return rf;
}

auto MUPAModel::signal(VaryingArray const &v, FixedArray const &) const -> QI_ARRAY(double) {
    using T = double;

    T const &PD = v[0];
    T const &R1 = 1. / v[1];
    T const &R2 = 1. / v[2];

    QI_DB(PD);
    QI_DB(R1);
    QI_DB(R2);

    AugMat R;
    R << -R2, 0, 0, 0,      //
        0, -R2, 0, 0,       //
        0, 0, -R1, PD * R1, //
        0, 0, 0, 0;

    // Setup readout segment matrices
    AugMat const Ard   = ((R + RF_1c(sequence.FA / sequence.Trf, 0)) * sequence.Trf).exp();
    AugMat const Rrd   = (R * sequence.TR).exp();
    AugMat const S     = Eigen::DiagonalMatrix<double, 4, 4>({0, 0, 1., 1.}).toDenseMatrix();
    AugMat const ramp  = (R * sequence.Tramp).exp();
    AugMat const RUFIS = S * Rrd * Ard;
    AugMat const seg   = RUFIS.pow(sequence.SPS);

    // Setup pulse matrices
    std::vector<AugMat> prep_mats(sequence.size());
    for (int is = 0; is < sequence.size(); is++) {
        auto const &name  = sequence.prep[is];
        auto const &pulse = sequence.prep_pulses[name];
        AugMat      C     = CalcPulse<AugMat>(pulse, R, &RF_1c);
        prep_mats[is]     = C;
    }

    // First calculate the system matrix and get SS
    QI_DBMAT(R);
    AugMat X = AugMat::Identity();
    for (int is = 0; is < sequence.size(); is++) {
        X = ramp * seg * ramp * prep_mats[is] * X;
    }
    AugVec m_aug = SolveSteadyState(X);
    QI_DBMAT(X);
    QI_DBVEC(m_aug);
    // Now loop through the segments and record the signal for each
    Eigen::ArrayXd sig(sequence.size());
    for (int is = 0; is < sequence.size(); is++) {
        m_aug     = ramp * prep_mats[is] * m_aug;
        auto m_gm = GeometricAvg(RUFIS, seg, m_aug, sequence.SPS);
        sig[is]   = m_gm[2] * sin(sequence.FA);
        m_aug     = ramp * seg * m_aug;
    }
    QI_DBVEC(sig);
    return sig;
}