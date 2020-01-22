// #define QI_DEBUG_BUILD 1
#include "Macro.h"
#include "mupa_model_mt.h"
#include "mupa_ss.hpp"

using AugMat = Eigen::Matrix<double, 5, 5>; // Short for Augmented Matrix
using AugVec = Eigen::Vector<double, 5>;

auto MUPAMTModel::signal(VaryingArray const &v, FixedArray const &) const -> QI_ARRAY(double) {
    using T = double;

    T const &PD   = v[0];
    T const &R1_f = 1. / v[1];
    T const &R1_b = 1.; // R1_f;
    T const &R2_f = 1. / v[2];
    T const &f_b  = v[3];
    T const &f_f  = 1. - f_b;
    T const &k    = 50.; // v[4];
    T const &k_bf = k * f_f;
    T const &k_fb = k * f_b;

    AugMat R;
    R << -R2_f, 0, 0, 0, 0,         //
        0, -R2_f, 0, 0, 0,          //
        0, 0, -R1_f, 0, f_f * R1_f, //
        0, 0, 0, -R1_b, f_b * R1_b, //
        0, 0, 0, 0, 0;

    AugMat K;
    K << 0, 0, 0, 0, 0,       //
        0, 0, 0, 0, 0,        //
        0, 0, -k_fb, k_bf, 0, //
        0, 0, k_fb, -k_bf, 0, //
        0, 0, 0, 0, 0;

    AugMat const RpK = R + K;

    auto RF_MT = [&](double const &B1x, double const &B1y) -> AugMat {
        double const W = M_PI * G0 * (B1x * B1x + B1y * B1y);
        AugMat       rf;
        rf << 0, 0, -B1y, 0, 0, //
            0, 0, B1x, 0, 0,    //
            B1y, -B1x, 0, 0, 0, //
            0, 0, 0, -W, 0,     //
            0, 0, 0, 0, 0;
        return rf;
    };

    QI_DBMAT(R);
    QI_DBMAT(K);

    // Setup readout segment matrices
    AugMat const Ard   = ((RpK + RF_MT(sequence.FA / sequence.Trf, 0)) * sequence.Trf).exp();
    AugMat const Rrd   = (RpK * (sequence.TR - sequence.Trf)).exp();
    AugMat const S     = Eigen::DiagonalMatrix<double, 5, 5>({0, 0, 1., 1., 1.}).toDenseMatrix();
    AugMat const ramp  = (RpK * sequence.Tramp).exp();
    AugMat const RUFIS = S * Rrd * Ard;
    AugMat const seg   = RUFIS.pow(sequence.SPS);

    // Setup pulse matrices
    std::vector<AugMat> prep_mats(sequence.size());
    for (int is = 0; is < sequence.size(); is++) {
        auto const &name  = sequence.prep[is];
        auto const &pulse = sequence.prep_pulses[name];
        AugMat      C     = CalcPulse<AugMat>(pulse, RpK, RF_MT);
        prep_mats[is]     = C;
    }

    // First calculate the system matrix
    AugMat X = AugMat::Identity();
    for (int is = 0; is < sequence.size(); is++) {
        X = ramp * seg * ramp * S * prep_mats[is] * X;
    }
    AugVec m_aug = SolveSteadyState(X);

    // Now loop through the segments and record the signal for each
    Eigen::ArrayXd sig(sequence.size());
    for (int is = 0; is < sequence.size(); is++) {
        m_aug             = ramp * S * prep_mats[is] * m_aug;
        auto       m_gm   = GeometricAvg(RUFIS, seg, m_aug, sequence.SPS);
        auto const signal = PD * m_gm[2] * sin(sequence.FA);
        sig[is]           = signal;
        m_aug             = ramp * seg * m_aug;
        QI_DB(signal);
    }
    return sig;
}