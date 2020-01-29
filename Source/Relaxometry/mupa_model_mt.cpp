// #define QI_DEBUG_BUILD 1
#include "mupa_model_mt.h"
#include "Macro.h"
#include "mupa_ss.hpp"

using AugMat = Eigen::Matrix<double, 5, 5>; // Short for Augmented Matrix
using AugVec = Eigen::Vector<double, 5>;

auto MUPAMTModel::signal(VaryingArray const &v, FixedArray const &) const -> QI_ARRAY(double) {
    using T = double;

    T const &PD   = v[0];
    T const &R1_f = 1. / v[1];
    T const &R1_b = R1_f;
    T const &R2_f = 1. / v[2];
    T const &f_b  = v[3];
    T const &f_f  = 1. - f_b;
    T const &k    = 4.3;
    T const &k_bf = k * f_f;
    T const &k_fb = k * f_b;
    T const &B1   = v[4];

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

    AugMat const S = Eigen::DiagonalMatrix<double, 5, 5>({0, 0, 1., 1., 1.}).toDenseMatrix();

    AugMat const RpK = R + K;

    auto RF_MT = [&](double const &B1x, double const &B1y) -> AugMat {
        double const W = M_PI * G0 * (B1x * B1x + B1y * B1y);
        AugMat       rf;
        rf << 0, 0, -B1y * B1, 0, 0,      //
            0, 0, B1x * B1, 0, 0,         //
            B1y * B1, -B1x * B1, 0, 0, 0, //
            0, 0, 0, -W * B1 * B1, 0,     //
            0, 0, 0, 0, 0;
        return rf;
    };

    // Setup readout segment matrices
    AugMat const        Rrd  = (RpK * (sequence.TR - sequence.Trf)).exp();
    AugMat const        ramp = (RpK * sequence.Tramp).exp();
    std::vector<AugMat> TR_mats(sequence.FA.rows());
    std::vector<AugMat> seg_mats(sequence.FA.rows());
    for (int is = 0; is < sequence.FA.rows(); is++) {
        AugMat const Ard = ((RpK + RF_MT(sequence.FA[is] / sequence.Trf, 0)) * sequence.Trf).exp();
        TR_mats[is]      = S * Rrd * Ard;
        seg_mats[is]     = TR_mats[is].pow(sequence.SPS);
    }

    // Setup pulse matrices
    AugVec              m0{0, 0, 1, 1, 1};
    std::vector<AugMat> prep_mats(sequence.size());
    for (int is = 0; is < sequence.size(); is++) {
        auto const & name = sequence.prep[is];
        auto const & p    = sequence.prep_pulses[name];
        double const W    = M_PI * G0 * B1 * B1 * p.B1_sq_mean;
        AugMat       C;
        C << 0, 0, 0, 0, 0,                                                          //
            0, 0, 0, 0, 0,                                                           //
            0, 0, exp(-R2_f * p.T_t) * cos(p.FA), 0, f_f * (1 - exp(-R1_f * p.T_l)), //
            0, 0, 0, exp(-W * p.T_act), f_b * (1. - exp(-R1_b * p.T_l)),             //
            0, 0, 0, 0, 1;
        prep_mats[is] = C;
    }

    // First calculate the system matrix
    AugMat X = AugMat::Identity();
    for (int is = 0; is < sequence.size(); is++) {
        X = ramp * seg_mats[is] * ramp * S * prep_mats[is] * X;
    }
    AugVec m_ss = SolveSteadyState(X);

    // Now loop through the segments and record the signal for each
    Eigen::ArrayXd sig(sequence.size());
    QI_DBVEC(m_ss);
    AugVec m_aug = m_ss;
    for (int is = 0; is < sequence.size(); is++) {
        m_aug             = ramp * S * prep_mats[is] * m_aug;
        auto       m_gm   = GeometricAvg(TR_mats[is], seg_mats[is], m_aug, sequence.SPS);
        auto const signal = PD * m_gm[2] * sin(sequence.FA[is]);
        sig[is]           = signal;
        m_aug             = ramp * seg_mats[is] * m_aug;
        QI_DBVEC(m_gm);
    }
    QI_DBVEC(v);
    QI_DBVEC(sig);
    return sig;
}