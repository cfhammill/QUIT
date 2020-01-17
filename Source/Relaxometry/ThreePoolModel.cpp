/*
 *  ThreePoolModel.hxx
 *
 *  Copyright (c) 2018 Tobias Wood.
 *
 *  Based on code by Sean Deoni
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "ThreePoolModel.h"
#include "Helpers.h"
#include "Log.h"

#include <Eigen/Dense>
// There's a bug in matrix log that leads to duplicate symbol definitions if
// this is included in the header file
#include <unsupported/Eigen/MatrixFunctions>

using namespace std::literals;

namespace {

Eigen::ArrayXd SSFP1(const double &          PD,
                     const double &          T1,
                     const double &          T2,
                     const double &          f0,
                     const double &          B1,
                     const QI::SSFPSequence &s) {
    {
        const double         E1     = exp(-s.TR / T1);
        const double         E2     = exp(-s.TR / T2);
        const Eigen::ArrayXd d      = (1 - E1 * cos(B1 * s.FA) - (E2 * E2) * (E1 - cos(B1 * s.FA)));
        const Eigen::ArrayXd G      = sin(B1 * s.FA) * (1 - E1) / d;
        const double         a      = E2;
        const Eigen::ArrayXd b      = E2 * (1 - E1) * (1 + cos(B1 * s.FA)) / d;
        const double         theta0 = 2.0 * M_PI * f0;
        const Eigen::ArrayXd theta  = theta0 - s.PhaseInc;
        const double         psi    = theta0 / 2.0;
        const Eigen::ArrayXd cos_th = cos(theta);
        const Eigen::ArrayXd sin_th = sin(theta);
        const double         cos_psi = cos(psi);
        const double         sin_psi = sin(psi);
        const Eigen::ArrayXd re_m =
            (cos_psi - a * (cos_th * cos_psi - sin_th * sin_psi)) * G / (1.0 - b * cos_th);
        const Eigen::ArrayXd im_m =
            (sin_psi - a * (cos_th * sin_psi + sin_th * cos_psi)) * G / (1.0 - b * cos_th);
        const Eigen::ArrayXd result = PD * sqrt(re_m.square() + im_m.square());
        // std::cout << PD << " " << T1 << " " << T2 << " " << f0 << " " << B1 << ":" <<
        // result.transpose() << std::endl;
        return result;
    }
}
} // namespace

namespace QI {

/* The inner two-pool model must not be scaled, as scaling will be done once signals are added */
ThreePoolModel::ThreePoolModel(SPGRSequence const &s1, SSFPSequence const &s2, const bool scale) :
    spgr{s1}, ssfp{s2}, scale_to_mean{scale}, two_pool{s1, s2, false} {
    bounds_lo << 1.0, 0.300, 0.010, 0.9, 0.040, 3.5, 1.0, 0.025, 0.001, 0.001;
    bounds_hi << 1.0, 0.800, 0.030, 1.5, 0.150, 5.0, 3.5, 0.600, 0.350, 0.999;
}

size_t ThreePoolModel::num_outputs() const {
    return 2;
}

int ThreePoolModel::output_size(int i) const {
    if (i == 0) {
        return spgr.size();
    } else if (i == 1) {
        return ssfp.size();
    } else {
        QI::Fail("Invalid output size: {}", i);
    }
}

bool ThreePoolModel::valid(const QI_ARRAYN(double, NV) & params) const {
    // Negative T1/T2 makes no sense
    if ((params[1] <= 0.) || (params[2] <= 0.))
        return false;
    else {
        if ((params[1] < params[3]) && (params[2] < params[4]) && (params[3] < params[5]) &&
            (params[4] < params[6]) && ((params[8] + params[9]) <= 1.0))
            return true;
        else
            return false;
    }
}

std::vector<Eigen::ArrayXd> ThreePoolModel::signals(const Eigen::ArrayXd &v,
                                                    const QI_ARRAYN(double, NF) & f) const {
    return {spgr_signal(v, f), ssfp_signal(v, f)};
}

Eigen::ArrayXd ThreePoolModel::signal(const Eigen::ArrayXd &v,
                                      const QI_ARRAYN(double, NF) & f) const {
    auto           sigs = signals(v, f);
    Eigen::ArrayXd sig(spgr.size() + ssfp.size());
    sig.head(spgr.size()) = sigs[0];
    sig.tail(ssfp.size()) = sigs[1];
    return sig;
}

Eigen::ArrayXd ThreePoolModel::spgr_signal(const Eigen::ArrayXd &v,
                                           const QI_ARRAYN(double, NF) & fixed) const {
    double f_ab = 1. - v[9];
    QI_ARRAYN(double, TwoPoolModel::NV) two_pool_varying;
    two_pool_varying << v[0] * f_ab, v[1], v[2], v[3], v[4], v[7], v[8] / f_ab;
    Eigen::VectorXd m_ab   = two_pool.spgr_signal(two_pool_varying, fixed);
    Eigen::VectorXd m_c    = SPGRSignal(v[0] * v[9], v[5], fixed[1], spgr);
    Eigen::VectorXd signal = m_ab + m_c;
    // std::cout << "SPGR\n" << m_ab.transpose() << "\n" << m_c.transpose() << "\n" <<
    // signal.transpose() << std::endl;
    if (scale_to_mean) {
        signal /= signal.mean();
    }
    return signal;
}

Eigen::ArrayXd ThreePoolModel::ssfp_signal(const Eigen::ArrayXd &v,
                                           const QI_ARRAYN(double, NF) & fixed) const {
    double f_ab = 1. - v[9];
    QI_ARRAYN(double, TwoPoolModel::NV) two_pool_varying;
    two_pool_varying << v[0] * f_ab, v[1], v[2], v[3], v[4], v[7], v[8] / f_ab;
    Eigen::VectorXd m_ab   = two_pool.ssfp_signal(two_pool_varying, fixed);
    Eigen::VectorXd m_c    = SSFP1(v[0] * v[9], v[5], v[6], fixed[0], fixed[1], ssfp);
    Eigen::VectorXd signal = m_ab + m_c;
    if (scale_to_mean) {
        signal /= signal.mean();
    }
    return signal;
}

} // End namespace QI
