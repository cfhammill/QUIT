/*
 *  SSFP.cpp
 *
 *  Copyright (c) 2018 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "SSFPSequence.h"

#define FA_PHASE_CHECK()\
    if (FA.rows() != PhaseInc.rows()) {\
        QI_FAIL("While reading " << this->name() << " number of phase increments " << PhaseInc.rows() <<\
                " did not match the number of flip angles " << FA.rows() );\
    }

namespace QI {

size_t SSFPBase::size() const {
    return FA.rows();
}

Eigen::ArrayXd SSFPSequence::weights(const double f0) const {
    Eigen::ArrayXd offset = PhaseInc + 2.*M_PI*f0*TR;
    Eigen::ArrayXd weights = 0.75 * (offset / 2).sin().square();
    return weights;
}

Eigen::ArrayXcd SSFPSequence::signal(std::shared_ptr<Model> m, const Eigen::VectorXd &p) const {
    return m->SSFP(p, FA, TR, PhaseInc);
}

void SSFPSequence::load(cereal::JSONInputArchive &ar) {
    ar(cereal::make_nvp("TR", TR));
    QI_SEQUENCE_LOAD_DEGREES( FA );
    QI_SEQUENCE_LOAD_DEGREES( PhaseInc );
    FA_PHASE_CHECK()
}

void SSFPSequence::save(cereal::JSONOutputArchive &ar) const {
    ar(cereal::make_nvp("TR", TR));
    QI_SEQUENCE_SAVE_DEGREES( FA );
    QI_SEQUENCE_SAVE_DEGREES( PhaseInc );
}

Eigen::ArrayXcd SSFPEchoSequence::signal(std::shared_ptr<Model> m, const Eigen::VectorXd &p) const {
    return m->SSFPEcho(p, FA, TR, PhaseInc);
}

Eigen::ArrayXd SSFPEchoSequence::signal_magnitude(std::shared_ptr<Model> m, const Eigen::VectorXd &p) const {
    return m->SSFPEchoMagnitude(p, FA, TR, PhaseInc);
}

void SSFPEchoSequence::load(cereal::JSONInputArchive &ar) {
    ar(cereal::make_nvp("TR", TR));
    QI_SEQUENCE_LOAD_DEGREES( FA );
    QI_SEQUENCE_LOAD_DEGREES( PhaseInc );
    FA_PHASE_CHECK()
}

void SSFPEchoSequence::save(cereal::JSONOutputArchive &ar) const {
    ar(cereal::make_nvp("TR", TR));
    QI_SEQUENCE_SAVE_DEGREES( FA );
    QI_SEQUENCE_SAVE_DEGREES( PhaseInc );
}

Eigen::ArrayXd SSFPFiniteSequence::weights(const double f0) const {
    Eigen::ArrayXd offset = PhaseInc + 2.*M_PI*f0*TR;
    Eigen::ArrayXd weights = 0.75 * (offset / 2).sin().square();
    return weights;
}

Eigen::ArrayXcd SSFPFiniteSequence::signal(std::shared_ptr<Model> m, const Eigen::VectorXd &p) const {
    return m->SSFPFinite(p, FA, TR, Trf, PhaseInc);
}

void SSFPFiniteSequence::load(cereal::JSONInputArchive &ar) {
    ar(cereal::make_nvp("TR", TR));
    ar(cereal::make_nvp("Trf", Trf));
    QI_SEQUENCE_LOAD_DEGREES( FA );
    QI_SEQUENCE_LOAD_DEGREES( PhaseInc );
    FA_PHASE_CHECK()
}

void SSFPFiniteSequence::save(cereal::JSONOutputArchive &ar) const {
    ar(cereal::make_nvp("TR", TR));
    ar(cereal::make_nvp("Trf", Trf));
    QI_SEQUENCE_SAVE_DEGREES( FA );
    QI_SEQUENCE_SAVE_DEGREES( PhaseInc );
}

Eigen::ArrayXd SSFPGSSequence::weights(const double f0) const {
    return Eigen::ArrayXd::Ones(size()); // Weight SPGR images higher than SSFP
}

Eigen::ArrayXcd SSFPGSSequence::signal(std::shared_ptr<Model> m, const Eigen::VectorXd &p) const {
    return m->SSFP_GS(p, FA, TR);
}

void SSFPGSSequence::load(cereal::JSONInputArchive &ar) {
    ar(cereal::make_nvp("TR", TR));
    QI_SEQUENCE_LOAD_DEGREES( FA );
}

void SSFPGSSequence::save(cereal::JSONOutputArchive &ar) const {
    ar(cereal::make_nvp("TR", TR));
    QI_SEQUENCE_SAVE_DEGREES( FA );
}

} // End namespace QI
