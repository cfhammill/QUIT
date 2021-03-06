/*
 *  ThreePoolModel.h
 *
 *  Copyright (c) 2018 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef QI_THREEPOOLMODEL_H
#define QI_THREEPOOLMODEL_H

#include "OnePoolSignals.h"
#include "SPGRSequence.h"
#include "SSFPSequence.h"
#include "SequenceGroup.h"
#include "TwoPoolModel.h"
#include <Eigen/Core>
#include <array>
#include <string>

namespace QI {

struct ThreePoolModel {
    using DataType      = double;
    using ParameterType = double;
    using SequenceType  = QI::SequenceGroup;

    static constexpr int NV = 10;
    static constexpr int ND = 0;
    static constexpr int NF = 2;
    static constexpr int NI = 2;

    static std::array<const std::string, NV> varying_names;
    static std::array<const std::string, NF> fixed_names;
    static const QI_ARRAYN(double, NF) fixed_defaults;

    SPGRSequence spgr;
    SSFPSequence ssfp;
    bool         scale_to_mean = false;
    TwoPoolModel two_pool;

    QI_ARRAYN(double, NV) bounds_lo;
    QI_ARRAYN(double, NV) bounds_hi;

    ThreePoolModel(SPGRSequence const &s1, SSFPSequence const &s2, const bool scale);
    bool   valid(const QI_ARRAYN(double, NV) & params) const; // For SRC
    size_t num_outputs() const;
    int    output_size(int i) const;

    Eigen::ArrayXd spgr_signal(const Eigen::ArrayXd &varying,
                               const QI_ARRAYN(double, NF) & fixed) const;

    Eigen::ArrayXd ssfp_signal(const Eigen::ArrayXd &varying,
                               const QI_ARRAYN(double, NF) & fixed) const;

    std::vector<Eigen::ArrayXd> signals(const Eigen::ArrayXd &varying,
                                        const QI_ARRAYN(double, NF) & fixed) const;
    Eigen::ArrayXd signal(const Eigen::ArrayXd &varying, const QI_ARRAYN(double, NF) & fixed) const;
};

} // End namespace QI

#endif // QI_THREEPOOLMODEL_H