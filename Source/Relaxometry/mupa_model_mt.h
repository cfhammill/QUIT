#pragma once

#include "Model.h"
#include "mupa_sequence.h"

struct MUPAMTModel : QI::Model<double, double, 4, 0> {
    MUPASequence &     sequence;
    double             G0;
    VaryingArray const start{50.0, 1.0, 0.1, 0.1};
    VaryingArray const lo{0.1, 0.5, 0.005, 1e-6};
    VaryingArray const hi{500.0, 5.0, 5.0, 0.9};

    std::array<std::string, NV> const varying_names{"PD", "T1_f", "T2_f", "f_b"};
    std::array<std::string, NF> const fixed_names{};

    auto signal(VaryingArray const &v, FixedArray const &) const -> QI_ARRAY(double);
};

template <> struct QI::NoiseFromModelType<MUPAMTModel> : QI::RealNoise {};
