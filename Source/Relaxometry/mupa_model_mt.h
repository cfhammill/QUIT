#pragma once

#include "ModelHelpers.h"
#include "mupa_sequence.h"

struct MUPAMTModel {
    using DataType      = double;
    using ParameterType = double;

    static constexpr int NV = 5; // Number of varying parameters
    // static constexpr int NV = 3;
    static constexpr int ND = 0; // Number of derived parameters
    static constexpr int NF = 0; // Number of fixed parameters
    static constexpr int NI = 1; // Number of inputs

    using VaryingArray = QI_ARRAYN(ParameterType, NV); // Type for the varying parameter array
    using FixedArray   = QI_ARRAYN(ParameterType, NF); // Type for the fixed parameter array

    // Sequence paramter structs
    MUPASequence &sequence;

    // Fitting start point and bounds
    // The PD will be scaled by the fitting function to keep parameters roughly the same magnitude
    VaryingArray const start{50.0, 1.0, 0.1, 0.1, 25.0};
    VaryingArray const lo{0.1, 0.5, 0.005, 1e-6, 0.1};
    VaryingArray const hi{500.0, 5.0, 5.0, 0.9, 100.0};
    // VaryingArray const start{30., 1., 0.1};
    // VaryingArray const lo{1, 0.01, 0.01};
    // VaryingArray const hi{150, 10.0, 10.0};

    // std::array<std::string, NV> const varying_names{"PD", "T1", "T2"};
    std::array<std::string, NV> const varying_names{"PD", "T1_f", "T2_f", "f_b", "k"};
    std::array<std::string, NF> const fixed_names{};
    // If fixed parameters not supplied, use these default values
    FixedArray const fixed_defaults{};

    auto signal(VaryingArray const &v, FixedArray const &) const -> QI_ARRAY(double);
};

template <> struct QI::NoiseFromModelType<MUPAMTModel> : QI::RealNoise {};
