/*
 *  qi_vfa_prep.cpp
 *
 *  Copyright (c) 2019 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <type_traits>

// #define QI_DEBUG_BUILD 1

#include "Args.h"
#include "ImageIO.h"
#include "Macro.h"
#include "Model.h"
#include "ModelFitFilter.h"
#include "SequenceBase.h"
#include "SimulateModel.h"
#include "Util.h"

#include "mupa_model_1c.h"
#include "mupa_model_mt.h"

template <typename ModelType> struct MUPACost {
    ModelType const &              model;
    typename ModelType::FixedArray fixed;
    QI_ARRAY(double) const data;

    template <typename T> bool operator()(T const *const vin, T *rin) const {
        Eigen::Map<typename ModelType::VaryingArray const> const varying(vin);
        Eigen::Map<QI_ARRAY(T)>                                  residuals(rin, data.rows());
        residuals = data - model.signal(varying, fixed);
        QI_DBVEC(residuals);
        return true;
    }
};

template <typename ModelType_> struct MUPAFit {
    // Boilerplate information required by ModelFitFilter
    static const bool Blocked = false; // = input is in blocks and outputs have multiple entries
    static const bool Indexed = false; // = the voxel index will be passed to the fit
    using RMSErrorType        = double;
    using FlagType            = int; // Almost always the number of iterations
    using ModelType           = ModelType_;
    ModelType model;

    // Have to tell the ModelFitFilter how many volumes we expect in each input
    int input_size(const int) const { return model.sequence.size(); }

    // This has to match the function signature that will be called in ModelFitFilter (which depends
    // on Blocked/Indexed. The return type is a simple struct indicating success, and on failure
    // also the reason for failure
    QI::FitReturnType
    fit(std::vector<Eigen::ArrayXd> const &inputs,  // Input: signal data
        Eigen::ArrayXd const &             fixed,   // Input: Fixed parameters
        typename ModelType::VaryingArray & varying, // Output: Varying parameters
        typename ModelType::RSDArray *     cov,
        RMSErrorType &                     rmse,      // Output: root-mean-square error
        std::vector<Eigen::ArrayXd> &      residuals, // Optional output: point residuals
        FlagType &                         iterations /* Usually iterations */) const {
        // First scale down the raw data so that PD will be roughly the same magnitude as other
        // parameters This is important for numerical stability in the optimiser

        QI_DBVEC(inputs[0]);

        double scale = inputs[0].maxCoeff();
        QI_DB(scale);
        if (scale < std::numeric_limits<double>::epsilon()) {
            varying = ModelType::VaryingArray::Zero();
            rmse    = 0.0;
            return {false, "Maximum data value was zero or less"};
        }
        Eigen::ArrayXd const data = inputs[0] / scale;

        // Setup Ceres
        ceres::Problem problem;
        using MUPADiff  = ceres::NumericDiffCostFunction<MUPACost<ModelType>,
                                                        ceres::CENTRAL,
                                                        ceres::DYNAMIC,
                                                        ModelType::NV>;
        auto *mupa_cost = new MUPADiff(new MUPACost<ModelType>{model, fixed, data},
                                       ceres::TAKE_OWNERSHIP,
                                       model.sequence.size());
        auto *loss      = new ceres::HuberLoss(1.0); // Don't know if this helps

        // This is where the parameters and cost functions actually get added to Ceres
        problem.AddResidualBlock(mupa_cost, loss, varying.data());

        // Set up parameter bounds
        for (int i = 0; i < ModelType::NV; i++) {
            problem.SetParameterLowerBound(varying.data(), i, model.lo[i]);
            problem.SetParameterUpperBound(varying.data(), i, model.hi[i]);
        }

        ceres::Solver::Options options;
        ceres::Solver::Summary summary;
        options.max_num_iterations  = 30;
        options.function_tolerance  = 1e-5;
        options.gradient_tolerance  = 1e-6;
        options.parameter_tolerance = 1e-4;
        options.logging_type        = ceres::SILENT;

        varying = model.start;
        ceres::Solve(options, &problem, &summary);
        if (!summary.IsSolutionUsable()) {
            return {false, summary.FullReport()};
        }
        iterations = summary.iterations.size();

        Eigen::ArrayXd const rs  = (data - model.signal(varying, fixed));
        double const         var = rs.square().sum();
        rmse                     = sqrt(var / data.rows()) * scale;
        if (residuals.size() > 0) {
            residuals[0] = rs * scale;
        }
        if (cov) {
            QI::GetRelativeStandardDeviation<ModelType>(
                problem, varying, var / (data.rows() - ModelType::NV), cov);
        }
        varying[0] *= scale; // Multiply signals/proton density back up

        return {true, ""};
    }
};

/*
 * Main
 */
int mupa_main(int argc, char **argv) {
    Eigen::initParallel();
    args::ArgumentParser parser("Calculates parametric maps from MUPA data "
                                "data.\nhttp://github.com/spinicist/QUIT");

    args::Positional<std::string> input_path(parser, "INPUT", "Input MUPA file");

    QI_COMMON_ARGS;

    args::Flag              mt(parser, "MT", "Use MT model", {"mt"});
    args::ValueFlag<double> g0(parser, "G0", "On-resonance absorption constant", {"G0"}, 1.75e-5);

    QI::ParseArgs(parser, argc, argv, verbose, threads);

    QI::CheckPos(input_path);

    QI::Log(verbose, "Reading sequence parameters");
    json doc = json_file ? QI::ReadJSON(json_file.Get()) : QI::ReadJSON(std::cin);

    MUPASequence sequence(doc["MUPA"]);

    auto process = [&](auto model, const std::string &model_name) {
        if (simulate) {
            QI::SimulateModel<decltype(model), false>(
                doc, model, {}, {input_path.Get()}, verbose, simulate.Get());
        } else {
            using FitType = MUPAFit<decltype(model)>;
            FitType fit{model};
            auto    fit_filter =
                QI::ModelFitFilter<FitType>::New(&fit, verbose, rsd, resids, subregion.Get());
            fit_filter->ReadInputs({input_path.Get()}, {}, mask.Get());
            fit_filter->Update();
            fit_filter->WriteOutputs(prefix.Get() + model_name);
        }
    };
    if (mt) {
        MUPAMTModel model{{}, sequence, g0.Get()};
        process(model, "MUPAMT_");
    } else {
        MUPAModel model{{}, sequence};
        process(model, "MUPA_");
    }
    QI::Log(verbose, "Finished.");
    return EXIT_SUCCESS;
}
