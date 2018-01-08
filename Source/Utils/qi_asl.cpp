/*
 *  qi_asl.cpp
 *
 *  Copyright (c) 2018 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <iostream>
#include <Eigen/Dense>

#include "Util.h"
#include "IO.h"
#include "Args.h"
#include "Models.h"
#include "Sequences.h"
#include "Types.h"

class CASL : public QI::ApplyF::Algorithm {
protected:
    const double m_T1, m_alpha, m_lambda, m_PLD, m_LD;
    const int m_inputsize;
public:
    CASL(const double T1, const double alpha, const double lambda,
         const double LD, const double PLD, const int inputsize) :
        m_T1(T1), m_alpha(alpha), m_lambda(lambda), m_LD(LD), m_PLD(PLD), m_inputsize(inputsize)
    {}

    size_t numInputs() const override  { return 1; }
    size_t numConsts() const override  { return 0; }
    size_t numOutputs() const override { return 1; }
    size_t dataSize() const override   { return m_inputsize; }
    float zero() const override { return 0.f; }
    std::vector<float> defaultConsts() const override {
        std::vector<float> def(0, 0);
        return def;
    }

    bool apply(const std::vector<TInput> &inputs, const std::vector<TConst> &consts,
               std::vector<TOutput> &outputs, TConst &residual,
               TInput &resids, TIters &its) const override
    {
        const int ts_size = m_inputsize / 2;
        const Eigen::Map<const Eigen::ArrayXf, 0, Eigen::InnerStride<>> even(inputs[0].GetDataPointer(), ts_size, Eigen::InnerStride<>(2));
        const Eigen::Map<const Eigen::ArrayXf, 0, Eigen::InnerStride<>> odd(inputs[0].GetDataPointer() + 1, ts_size, Eigen::InnerStride<>(2));

        const double diff = (odd.cast<double>() - even.cast<double>()).mean();
        const double SI_PD = odd.cast<double>().mean();
        const double CBF = (6000 * m_lambda * diff * exp(m_PLD / m_T1)) / 
                           (2. * m_alpha * m_T1 * SI_PD * (1. - exp(-m_LD / m_T1)));
        // std::cout << "l " << m_lambda << " diff " << diff << " PLD " << m_PLD << " T1 " << m_T1 << " e(PLD) " << exp(m_PLD / m_T1) << std::endl;
        // std::cout << "a " << m_alpha << " PD " << SI_PD << " LD " << m_LD << " (1 - exp()) " << (1. - exp(-m_LD / m_T1)) << std::endl;
        // std::cout << "CBF: " << CBF << std::endl;
        outputs[0] = CBF;
        residual = 0;
        resids.Fill(0.);
        its = 0;
        return true;
    }
};

/*
 * Main
 */
int main(int argc, char **argv) {
    Eigen::initParallel();
    args::ArgumentParser parser("Calculates CBF from ASL data.\nhttp://github.com/spinicist/QUIT");
    
    args::Positional<std::string> input_path(parser, "ASL_FILE", "Input ASL file");
    
    args::HelpFlag help(parser, "HELP", "Show this help message", {'h', "help"});
    args::Flag     verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::Flag     noprompt(parser, "NOPROMPT", "Suppress input prompts", {'n', "no-prompt"});
    args::ValueFlag<int> threads(parser, "THREADS", "Use N threads (default=4, 0=hardware limit)", {'T', "threads"}, 4);
    args::ValueFlag<std::string> outarg(parser, "OUTPREFIX", "Add a prefix to output filename", {'o', "out"});
    args::ValueFlag<std::string> mask(parser, "MASK", "Only process voxels within the mask", {'m', "mask"});
    args::ValueFlag<double> T1_blood(parser, "BLOOD T1", "Value of blood T1 to use (seconds), default 2.429", {'t', "T1"}, 2.429);
    args::ValueFlag<double> alpha(parser, "ALPHA", "Labelling efficiency, default 0.9", {'a', "alpha"}, 0.9);
    args::ValueFlag<double> lambda(parser, "LAMBDA", "Blood-brain partition co-efficent, default 0.9 mL/g", {'l', "lambda"}, 0.9);
    args::ValueFlag<std::string> subregion(parser, "SUBREGION", "Process subregion starting at voxel I,J,K with size SI,SJ,SK", {'s', "subregion"});
    QI::ParseArgs(parser, argc, argv);
    bool prompt = !noprompt;

    if (verbose) std::cout << "Reading ASL data from: " << QI::CheckPos(input_path) << std::endl;
    auto input = QI::ReadVectorImage(QI::CheckPos(input_path));

    double LD, PLD;
    if (prompt) std::cout << "Enter label delay (seconds): " << std::flush;
    std::cin >> LD;
    if (prompt) std::cout << "Enter post-label delay (seconds): " << std::flush;
    std::cin >> PLD;

    auto apply = QI::ApplyF::New();
    if (verbose) {
        std::cout << "T1 blood: " << T1_blood.Get() << " Alpha: " << alpha.Get() << " Lambda: " << lambda.Get() << "\n";
        std::cout << "Label time: " << LD << " Post-label delay: " << PLD << std::endl;
    }
    std::shared_ptr<CASL> algo = std::make_shared<CASL>(T1_blood.Get(), alpha.Get(), lambda.Get(),
                                                        LD, PLD, input->GetNumberOfComponentsPerPixel());
    apply->SetVerbose(verbose);
    apply->SetAlgorithm(algo);
    apply->SetOutputAllResiduals(false);
    if (verbose) std::cout << "Using " << threads.Get() << " threads" << std::endl;
    apply->SetPoolsize(threads.Get());
    apply->SetSplitsPerThread(threads.Get());
    apply->SetInput(0, input);
    if (mask) apply->SetMask(QI::ReadImage(mask.Get()));
    if (subregion) {
        apply->SetSubregion(QI::RegionArg(args::get(subregion)));
    }
    if (verbose) {
        std::cout << "Processing" << std::endl;
        auto monitor = QI::GenericMonitor::New();
        apply->AddObserver(itk::ProgressEvent(), monitor);
    }
    apply->Update();
    if (verbose) {
        std::cout << "Elapsed time was " << apply->GetTotalTime() << "s" << std::endl;
        std::cout << "Writing results files." << std::endl;
    }
    const std::string outPrefix = outarg ? outarg.Get() : QI::Basename(input_path.Get());
    QI::WriteImage(apply->GetOutput(0), outPrefix + "CBF" + QI::OutExt());
    return EXIT_SUCCESS;
}