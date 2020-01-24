/*
 *  rf_integral.cpp
 *
 *  Copyright (c) 2020 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

// #define QI_DEBUG_BUILD 1

#include "Args.h"
#include "JSON.h"
#include "Macro.h"
#include "Util.h"

#include "mupa_model_1c.h"
#include "mupa_model_mt.h"
#include "mupa_sequence.h"
#include "mupa_ss.hpp"

/*
 * Main
 */
int mupa_rf_main(int argc, char **argv) {
    Eigen::initParallel();
    args::ArgumentParser parser("Calculates the time integral of a pulse in JSON format "
                                "\nhttp://github.com/spinicist/QUIT");
    args::HelpFlag       help(parser, "HELP", "Show this help message", {'h', "help"});
    args::Flag           verbose(parser, "VERBOSE", "Print more messages", {'v', "verbose"});
    args::ValueFlag<int> threads(parser,
                                 "THREADS",
                                 "Use N threads (default=hardware limit or $QUIT_THREADS)",
                                 {'T', "threads"},
                                 QI::GetDefaultThreads());

    args::Positional<std::string> json_file(parser, "INPUT", "Input JSON file");
    // args::Positional<std::string> output_path(parser, "OUTPUT", "Output JSON file");

    args::Flag mt(parser, "MT", "Use MT model", {"mt"});
    // args::ValueFlag<double> g0(parser, "G0", "On-resonance absorption constant",
    // {"G0"}, 1.75e-5);

    QI::ParseArgs(parser, argc, argv, verbose, threads);
    QI::CheckPos(json_file);

    QI::Log(verbose, "Reading sequence parameters");
    json         doc = json_file ? QI::ReadJSON(json_file.Get()) : QI::ReadJSON(std::cin);
    MUPASequence sequence(doc["MUPA"]);

    double const R1 = 1. / 1.0;
    double const R2 = 1. / 0.1;
    double const PD = 100;

    if (mt) {

        using AugMat      = Eigen::Matrix<double, 5, 5>;
        using AugVec      = Eigen::Vector<double, 5>;
        double const f_b  = 0.2;
        double const f_f  = 1. - f_b;
        double const k    = 5.;
        double const k_bf = k * f_f;
        double const k_fb = k * f_b;

        AugMat R;
        R << -R2, 0, 0, 0, 0,       //
            0, -R2, 0, 0, 0,        //
            0, 0, -R1, 0, R1 * f_f, //
            0, 0, 0, -R1, R1 * f_b, //
            0, 0, 0, 0, 0;

        AugMat K;
        K << 0, 0, 0, 0, 0,       //
            0, 0, 0, 0, 0,        //
            0, 0, -k_fb, k_bf, 0, //
            0, 0, k_fb, -k_bf, 0, //
            0, 0, 0, 0, 0;

        AugMat const RpK   = R + K;
        auto         RF_MT = [&](double const &B1x, double const &B1y) -> AugMat {
            double const W = M_PI * 1.5e-6 * (B1x * B1x + B1y * B1y);
            AugMat       rf;
            rf << 0, 0, -B1y, 0, 0, //
                0, 0, B1x, 0, 0,    //
                B1y, -B1x, 0, 0, 0, //
                0, 0, 0, -W, 0,     //
                0, 0, 0, 0, 0;
            return rf;
        };

        AugVec m0{0, 0, f_f, f_b, 1.};

        AugMat const equ_check = (K * 10e-3).exp();
        AugVec const equ_m     = equ_check * m0;
        fmt::print("K\n{}\nRpK.exp()\n{}\n", K, equ_check);
        fmt::print("m0:   {}\nequ_m: {}\n", m0.transpose(), equ_m.transpose());

        for (auto p : sequence.prep_pulses) {
            auto const &name  = p.first;
            auto const &pulse = p.second;
            AugMat      C     = CalcPulse<AugMat>(pulse, RpK, RF_MT);
            AugVec      m     = PD * C * m0;
            fmt::print("Pulse Name: {}\t Eff_f: {}\t Eff_b: {}\n ", name, m[2] / f_f, m[3] / f_b);
        }

    } else {

        using AugMat = Eigen::Matrix<double, 4, 4>;
        using AugVec = Eigen::Vector<double, 4>;
        AugMat R;
        R << -R2, 0, 0, 0,      //
            0, -R2, 0, 0,       //
            0, 0, -R1, PD * R1, //
            0, 0, 0, 0;

        AugVec m0{0, 0, PD, 1.};
        for (int is = 0; is < sequence.size(); is++) {
            auto const &name  = sequence.prep[is];
            auto const &pulse = sequence.prep_pulses[name];
            AugMat      C     = CalcPulse<AugMat>(pulse, R, &RF_1c);
            AugVec      m     = C * m0;
            fmt::print("Pulse Name: {}\t m: {}\t Eff: {}\n", name, m.transpose(), m[2]);
        }
    }
    QI::Log(verbose, "Finished.");
    return EXIT_SUCCESS;
}
