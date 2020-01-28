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

    QI::ParseArgs(parser, argc, argv, verbose, threads);
    QI::CheckPos(json_file);

    QI::Log(verbose, "Reading sequence parameters");
    json         doc = json_file ? QI::ReadJSON(json_file.Get()) : QI::ReadJSON(std::cin);
    MUPASequence sequence(doc["MUPA"]);

    double const R1 = 1. / 1.0;
    double const R2 = 1. / 0.1;
    double const PD = 1.;

    using AugMat = Eigen::Matrix<double, 4, 4>;
    using AugVec = Eigen::Vector<double, 4>;
    AugMat R;
    R << -R2, 0, 0, 0,      //
        0, -R2, 0, 0,       //
        0, 0, -R1, PD * R1, //
        0, 0, 0, 0;

    AugVec m0{0, 0, PD, 1.};
    for (auto const &p : sequence.prep_pulses) {
        auto const &name  = p.first;
        auto const &pulse = p.second;

        double int_b1     = 0;
        double int_b1_sq  = 0;
        double eff_tv     = 0;
        double eff_long   = 0;
        double t_total    = 0;
        double tact_total = 0;
        AugMat C_rf       = AugMat::Identity();
        AugMat C_both     = AugMat::Identity();
        AugVec m_rf       = m0;
        for (long ii = 0; ii < pulse.B1x.rows(); ii++) {
            double const &t     = pulse.timestep[ii];
            double        b1_sq = (pulse.B1x[ii] * pulse.B1x[ii] + pulse.B1y[ii] * pulse.B1y[ii]);
            int_b1 += sqrt(b1_sq) * t;
            int_b1_sq += b1_sq * t;

            AugMat const rf     = RF_1c(pulse.B1x[ii], pulse.B1y[ii]);
            AugMat const A_rf   = (rf * t).exp();
            AugMat const both   = R + rf;
            AugMat const A_both = (both * t).exp();

            C_rf   = A_rf * C_rf;
            C_both = A_both * C_both;
            m_rf   = A_rf * m_rf;
            t_total += t;
            if (b1_sq > 0.) {
                tact_total += t;
            }
            eff_tv += sqrt(m_rf[0] * m_rf[0] + m_rf[1] * m_rf[1]) * t;
            eff_long += std::abs(m_rf[2]) * t;
        }

        double       eff_flip = atan2(m_rf[2], m_rf.head(2).norm()) * 180 / M_PI - 90;
        double const W        = M_PI * 1.4e-5 * int_b1_sq / tact_total;
        fmt::print("Pulse Name: {}\n\ttime {}\n\tactive {}\n\tint_b1 {}\n\tint_b1_sq {}\n\tW {} "
                   "Sat {}\n\tEffective Transverse Time: {} "
                   "\n\tEffective Longitudinal Time: {}\n\tEffective flip-angle: {}\n\tfinal: {}\n",
                   name,
                   t_total,
                   tact_total,
                   int_b1,
                   int_b1_sq,
                   W,
                   exp(-W * tact_total),
                   eff_tv,
                   eff_long,
                   eff_flip,
                   m_rf.transpose());
    }

    double gB1       = sequence.FA[0] / sequence.Trf;
    double int_b1_sq = (gB1 * gB1) * sequence.Trf;
    double W         = M_PI * 1.4e-5 * int_b1_sq / sequence.Trf;
    fmt::print("Excitation:\n\tFA: {} Trf: {}us\n\tint_b1_sq {}\n\tW {} Sat {}\n",
               sequence.FA.transpose() * 180 / M_PI,
               sequence.Trf * 1e6,
               int_b1_sq,
               W,
               exp(-W * sequence.Trf * sequence.SPS));

    QI::Log(verbose, "Finished.");
    return EXIT_SUCCESS;
}
