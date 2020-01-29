#include "mupa_pulse.h"
#include "JSON.h"
#include "Log.h"

void from_json(const json &j, RFPulse &p) {
    p.B1x      = QI::ArrayFromJSON(j, "B1x");
    p.B1y      = QI::ArrayFromJSON(j, "B1y");
    p.timestep = QI::ArrayFromJSON(j, "timestep", 1.e-6);

    if (p.B1x.rows() != p.B1y.rows()) {
        QI::Fail("B1x and B1y array lengths did not match {} vs {}", p.B1x.rows(), p.B1y.rows());
    }
    if (p.timestep.rows() != p.B1x.rows()) {
        QI::Fail(
            "B1 timestep array lengths did not match {} vs {}", p.timestep.rows(), p.B1x.rows());
    }
}

void to_json(json &j, RFPulse const &p) {
    j = json{{"B1x", p.B1x}, {"B1y", p.B1y}, {"timestep", p.timestep * 1.e6}};
}

void from_json(const json &j, PrepPulse &p) {
    p.FA = j.at("FA").get<double>() * M_PI / 180.0;
    j.at("T_l").get_to(p.T_l);
    j.at("T_t").get_to(p.T_t);
    j.at("T_act").get_to(p.T_act);
    j.at("B1_sq_mean").get_to(p.B1_sq_mean);
}

void to_json(json &j, PrepPulse const &p) {
    j = json{{"FA", p.FA * 180 / M_PI},
             {"T_l", p.T_l},
             {"T_t", p.T_t},
             {"T_act", p.T_act},
             {"B1_sq_mean", p.B1_sq_mean}};
}