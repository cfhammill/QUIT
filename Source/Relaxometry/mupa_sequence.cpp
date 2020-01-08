#include "mupa_sequence.h"

void from_json(const json &j, MUPASequence &s) {
    j.at("TR").get_to(s.TR);
    j.at("Trf").get_to(s.Trf);
    j.at("Tramp").get_to(s.Tramp);
    s.FA = j.at("FA").get<double>() * M_PI / 180.0;
    j.at("SPS").get_to(s.SPS);
    j.at("prep_pulses").get_to(s.prep_pulses);
    j.at("prep").get_to(s.prep);
}