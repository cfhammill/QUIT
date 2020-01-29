#pragma once

#include "SequenceBase.h"
#include "mupa_pulse.h"
#include <unordered_map>

struct MUPASequence : QI::SequenceBase {
    double                                   TR, Tramp, Trf;
    Eigen::ArrayXd                           FA;
    int                                      SPS;
    std::unordered_map<std::string, PrepPulse> prep_pulses;
    std::vector<std::string>                 prep;
    QI_SEQUENCE_DECLARE(MUPA);
    Eigen::Index size() const override { return prep.size(); };
};
void from_json(const json &j, MUPASequence &s);
