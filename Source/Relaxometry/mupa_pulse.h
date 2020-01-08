#pragma once

struct RFPulse {
    Eigen::ArrayXd B1x, B1y, timestep;
};
void from_json(const json &j, RFPulse &s);
