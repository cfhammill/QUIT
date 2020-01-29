#pragma once

struct RFPulse {
    Eigen::ArrayXd B1x, B1y, timestep;
};
void from_json(const json &j, RFPulse &s);

struct PrepPulse {
    double FA, T_l, T_t, B1_sq_mean;
};
void from_json(const json &j, PrepPulse &p);
void to_json(json &j, const PrepPulse &s);