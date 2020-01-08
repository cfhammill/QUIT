#include "JSON.h"
#include "Log.h"
#include "mupa_pulse.h"

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