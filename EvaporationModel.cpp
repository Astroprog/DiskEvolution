#include "EvaporationModel.h"


EvaporationClarke::EvaporationClarke() {
	luminosity = normLuminosity;
	photoRadius = 1e14;
}

EvaporationClarke::EvaporationClarke(const double lum, const double r) {
	luminosity = lum;
	photoRadius = r;
}

double EvaporationClarke::getLossAtRadius(const double r, const double data, const double dt) {
	if (r * au >= photoRadius) {
        double soundspeed = sqrt(kb * 1e4 / (2.3 * mp));
        double numberDensity = 5.7e4 * sqrt(luminosity / normLuminosity) * pow(photoRadius * 1e-14, -1.5) * pow(r * au / photoRadius, -2.5);
        double densityLoss = 2 * soundspeed * numberDensity * mp;
        if (densityLoss * r * dt >= data) {
            densityLoss = data / (r * dt);
        }
        return densityLoss;
    } else {
        return 0.0;
    }
}
