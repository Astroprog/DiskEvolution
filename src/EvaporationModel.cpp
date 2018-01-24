#include "EvaporationModel.h"
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>
#include <iostream>
#include <fstream>


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


EvaporationIon100::EvaporationIon100() {}

EvaporationIon100::EvaporationIon100(const std::string filename) {
    std::ifstream instream;
    std::string line;

    instream.open(filename);
    bool first = true;

    while (std::getline(instream, line)) {
        boost::trim(line);
        boost::tokenizer<boost::char_separator<char>> tk(line, boost::char_separator<char>(" "));
        for (auto i(tk.begin()); i != tk.end(); ++i) {
            if (first) {
                radii.push_back(std::stod(*i));
            } else {
                losses.push_back(std::stod(*i));
            }
        }
        first = false;
    }

    double step = radii[1] - radii[0];
    spline = boost::math::cubic_b_spline<double>(losses.begin(), losses.end(), radii[0], step);
}

double EvaporationIon100::getLossAtRadius(const double r, const double data, const double dt) {
    double densityLoss = spline(r);
    if (densityLoss * r * dt >= data) {
        densityLoss = data / (r * dt);
    }
    return densityLoss;
}





