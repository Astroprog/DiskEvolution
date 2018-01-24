#pragma once

#include <string>
#include <vector>
#include <cmath>
#include <boost/math/interpolators/cubic_b_spline.hpp>

class EvaporationModel
{
public:
	virtual ~EvaporationModel() {}
	virtual double getLossAtRadius(const double r, const double data, const double dt) = 0;
};


class EvaporationClarke : public EvaporationModel
{
private:
    const double au = 1.495978707e13;
    const double G = 6.67408e-8;
    const double kb = 1.38064852e-16;
    const double mp = 1.672621898e-24;
    const double T0 = 280;
    const double normLuminosity = 1e41;

    double luminosity;
    double photoRadius;

public:
	EvaporationClarke();
	EvaporationClarke(const double lum, const double r);

	void setLuminosity(const double lum);
	void setPhotoRadius(const double r);

	virtual double getLossAtRadius(const double r, const double data, const double dt);
};

class EvaporationIon100 : public EvaporationModel
{
private:
	std::vector<double> radii;
	std::vector<double> losses;
	boost::math::cubic_b_spline<double> spline;

public:
	EvaporationIon100();
	EvaporationIon100(const std::string filename);

	virtual double getLossAtRadius(const double r, const double data, const double dt);
};

