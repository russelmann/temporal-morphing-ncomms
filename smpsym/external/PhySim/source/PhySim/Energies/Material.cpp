#pragma once


#include <PhySim/Energies/Material.h>

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	const string Material::Property::CollT = "COLLT";
	const string Material::Property::CollK = "COLLK";
	const string Material::Property::Pressure = "PRESSURE";
	const string Material::Property::Density = "DENSITY";
	const string Material::Property::Poisson = "POISSON";
	const string Material::Property::StretchK = "STRETCH";
	const string Material::Property::BendingK = "BENDING";
	const string Material::Property::ShearK = "SHEAR";
	const string Material::Property::TwistK = "TWIST";
	const string Material::Property::Young = "YOUNG";
	const string Material::Property::Bulk = "BULK";
	const string Material::Property::Lame1 = "LAME1";
	const string Material::Property::Lame2 = "LAME2";
	const string Material::Property::ShearMod = "LAME2";
	const string Material::Property::Mooney01 = "MR01";
	const string Material::Property::Mooney10 = "MR10";
	const string Material::Property::Thickness = "THICK";

	void Material::InitFromDensity(Real density)
	{
		this->AddProperty(Property::Density, density);
	}

	void Material::InitRealisticFromYoungPoisson(Real young, Real poisson, Real density)
	{
		this->AddProperty(Property::Young, young);
		this->AddProperty(Property::Poisson, poisson);
		this->AddProperty(Property::Density, density);

		Real t0 = 1 + poisson;
		Real t1 = 1 - 2 * poisson;
		Real lame1 = (young*poisson) / (t0*t1);
		Real lame2 = young / (2 * t0);
		Real bulk = young / (3 * t1);

		this->AddProperty(Property::Lame1, lame1);
		this->AddProperty(Property::Lame2, lame2);
		this->AddProperty(Property::Bulk, bulk);
	}

	void Material::InitRealisticFromLameParameter(Real lame1, Real lame2, Real density)
	{
		this->AddProperty(Property::Lame1, lame1);
		this->AddProperty(Property::Lame2, lame2);
		this->AddProperty(Property::Density, density);

		Real l = lame1;
		Real G = lame2;
		Real bulk = l + (2 / 3)*G;
		Real poisson = l / (2 * (l + G));
		Real young = G*(3 * l + 2 * G) / (l + G);

		this->AddProperty(Property::Young, young);
		this->AddProperty(Property::Poisson, poisson);
		this->AddProperty(Property::Bulk, bulk);
	}
}

