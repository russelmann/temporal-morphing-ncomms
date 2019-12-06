//==========================================================
//
//	PhySim library. Generic library for physical simulation.
//
//	Authors:
//			Jesus Perez Rodriguez, IST Austria
//
//==========================================================

#pragma once

#include <PhySim/CommonIncludes.h>
#include <PhySim/PhySimInterface.h>

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	class Material
	{
	public:
		struct Property
		{
			static const string Thickness;
			static const string CollT;
			static const string CollK;
			static const string Pressure;
			static const string Density;
			static const string Poisson;
			static const string StretchK;
			static const string BendingK;
			static const string TwistK;
			static const string ShearK;
			static const string Young;
			static const string Bulk;
			static const string Lame1;
			static const string Lame2;
			static const string Mooney10;
			static const string Mooney01;
			static const string ShearMod;
		};

	public:

		Material()
		{
			// Nothing to do here...
		}

		virtual ~Material()
		{
			this->m_matProperties.clear();
		}

		inline bool HasProperty(const string& name) const
		{
			return this->m_matProperties.find(name) != this->m_matProperties.end();
		}

		inline void AddProperty(const string& name, const Real& value)
		{
			this->m_matProperties.insert(pair<string, Real>(name, value));

			// Set actual value
			(*this)[name] = value;
		}

		inline void ClearProperties()
		{
			this->m_matProperties.clear();
		}

		inline Real& operator[](const string& name)
		{
			if (!this->HasProperty(name))
            assert(false); ;
			return this->m_matProperties.at(name);
		}

		inline const Real& operator[](const string& name) const
		{
			if (!this->HasProperty(name))
            assert(false); ;
			return this->m_matProperties.at(name);
		}

		void InitFromDensity(Real density);
		void InitRealisticFromYoungPoisson(Real young, Real poisson, Real density);
		void InitRealisticFromLameParameter(Real lame1, Real lame2, Real density);

	protected:

		map<string, Real> m_matProperties;

	};

}
