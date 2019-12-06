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

	class NodeEBD 
	{
		bool				m_valid;
		VectorXd			m_visopar;
		Polytope*			m_pMaster;

	public:
		NodeEBD();
		NodeEBD(const NodeEBD& toCopy);
		NodeEBD(Polytope* pMaster, const VectorXd& vpar);
		virtual ~NodeEBD(void);

		inline bool Valid() const { return this->m_valid; }
		inline const Polytope* Master() const { return this->m_pMaster; }
		inline const VectorXd Parameters() const { return this->m_visopar; }

		Vector3d GetPosition(Space s = Space::DEF) const;

	};
}
