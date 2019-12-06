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

	class DoFSet : public IDoFSet
	{
	protected:

		int m_numDim;
		int m_idFull;
		int m_idFree;
		int m_offsetFull;
		int m_offsetFree;
		bool m_fixed;

		//Real mass;
		vector<VectorXd> m_vpos;
		vector<VectorXd> m_vvel;

	public:

		inline DoFSet(int numDim)
		{
			this->m_numDim = numDim;
			m_vpos.resize(NUMSPACES);
			m_vvel.resize(NUMSPACES);
			m_offsetFull = -1;
			m_offsetFree = -1;
			m_idFull = -1;
			m_idFree = -1;
			m_fixed = false;
			m_vpos[0].setZero();
			m_vpos[1].setZero();
			m_vvel[0].setZero();
			m_vvel[1].setZero();
		}

		inline virtual ~DoFSet() {}

		inline int GetNumDim(void) const { return m_numDim; }

		inline int GetId_Full(void) const { return m_idFull; }
		inline void SetId_Full(int _id) { m_idFull = _id; }

		inline int GetId_Free(void) const { return m_idFree; }
		inline void SetId_Free(int _id) { m_idFree = _id; }

		inline int GetOffset_Full() const { return m_offsetFull; }
		inline void SetOffset_Full(int os) { m_offsetFull = os; }

		inline int GetOffset_Free() const { return m_offsetFree; }
		inline void SetOffset_Free(int os) { m_offsetFree = os; }

		inline bool IsFixed(void) const { return m_fixed; }

		inline void Fix(void) { m_fixed = true; }
		inline void Unfix(void) { m_fixed = false; }

		inline const VectorXd& GetPosition(Space s = Space::DEF) const { return m_vpos[s]; }
		inline void SetPosition(const VectorXd& x, Space s = Space::DEF) { m_vpos[s] = x; }

		inline const VectorXd& GetVelocity(Space s = Space::DEF) const { return m_vvel[s]; }
		inline void SetVelocity(const VectorXd& x, Space s = Space::DEF) { m_vpos[s] = x; }

		inline const VectorXd& GetVelocity_x() const { return m_vvel[Space::DEF]; }
		inline void SetVelocity_x(const VectorXd& v) { m_vvel[Space::DEF] = v; }

		inline const VectorXd& GetVelocity_0() const { return m_vvel[Space::MAT]; }
		inline void SetVelocity_0(const VectorXd& v) { m_vvel[Space::MAT] = v; }

		inline const VectorXd& GetPosition_x() const { return m_vpos[Space::DEF]; }
		inline void SetPosition_x(const VectorXd& x) { m_vpos[Space::DEF] = x; }

		inline const VectorXd& GetPosition_0() const { return m_vpos[Space::MAT]; }
		inline void SetPosition_0(const VectorXd& x) { m_vpos[Space::MAT] = x; }

		virtual void SetPositions(const VectorXd& vpos)
		{
			for (int i = 0; i < NUMSPACES; ++i)
				this->m_vpos[i] = vpos;
		}

		virtual void SetVelocities(const VectorXd& vvel)
		{
			for (int i = 0; i < NUMSPACES; ++i)
				this->m_vvel[i] = vvel;
		}

	};

}
