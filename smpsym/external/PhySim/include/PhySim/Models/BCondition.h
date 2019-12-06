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

	enum BCType
	{
		Fixed,
		Force,
		Gravity
	};

	struct BCSetup
	{
		int						m_maxStep;
		int						m_incStep;
		Real					m_maxError;
		BCType					m_type;
		vector<DoFSet*>			m_vDoF;
		vector<VectorXd>		m_vini;
		VectorXd				m_vendT;
		VectorXd				m_vendR;


		BCSetup() 
		{
			this->m_maxStep = 1;
			this->m_incStep = 1;
			this->m_vDoF.clear();
			this->m_vini.clear();
			this->m_vendR.setZero(3);
			this->m_vendT.setZero(3);
			this->m_type = BCType::Fixed;
		}
	};

	class BCondition
	{
	public:
		int					m_curStep;
		BCSetup				m_setup;

		Model*				m_pModel;
		IEnergyElement*		m_pEle;

	public:
		BCondition(Model* pModel, const BCSetup& bcStup)
		{
			assert(pModel != NULL);

			this->m_pModel = pModel;
			this->m_setup = bcStup;
			
			this->m_curStep = 0;
		}

		virtual ~BCondition(void)
		{
			// Nothing to do here...
		}

		inline IEnergyElement*& Element() { return this->m_pEle; }

		inline virtual BCSetup& Setup() { return this->m_setup; }
	
		inline virtual int& CurStep() { return this->m_curStep; }
		inline virtual int& MaxStep() { return m_setup.m_maxStep; }

		virtual bool IsFullyLoaded() const { return this->m_curStep == this->m_setup.m_maxStep; }

		virtual bool StepLoading();

		virtual bool ResetLoading();

		virtual bool FullLoading();

		virtual void GetIniValues(vector<VectorXd>& vval)
		{
			vval.resize(this->m_setup.m_vini.size());
			for (size_t i = 0; i < vval.size(); ++i)
				vval[i] = this->m_setup.m_vini[i];
		}

		virtual void GetCurValues(vector<VectorXd>& vval)
		{
			Real a = ((Real)m_curStep / (Real)m_setup.m_maxStep);

			this->GetValues(a, vval);
		}

		virtual void GetEndValues(vector<VectorXd>& vval)
		{
			this->GetValues(1, vval);
		}

		virtual void GetValues(Real alpha, vector<VectorXd>& vval)
		{
			vval.resize(this->m_setup.m_vini.size());

			VectorXd vcurT = alpha*this->m_setup.m_vendT;
			VectorXd vcurR = alpha*this->m_setup.m_vendR;

			for (size_t i = 0; i < vval.size(); ++i)
			{
				Vector3d vposIni, vposRot;
				vposIni = m_setup.m_vini[i];
				eulerRotation(vcurR.data(), vposIni.data(), vposRot.data());

				vval[i] = vcurT + vposRot;
			}
		}

		virtual void Init() = 0;
		virtual void Apply() = 0;

	};

}
