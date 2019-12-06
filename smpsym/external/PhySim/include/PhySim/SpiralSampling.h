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

#include <PhySim/Geometry/DoFSet.h>
#include <PhySim/Models/Model_FEM.h>
#include <PhySim/Models/BC_Fixed.h>
#include <PhySim/Models/BC_Force.h>
#include <PhySim/Models/BC_Gravity.h>
#include <PhySim/Solvers/OptimSolver_USQP_LS.h>
#include <PhySim/Solvers/OptimSolver_KnitroBridge.h>
#include <PhySim/Solvers/OptimProblem_BasicStatic.h>
#include <PhySim/Geometry/TetGenReader.h>

#include <PhySim/Energies/EnergyElement_FEM.h>

#include <PhySim/Utils/json.hpp>

#include <fstream>
#include <streambuf>

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	using Json = nlohmann::json;

	class SpiralSampling 
	{
	public:
		struct Sample
		{
			string				m_type;
			float				m_stretchX;
			float				m_stretchY;
			float				m_shearX;
			float				m_shearY;
			float				m_bendingInX;
			float				m_bendingInY;
			float				m_bendingOutX1;
			float				m_bendingOutY1;
			float				m_bendingOutX2;
			float				m_bendingOutY2;
			float				m_twistX1;
			float				m_twistY1;
			float				m_twistX2;
			float				m_twistY2;
			Sample()
			{
				m_type = "RANDOM";
				m_stretchX = 0;
				m_stretchY = 0;
				m_shearX = 0;
				m_shearY = 0;
				m_bendingInX = 0;
				m_bendingInY = 0;
				m_bendingOutX1 = 0;
				m_bendingOutY1 = 0;
				m_bendingOutX2 = 0;
				m_bendingOutY2 = 0;
				m_twistX1 = 0;
				m_twistY1 = 0;
				m_twistX2 = 0;
				m_twistY2 = 0;
			}
		};

	private:

		// Current file
		
		int											m_countSpiral;
		shared_ptr<IOptimProblem>					m_pProblem;
		shared_ptr<IOptimSolver>					m_pSolver;
		shared_ptr<Model_FEM>						m_pModel;
		MatrixXi									m_mSurface;
		MatrixXd									m_mVertices;
		vector<BCondition*>							m_vpBC;
		vector<NodeEBD>								m_vTP;
		Vector3d									m_vBBmax;
		Vector3d									m_vBBmin;
		Vector3d									m_vBBran;
		vector<Frame3d>								m_vF0;
		vector<Frame3d>								m_vFx;
		MatrixXd									m_mV;
		MatrixXi									m_mE;
		MatrixXi									m_mT;
		vector<MatrixXi>							m_vmF;

		// Current sample

		int					m_countSample;
		vector<Sample>		m_vsampleList;
		vector<string>		m_vspiralList;
		int					m_numSamples;
		int					m_numSpirals;

		// Sampling setup

		Json m_jsetup;
		string m_jinputFolder;
		string m_joutputFolder;
		string m_jpreSolveFolder;
		int m_jiniSpiralIndex;
		int m_jiniSampleIndex;
		int m_jendSpiralIndex;
		int m_jendSampleIndex;
		double m_jyoungModulus;
		double m_jpoissonRatio;
		double m_jdensity;
		int m_jsampStretchNum;
		int m_jsampShearNum;
		int m_jsampBendingInNum;
		int m_jsampBendingOutNum;
		int m_jsampTwistNum;
		int m_jsampRandomNum00;
		int m_jsampRandomNum20;
		int m_jsampRandomNum40;
		int m_jsampRandomNum60;
		int m_jsampRandomNum80;
		double m_jsampStretchMin;
		double m_jsampStretchMax;
		double m_jsampShearMin;
		double m_jsampShearMax;
		double m_jsampBendingInMin;
		double m_jsampBendingInMax;
		double m_jsampBendingOutMin;
		double m_jsampBendingOutMax;
		double m_jsampTwistMin;
		double m_jsampTwistMax;
		double m_jsampRandomMod;
		bool m_jquadratic;
		bool m_jpreSolved;
		int m_jmaxIters;
		double m_jmaxError;
		double m_jmaxStep;
		int m_jmaxBCStep;
		int m_jincBCStep;
		double m_jmaxBCError;

	public:
		SpiralSampling()
		{
			this->m_pModel.reset();
			this->m_pProblem.reset();
			this->m_pSolver.reset();
			m_vpBC.resize(5);
			for (int i = 0; i < 5; ++i)
				m_vpBC[i] = NULL; // Init.
		}

		virtual ~SpiralSampling(void)
		{
			// Nothing to do here...
		}

		shared_ptr<IOptimSolver> Solver() { return this->m_pSolver; }
		shared_ptr<IOptimProblem> Problem() { return this->m_pProblem; }
		shared_ptr<Model_FEM> Model() { return this->m_pModel; }
		const MatrixXd& Vertices() { return this->m_mVertices; }
		const MatrixXi& Surface() { return this->m_mSurface; }

		vector<BCondition*>& BC() { return this->m_vpBC; }
		vector<NodeEBD>& TP() { return this->m_vTP; }

		virtual float& StretchX() { return m_vsampleList[m_countSample].m_stretchX; }
		virtual float& StretchY() { return m_vsampleList[m_countSample].m_stretchY; }
		virtual float& ShearX() { return m_vsampleList[m_countSample].m_shearX; }
		virtual float& ShearY() { return m_vsampleList[m_countSample].m_shearY; }
		virtual float& BendingInX() { return m_vsampleList[m_countSample].m_bendingInX; }
		virtual float& BendingInY() { return m_vsampleList[m_countSample].m_bendingInY; }
		virtual float& BendingOutX1() { return m_vsampleList[m_countSample].m_bendingOutX1; }
		virtual float& BendingOutY1() { return m_vsampleList[m_countSample].m_bendingOutY1; }
		virtual float& BendingOutX2() { return m_vsampleList[m_countSample].m_bendingOutX2; }
		virtual float& BendingOutY2() { return m_vsampleList[m_countSample].m_bendingOutY2; }
		virtual float& TwistX1() { return m_vsampleList[m_countSample].m_twistX1; }
		virtual float& TwistY1() { return m_vsampleList[m_countSample].m_twistY1; }
		virtual float& TwistX2() { return m_vsampleList[m_countSample].m_twistX2; }
		virtual float& TwistY2() { return m_vsampleList[m_countSample].m_twistY2; }

		virtual Vector3d& BBMin() { return this->m_vBBmin; }
		virtual Vector3d& BBMax() { return this->m_vBBmax; }

		virtual bool InitSampling(string& setupFile)
		{
			string setupFileContent;

			ifstream setupFileIn(setupFile);
			setupFileIn.seekg(0, ios::end);

			setupFileContent.reserve(setupFileIn.tellg());

			setupFileIn.seekg(0, ios::beg);

			setupFileContent.assign((istreambuf_iterator<char>(setupFileIn)), istreambuf_iterator<char>());

			// Read setup JSON

			m_jsetup = Json::parse(setupFileContent);
			
			// Files

			m_jinputFolder = m_jsetup["folderInput"].get<string>();
			m_joutputFolder = m_jsetup["folderOutput"].get<string>();
			m_jpreSolveFolder = m_jsetup["folderPreSolve"].get<string>();
			m_jiniSpiralIndex = m_jsetup["iniSpiralIndex"].get<int>();
			m_jendSpiralIndex = m_jsetup["endSpiralIndex"].get<int>();
			m_jiniSampleIndex = m_jsetup["iniSampleIndex"].get<int>();
			m_jendSampleIndex = m_jsetup["endSampleIndex"].get<int>();

			// Material

			m_jyoungModulus = m_jsetup["material"]["youngModulus"].get<double>();
			m_jpoissonRatio = m_jsetup["material"]["poissonRatio"].get<double>();
			m_jdensity = m_jsetup["material"]["density"].get<double>();

			// Sampling

			m_jsampStretchNum = m_jsetup["samples"]["sampStretchNum"].get<int>();
			m_jsampShearNum = m_jsetup["samples"]["sampShearNum"].get<int>();
			m_jsampBendingInNum = m_jsetup["samples"]["sampBendingInNum"].get<int>();
			m_jsampBendingOutNum = m_jsetup["samples"]["sampBendingOutNum"].get<int>();
			m_jsampTwistNum = m_jsetup["samples"]["sampTwistNum"].get<int>();
			m_jsampRandomNum20 = m_jsetup["samples"]["sampRandomNum20"].get<int>();
			m_jsampRandomNum40 = m_jsetup["samples"]["sampRandomNum40"].get<int>();
			m_jsampRandomNum60 = m_jsetup["samples"]["sampRandomNum60"].get<int>();
			m_jsampRandomNum80 = m_jsetup["samples"]["sampRandomNum80"].get<int>();
			m_jsampRandomNum00 = m_jsetup["samples"]["sampRandomNum00"].get<int>();

			m_jsampStretchMin = m_jsetup["ranges"]["sampStretchMin"].get<double>();
			m_jsampStretchMax = m_jsetup["ranges"]["sampStretchMax"].get<double>();
			m_jsampShearMin = m_jsetup["ranges"]["sampShearMin"].get<double>();
			m_jsampShearMax = m_jsetup["ranges"]["sampShearMax"].get<double>();
			m_jsampBendingInMin = m_jsetup["ranges"]["sampBendingInMin"].get<double>();
			m_jsampBendingInMax = m_jsetup["ranges"]["sampBendingInMax"].get<double>();
			m_jsampBendingOutMin = m_jsetup["ranges"]["sampBendingOutMin"].get<double>();
			m_jsampBendingOutMax = m_jsetup["ranges"]["sampBendingOutMax"].get<double>();
			m_jsampTwistMin = m_jsetup["ranges"]["sampTwistMin"].get<double>();
			m_jsampTwistMax = m_jsetup["ranges"]["sampTwistMax"].get<double>();

			// Options

			m_jquadratic = (bool) m_jsetup["options"]["quadratic"].get<int>();
			m_jpreSolved = (bool) m_jsetup["options"]["preSolved"].get<int>();
			m_jmaxIters = m_jsetup["options"]["maxIters"].get<int>();
			m_jmaxError = m_jsetup["options"]["maxError"].get<double>();
			m_jmaxStep = m_jsetup["options"]["maxStep"].get<double>();
			m_jmaxBCStep = m_jsetup["options"]["maxBCStep"].get<int>();
			m_jincBCStep = m_jsetup["options"]["incBCStep"].get<int>();
			m_jmaxBCError = m_jsetup["options"]["maxBCError"].get<double>();

			// List input spirals

			vector<string> vfileList;

			if (!listDirectoryFiles(m_jinputFolder, vfileList))
				return false;

			// One file per-spiral

			m_vspiralList.clear();

			for (int i = 0; i < (int) vfileList.size(); ++i)
			{
				vector<string> vfileParts;
				splitString(vfileList[i], ".0.stl", vfileParts);
				if (vfileParts.size() == 2)
					this->m_vspiralList.push_back(vfileParts[0]);
			}

			m_numSpirals = (int) m_vspiralList.size();

			// Initialize samples

			//m_numSamples =
			//	(m_jsampShearNum != 0) ? pow(m_jsampShearNum, 2) : 0 + // 2^k shear samples
			//	(m_jsampStretchNum != 0) ? pow(m_jsampStretchNum, 2) : 0 + // 2^k stretch samples
			//	(m_jsampBendingInNum != 0) ? pow(m_jsampBendingInNum, 2) : 0 + // 2^k bending in samples
			//	(m_jsampBendingOutNum != 0) ? pow(m_jsampBendingOutNum, 4) : 0 + // 4^k bending out samples
			//	(m_jsampTwistNum != 0) ? pow(m_jsampTwistNum, 4) : 0 + // 4^k twist samples
			//	m_jsampRandomNum00;

			this->m_vsampleList.clear();

			// Add stretch

			if (m_jsampStretchNum > 0)
			{
				Real ranStretch = (m_jsampStretchMax - m_jsampStretchMin);
				Real incStretch = ranStretch / (m_jsampStretchNum - 1);
				for (int i = 0; i < m_jsampStretchNum; ++i)
				{
					for (int j = 0; j <= i; ++j)
					{
						Sample s;
						s.m_stretchX = m_jsampStretchMin + incStretch*i;
						s.m_stretchY = m_jsampStretchMin + incStretch*j;
						s.m_type = "STRETCH";
						this->m_vsampleList.push_back(s);
					}
				}
			}

			if (m_jsampShearNum > 0)
			{
				Real ranShear = (m_jsampShearMax - m_jsampShearMin);
				Real incShear = ranShear / (m_jsampShearNum - 1);
				for (int i = 0; i < m_jsampShearNum; ++i)
				{
					for (int j = 0; j <= i; ++j)
					{
						Sample s;
						s.m_shearX = m_jsampShearMin + incShear*i;
						s.m_shearY = m_jsampShearMin + incShear*j;
						s.m_type = "SHEAR";
						this->m_vsampleList.push_back(s);
					}
				}
			}

			if (m_jsampBendingInNum > 0)
			{
				Real ranBend = (m_jsampBendingInMax - m_jsampBendingInMin);
				Real incBend = ranBend / (m_jsampBendingInNum - 1);
				for (int i = 0; i < m_jsampBendingInNum; ++i)
				{
					for (int j = 0; j <= i; ++j)
					{
						Sample s;
						s.m_bendingInX = m_jsampBendingInMin + incBend*i;
						s.m_bendingInY = m_jsampBendingInMin + incBend*j;
						s.m_type = "BENDINGIN";
						this->m_vsampleList.push_back(s);
					}
				}
			}

			if (m_jsampBendingOutNum > 0)
			{
				Real ranBend = (m_jsampBendingOutMax - m_jsampBendingOutMin);
				Real incBend = ranBend / (m_jsampBendingOutNum - 1);
				for (int i = 0; i < m_jsampBendingOutNum; ++i)
				{
					for (int j = 0; j <= i; ++j)
					{
						for (int k = 0; k <= j; ++k)
						{
							for (int l = 0; l <= k; ++l)
							{
								Sample s;
								s.m_bendingOutX1 = m_jsampBendingOutMin + incBend*i;
								s.m_bendingOutX2 = m_jsampBendingOutMin + incBend*j;
								s.m_bendingOutY1 = m_jsampBendingOutMin + incBend*k;
								s.m_bendingOutY2 = m_jsampBendingOutMin + incBend*l;
								s.m_type = "BENDINGOUT";
								this->m_vsampleList.push_back(s);
							}
						}
					}
				}
			}

			if (m_jsampTwistNum > 0)
			{
				Real ranTwist = (m_jsampTwistMax - m_jsampTwistMin);
				Real incTwist = ranTwist / (m_jsampTwistNum - 1);
				for (int i = 0; i < m_jsampTwistNum; ++i)
				{
					for (int j = 0; j <= i; ++j)
					{
						for (int k = 0; k <= j; ++k)
						{
							for (int l = 0; l <= k; ++l)
							{
								Sample s;
								s.m_twistX1 = m_jsampTwistMin + incTwist*i;
								s.m_twistX2 = m_jsampTwistMin + incTwist*j;
								s.m_twistY1 = m_jsampTwistMin + incTwist*k;
								s.m_twistY2 = m_jsampTwistMin + incTwist*l;
								s.m_type = "TWIST";
								this->m_vsampleList.push_back(s);
							}
						}
					}
				}
			}

			srand(23654481256987);
			this->AddRandomSamples(m_jsampRandomNum20, 0.2);
			this->AddRandomSamples(m_jsampRandomNum40, 0.4);
			this->AddRandomSamples(m_jsampRandomNum60, 0.6);
			this->AddRandomSamples(m_jsampRandomNum80, 0.8);
			this->AddRandomSamples(m_jsampRandomNum00, 1.0);

			this->m_numSamples = (int) this->m_vsampleList.size();

			// Sampling log

			ostringstream logFile;
			logFile << this->m_joutputFolder << "/samplingLog.txt";
			ofstream logOut(logFile.str().c_str(), fstream::trunc);

			ostringstream fiLogFile;
			ostringstream siLogFile;
			ostringstream tiLogFile;
			fiLogFile << this->m_joutputFolder << "/firstInvariant.csv";
			siLogFile << this->m_joutputFolder << "/secondInvariant.csv";
			tiLogFile << this->m_joutputFolder << "/thirdInvariant.csv";
			ofstream filogOut(fiLogFile.str().c_str(), fstream::trunc);
			ofstream silogOut(siLogFile.str().c_str(), fstream::trunc);
			ofstream tilogOut(tiLogFile.str().c_str(), fstream::trunc);

			// Start first spiral

			this->m_countSpiral = m_jiniSpiralIndex - 1;
			this->m_countSample = m_jiniSampleIndex - 1;

			if (this->m_jendSpiralIndex == -1) this->m_jendSpiralIndex = this->m_numSpirals - 1;
			if (this->m_jendSampleIndex == -1) this->m_jendSampleIndex = this->m_numSamples - 1;

			this->NextSpiral();
			this->NextSample();

			// Create output log file

			InitSpiral();

			return true;		
		}

		virtual void AddRandomSamples(int num, Real randomMod)
		{
			Sample s;

			for (int i = 0; i < num; ++i)
			{
				s.m_stretchX = (m_jsampStretchMin + ((float)rand() / (float)RAND_MAX)*(m_jsampStretchMax - m_jsampStretchMin))*randomMod;
				s.m_stretchY = (m_jsampStretchMin + ((float)rand() / (float)RAND_MAX)*(m_jsampStretchMax - m_jsampStretchMin))*randomMod;
				s.m_shearX = (m_jsampShearMin + ((float)rand() / (float)RAND_MAX)*(m_jsampShearMax - m_jsampShearMin))*randomMod;
				s.m_shearY = (m_jsampShearMin + ((float)rand() / (float)RAND_MAX)*(m_jsampShearMax - m_jsampShearMin))*randomMod;
				s.m_bendingInX = (m_jsampBendingInMin + ((float)rand() / (float)RAND_MAX)*(m_jsampBendingInMax - m_jsampBendingInMin))*randomMod;
				s.m_bendingInY = (m_jsampBendingInMin + ((float)rand() / (float)RAND_MAX)*(m_jsampBendingInMax - m_jsampBendingInMin))*randomMod;
				s.m_bendingOutX1 = (m_jsampBendingOutMin + ((float)rand() / (float)RAND_MAX)*(m_jsampBendingOutMax - m_jsampBendingOutMin))*randomMod;
				s.m_bendingOutX2 = (m_jsampBendingOutMin + ((float)rand() / (float)RAND_MAX)*(m_jsampBendingOutMax - m_jsampBendingOutMin))*randomMod;
				s.m_bendingOutY1 = (m_jsampBendingOutMin + ((float)rand() / (float)RAND_MAX)*(m_jsampBendingOutMax - m_jsampBendingOutMin))*randomMod;
				s.m_bendingOutY2 = (m_jsampBendingOutMin + ((float)rand() / (float)RAND_MAX)*(m_jsampBendingOutMax - m_jsampBendingOutMin))*randomMod;
				s.m_twistX1 = (m_jsampTwistMin + ((float)rand() / (float)RAND_MAX)*(m_jsampTwistMax - m_jsampTwistMin))*randomMod;
				s.m_twistX2 = (m_jsampTwistMin + ((float)rand() / (float)RAND_MAX)*(m_jsampTwistMax - m_jsampTwistMin))*randomMod;
				s.m_twistY1 = (m_jsampTwistMin + ((float)rand() / (float)RAND_MAX)*(m_jsampTwistMax - m_jsampTwistMin))*randomMod;
				s.m_twistY2 = (m_jsampTwistMin + ((float)rand() / (float)RAND_MAX)*(m_jsampTwistMax - m_jsampTwistMin))*randomMod;
				s.m_type = "RANDOM";
				this->m_vsampleList.push_back(s);
			}
		}

		virtual bool InitSpiral()
		{
			int idx = m_countSpiral;

			vector<string> vdescriptor;
			splitString(m_vspiralList[idx], "_", vdescriptor);
			float thick = stof(vdescriptor[1]);
			float scale = stof(vdescriptor[3]);
			float twist = stof(vdescriptor[5]);

			// Load tetrahedral model

			ostringstream meshInputFile;
			meshInputFile << m_jinputFolder << "/";
			meshInputFile << m_vspiralList[idx];

			if (!TetGenReader::readPoly(meshInputFile.str(), m_mV, m_mT, m_mE, m_vmF))
			{
				return false;
			}

			// Create FEM model

			m_mVertices = m_mV;
			m_mSurface = m_vmF[0];

			this->m_pModel = shared_ptr<Model_FEM>(new Model_FEM());

			this->m_pModel->SetPresolveFixed(!this->m_jpreSolved);
			this->m_pModel->GetOptions().m_mNodes = m_mV;
			this->m_pModel->GetOptions().m_mElems = m_mT;

			if (!this->m_jquadratic)
			{
				logSimu("\n[TRACE] Creating basic linear tetrahedron simulation");

				this->m_pModel->GetOptions().m_numQuadrature = 1;
				this->m_pModel->GetOptions().m_discretization = PhySim::Discretization::Tetrahedra4;
			}
			else
			{
				logSimu("\n[TRACE] Creating basic quadratic tetrahedron simulation");

				this->m_pModel->GetOptions().m_numQuadrature = 8;
				this->m_pModel->GetOptions().m_discretization = PhySim::Discretization::Tetrahedra10;
			}
			this->m_pModel->GetOptions().m_materialModel = PhySim::MaterialModel::CoNH;
			this->m_pModel->GetOptions().m_material.InitRealisticFromYoungPoisson(m_jyoungModulus, 
																				  m_jpoissonRatio, 
																				  m_jdensity);
			this->m_pModel->Init();

			// Create boundary conditions

			for (int i = 0; i < 4; ++i)
				if (m_vpBC[i] != NULL)
					delete m_vpBC[i];

			computeBoundingBox(m_mV, m_vBBmin, m_vBBmax);
			m_vBBran = m_vBBmax- m_vBBmin;

			Vector3d vminBC0(m_vBBmin.x() - 1e-6, -2 * thick, m_vBBmin.z() - 1e-6);
			Vector3d vmaxBC0(m_vBBmin.x() + 1e-6, +2 * thick, m_vBBmax.z() + 1e-6);
			Vector3d vminBC1(m_vBBmax.x() - 1e-6, -2 * thick, m_vBBmin.z() - 1e-6);
			Vector3d vmaxBC1(m_vBBmax.x() + 1e-6, +2 * thick + 1e-3, m_vBBmax.z() + 1e-6);
			Vector3d vminBC2(-2 * thick, m_vBBmax.y() - 1e-6, m_vBBmin.z() - 1e-6);
			Vector3d vmaxBC2(+2 * thick, m_vBBmax.y() + 1e-6, m_vBBmax.z() + 1e-6);
			Vector3d vminBC3(-2 * thick, m_vBBmin.y() - 1e-6, m_vBBmin.z() - 1e-6);
			Vector3d vmaxBC3(+2 * thick, m_vBBmin.y() + 1e-6, m_vBBmax.z() + 1e-6);
			m_vpBC[0] = this->CreateBCFromSelBoxBounds(vminBC0, vmaxBC0);
			m_vpBC[1] = this->CreateBCFromSelBoxBounds(vminBC1, vmaxBC1);
			m_vpBC[2] = this->CreateBCFromSelBoxBounds(vminBC2, vmaxBC2);
			m_vpBC[3] = this->CreateBCFromSelBoxBounds(vminBC3, vmaxBC3);
			this->m_pModel->AddBoundaryCondition(m_vpBC[0]);
			this->m_pModel->AddBoundaryCondition(m_vpBC[1]);
			this->m_pModel->AddBoundaryCondition(m_vpBC[2]);
			this->m_pModel->AddBoundaryCondition(m_vpBC[3]);

			BCSetup bcForceSetup;

			bcForceSetup.m_type = BCType::Force;
			bcForceSetup.m_vDoF = m_pModel->GetDoFSets();
			int numDoF = (int)bcForceSetup.m_vDoF.size();

			bcForceSetup.m_maxStep = m_jmaxBCStep;
			bcForceSetup.m_incStep = m_jincBCStep;
			bcForceSetup.m_maxError = m_jmaxBCError;
			bcForceSetup.m_vendT = Vector3d(0.0, 0.0, 0);
			bcForceSetup.m_vendR.setZero();

			bcForceSetup.m_vini.resize(numDoF);
			for (size_t i = 0; i < numDoF; ++i)
				bcForceSetup.m_vini[i] = Vector3d(0.0, 0.0, 0);

			m_vpBC[4] = new BC_Force(m_pModel.get(), bcForceSetup);
			this->m_pModel->AddBoundaryCondition(m_vpBC[4]);

			// Tracked points

			vector<Vector3d> vpos(15);
			vpos[0] = Vector3d(m_vBBmin.x(), 0, 0);
			vpos[1] = Vector3d(m_vBBmin.x(), 0, 1e-6);
			vpos[2] = Vector3d(m_vBBmin.x(), 1e-6, 0);

			vpos[3] = Vector3d(0, m_vBBmin.y(), 0);
			vpos[4] = Vector3d(0, m_vBBmin.y(), 1e-6);
			vpos[5] = Vector3d(-1e-6, m_vBBmin.y(), 0);

			vpos[6] = Vector3d(m_vBBmax.x(), 0, 0);
			vpos[7] = Vector3d(m_vBBmax.x(), 0, 1e-6);
			vpos[8] = Vector3d(m_vBBmax.x(), -1e-6, 0);

			vpos[9] = Vector3d(0, m_vBBmax.y(), 0);
			vpos[10] = Vector3d(0, m_vBBmax.y(), 1e-6);
			vpos[11] = Vector3d(1e-6, m_vBBmax.y(), 0);

			vpos[12] = Vector3d(0, 0, 0);
			vpos[13] = Vector3d(0, 0, 1e-6);
			vpos[14] = Vector3d(0, -1e-6, 0);

			this->m_vTP = this->m_pModel->ComputeEmbedding(vpos);

			m_vF0.resize(5);

			for (size_t i = 0; i < 5; ++i)
			{
				m_vF0[i].nor = (this->m_vTP[3 * i + 1].GetPosition() - this->m_vTP[3 * i + 0].GetPosition()).normalized();
				m_vF0[i].bin = (this->m_vTP[3 * i + 2].GetPosition() - this->m_vTP[3 * i + 0].GetPosition()).normalized();
				m_vF0[i].tan = m_vF0[i].nor.cross(m_vF0[i].bin);
			}

			// Prepare for simulation

			this->m_pModel->PrepareForSimulation();

			// Create optimization problem

			this->m_pProblem = shared_ptr<OptimProblem_BasicStatic>(new PhySim::OptimProblem_BasicStatic(this->m_pModel.get()));

			// Create optimization solver

			PhySim::OptimSolverOptions options = PhySim::OptimSolver_SQP::CreateDefaultOptions();
			options.lineSearchIters = 25;
			options.maxError = m_jmaxError;
			options.maxIters = m_jmaxIters;
			options.maxStepSize = m_jmaxStep;
			options.lsSolverType = LSSolverType::LS_EigenLDLT;

			//if (this->m_jpreSolved)
			//{
			//	options.maxIters = m_jmaxIters *= 10;
			//}

			this->m_pSolver = shared_ptr<OptimSolver_USQP_LS>(new PhySim::OptimSolver_USQP_LS(this->m_pProblem.get(), options));

			//this->m_pSolver = shared_ptr<OptimSolver_KnitroBridge>(new PhySim::OptimSolver_KnitroBridge(this->m_pProblem.get()));


			// Start first sample

			m_countSample = m_jiniSampleIndex;

			if (!this->m_jpreSolved)
			{
				// Output default

				ostringstream directFile;
				ostringstream objectFile;
				directFile << m_joutputFolder << "/" << m_vspiralList[m_countSpiral];
				objectFile << directFile.str() << "/" << "default" << ".obj";

				mkdir(directFile.str().c_str());
				writeTriMesh_Obj(objectFile.str(),
								 this->Vertices(),
								 this->Surface());
			}

			return true;
		}

		bool InitSample()
		{
			int idx = this->m_countSample;

			Sample& S = this->m_vsampleList[idx];

			m_vpBC[0]->Setup().m_vendT = Vector3d(-S.m_stretchX*m_vBBran.x(), -S.m_shearX*m_vBBran.y(), 0);

			Vector3d t0 = Vector3d(S.m_twistX1 / 180 * M_PI, 0, 0);
			Vector3d b0 = Vector3d(0, S.m_bendingOutX1 / 180 * M_PI, S.m_bendingInX / 180 * M_PI);
			Matrix3d mT0 = transformEulerToMatrix(t0)*transformEulerToMatrix(b0);
			m_vpBC[0]->Setup().m_vendR = transformMatrixToEuler(mT0);

			//m_vpBC[0]->Setup().m_vendR = Vector3d(S.m_twistX1 / 180 * M_PI, S.m_bendingOutX1 / 180 * M_PI, S.m_bendingInX / 180 * M_PI);
		
			m_vpBC[1]->Setup().m_vendT = Vector3d(S.m_stretchX*m_vBBran.x(), S.m_shearX*m_vBBran.y(), 0);

			Vector3d t1 = Vector3d(-S.m_twistX2 / 180 * M_PI, 0, 0);
			Vector3d b1 = Vector3d(0, -S.m_bendingOutX2 / 180 * M_PI, S.m_bendingInX / 180 * M_PI);
			Matrix3d mT1 = transformEulerToMatrix(t1)*transformEulerToMatrix(b1);
			m_vpBC[1]->Setup().m_vendR = transformMatrixToEuler(mT1);

			//m_vpBC[1]->Setup().m_vendR = Vector3d(-S.m_twistX2 / 180 * M_PI, -S.m_bendingOutX2 / 180 * M_PI, S.m_bendingInX / 180 * M_PI);

			m_vpBC[2]->Setup().m_vendT = Vector3d(S.m_shearY*m_vBBran.x(), S.m_stretchY*m_vBBran.y(), 0);

			Vector3d t2 = Vector3d(0, S.m_twistY1 / 180 * M_PI, 0);
			Vector3d b2 = Vector3d(S.m_bendingOutY1 / 180 * M_PI, 0, -S.m_bendingInY / 180 * M_PI);
			Matrix3d mT2 = transformEulerToMatrix(t2)*transformEulerToMatrix(b2);
			m_vpBC[2]->Setup().m_vendR = transformMatrixToEuler(mT2);

			//m_vpBC[2]->Setup().m_vendR = Vector3d(S.m_bendingOutY1 / 180 * M_PI, S.m_twistY1 / 180 * M_PI, -S.m_bendingInY / 180 * M_PI);

			m_vpBC[3]->Setup().m_vendT = Vector3d(-S.m_shearY*m_vBBran.x(), -S.m_stretchY*m_vBBran.y(), 0);

			Vector3d t3 = Vector3d(0, -S.m_twistY2 / 180 * M_PI, 0);
			Vector3d b3 = Vector3d(-S.m_bendingOutY2 / 180 * M_PI, 0, -S.m_bendingInY / 180 * M_PI);
			Matrix3d mT3 = transformEulerToMatrix(t3)*transformEulerToMatrix(b3);
			m_vpBC[3]->Setup().m_vendR = transformMatrixToEuler(mT3);

			//m_vpBC[3]->Setup().m_vendR = Vector3d(-S.m_bendingOutY2 / 180 * M_PI, -S.m_twistY2 / 180 * M_PI, -S.m_bendingInY / 180 * M_PI);

			IModel::StateP state = this->m_pModel->CreateState();
			this->m_pModel->GetState(state, Space::MAT);
			this->m_pModel->SetState(state, Space::DEF);

			if (!this->m_jpreSolved)
			{
				this->m_pModel->ResetBoundaryConditions();
			}
			else
			{
				ostringstream presolvedFile;

				presolvedFile << this->m_jpreSolveFolder << "/" << m_vspiralList[this->m_countSpiral] << "/";
				presolvedFile << "sample_" << std::setw(6) << std::setfill('0') << m_countSample << ".csv";

				MatrixXd mVSolved;
				if (readMatrix_CSV(presolvedFile.str(), mVSolved))
				{
					m_pModel->SetSubelementPositions(mVSolved);
					this->m_pModel->FullBoundaryConditions();
				}

			}

			return true;
		}

		BC_Fixed* CreateBCFromSelBoxBounds(const Vector3d& vmin, const Vector3d& vmax)
		{
			BCSetup bcSetup;

			bcSetup.m_type = BCType::Fixed;

			bcSetup.m_vDoF = m_pModel->SelectDoF(vmin, vmax);
			size_t numDoF = bcSetup.m_vDoF.size();

			bcSetup.m_maxStep = m_jmaxBCStep;
			bcSetup.m_incStep = m_jincBCStep;
			bcSetup.m_maxError = m_jmaxBCError;
			bcSetup.m_vendT.setZero();
			bcSetup.m_vendR.setZero();

			bcSetup.m_vini.resize(numDoF);
			for (size_t i = 0; i < numDoF; ++i)
				bcSetup.m_vini[i] = bcSetup.m_vDoF[i]->GetPosition_0();
				
			return new BC_Fixed(m_pModel.get(), bcSetup);
		}

		bool SampleSpiral()
		{
			if (!HasSpiral())
				return false;

			int curTrial = 1;
			int maxTrial = 10;
			int numSucc = 0;
			int maxSucc = 3;

			// Solver max step

			OptimSolver_SQP* pSolver = static_cast<OptimSolver_SQP*>(this->m_pSolver.get());

			Real maxStepSize = pSolver->Options().maxStepSize;

			// BC max step

			iVector vmaxStep(5);
			for (int i = 0; i < 5; ++i)
				vmaxStep[i] = m_vpBC[i]->Setup().m_maxStep;

			while (HasSample())
			{
				this->InitSample();

				SolveResult result;

				if (!this->m_jpreSolved)
				{
					for (int i = 0; i < (int) this->m_vpBC.size(); ++i)
					{
						m_vpBC[i]->Setup().m_maxStep = vmaxStep[i]*curTrial;
					}
				}
				else
				{
					pSolver->Options().maxStepSize = maxStepSize/curTrial;
				}

				result = this->m_pSolver->SolveFull();
				if (result != SolveResult::Success)
				{
					numSucc = 0;

					if (curTrial <= maxTrial)
					{
						curTrial++;
						continue;
					}
				}
				else
				{
					if (numSucc < maxSucc)
					{
						numSucc++;
					}
					else
					{
						curTrial = max(curTrial-1, 1);
					}
				}
				
				// Check strain here

				const vector<EnergyElement_FEM*>& velements = this->m_pModel->GetEnergyElements_FEM();
				int numEleDef = (int) velements.size();
				dVector vFI(numEleDef);
				dVector vSI(numEleDef);
				dVector vTI(numEleDef);
				for (int k = 0; k < numEleDef; ++k)
				{
					MatrixXd F = velements[k]->ComputeDeformationGradient();
					MatrixXd C = F.transpose()*F;
					Real fi = C.trace();
					Real ti = C.determinant();
					Real si = 0.5*(fi*fi - (C*C).trace());
					vFI[k] = fi;
					vSI[k] = si;
					vTI[k] = ti;
				}

				// Sort
				int numEleFilMin = 0.01*numEleDef;
				int numEleFilMax = 0.99*numEleDef;
				int numEleFil = numEleFilMax - numEleFilMin;
				std::sort(vFI.begin(), vFI.end());
				std::sort(vSI.begin(), vSI.end());
				std::sort(vTI.begin(), vTI.end());
				VectorXd vFIFil(numEleFil);
				VectorXd vSIFil(numEleFil);
				VectorXd vTIFil(numEleFil);
				int count = 0;
				for (int k = numEleFilMin; k < numEleFilMax; ++k)
				{
					vFIFil(count) = vFI[k];
					vSIFil(count) = vSI[k];
					vTIFil(count) = vTI[k];
					count++;
				}

				string fiLogFile = this->m_joutputFolder + "/firstInvariant.csv";
				string siLogFile = this->m_joutputFolder + "/secondInvariant.csv";
				string tiLogFile = this->m_joutputFolder + "/thirdInvariant.csv";

				fstream fiLogOut;
				fstream siLogOut;
				fstream tiLogOut;
				fiLogOut.open(fiLogFile.c_str(), fstream::out | fstream::app);
				siLogOut.open(siLogFile.c_str(), fstream::out | fstream::app);
				tiLogOut.open(tiLogFile.c_str(), fstream::out | fstream::app);
				fiLogOut << vFIFil.minCoeff() << ", " << vFIFil.maxCoeff() << ", " << vFIFil.mean() << ", " << vFIFil(numEleFil/2) << endl;
				siLogOut << vSIFil.minCoeff() << ", " << vSIFil.maxCoeff() << ", " << vSIFil.mean() << ", " << vSIFil(numEleFil/2) << endl;
				tiLogOut << vTIFil.minCoeff() << ", " << vTIFil.maxCoeff() << ", " << vTIFil.mean() << ", " << vTIFil(numEleFil/2) << endl;
				fiLogOut.flush();
				fiLogOut.close();
				siLogOut.flush();
				siLogOut.close();
				tiLogOut.flush();
				tiLogOut.close();

				// <Write results here>

				Json joutput;

				// Original setup

				joutput["setup"] = m_jsetup;

				// Boundary

				for (int i = 0; i < 4; ++i)
				{
					const VectorXd& vendT = this->m_vpBC[i]->Setup().m_vendT;
					const VectorXd& vendR = this->m_vpBC[i]->Setup().m_vendR;
					Json jbc;
					for (int j = 0; j < 3; ++j)
					{
						jbc["T"].emplace_back(vendT(j));
						jbc["R"].emplace_back(vendR(j));
					}
					joutput["boundary"]["transform"].emplace_back(jbc);
				}

				const Sample& s = this->m_vsampleList[this->m_countSample];
				joutput["boundary"]["parameters"]["type"] = s.m_type;
				joutput["boundary"]["parameters"]["stretchX"] = s.m_stretchX;
				joutput["boundary"]["parameters"]["stretchY"] = s.m_stretchY;
				joutput["boundary"]["parameters"]["shearX"] = s.m_shearX;
				joutput["boundary"]["parameters"]["shearY"] = s.m_shearY;
				joutput["boundary"]["parameters"]["bendingInX"] = s.m_bendingInX;
				joutput["boundary"]["parameters"]["bendingInY"] = s.m_bendingInY;
				joutput["boundary"]["parameters"]["bendingOutX1"] = s.m_bendingOutX1;
				joutput["boundary"]["parameters"]["bendingOutX2"] = s.m_bendingOutX2;
				joutput["boundary"]["parameters"]["bendingOutY1"] = s.m_bendingOutY1;
				joutput["boundary"]["parameters"]["bendingOutY2"] = s.m_bendingOutY2;
				joutput["boundary"]["parameters"]["twistX1"] = s.m_twistX1;
				joutput["boundary"]["parameters"]["twistX2"] = s.m_twistX2;
				joutput["boundary"]["parameters"]["twistY1"] = s.m_twistY1;
				joutput["boundary"]["parameters"]["twistY2"] = s.m_twistY2;

				// Features 0

				for (int i = 0; i < 5; ++i)
				{
					Json jpos;
					Vector3d vp = this->m_vTP[3 * i + 0].GetPosition(Space::MAT);
					for (int j = 0; j < 3; ++j)
					{
						jpos["pos"].emplace_back(vp(j));
					}
					Json jtan;
					Json jnor;
					Json jbin;
					for (int k = 0; k < 3; ++k)
					{
						jtan["tan"].emplace_back(m_vF0[i].tan(k));
						jnor["nor"].emplace_back(m_vF0[i].nor(k));
						jbin["bin"].emplace_back(m_vF0[i].bin(k));
					}

					Json jfea;
					jfea["feature"].emplace_back(jpos);
					jfea["feature"].emplace_back(jtan);
					jfea["feature"].emplace_back(jnor);
					jfea["feature"].emplace_back(jbin);
					joutput["features"].push_back(jfea);
				}

				// Features x

				m_vFx.resize(5);

				for (size_t i = 0; i < 5; ++i)
				{
					m_vFx[i].nor = (this->m_vTP[3 * i + 1].GetPosition() - this->m_vTP[3 * i + 0].GetPosition()).normalized();
					m_vFx[i].bin = (this->m_vTP[3 * i + 2].GetPosition() - this->m_vTP[3 * i + 0].GetPosition()).normalized();
					m_vFx[i].tan = m_vFx[i].nor.cross(m_vFx[i].bin);
				}

				for (int i = 0; i < 5; ++i)
				{
					Json jpos;
					Vector3d vp = this->m_vTP[3 * i + 0].GetPosition(Space::DEF);
					for (int j = 0; j < 3; ++j)
					{
						jpos["pos"].emplace_back(vp(j));
					}
					Json jtan;
					Json jnor;
					Json jbin;
					for (int k = 0; k < 3; ++k)
					{
						jtan["tan"].emplace_back(m_vFx[i].tan(k));
						jnor["nor"].emplace_back(m_vFx[i].nor(k));
						jbin["bin"].emplace_back(m_vFx[i].bin(k));
					}

					Json jfea;
					jfea["feature"].emplace_back(jpos);
					jfea["feature"].emplace_back(jtan);
					jfea["feature"].emplace_back(jnor);
					jfea["feature"].emplace_back(jbin);
					joutput["results"]["features"].push_back(jfea);
				}

				// Energy

				if (result == SolveResult::Success)
				{
					joutput["results"]["state"] = "SUCCESS";
				}
				else if (result == SolveResult::MaxIter)
				{
					joutput["results"]["state"] = "MAXITER";
				}
				else if (result == SolveResult::NonDesc)
				{
					joutput["results"]["state"] = "NONDESC";
				}
				else if (result == SolveResult::Failure)
				{
					joutput["results"]["state"] = "FAILURE";
				}
				else
				{
					joutput["results"]["state"] = "OTHER";
				}

				joutput["results"]["energy"] = this->m_pModel->GetEnergy();
				joutput["results"]["gradient"] = this->m_pModel->GetGradient().norm();

				string summary = joutput.dump(4);

				// Write sample summary

				ostringstream directFile;
				ostringstream resultFile;
				ostringstream summaryFile;
				ostringstream objectFile;

				directFile << this->m_joutputFolder << "/" << m_vspiralList[this->m_countSpiral];
				summaryFile << directFile.str() << "/" << "sample_" << std::setw(6) << std::setfill('0') << this->m_countSample << ".json";
				objectFile << directFile.str() << "/" << "sample_" << std::setw(6) << std::setfill('0') << this->m_countSample << ".obj";
				resultFile << directFile.str() << "/" << "sample_" << std::setw(6) << std::setfill('0') << this->m_countSample << ".csv";
				
				// Write result summary
				
				ofstream summaryWriter;
				summaryWriter.open(summaryFile.str().c_str(), ios::out);
				summaryWriter << summary;
				summaryWriter.close();

				// Write result object

				int numD = 3;
				int numV = this->m_mVertices.rows();
				VectorXd vx;
				this->m_pModel->GetFullDOFPosition(vx);
				VectorXd vxSub = vx.segment(0, numV*numD);
				m_mVertices = MatrixXd(vxSub);
				m_mVertices.resize(numD, numV);
				m_mVertices.transposeInPlace();

				writeTriMesh_Obj(objectFile.str(), this->Vertices(), this->Surface());

				// Write result matrix

				writeMatrix_CSV(resultFile.str(), m_mVertices);

				// Log result summary

				ostringstream logLine;
				logLine << this->m_countSpiral << "_";
				logLine << this->m_countSample << "_";
				logLine << s.m_type << "_";
				logLine << curTrial << "_";
				logLine << this->m_pModel->GetEnergy() << "_";
				logLine << this->m_pModel->GetGradient().norm() << "_";
				logLine << joutput["results"]["state"].get<string>() << "\n";
				
				string logFile = this->m_joutputFolder + "/samplingLog.txt";

				fstream logOut;
				logOut.open(logFile.c_str(),
							fstream::out | 
							fstream::app);
				logOut << logLine.str();
				logOut.flush();
				logOut.close();

				// </Write results here>

				this->NextSample();
			}

			return HasSpiral();
		}

		bool NextSpiral()
		{
			float thick = 0;
			float scale = 0;
			float twist = 0;

			do
			{
				this->m_countSpiral++;

				if (!HasSpiral())
					return false;

				int idx = m_countSpiral;

				vector<string> vdescriptor;
				splitString(m_vspiralList[idx], "_", vdescriptor);
				thick = stof(vdescriptor[1]);
				scale = stof(vdescriptor[3]);
				twist = stof(vdescriptor[5]);
			} while (twist == 0);

			if (HasSpiral())
			{
				return true;
			}
			else
			{
				return false;
			}
		}

		bool NextSample()
		{
			m_countSample++;

			if (HasSample())
			{
				return true;
			}
			else
			{
				return false;
			}	
		}

		bool HasSpiral()
		{
			return this->m_countSpiral <= this->m_jendSpiralIndex;
		}

		bool HasSample()
		{
			return this->m_countSample <= this->m_jendSampleIndex;
		}

	};
}
