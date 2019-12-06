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
#include <PhySim/Geometry/TetGenReader.h>

#include <PhySim/Solvers/OptimProblem_BasicStatic.h>
#include <PhySim/Solvers/OptimProblem_BasicDynamic.h>

#include <PhySim/Solvers/OptimSolver_USQP_LS.h>

#include <PhySim/Utils/json.hpp>
#include <PhySim/Utils/objload.h>

#include <fstream>
#include <streambuf>

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	using Json = nlohmann::json;

	class Simulator_Volume
	{

	private:

		vector<BCondition*>						m_vpBC;
		shared_ptr<OptimProblem>				m_pProblem;
		shared_ptr<OptimSolver_SQP>					m_pSolver;
		shared_ptr<Model_FEM>					m_pModel;
		MatrixXi								m_mSurface;
		MatrixXd								m_mVertices;
		Json									m_jsetup;


	public:
		Simulator_Volume()
		{
			this->m_pModel.reset();
			this->m_pProblem.reset();
			this->m_pSolver.reset();
			this->m_vpBC.clear();
		}

		virtual ~Simulator_Volume(void)
		{
			// Nothing to do here...
		}

		shared_ptr<OptimProblem> Problem() { return this->m_pProblem; }
		shared_ptr<OptimSolver_SQP> Solver() { return this->m_pSolver; }
		shared_ptr<Model_FEM> Model() { return this->m_pModel; }
		const MatrixXd& Vertices() { return this->m_mVertices; }
		const MatrixXi& Surface() { return this->m_mSurface; }
		vector<BCondition*>& BC() { return this->m_vpBC; }

		virtual bool LoadMaterialParam(Json& ele, Material& mat)
		{
			if (ele.find("materialParam") == ele.end())
			{
				return false;
			}
			else
			{
				if (ele["materialParam"].find("young") != ele["materialParam"].end()) mat.AddProperty(Material::Property::Young, ele["materialParam"]["young"]);
				if (ele["materialParam"].find("poisson") != ele["materialParam"].end()) mat.AddProperty(Material::Property::Poisson, ele["materialParam"]["poisson"]);
				if (ele["materialParam"].find("density") != ele["materialParam"].end()) mat.AddProperty(Material::Property::Density, ele["materialParam"]["density"]);
				if (ele["materialParam"].find("moone10") != ele["materialParam"].end()) mat.AddProperty(Material::Property::Mooney10, ele["materialParam"]["moone10"]);
				if (ele["materialParam"].find("moone01") != ele["materialParam"].end()) mat.AddProperty(Material::Property::Mooney01, ele["materialParam"]["moone01"]);
				if (ele["materialParam"].find("lame1") != ele["materialParam"].end()) mat.AddProperty(Material::Property::Lame1, ele["materialParam"]["lame1"]);
				if (ele["materialParam"].find("lame2") != ele["materialParam"].end()) mat.AddProperty(Material::Property::Lame2, ele["materialParam"]["lame2"]);
				
				if (mat.HasProperty(Material::Property::Young) &&
					mat.HasProperty(Material::Property::Poisson) &&
					mat.HasProperty(Material::Property::Density))
				{
					mat.InitRealisticFromYoungPoisson(mat[Material::Property::Young],
													  mat[Material::Property::Poisson],
													  mat[Material::Property::Density]);
					return true;
				}

				if (mat.HasProperty(Material::Property::Lame1) &&
					mat.HasProperty(Material::Property::Lame2) &&
					mat.HasProperty(Material::Property::Density))
				{
					mat.InitRealisticFromLameParameter(mat[Material::Property::Lame1],
													   mat[Material::Property::Lame2],
													   mat[Material::Property::Density]);
					return true;
				}
			}

			return true;
		}

		virtual bool LoadMaterialModel(Json& ele, MaterialModel& mat)
		{
			if (ele.find("materialModel") == ele.end())
			{
				return false;
			}
			else
			{
				if (ele["materialModel"].get<string>().compare("StVK") == 0)
				{
					mat = PhySim::MaterialModel::StVK;
					return true;
				}

				if (ele["materialModel"].get<string>().compare("CoNH") == 0)
				{
					mat = PhySim::MaterialModel::CoNH;
					return true;
				}

				if (ele["materialModel"].get<string>().compare("InNH") == 0)
				{
					mat = PhySim::MaterialModel::InNH;
					return true;
				}

				if (ele["materialModel"].get<string>().compare("CoMR") == 0)
				{
					mat = PhySim::MaterialModel::CoMR;
					return true;
				}

				if (ele["materialModel"].get<string>().compare("InMR") == 0)
				{
					mat = PhySim::MaterialModel::InMR;
					return true;
				}
			}

			return false;
		}

		virtual bool LoadModel(Json& ele, Model_FEM& model)
		{
			if (ele.find("model") == ele.end())
			{
				return false;
			}
			else
			{
				// Load tetrahedral mesh

				MatrixXd mV;
				MatrixXi mT;
				MatrixXi mE;
				vector<MatrixXi> vmF;

				string inputMesh = ele["model"]["input"].get<string>();
				if (!TetGenReader::readPoly(inputMesh, mV, mT, mE, vmF))
				{
					logSimu("\n[ERROR] Invalid input tetrahedral mesh");
					return false;
				}

				m_mVertices = mV;
				m_mSurface = vmF[0];

				// Create FEM model

				model.SetPresolveFixed(ele["model"]["preSolveFix"].get<bool>());
				model.GetOptions().m_mNodes = mV;
				model.GetOptions().m_mElems = mT;

				Material matParam;
				MaterialModel matModel;
				this->LoadMaterialParam(ele["model"], matParam);
				this->LoadMaterialModel(ele["model"], matModel);

				if (ele["model"]["element"].get<string>().compare("Tet4") == 0)
					model.GetOptions().m_discretization = PhySim::Discretization::Tetrahedra4;

				if (ele["model"]["element"].get<string>().compare("Tet10") == 0)
					model.GetOptions().m_discretization = PhySim::Discretization::Tetrahedra10;

				model.GetOptions().m_numQuadrature = ele["model"]["quadrature"].get<int>();

				this->m_pModel->GetOptions().m_materialModel = matModel;
				this->m_pModel->GetOptions().m_material = matParam;
				this->m_pModel->Init();
			}

			return true;
		}

		virtual bool LoadProblem(Json& ele, shared_ptr<OptimProblem>& pProblem)
		{
			if (ele.find("problem") == ele.end())
			{
				return false;
			}
			else
			{
				string typestr = ele["problem"]["type"].get<string>();
				
				if (typestr.compare("Static") == 0)
				{
					pProblem = shared_ptr<OptimProblem>(new OptimProblem_BasicStatic(m_pModel.get()));
					return true;
				}

				if (typestr.compare("Dynamic") == 0)
				{
					pProblem = shared_ptr<OptimProblem>(new OptimProblem_BasicDynamic(m_pModel.get(), 0.01));
					return true;
				}
			}

			return true;
		}

		virtual bool LoadSolver(Json& ele, shared_ptr<OptimSolver_SQP>& pSolver)
		{
			if (ele.find("solver") == ele.end())
			{
				return false;
			}
			else
			{
				OptimSolverOptions options;

				if (ele["solver"].find("maxIters") != ele["solver"].end())
					options.maxIters = ele["solver"]["maxIters"].get<int>();

				if (ele["solver"].find("maxError") != ele["solver"].end())
					options.maxError = ele["solver"]["maxError"].get<Real>();
				
				if (ele["solver"].find("linear") != ele["solver"].end())
				{
					string typestr = ele["solver"]["linear"]["type"].get<string>();

					if (typestr.compare("CUDASC") == 0) options.lsSolverType = LS_CUDASC;
					if (typestr.compare("EigenCG") == 0) options.lsSolverType = LS_EigenCG;
					if (typestr.compare("EigenLU") == 0) options.lsSolverType = LS_EigenLU;
					if (typestr.compare("BiCGSTAB") == 0) options.lsSolverType = LS_BiCGSTAB;
					if (typestr.compare("EigenLDLT") == 0) options.lsSolverType = LS_EigenLDLT;
					if (typestr.compare("SSparseSPQR") == 0) options.lsSolverType = LS_SSparseSPQR;
					if (typestr.compare("CholmodLDLT") == 0) options.lsSolverType = LS_CholmodLDLT;
				}

				pSolver = shared_ptr<OptimSolver_SQP>(static_cast<OptimSolver_SQP*>(new OptimSolver_USQP_LS(m_pProblem.get(), options)));
			}

			return true;
		}

		virtual bool LoadVector(Json& ele, VectorXd& vec)
		{
			int numEle = (int)ele.size();

			vec.resize(numEle);
			for (int i = 0; i < numEle; ++i)
				vec(i) = ele[i].get<Real>();

			return true;
		}

		virtual bool LoadBoundary(Json& ele, vector<BCondition*>& vpBC)
		{
			if (ele.find("boundary") == ele.end())
			{
				return false;
			}
			else
			{
				int numBC = ele["boundary"].size();

				// Load all boundary conditions

				for (int i = 0; i < numBC; ++i)
				{
					BCondition* pBC;
					if (this->LoadBC(ele["boundary"][i], &pBC))
						vpBC.push_back(pBC); // Add BC to list
				}
			}

			return true;
		}

		virtual bool LoadBC(Json& ele, BCondition** ppBC)
		{
			BCSetup bcSetup;

			// Types

			string bctypestr = ele["type"].get<string>();

			if (bctypestr.compare("FIXED") == 0)
				bcSetup.m_type = BCType::Fixed;

			if (bctypestr.compare("FORCE") == 0)
				bcSetup.m_type = BCType::Force;

			if (bctypestr.compare("GRAVITY") == 0)
				bcSetup.m_type = BCType::Gravity;

			// Selection

			if (bcSetup.m_type != BCType::Gravity)
			{
				string seltypestr = ele["selection"]["type"].get<string>();

				if (seltypestr.compare("BOXFILE") == 0 ||
					seltypestr.compare("BOXRANGE") == 0)
				{
					Vector3d vmin;
					Vector3d vmax;
					if (seltypestr.compare("BOXFILE") == 0)
					{
						string filestr = ele["selection"]["file"];

						MatrixXd mV; 
						MatrixXi mF;
						if (readTriMesh_Obj(filestr, mV, mF))
							computeBoundingBox(mV, vmin, vmax);
					}

					if (seltypestr.compare("BOXRANGE") == 0)
					{
						VectorXd vminBox;
						VectorXd vmaxBox;
						this->LoadVector(ele["selection"]["boxMin"], vminBox);
						this->LoadVector(ele["selection"]["boxMax"], vmaxBox);
						vmin = vminBox;
						vmax = vmaxBox;
					}

					bcSetup.m_vDoF = m_pModel->SelectDoF(vmin, vmax);
				}

				if (seltypestr.compare("INDICES") == 0)
				{
					// TODO
				}
			}
			else
			{
				bcSetup.m_vDoF = this->m_pModel->GetDoFSets();
			}

			size_t numDoF = bcSetup.m_vDoF.size();

			// Loading

			bcSetup.m_incStep = ele["incBCStep"].get<int>();
			bcSetup.m_maxStep = ele["maxBCStep"].get<int>();
			bcSetup.m_maxError = ele["maxToStep"].get<double>();


			// Values

			VectorXd vini;
			this->LoadVector(ele["iniV"], vini);
			this->LoadVector(ele["endT"], bcSetup.m_vendT);
			this->LoadVector(ele["endR"], bcSetup.m_vendR);

			if (bcSetup.m_type == BCType::Fixed)
			{
				bcSetup.m_vini.resize(numDoF);
				for (int i = 0; i < numDoF; ++i)
					bcSetup.m_vini[i] = bcSetup.m_vDoF[i]->GetPosition_0();
			}

			if (bcSetup.m_type == BCType::Force)
			{
				bcSetup.m_vini.resize(numDoF);
				for (int i = 0; i < numDoF; ++i)
					bcSetup.m_vini[i] = vini / numDoF;
			}

			if (bcSetup.m_type == BCType::Gravity)
			{
				bcSetup.m_vini.resize(numDoF);
				for (int i = 0; i < numDoF; ++i)
					bcSetup.m_vini[i] = vini;
			}
			
			// Create

			if (bcSetup.m_type == BCType::Fixed) *ppBC = new BC_Fixed(m_pModel.get(), bcSetup);
			if (bcSetup.m_type == BCType::Force) *ppBC = new BC_Force(m_pModel.get(), bcSetup);
			if (bcSetup.m_type == BCType::Gravity) *ppBC = new BC_Gravity(m_pModel.get(), bcSetup);
			
			return true;
		}

		virtual bool Init(string& setupFile)
		{
			// Read setup JSON

			string setupFileContent;

			ifstream setupFileIn(setupFile);
			setupFileIn.seekg(0, ios::end);

			setupFileContent.reserve(setupFileIn.tellg());

			setupFileIn.seekg(0, ios::beg);

			setupFileContent.assign((istreambuf_iterator<char>(setupFileIn)), istreambuf_iterator<char>());

			try
			{
				m_jsetup = Json::parse(setupFileContent);
			}
			catch (exception e)
			{
				logSimu("\n[ERROR] Invalid setup file: %s", e.what());
				return false;
			}

			// Load model

			this->m_pModel = shared_ptr<Model_FEM>(new Model_FEM);
			if (!this->LoadModel(m_jsetup, *this->m_pModel.get()))
			{
				logSimu("\n[ERROR] Impossible to load model");
				return false;
			}

			// Load boundary

			this->m_vpBC.clear();
			LoadBoundary(m_jsetup, this->m_vpBC);
			for (int i = 0; i < (int) m_vpBC.size(); ++i)
				this->m_pModel->AddBoundaryCondition(this->m_vpBC[i]);

			// Prepare for simulation

			this->m_pModel->PrepareForSimulation();

			// Load problem/solver

			this->LoadProblem(m_jsetup, this->m_pProblem);
			this->LoadSolver(m_jsetup, this->m_pSolver);
		}

		virtual SolveResult SolveStep()
		{
			return this->SolveStep();
		}

		virtual SolveResult SolveFull()
		{
			return this->SolveFull();
		}

	};
}