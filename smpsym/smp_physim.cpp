// This file is part of smpsym, an inverse design and simulation tool for
// research paper Guseinov R. et al "Programming temporal morphing of
// self-actuated shells"
//
// Copyright (C) 2019 Ruslan Guseinov <guseynov.ruslan@gmail.com>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

#include "smp_physim.h"

#include "PhySim/Solvers/OptimSolver_Staggered.h"
#include "PhySim/Solvers/OptimProblem_Subproblem.h"

#include "PhySim/Geometry/RigidBody.h"

namespace smpup {

	void SmpPhySim::init(const StaticSimulatorStruct& sss)
	{
		*(StaticSimulatorStruct*)this = sss;

		PhySim::Material material;
		if (mem_model == MEM_COARSE || mem_model == MEM_FINE)
		{
			material.InitRealisticFromLameParameter(0, get_smp_options().lame2, 1);
			material.AddProperty(PhySim::Material::Property::Thickness, get_smp_options().mem_h);
		}
		else
			material.InitFromDensity(1);
		material.AddProperty(PhySim::Material::Property::CollK, 0);
		material.AddProperty(PhySim::Material::Property::CollT, 0);
		material.AddProperty(PhySim::Material::Property::StretchK, 0);
		material.AddProperty(PhySim::Material::Property::BendingK, 0);
		material.AddProperty(PhySim::Material::Property::ShearK, get_smp_options().scissor);

		myModel.GetOptions().m_springs.m_vlinear = sss.linear;
		myModel.GetOptions().m_springs.m_mcrosses = sss.cross;
		myModel.GetOptions().m_material = material;
		myModel.GetOptions().m_mNodes = sss.nodes;
		myModel.GetOptions().m_mElems = sss.edges;
		myModel.GetOptions().m_mCollPairs = sss.flip;
		myModel.GetOptions().m_vCollDists.setZero(sss.flip.rows());

		myModel.GetOptions().m_linearType = PhySim::Model_MassSpring::SpringType::UPoly8;

		if (mem_model == MEM_COARSE || mem_model == MEM_FINE)
		{
			myModel.GetOptions().m_memDiscretization = PhySim::Discretization::Triangles;
			myModel.GetOptions().m_memMaterialModel = PhySim::MaterialModel::InNH;
			myModel.GetOptions().m_mmemElems = sss.memtri;
		}

		myModel.Init();

		if (mem_model == MEM_COARSE || mem_model == MEM_FINE)
		{
			myModel.SetMembraneRestStrains(sss.memrest);
		}

		std::vector<PhySim::Polytope*> vpolys;
		for (int i = 0; i < sss.bodies.size(); ++i)
		{
			int nodes = sss.bodies[i].size();
			std::vector<PhySim::Node*> vpolyNodes(nodes);
			for (int j = 0; j < nodes; ++j)
				vpolyNodes[j] = myModel.GetNodes()[sss.bodies[i][j]];
			vpolys.push_back(new PhySim::Cell_Poly(vpolyNodes));
		}

		myReduc.Init(&myModel, vpolys);
	}

	void SmpPhySim::precompute()
	{
		//myReduc.GetRigidBodies()[0]->CenterDoF()->Fix();
		//myReduc.GetRigidBodies()[0]->QuaterDoF()->Fix();

		myModel.FiniteDifferent_Gradient() = false;
		myModel.FiniteDifferent_Hessian() = false;
		myReduc.FiniteDifferent_Gradient() = false;
		myReduc.FiniteDifferent_Hessian() = false;

		myModel.PrepareForSimulation();
		myReduc.PrepareForSimulation();
	
		//myModel.TestGlobalGradient();
		//myModel.TestGlobalHessian();
		//myReduc.TestGlobalGradient();
		//myReduc.TestGlobalHessian();

		//Eigen::VectorXd vg = myReduc.GetGradient();
		//std::string vgstr = PhySim::matrixToString_CSV(vg);
		//PhySim::logSimu("\n[DEBUG] Gradient:\n%s", vgstr.c_str());

		myProblem = std::make_shared<PhySim::OptimProblem_BasicStatic>(&myReduc);
		auto pSolverSQP = std::make_shared<PhySim::OptimSolver_USQP_LS>();
		pSolverSQP->GetOptions().m_pProblem = myProblem.get();
		pSolverSQP->GetOptions().m_maxError = get_settings().tol;
		pSolverSQP->GetOptions().m_maxStepSize = get_settings().max_step_size;
		pSolverSQP->GetOptions().m_maxIters = get_settings().max_iter_dry;
		pSolverSQP->RegisterCallback_Step(physim_callback, this);
		mySolver = pSolverSQP;

		//myProblem = new PhySim::OptimProblem_BasicStatic(&myReduc);
		//PhySim::OptimSolver_SQP* pSolverSQP = new PhySim::OptimSolver_USQP_LS();
		//pSolverSQP->GetOptions().m_pProblem = myProblem;
		//pSolverSQP->GetOptions().m_maxError = get_settings().tol;
		//pSolverSQP->GetOptions().m_maxStepSize = get_settings().max_step_size;
		//pSolverSQP->GetOptions().m_maxIters = get_settings().max_iter_dry;
		//pSolverSQP->RegisterCallback_Step(physim_callback, this);
		//mySolver = pSolverSQP;

		// Using staggered solver (TODO: Debug this part)

		//myProblem = new PhySim::OptimProblem_BasicStatic(&myReduc);

		//PhySim::OptimSolver_Staggered* pSolverStaggered = new PhySim::OptimSolver_Staggered();
		//pSolverStaggered->GetOptions().m_pProblem = myProblem;
		//pSolverStaggered->GetOptions().m_maxError = get_settings().tol;
		//pSolverStaggered->GetOptions().m_maxIters = get_settings().max_iter_dry;

		//PhySim::iVector vdofRigidBody(myReduc.NumRigidDoF());
		//PhySim::iVector vdofMembrane(myReduc.NumOtherDoF());
		//for (int i = 0; i < myReduc.NumRigidDoF(); ++i)
		//	vdofRigidBody[i] = i;
		//for (int i = 0; i < myReduc.NumOtherDoF(); ++i)
		//	vdofMembrane[i] = myReduc.NumRigidDoF() + i;

		//PhySim::OptimSolver_SQP* pSolverRigidBody = new PhySim::OptimSolver_USQP_LS();
		//PhySim::OptimSolver_SQP* pSolverMembrane = new PhySim::OptimSolver_USQP_LS();
		//PhySim::OptimProblem_Subproblem* pProblemRigidBody = new PhySim::OptimProblem_Subproblem(myProblem, vdofRigidBody);
		//PhySim::OptimProblem_Subproblem* pProblemMembrane = new PhySim::OptimProblem_Subproblem(myProblem, vdofMembrane);

		//pSolverRigidBody->GetOptions().m_pProblem = pProblemRigidBody;
		//pSolverRigidBody->GetOptions().m_maxError = get_settings().tol;
		//pSolverRigidBody->GetOptions().m_maxStepSize = get_settings().max_step_size;
		//pSolverRigidBody->Init();

		//pSolverMembrane->GetOptions().m_pProblem = pProblemMembrane;
		//pSolverMembrane->GetOptions().m_maxError = get_settings().tol;
		//pSolverMembrane->GetOptions().m_maxStepSize = get_settings().max_step_size;
		//pSolverMembrane->Init();

		//pSolverStaggered->GetOptions().m_vpProblem.resize(2);
		//pSolverStaggered->GetOptions().m_vpSolvers.resize(2);
		//pSolverStaggered->GetOptions().m_vpProblem[0].reset(pProblemRigidBody);
		//pSolverStaggered->GetOptions().m_vpProblem[1].reset(pProblemMembrane);
		//pSolverStaggered->GetOptions().m_vpSolvers[0].reset(pSolverRigidBody);
		//pSolverStaggered->GetOptions().m_vpSolvers[1].reset(pSolverMembrane);

		//mySolver = pSolverStaggered;

		mySolver->Init();
		mySolver->RegisterCallback_Step(physim_callback, this);
	}

	bool SmpPhySim::step(double dt, bool dry)
	{
		SolverState solver_state(myModel, myReduc);

		current_solver_info.reset(new SolverInfo());
		if (dry) mySolver->GetOptions().m_maxIters = get_settings().max_iter_dry;
		else mySolver->GetOptions().m_maxIters = get_settings().max_iter_wet;
		std::clock_t start = std::clock();
		PhySim::SolveResult status = mySolver->SolveFull();
		current_solver_info->time = (std::clock() - start) / (double)CLOCKS_PER_SEC;
		
		std::string string_status;
		switch (status) {
		case PhySim::SolveResult::Failure:  string_status = "failure"; break;
		case PhySim::SolveResult::MaxIter:  string_status = "max iteration stop"; break;
		case PhySim::SolveResult::NonDesc:  string_status = "non-descending stop"; break;
		case PhySim::SolveResult::NonSPD:   string_status = "non-SPD"; break;
		case PhySim::SolveResult::Singular: string_status = "singular matrix"; break;
		case PhySim::SolveResult::Stopped:  string_status = "stopped externally"; break;
		case PhySim::SolveResult::Success:  string_status = "success"; break;
		}
		std::clog << "SolveFull " << string_status << "\n";
		current_solver_info->status = string_status;

		if (get_settings().compute_condition)
		{
			int N = myReduc.GetNumFreeDOF();
			Eigen::SparseMatrix<double> Ht = myReduc.GetHessian().selfadjointView<Eigen::Lower>();
			Eigen::MatrixXd H = Ht.toDense();
			Eigen::BDCSVD<Eigen::MatrixXd> svdSolver(H);
			const Eigen::VectorXd& vs = svdSolver.singularValues();
			double maxEigen = -HUGE_VAL;
			double minEigen = HUGE_VAL;
			for (int i = 0; i < N; ++i)
			{
				if (vs(i) != 0 && vs(i) > maxEigen)
					maxEigen = vs(i);
				if (vs(i) != 0 && vs(i) < minEigen)
					minEigen = vs(i);
			}
			current_solver_info->condition = maxEigen / minEigen;
		}

		if (status == PhySim::SolveResult::Success || status == PhySim::SolveResult::MaxIter || status == PhySim::SolveResult::NonDesc)
		{
			if (get_settings().freeze_centroid)
			{
				Eigen::Vector3d cm = Eigen::Vector3d::Zero();
				for (int i = 0; i < bodies.size(); ++i)
					cm += myReduc.GetRigidBodies()[i]->GetCentroid();
				cm /= bodies.size();
				Eigen::Vector3d vt = cm;

				//PhySim::RigidBody* pRB = myReduc.GetRigidBodies()[0];
				//Eigen::Matrix3d mR = pRB->GetRotation();
				//Eigen::Vector3d vt = pRB->GetCentroid();

				for (int i = 0; i < (int)myModel.GetNodes().size(); ++i)
				{
					Eigen::VectorXd vx = myModel.GetNodes()[i]->DoF()->GetPosition_x();
					//myModel.GetNodes()[i]->DoF()->SetPosition_x(mR.transpose()*(vx - vt));
					myModel.GetNodes()[i]->DoF()->SetPosition_x(vx - vt);
				}

				for (int i = 0; i < (int)myReduc.GetRigidBodies().size(); ++i)
				{
					myReduc.GetRigidBodies()[i]->RecomputeRotation();
					myReduc.GetRigidBodies()[i]->RecomputeCentroid();
				}

				Eigen::VectorXd pos;
				myReduc.GetFreeDOFPosition(pos);
				int numMemNodes = myReduc.NumOtherDoF()/3;
				for (int i = 0; i < numMemNodes; ++i)
				{
					Eigen::VectorXd vx = pos.segment(myReduc.NumRigidDoF() + i * 3, 3);
					//pos.segment(myReduc.NumRigidDoF() + i * 3, 3) = mR.transpose()*(vx - vt);
					pos.segment(myReduc.NumRigidDoF() + i * 3, 3) = vx - vt;
				}

				myReduc.SetFreeDOFPosition(pos);
			}

			return true;
		}

		solver_state.apply_state(myModel, myReduc);
		return false;
	}

	void SmpPhySim::get_frame_info(FrameInfo& finfo)
	{
		bool is_loading = (0 < finfo.body_dofs.size());

		if (is_loading)
		{
			finfo.apply_state(myModel, myReduc);
			myModel.SetRestLinearLengths(finfo.restlen);
			update_brackets(finfo.time, finfo.plastic_fraction); // no plasticity for the dry simulation step
		}
		else
		{
			finfo.extract_state(myModel, myReduc);
			finfo.restlen = get_restlen();
		}
		
		myReduc.ComputeAndStore_Energy();
		myReduc.ComputeAndStore_Gradient();

		finfo.nodes = get_positions();

		const auto& rbodies = myReduc.GetRigidBodies();
		finfo.body_centroids.resize(rbodies.size(), 3);
		for (int i = 0; i < rbodies.size(); ++i)
			finfo.body_centroids.row(i) = rbodies[i]->GetCentroid().transpose();

		auto segments = myModel.GetEnergyElements_SpringsLinear();

		finfo.linear_spring_energy.resize(myModel.GetEnergyElements_SpringsLinear().size());
		for (int i = 0; i < myModel.GetEnergyElements_SpringsLinear().size(); ++i)
			finfo.linear_spring_energy(i) = segments[i]->GetElementEnergy();

		finfo.total_segment_energy = finfo.linear_spring_energy.sum();

		finfo.collided.resize(myModel.GetEnergyElements_SpringsLinear().size());
		for (int i = 0; i < myModel.GetEnergyElements_SpringsLinear().size(); ++i)
		{
			if (segments[i]->GetRestLength() == bumpers(i)) finfo.collided(i) = 1;
			else finfo.collided(i) = 0;
		}

		finfo.segment_forces.setZero(linear.size());
		for (int i = 0; i < finfo.segment_forces.size(); ++i)
			finfo.segment_forces(i) = segments[i]->GetElementGradient().head(3).norm();

		auto scissors = myModel.GetEnergyElements_SpringsCross();

		double total_cross_strain = 0;
		double total_cross_energy = 0;
		for (int i = 0; i < scissors.size(); ++i)
		{
			total_cross_strain += scissors[i]->ComputeStrain();
			total_cross_energy += scissors[i]->GetElementEnergy();
		}
		finfo.total_cross_strain = total_cross_strain;
		finfo.total_cross_energy = total_cross_energy;

		auto membrane = myModel.GetEnergyElements_Membrane(); // TODO: depends on membrane model: NH or segment

		for (int i = 0; i < membrane.size(); ++i)
			finfo.total_membrane_energy += membrane[i]->GetElementEnergy();

		// current length
		//finfo.restlen.setZero(linear.size());
		//for (int i = 0; i < finfo.restlen.size(); ++i)
		//	finfo.restlen(i) = (finfo.nodes.row(edges(linear(i), 0)) - finfo.nodes.row(edges(linear(i), 1))).norm();

		finfo.body_origins.resize(bodies.size(), 3);
		for (int i = 0; i < bodies.size(); ++i)
			finfo.body_origins.row(i) = finfo.nodes.row(bodies[i][0]);

		if (!is_loading)
		{
			finfo.solver_info = current_solver_info;
			if (current_solver_info) finfo.solver_info->num_iter = finfo.solver_info->iter.size();
			current_solver_info.reset();
		}
	}

	void SmpPhySim::update_brackets(double t, double plastic_fraction)
	{
		std::vector<Eigen::VectorXd> stiffness_poly = get_stiffness_polynomial();

		Eigen::VectorXd restlen = get_restlen();

		Eigen::MatrixXd poly;
		get_bracket_model()->make_energy_poly_coeff(thickness, segment_initial_lengths, t, poly);
		for (int i = 0; i < stiffness_poly.size(); ++i)
		{
			if (bumpers(i) < 0) continue; // membrane mimicking segment spring
			stiffness_poly[i] = poly.row(i).transpose();
			if (plastic_fraction != 0)
			{
				for (int j = 1; j < 9; ++j)
					stiffness_poly[i](j) /= std::pow(1 - plastic_fraction, j - 1);
			}
		}

		auto pos = get_positions();

		for (int i = 0; i < bumpers.size(); ++i)
		{
			if (bumpers(i) < 0) continue; // membrane mimicking segment spring

			double ps = (pos.row(edges(linear(i), 0)) - pos.row(edges(linear(i), 1))).norm();

			if (plastic_fraction != 0)
				restlen(i) = std::min(restlen(i), segment_initial_lengths(i) + (ps - segment_initial_lengths(i)) * std::min(1.0, plastic_fraction));

			if (ps <= bumpers(i))
			{
				stiffness_poly[i].setZero();
				stiffness_poly[i](2) = 99.9; //TODO: magic number, it was for stress, but now we have energy here, make sure it works
				restlen(i) = bumpers(i);
			}
		}

		set_stiffness_polynomial(stiffness_poly);
		set_rest_lengths(restlen);
	}

}
