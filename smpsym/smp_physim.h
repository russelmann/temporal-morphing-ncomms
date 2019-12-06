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

#ifndef SMP_PHYSIM_H
#define SMP_PHYSIM_H

#include "interface.h"

#include "PhySim/Geometry/Node.h"
#include "PhySim/Geometry/Edge.h"
#include "PhySim/Geometry/Cell_Poly.h"
#include "PhySim/Energies/EnergyElement_SpringLinear.h"
#include "PhySim/Energies/EnergyElement_SpringCross.h"
#include "PhySim/Energies/EnergyElement_FEM_MembraneInNH.h"
#include "PhySim/Models/Model_MassSpring_SMP.h"
#include "PhySim/Models/Model_Reduction_RBCloud.h"
#include "PhySim/Solvers/OptimSolver_USQP_LS.h"
#include "PhySim/Solvers/OptimProblem_BasicStatic.h"

#include "igl/readOBJ.h"


namespace smpup {

	struct SolverIterInfo
	{
		int num;
		double energy;
		double residual;

		template<typename Archive>
		void serialize(Archive& ar)
		{
			ar(cereal::make_nvp("num", num));
			ar(cereal::make_nvp("energy", energy));
			ar(cereal::make_nvp("residual", residual));
		}
	};

	struct SolverInfo
	{
		std::string status;
		float condition;
		double time; // full solver time, sec
		int num_iter;
		std::vector<SolverIterInfo> iter;

		SolverInfo()
		{
			condition = NAN;
		}

		template<typename Archive>
		void serialize(Archive& ar)
		{
			ar(cereal::make_nvp("status", status));
			ar(cereal::make_nvp("condition", condition));
			ar(cereal::make_nvp("num_iter", num_iter));
			ar(cereal::make_nvp("time", time));
			if (Archive::is_loading::value) ar(cereal::make_nvp("iter", iter));
			else
			{
				// store only first and last iterations to reduce storage size
				ar(cereal::make_nvp("iter", std::vector<SolverIterInfo>{ iter.front(), iter.back() }));
			}
		}
	};

	struct SolverState
	{
		Eigen::VectorXd body_dofs;            // all free dofs of myReduc (rigid body positions and current rotations)
		std::vector<Eigen::Matrix3d> body_rotation;
		Eigen::VectorXd membrane_dofs;

		SolverState() { }

		SolverState(const PhySim::Model_MassSpring_SMP& model, const PhySim::Model_Reduction_RBCloud& reduc)
		{
			extract_state(model, reduc);
		}

		void extract_state(const PhySim::Model_MassSpring_SMP& model, const PhySim::Model_Reduction_RBCloud& reduc)
		{
			PhySim::IModel::StateP pModel = model.CreateState(PhySim::Space::DEF);
			PhySim::IModel::StateP pReduc = reduc.CreateState(PhySim::Space::DEF);
			PhySim::Model_MassSpring_SMP::State* pModelState = static_cast<PhySim::Model_MassSpring_SMP::State*>(pModel.get());
			PhySim::Model_Reduction_RBCloud::State* pReducState = static_cast<PhySim::Model_Reduction_RBCloud::State*>(pReduc.get());

			model.GetState(pModel, PhySim::Space::DEF);
			reduc.GetState(pReduc, PhySim::Space::DEF);

			body_dofs = pReducState->m_vx;
			body_rotation = pReducState->m_vmR;
			membrane_dofs = pModelState->m_vx;
		}

		void apply_state(PhySim::Model_MassSpring_SMP& model, PhySim::Model_Reduction_RBCloud& reduc) const
		{
			PhySim::IModel::StateP pModel = model.CreateState(PhySim::Space::DEF);
			PhySim::IModel::StateP pReduc = reduc.CreateState(PhySim::Space::DEF);
			PhySim::Model_MassSpring_SMP::State* pModelState = static_cast<PhySim::Model_MassSpring_SMP::State*>(pModel.get());
			PhySim::Model_Reduction_RBCloud::State* pReducState = static_cast<PhySim::Model_Reduction_RBCloud::State*>(pReduc.get());

			pModelState->m_vx = membrane_dofs;
			pReducState->m_vx = body_dofs;
			pReducState->m_vmR = body_rotation;

			reduc.SetState(pReduc, PhySim::Space::DEF);
			model.SetState(pModel, PhySim::Space::DEF);
		}
	};

	// Simulation frame information tracking data
	struct FrameInfo : public SolverState
	{
		double time;
		double plastic_fraction;

		Eigen::VectorXd restlen;         // rest lengths of all segment springs of myModel

		Eigen::MatrixXd nodes;           // #nodes by 3 : node positions
		Eigen::MatrixXd body_centroids;
		
		Eigen::Matrix<float, -1, -1, Eigen::RowMajor> body_vertices;
		Eigen::Matrix<float, -1, -1, Eigen::RowMajor> body_normals;
		Eigen::Matrix<float, -1, -1, Eigen::RowMajor> bracket_vertices;
		Eigen::Matrix<float, -1, -1, Eigen::RowMajor> bracket_normals;
		Eigen::Matrix<float, -1, -1, Eigen::RowMajor> mem_normals;

		Eigen::VectorXd linear_spring_energy;
		float total_cross_strain;
		float total_cross_energy;
		float total_segment_energy;
		float total_membrane_energy;
		Eigen::VectorXd segment_forces;  // force magnitudes per segment springs
		Eigen::VectorXi collided;        // 1 if bumper collision has been detected
		Eigen::MatrixXd body_origins;

		std::shared_ptr<SolverInfo> solver_info;

		template<typename Archive>
		void serialize(Archive& ar)
		{
			ar(cereal::make_nvp("time", time));
			ar(cereal::make_nvp("plastic_fraction", plastic_fraction));
			ar(cereal::make_nvp("body_dofs", body_dofs)); // TODO: body rotations are not stored properly (effect of keeping it near zero?)
			ar(cereal::make_nvp("body_rotation", body_rotation));
			ar(cereal::make_nvp("membrane_dofs", membrane_dofs));
			ar(cereal::make_nvp("restlen", restlen));
			ar(cereal::make_nvp("solver_info", solver_info));
		}
	};


	class SmpPhySim : public StaticSimulatorStruct
	{
	private:
		PhySim::Model_MassSpring_SMP    myModel; // PhySim model 
		PhySim::Model_Reduction_RBCloud myReduc; // PhySim reduction

		std::shared_ptr<PhySim::OptimProblem> myProblem;  // PhySim problem
		std::shared_ptr<PhySim::OptimSolver>  mySolver;   // PhySim solver

		void(*external_physim_callback)(PhySim::IOptimSolver*, void*) = nullptr;
		void* external_callback_receiver;
		std::shared_ptr<SolverInfo> current_solver_info;

	public:
		SmpPhySim(std::shared_ptr<LinearMemModel> linear_mem_model = nullptr, std::shared_ptr<BracketModel>bracket_model = nullptr)
			: StaticSimulatorStruct(linear_mem_model, bracket_model) { }

		SmpPhySim(SmpOptions options, SimulatorSettings settings, std::shared_ptr<LinearMemModel> linear_mem_model, std::shared_ptr<BracketModel> bracket_model)
			: StaticSimulatorStruct(options, settings, linear_mem_model, bracket_model) { }


		// Initialize PhySim interface
		void init(const StaticSimulatorStruct& sss);

		// Precompute PhySim interface
		void precompute();

		static void physim_callback(PhySim::IOptimSolver* solver, void* receiver)
		{
			auto simulator = (SmpPhySim*)receiver;
			SolverIterInfo iter_info;
			iter_info.num = simulator->current_solver_info->iter.size() + 1;
			iter_info.energy = simulator->myReduc.GetEnergy();
			iter_info.residual = simulator->myReduc.GetGradient().norm();
			simulator->current_solver_info->iter.emplace_back(iter_info);
			if (simulator->external_physim_callback) simulator->external_physim_callback(solver, simulator->external_callback_receiver);
		}

		// Register callback
		void register_callback(void(*callback_fun)(PhySim::IOptimSolver*, void*), void *receiver = nullptr)
		{
			external_physim_callback = callback_fun;
			external_callback_receiver = receiver;
		}

		// Run PhySim solver
		bool step(double dt, bool dry);

		// Update bracket parameters
		void update_brackets(double t, double plastic_fraction);


		// Set rest cross spring ratios
		inline void set_rest_crossratios(const Eigen::VectorXd& rest)
		{
			myModel.SetRestCrossStrain(rest);
		}

		// Set bracket rest lengths
		inline void set_rest_lengths(const Eigen::VectorXd& rest)
		{
			myModel.SetRestLinearLengths(rest);
		}

		// Set bracket thicknesses
		inline void set_thickness(const Eigen::VectorXd& thick)
		{
			thickness = thick;
		}

		// Set segment springs stress-strain polynomials
		inline void set_stiffness_polynomial(const std::vector<Eigen::VectorXd>& polyk)
		{
			myModel.SetStretchPolyCoeff(polyk);
		}

		
		// Get current PhySim node positions
		Eigen::MatrixXd get_positions() const
		{
			Eigen::MatrixXd pos(myModel.GetNodes().size(), 3);
			for (int i = 0; i < myModel.GetNodes().size(); ++i)
				pos.row(i) = myModel.GetNodes()[i]->DoF()->GetPosition();
			return pos;
		}

		// Get current PhySim segment springs rest lengths
		inline Eigen::VectorXd get_restlen() const
		{
			return myModel.GetRestLinearLengths();
		}

		// Get segment springs stress-strain polynomials
		inline std::vector<Eigen::VectorXd> get_stiffness_polynomial() const
		{
			std::vector<Eigen::VectorXd> polyk;
			myModel.GetStretchPolyCoeff(polyk);
			return polyk;
		}

		// Record frame information
		void get_frame_info(FrameInfo& finfo);

		inline const std::shared_ptr<const SolverInfo> get_current_solver_info() const { return current_solver_info; }

		std::shared_ptr<const PhySim::OptimProblem> get_problem() const { return myProblem; }

		~SmpPhySim()
		{
		}

	};

}

#endif
