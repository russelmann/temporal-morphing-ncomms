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

#ifndef BRACKET_MODEL_H
#define BRACKET_MODEL_H

#include "utils/self_serializer.h"
#include "Eigen/Core"

namespace smpup {

	// Assumes rest length is zero
	class LinearMemModel : public Eigen::MatrixXd
	{
	public:
		// Compute energy polynomial based on initial length of segment membrane
		Eigen::VectorXd make_energy_poly_coeff(double linear_mem_init_length) const 
		{
			Eigen::VectorXd poly;
			Eigen::ArrayXd naturals(cols());
			for (int i = 0; i < cols(); ++i)
				naturals(i) = i;
			Eigen::VectorXd len_pow = Eigen::ArrayXd::Constant(cols(), linear_mem_init_length).pow(naturals);
			poly = *this * len_pow;
			return poly;
		}

		// Compute force polynomial based on initial length of segment membrane
		Eigen::VectorXd make_force_poly_coeff(double linear_mem_init_length) const
		{
			Eigen::VectorXd poly;
			Eigen::ArrayXd naturals(cols());
			for (int i = 0; i < cols(); ++i)
				naturals(i) = i;
			Eigen::VectorXd len_pow = Eigen::ArrayXd::Constant(cols(), linear_mem_init_length).pow(naturals);
			poly = this->bottomRows(rows() - 1) * len_pow;
			for (int i = 0; i < poly.size(); ++i)
				poly(i) *= i + 1;
			return poly;
		}

	};

	class BracketModel
	{
	private:
		// grid assumes order: thickness, length, time
		Eigen::Array3d grid_min;
		Eigen::Array3d grid_step;
		Eigen::Array3i grid_num;
		Eigen::MatrixXd beta; // polynomial coefficients, per column, starting from power 2
	public:
		enum GRID { GRID_THICK, GRID_LENGTH, GRID_TIME };

		BracketModel() { }

		BracketModel(Eigen::VectorXd beta_single)
		{
			grid_min = Eigen::Array3d{ 0, 0, 0 };
			grid_step = Eigen::Array3d{ 1, 10, 2 };
			grid_num = Eigen::Array3i{ 2, 2, 2 };
			beta.resize(beta_single.size(), 8);
			for (int i = 0; i < 8; ++i)
				beta.col(i) = beta_single;
		}

		Eigen::VectorXd make_energy_poly_coeff(double thickness, double bracket_length, double time) const;
		Eigen::VectorXd make_force_poly_coeff(double thickness, double bracket_length, double time) const;

		void make_energy_poly_coeff(const Eigen::VectorXd& thickness, const Eigen::VectorXd& bracket_length, double time, Eigen::MatrixXd& poly) const;

		double configure_thickness(double bracket_length, double target_length, double time, double target_force);

		double compute_time(double bracket_length, double target_length, double target_force, double thick);

		Eigen::Vector2d compute_time_range(double bracket_length, double target_length, double target_force);

		// Compute equilibrium with membrane: time, seconds
		double compute_mem_deformation(double thickness, double bracket_length, double time, Eigen::VectorXd& mem_poly, double init_mem_length);

		inline Eigen::Array3d& get_grid_min() { return grid_min; };
		inline Eigen::Array3d& get_grid_step() { return grid_step; };
		inline Eigen::Array3i& get_grid_num() { return grid_num; };
		//inline Eigen::Array3d get_grid_max() { return grid_min + grid_step * (grid_num - 1); };
		inline Eigen::MatrixXd& get_beta() { return beta; }
	};

	template<class Archive>
	void serialize(Archive& ar, smpup::BracketModel& bracket_model)
	{
		ar(cereal::make_nvp("grid_min", bracket_model.get_grid_min()));
		ar(cereal::make_nvp("grid_step", bracket_model.get_grid_step()));
		ar(cereal::make_nvp("grid_num", bracket_model.get_grid_num()));
		ar(cereal::make_nvp("beta", bracket_model.get_beta()));
	}

}

#endif
