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

#include "bracket_model.h"
#include "utils/umath.h"

namespace smpup {

	// Compute bracket energy polynomial (time, sec)
	Eigen::VectorXd BracketModel::make_energy_poly_coeff(double thickness, double bracket_length, double time) const
	{
		Eigen::VectorXd poly = Eigen::VectorXd::Zero(9);
		time /= 60;

		Eigen::Array3d point;
		point(GRID_THICK ) = thickness;
		point(GRID_LENGTH) = bracket_length;
		point(GRID_TIME  ) = time;

		Eigen::Array3i corner = ((point - grid_min) / grid_step).unaryExpr([](double x) { return std::floor(x); }).cast<int>(); //TODO: use floor() from newer Eigen
		//if ((corner < -1).any() || (corner >= grid_num ).any()) continue;

		// Stick to closest boundary of sampled space
		for (int i = 0; i < 3; ++i)
			if (corner(i) < 0)
			{
				point(i) = grid_min(i);
				corner(i) = 0;
			}
			else if (grid_num(i) <= corner(i))
			{
				point(i) = grid_min(i) + grid_step(i) * grid_num(i);
				corner(i) = grid_num(i) - 1;
			}

		Eigen::Array3d alpha = corner.cast<double>() + 1 - (point - grid_min) / grid_step;

		for (int i = 0; i < 3; ++i)
		{
			if (corner(i) < 0)
			{
				corner(i) += 1;
				alpha(i) += 1;
			}
			else if (corner(i) >= grid_num(i) - 1)
			{
				corner(i) -= 1;
				alpha(i) -= 1;
			}
		}

		Eigen::VectorXd betax = Eigen::VectorXd::Zero(beta.rows());
		for (int c0 = 0; c0 < 2; ++c0)
			for (int c1 = 0; c1 < 2; ++c1)
				for (int c2 = 0; c2 < 2; ++c2)
				{
					Eigen::Array3i c{ c0, c1, c2 };
					Eigen::Array3i cx = corner + c;
					int ind = (cx(0) * grid_num(1) + cx(1)) * grid_num(2) + cx(2);
					betax += (c.cast<double>() + (1 - 2 * c).cast<double>() * alpha).prod() * beta.col(ind);
				}
		poly.segment(2, betax.size()) = betax;

		// Correct stress-strain curve to be fixed in terms of strain (not displacement)
		if (point(1) != bracket_length)
		{
			double sc = point(1) / bracket_length;
			double scp = sc;
			for (int i = 2; i < poly.cols(); ++i)
			{
				poly(i) *= scp;
				scp *= sc;
			}
		}
		return poly;
	}

	Eigen::VectorXd BracketModel::make_force_poly_coeff(double thickness, double bracket_length, double time) const
	{
		Eigen::VectorXd poly;
		Eigen::VectorXd energy_poly = make_energy_poly_coeff(thickness, bracket_length, time);
		poly = energy_poly.array().tail(energy_poly.size() - 1) * natural_seriesd(1, energy_poly.size() - 1);
		return poly;
	}

	void BracketModel::make_energy_poly_coeff(const Eigen::VectorXd& thickness, const Eigen::VectorXd& bracket_length, double time, Eigen::MatrixXd& poly) const
	{
		int m = thickness.rows();
		poly = Eigen::MatrixXd::Zero(m, 9);
		for (int k = 0; k < thickness.size(); ++k)
			poly.row(k) = make_energy_poly_coeff(thickness(k), bracket_length(k), time);
	}

	double BracketModel::configure_thickness(double bracket_length, double target_length, double time, double target_force)
	{
		Eigen::VectorXd br_coef = Eigen::ArrayXd::Constant(8, bracket_length - target_length).pow(natural_seriesd(8));
		double force;

		double th_a = grid_min(GRID_THICK);
		force = br_coef.dot(make_force_poly_coeff(th_a, bracket_length, time));
		if (target_force < force) return th_a;

		double th_b = grid_min(GRID_THICK) + grid_step(GRID_THICK) * (grid_num(GRID_THICK) - 1);
		force = br_coef.dot(make_force_poly_coeff(th_b, bracket_length, time));
		if (force < target_force) return th_b;

		double th;
		// TODO: make a better estimation of precision, e.g. based on max thickness error or force error
		for (int i = 0; i < 20; ++i)
		{
			th = (th_a + th_b) / 2.0;
			force = br_coef.dot(make_force_poly_coeff(th, bracket_length, time));
			if (force < target_force) th_a = th;
			else th_b = th;
		}
		return th;
	}

	double BracketModel::compute_time(double bracket_length, double target_length, double target_force, double thick)
	{
		Eigen::VectorXd br_coef = Eigen::ArrayXd::Constant(8, bracket_length - target_length).pow(natural_seriesd(8));
		double force;
		double time_a = grid_min(GRID_TIME);
		force = br_coef.dot(make_force_poly_coeff(thick, bracket_length, time_a));
		if (force < target_force) return time_a;
		double time_b = (grid_min(GRID_TIME) + grid_step(GRID_TIME) * (grid_num(GRID_TIME) - 1)) * 60;
		force = br_coef.dot(make_force_poly_coeff(thick, bracket_length, time_b));
		if (target_force < force) return time_b;
		double time;
		// TODO: make a better estimation of precision, e.g. based on max thickness error or force error
		for (int i = 0; i < 20; ++i)
		{
			time = (time_a + time_b) / 2.0;
			force = br_coef.dot(make_force_poly_coeff(thick, bracket_length, time));
			if (target_force < force) time_a = time;
			else time_b = time;
		}
		return time;
	}

	Eigen::Vector2d BracketModel::compute_time_range(double bracket_length, double target_length, double target_force)
	{
		Eigen::Vector2d time_range;
		Eigen::VectorXd br_coef = Eigen::ArrayXd::Constant(8, bracket_length - target_length).pow(natural_seriesd(8));
		Eigen::Vector2d th_range = { grid_min(GRID_THICK), grid_min(GRID_THICK) + grid_step(GRID_THICK) * (grid_num(GRID_THICK) - 1) };
			
		for (int k = 0; k < 2; ++k)
		{
			double force;
			double time_a = grid_min(GRID_TIME);
			force = br_coef.dot(make_force_poly_coeff(th_range(k), bracket_length, time_a));
			if (force < target_force) time_range(k) = time_a;
			else {
				double time_b = (grid_min(GRID_TIME) + grid_step(GRID_TIME) * (grid_num(GRID_TIME) - 1)) * 60;
				force = br_coef.dot(make_force_poly_coeff(th_range(k), bracket_length, time_b));
				if (target_force < force) time_range(k) = time_b;
				else
				{
					double time;
					// TODO: make a better estimation of precision, e.g. based on max thickness error or force error
					for (int i = 0; i < 20; ++i)
					{
						time = (time_a + time_b) / 2.0;
						force = br_coef.dot(make_force_poly_coeff(th_range(k), bracket_length, time));
						if (target_force < force) time_a = time;
						else time_b = time;
					}
					time_range(k) = time;
				}
			}
		}

		return time_range;
	}

	double BracketModel::compute_mem_deformation(double thickness, double bracket_length, double time, Eigen::VectorXd& mem_poly, double init_mem_length)
	{
		auto br_poly = make_force_poly_coeff(thickness, bracket_length, time);

		double disp_a = 0;
		Eigen::VectorXd br_coef = Eigen::ArrayXd::Constant(8, disp_a).pow(natural_seriesd(8));
		double br_force = br_coef.dot(br_poly) * 4.0;
		Eigen::VectorXd mem_length_pow = Eigen::ArrayXd::Constant(2, disp_a - init_mem_length).pow(natural_seriesd(2));
		double mem_force = -mem_poly.dot(mem_length_pow); // negative since the force is along negative direction
		if (mem_force < br_force) return disp_a;

		double disp_b = init_mem_length;
		br_coef = Eigen::ArrayXd::Constant(8, disp_b).pow(natural_seriesd(8));
		br_force = br_coef.dot(br_poly) * 4.0;
		mem_length_pow = Eigen::ArrayXd::Constant(2, disp_b - init_mem_length).pow(natural_seriesd(2));
		mem_force = -mem_poly.dot(mem_length_pow); // negative since the force is along negative direction
		if (br_force < mem_force) return disp_b;

		double disp;
		// TODO: make a better estimation of precision, e.g. based on max thickness error or force error
		for (int i = 0; i < 20; ++i)
		{
			disp = (disp_a + disp_b) / 2.0;

			br_coef = Eigen::ArrayXd::Constant(8, disp).pow(natural_seriesd(8));
			br_force = br_coef.dot(br_poly) * 4.0;
			mem_length_pow = Eigen::ArrayXd::Constant(2, disp - init_mem_length).pow(natural_seriesd(2));
			mem_force = -mem_poly.dot(mem_length_pow); // negative since the force is along negative direction
			if (br_force < mem_force) disp_a = disp;
			else disp_b = disp;
		}

		return disp;
	}

}
