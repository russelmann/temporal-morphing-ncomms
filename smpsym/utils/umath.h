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

#ifndef UMATH_H
#define UMATH_H

#define _USE_MATH_DEFINES
#include <cmath>
#include "smp_addons.h"
#include "Eigen/Eigen"

inline int next3(int i)
{
	if (i == 0) return 1;
	if (i == 1) return 2;
	return 0;
}

inline int prev3(int i)
{
	if (i == 0) return 2;
	if (i == 1) return 0;
	return 1;
}

inline double angle_vectors(Eigen::Vector3d a, Eigen::Vector3d b)
{
	return atan2(a.cross(b).norm(), a.dot(b));
}

Eigen::Matrix3d vectors2rotation(Eigen::Vector3d dirA, Eigen::Vector3d sideA, Eigen::Vector3d dirB, Eigen::Vector3d sideB);

// Vector of numbers in range 0 .. size-1
template<typename _Scalar>
Eigen::Array<_Scalar, -1, 1> natural_series(int size);

// Vector of numbers in range a .. b
template<typename _Scalar>
Eigen::Array<_Scalar, -1, 1> natural_series(int a, int b);

Eigen::ArrayXi natural_seriesi(int size);
Eigen::ArrayXi natural_seriesi(int a, int b);
Eigen::ArrayXd natural_seriesd(int size);
Eigen::ArrayXd natural_seriesd(int a, int b);

// Assumes conformality, no conformality check
Eigen::VectorXd compute_scaling_factors(const Eigen::MatrixXd& V0, const Eigen::MatrixXd& V1, const Eigen::MatrixXi& F);

#endif