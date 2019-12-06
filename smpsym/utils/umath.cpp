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

#include "utils/umath.h"

Eigen::Matrix3d vectors2rotation(Eigen::Vector3d dirA, Eigen::Vector3d sideA, Eigen::Vector3d dirB, Eigen::Vector3d sideB)
{
	Eigen::Matrix3d R;

	dirA.normalize();
	sideA -= dirA.dot(sideA) * dirA;
	sideA.normalize();
	dirB.normalize();
	sideB -= dirB.dot(sideB) * dirB;
	sideB.normalize();

	Eigen::Matrix3d A, B;
	A.col(0) = dirA;
	A.col(1) = sideA;
	A.col(2) = dirA.cross(sideA);
	B.col(0) = dirB;
	B.col(1) = sideB;
	B.col(2) = dirB.cross(sideB);

	R = A * B.inverse();

	return R;
}

// Vector of numbers in range 0 .. size-1
template<typename _Scalar>
Eigen::Array<_Scalar, -1, 1> natural_series(int size)
{
	Eigen::Array<_Scalar, -1, 1> vec(size);
	for (int i = 0; i < vec.size(); ++i)
		vec(i) = i;
	return vec;
}

// Vector of numbers in range a .. b
template<typename _Scalar>
Eigen::Array<_Scalar, -1, 1> natural_series(int a, int b)
{
	Eigen::Array<_Scalar, -1, 1> vec;
	vec.setConstant(std::abs(b - a) + 1, a);
	for (int i = 0; i < vec.size(); ++i)
		vec(i) += (b < a) ? -i : i;
	return vec;
}

Eigen::ArrayXi natural_seriesi(int size) { return natural_series<int>(size); }
Eigen::ArrayXi natural_seriesi(int a, int b) { return natural_series<int>(a, b); }
Eigen::ArrayXd natural_seriesd(int size) { return natural_series<double>(size); }
Eigen::ArrayXd natural_seriesd(int a, int b) { return natural_series<double>(a, b); }

// Assumes conformality, no conformality check
Eigen::VectorXd compute_scaling_factors(const Eigen::MatrixXd& V0, const Eigen::MatrixXd& V1, const Eigen::MatrixXi& F)
{
	//TODO: test this function
	Eigen::VectorXd u(V0.rows());
	for (int i = 0; i < F.rows(); ++i)
	{
		Eigen::Vector3d len0, len1;
		for (int j = 0; j < 3; ++j)
		{
			len0(j) = (V0.row(F(i, (j + 1) % 3)) - V0.row(F(i, (j + 2) % 3))).norm();
			len1(j) = (V1.row(F(i, (j + 1) % 3)) - V1.row(F(i, (j + 2) % 3))).norm();
		}
		for (int j = 0; j < 3; ++j)
		{
			u(F(i, j)) = len0(j) / len1(j) * len1((j + 1) % 3) / len0((j + 1) % 3) * len1((j + 2) % 3) / len0((j + 2) % 3);
		}
	}
	u = u.array().log();
	return u;
}