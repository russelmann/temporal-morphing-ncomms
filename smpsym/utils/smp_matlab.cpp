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

#include "smp_matlab.h"

typedef unsigned long long ullong;

matdata::TypedArray<double> const Matlab::createArray(const Eigen::MatrixXd& matrix)
{
	matdata::TypedArray<double> const arr =
		factory.createArray({ ullong(matrix.rows()), ullong(matrix.cols()) }, matrix.data(), matrix.data() + matrix.size());
	return arr;
}

void Matlab::figure(const std::string& name, bool hold)
{
	const std::shared_ptr<matlab::engine::StreamBuffer> output;
	const std::shared_ptr<matlab::engine::StreamBuffer> error;

	try {
		matlab_ptr->eval(mateng::convertUTF8StringToUTF16String("figure(" + name + ")"), output, error);
	}
	catch (const std::exception& ex) {
		matlab_ptr->eval(mateng::convertUTF8StringToUTF16String(name + " = figure()"), output, error);
	}

	if (hold) matlab_ptr->eval(mateng::convertUTF8StringToUTF16String("hold on"), output, error);
}

void Matlab::figure(const std::string& name, const std::string& title, const std::string& xlabel, const std::string& ylabel, bool hold)
{
	figure(name, hold);

	const std::shared_ptr<matlab::engine::StreamBuffer> output;
	const std::shared_ptr<matlab::engine::StreamBuffer> error;

	if (!name.empty()) matlab_ptr->eval(mateng::convertUTF8StringToUTF16String(name + ".Name = '" + title + "'"), output, error);
	matlab_ptr->eval(mateng::convertUTF8StringToUTF16String("title('" + title + "')"), output, error);
	matlab_ptr->eval(mateng::convertUTF8StringToUTF16String("xlabel('" + xlabel + "')"), output, error);
	matlab_ptr->eval(mateng::convertUTF8StringToUTF16String("ylabel('" + ylabel + "')"), output, error);
	matlab_ptr->eval(mateng::convertUTF8StringToUTF16String("hold on"), output, error);
}

void Matlab::plot(const Eigen::Ref<const Eigen::VectorXd>& X, const Eigen::Ref<const Eigen::VectorXd> Y, std::string LineSpec, Eigen::VectorXd color, double line_width)
{
	std::vector<matlab::data::Array> args;
	if (X.size() > 0) args.emplace_back(createArray(X));
	args.emplace_back(createArray(Y));
	args.emplace_back(factory.createCharArray(LineSpec));
	if (color.size() > 0) args.emplace_back(factory.createCharArray("Color"));
	if (color.size() > 0) args.emplace_back(createArray(color));
	args.emplace_back(factory.createCharArray("LineWidth"));
	args.emplace_back(createArray(Eigen::Matrix<double, 1, 1>{line_width}));
	matlab_ptr->feval("plot", args);
}

void Matlab::plot(const Eigen::Ref<const Eigen::VectorXd>& Y, std::string LineSpec, Eigen::VectorXd color)
{
	plot(Eigen::VectorXd(), Y, LineSpec);
}

void Matlab::histogram(const Eigen::Ref<const Eigen::VectorXd>& X, int nbins, Eigen::VectorXd color)
{
	std::vector<matlab::data::Array> args;
	args.emplace_back(createArray(X));
	if (nbins > 0) args.emplace_back(factory.createScalar(nbins));
	if (color.size() > 0)
	{
		args.emplace_back(factory.createCharArray("FaceColor"));
		args.emplace_back(createArray(color));
	}
	matlab_ptr->feval("histogram", args);
}