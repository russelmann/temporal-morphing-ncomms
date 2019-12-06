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

#ifndef SMP_MATLAB_H
#define SMP_MATLAB_H

#include "smp_addons.h"
#include "Eigen/Eigen"

#include "MatlabDataArray.hpp"
#include "MatlabEngine.hpp"

namespace mateng = matlab::engine;
namespace matdata = matlab::data;

class Matlab
{
private:
	std::unique_ptr<mateng::MATLABEngine> matlab_ptr;
	matdata::ArrayFactory factory;

public:

	Matlab() { }

	// Start Matlab engine
	inline bool start()
    {
        try {
            matlab_ptr = mateng::startMATLAB();
        } catch (const std::exception& e) {
            return false;
        };
        return true;
    }

	// Get reference to Matlab engine
	inline mateng::MATLABEngine& engine() { return *matlab_ptr.get(); }

	// Create matlab array
	matdata::TypedArray<double> const createArray(const Eigen::MatrixXd& matrix);

	// Call command
	inline void eval(std::string command)
	{
		matlab_ptr->eval(mateng::convertUTF8StringToUTF16String(command));
	}

	// Put variable
	inline void put_variable(std::string name, const Eigen::MatrixXd& X)
	{
		auto argX = createArray(X);
		matlab_ptr->setVariable(name, argX);
	}

	// Call matlab figure
	void figure(const std::string& name, bool hold = true);

	// Call matlab figure
	void figure(const std::string& name, const std::string& title, const std::string& xlabel, const std::string& ylabel, bool hold = true);

	// Call matlab plot(X, Y) for Eigen vectors
	void plot(const Eigen::Ref<const Eigen::VectorXd>& X, const Eigen::Ref<const Eigen::VectorXd> Y, std::string LineSpec = "", Eigen::VectorXd color = Eigen::VectorXd(), double line_width = 0.5);

	// Call matlab plot(X) for Eigen vectors
	void plot(const Eigen::Ref<const Eigen::VectorXd>& X, std::string LineSpec = "", Eigen::VectorXd color = Eigen::VectorXd());

	// Call matlab hist(X) for Eigen vector
	void histogram(const Eigen::Ref<const Eigen::VectorXd>& X, int nbins = 0, Eigen::VectorXd color = Eigen::VectorXd());
};

#endif
