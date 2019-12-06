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

#ifndef SMP_ADDONS_H
#define SMP_ADDONS_H

#include "cereal/cereal.hpp"
#include "cereal/archives/json.hpp"
#include "cereal/types/vector.hpp"
#include "cereal/types/polymorphic.hpp"
#include "cereal/archives/portable_binary.hpp"
#include "cinder/gl/gl.h"
#include "Eigen/Eigen"

namespace cereal
{
	template<class Archive>
	void save(Archive& ar, std::set<int> const& setT)
	{
		ar(cereal::make_nvp("set", std::vector<int>(setT.begin(), setT.end())));
	}

	template<class Archive>
	void load(Archive& ar, std::set<int>& setT)
	{
		std::vector<int> vecT;
		ar(cereal::make_nvp("set", vecT));
		setT = std::set<int>(vecT.begin(), vecT.end());
	}
}

template<typename Derived>
cinder::vec3& operator <<(cinder::vec3& v, const Eigen::MatrixBase<Derived>& from)
{
	v.x = from(0);
	v.y = from(1);
	v.z = from(2);
	return v;
}

#endif