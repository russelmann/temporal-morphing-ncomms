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

#ifndef SELF_SERIALIZER_H
#define SELF_SERIALIZER_H

#include "smp_addons.h"

#ifdef USE_BOOST
#include "boost/filesystem.hpp"
namespace filesystem = boost::filesystem;
#else
#include <filesystem>
namespace filesystem = std::experimental::filesystem;
#endif

namespace cereal
{
	template<typename T>
	class SelfSerializer
	{
	private:
		std::string name;
		filesystem::path default_dir;

		inline T& child() { return *static_cast<T*>(this); }

	public:
		SelfSerializer(std::string name, filesystem::path default_dir = filesystem::path())
			: name(name), default_dir(default_dir) { }

		void serialization(filesystem::path path = "")
		{
			if (path.empty()) path = default_dir;
			//T* obj = static_cast<T*>(this);
			if (!load(path)) save(path);
		}

		bool load(filesystem::path path = "")
		{
			if (path.empty()) path = default_dir;
			std::ifstream ss(path.append(name + ".json").c_str());
            if (ss.good())
			{
				cereal::JSONInputArchive archive(ss);
				try {
					archive(cereal::make_nvp(name, child()));
				}
				catch (const std::exception&) {
					return false;
				}
                return true;
			}
            return false;
		}

		bool save(filesystem::path path = "")
		{
			if (path.empty()) path = default_dir;
			std::ofstream ss(path.append(name + ".json").c_str());
            if (ss.good())
			{
				cereal::JSONOutputArchive archive(ss);
				archive(cereal::make_nvp(name, child()));
                return true;
			}
            return false;
		}

		inline void set_default_directory(filesystem::path path) { default_dir = path; }
		inline void set_default_directory(std::string path) { default_dir = path; }

		inline filesystem::path get_default_directory() { return default_dir; }

	};

}

#endif
