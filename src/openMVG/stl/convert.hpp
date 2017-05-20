// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_STL_CONVERT_HPP
#define OPENMVG_STL_CONVERT_HPP

#include <sstream>
#include <string>
#include <vector>

namespace stl
{
	template <class T>
	std::string convert_to_string(T value)
	{
		  std::stringstream ss;
		  ss << value;
		  return ss.str();
	}
} // namespace stl
#endif // OPENMVG_STL_SPLIT_HPP
