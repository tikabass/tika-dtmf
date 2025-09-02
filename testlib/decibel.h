// *****************************************************************************
// MIT License
//
// Copyright (c) 2025 tikabass <tika.devel@gmail.com>
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
// *****************************************************************************
/**
 * @file	decibel.h
 * 			Decibel calculus support 
 */

#ifndef TIKA_DECIBEL_H
#define TIKA_DECIBEL_H
#pragma once

#include <cmath>
#include <concepts>
#include <type_traits>

// *****************************************************************************
/**
 * @brief Magnitude decibel to linear factor. 
 * 
 * @tparam Float 	Any floating po,int type.
 * @param x 		Value in decibels
 * @return Float 	The corresponding magnitude gain factor.
 */

template <typename Float, typename R = std::enable_if<std::is_floating_point_v<Float>, Float>>
inline typename R::type dBm(Float x) noexcept
{
	return std::pow(Float(10), x / Float(10));
}

// *****************************************************************************

template <typename Float, typename R = std::enable_if<std::is_floating_point_v<Float>, Float>>
inline typename R::type to_dBm(Float x) noexcept
{
	return Float(10) * std::log10(x);
}

#endif // ****************************************************************** EOF
