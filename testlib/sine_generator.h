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
 * @file 	sine_generator.h
 * @brief 	Sinusoid signal generator
 * @internal
 */

#ifndef TIKA_SINE_GENERATOR_H
#define TIKA_SINE_GENERATOR_H
#pragma once

#include <cassert>
#include <cmath>
#include <concepts>
#include <numbers>
#include <type_traits>

// *****************************************************************************

template <typename Float>
class sine_generator
{
  public:
	using float_t = typename std::enable_if<std::is_floating_point_v<Float>, Float>::type;

	// -------------------------------------------------------------------------
  public:
	sine_generator() noexcept = default;
	sine_generator(float_t freq) noexcept
	{
		assert(0 <= freq && freq <= .5);
		set_frequency(freq);
	}

	// -------------------------------------------------------------------------

	float_t operator()() noexcept
	{
		auto sample = std::sin(phase_);
		phase_ += (step_ * (1 + detune_));
		limit_phase();
		return sample;
	}

	// -------------------------------------------------------------------------

	void detune(float_t ratio) noexcept
	{
		detune_ = ratio;
	}

	// -------------------------------------------------------------------------

	void phase_offset(float_t ofst) noexcept
	{
		phase_ += ofst;
		limit_phase();
	}

	// -------------------------------------------------------------------------

	void set_frequency(float_t freq) noexcept
	{
		assert(0 <= freq && freq <= .5);
		step_ = two_pi * freq;
	}

	// -------------------------------------------------------------------------

	void set_phase(float_t phase) noexcept
	{
		assert(-two_pi <= phase && phase <= two_pi);
		phase_ = phase;
		limit_phase();
	}

	// -------------------------------------------------------------------------
  private:
	void limit_phase() noexcept
	{
		auto turns = std::floor(phase_ * one_over_two_pi);
		phase_ -= turns * two_pi;
	}

	// -------------------------------------------------------------------------
  private:
	float_t phase_  = {};
	float_t step_   = {};
	float_t detune_ = {}; // proportional detune offset.

	static constexpr auto pi              = float_t(3.141592653589793238462643383279L);
	static constexpr auto two_pi          = pi + pi;
	static constexpr auto one_over_two_pi = float_t(1) / two_pi;
};

#endif // ****************************************************************** EOF
