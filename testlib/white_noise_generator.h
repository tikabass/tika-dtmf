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
 * @file 	white_noise_generator.h
 * @brief 	White noise generator
 * @internal
 */

#ifndef TIKA_WHITE_NOISE_GENERATOR_H
#define TIKA_WHITE_NOISE_GENERATOR_H
#pragma once

#include "decibel.h"
#include <random>
#include <type_traits>

// *****************************************************************************
/**
 * 	Normal distribution gaussian white noise generator.
 * 
 * 	@note  Defaults to -20dB power output.
 * 
 *  @tparam Float  A floating point type.
 */

template <typename Float>
class white_noise_generator
{
  public:
	using float_t = typename std::enable_if<std::is_floating_point_v<Float>, Float>::type;

  public:
	// -------------------------------------------------------------------------
	/**
	 * @brief 	Betrieves the next sample from the generator
	 * 
	 * @return A single float_t sample
	 */

	float_t operator()() noexcept
	{
		return gain_ * rnd_(gen_);
	}

	// -------------------------------------------------------------------------
	/**
	 * @brief Set the power output ofthe generator in dB.
	 * 
	 * 	The amplitude gain is adjusted to @p pwr / 2.
	 * 
	 * @param [in] pwr 	The desired power level of the generated noise signal
	 */

	void set_power_dB(Float pwr) noexcept
	{
		gain_ = dBm(pwr / 2);
	}

	// -------------------------------------------------------------------------
	/**
	 * @brief Sets the magnitude of the output signal.
	 * 
	 * @param amplitude Desired amplitude.
	 */

	void set_gain(Float amplitude) noexcept
	{
		gain_ = amplitude;
	}

	// -------------------------------------------------------------------------
	/**
	 * @brief Get the power output of the generated signal 
	 * 
	 * @return Power output in dB
	 */

	float get_power_dB() const noexcept
	{
		return to_dBm(gain_ * gain_);
	}

	// -------------------------------------------------------------------------
	/**
	 * @brief Get the amplitude of the generated signal 
	 * 
	 * @return Amplitude.
	 */

	float get_gain() const noexcept
	{
		return gain_;
	}

	// *************************************************************************
  private:
	float_t gain_ = dBm(float_t(-20 / 2));
	std::random_device dev_;
	std::mt19937 gen_{dev_()};
	std::normal_distribution<float_t> rnd_{float_t(0), float_t(1)};
};

#endif // ****************************************************************** EOF
