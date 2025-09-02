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
 * 	@file dtmf_generator.h
 *
 * 			Definition of class dtmf_generator.
 */
// *****************************************************************************

#ifndef TIKA_DTMF_GENERATOR_H
#define TIKA_DTMF_GENERATOR_H
#pragma once

#include "decibel.h"
#include "sine_generator.h"
#include <algorithm>

// *****************************************************************************
/**
 * 	Generates DTMF tones.
 */
class dtmf_generator
{
  public:
	// -------------------------------------------------------------------------
	/**
	 *
	 */

	void press_key(char key) noexcept
	{
		auto [row, col] = key_to_grid_pos(key);
		row_            = row;
		col_            = col;
		if (row_ >= 0 && col_ >= 0)
		{
			row_sine_.set_frequency(row_frequency_[row] / sample_rate_);
			col_sine_.set_frequency(col_frequency_[col] / sample_rate_);
			current_key_ = key;
		}
		else
		{
			row_sine_.set_frequency(0);
			col_sine_.set_frequency(0);
			row_sine_.set_phase(0);
			col_sine_.set_phase(0);
			current_key_ = -1;
		}
	}

	// -------------------------------------------------------------------------

	void release_key() noexcept
	{
		press_key(-1);
	}

	// -------------------------------------------------------------------------
	/**
	 * 	Generates audio output.
	 *
	 * 	@returns	One sample of audio in normalized floating point format.
	 */

	float operator()() noexcept
	{
		if (current_key_ > 0)
		{
			const auto x = row_sine_();
			const auto y = col_sine_();
			return (x + (row_col_diff_ * y)) * gain_;
		}
		else
		{
			return 0;
		}
	}

	// -------------------------------------------------------------------------

	void set_gain(float gain) noexcept
	{
		gain_ = gain;
	}

	// -------------------------------------------------------------------------

	void set_gain_dB(float gain) noexcept
	{
		set_gain(dBm(gain_ * 2));
	}

	// -------------------------------------------------------------------------

	float get_gain() const noexcept
	{
		return gain_;
	}

	// -------------------------------------------------------------------------

	float get_gain_dB() const noexcept
	{
		return to_dBm(std::sqrt(gain_));
	}

	// -------------------------------------------------------------------------

	void set_power_difference(float gain) noexcept
	{
		row_col_diff_ = gain;
	}

	// -------------------------------------------------------------------------

	void set_power_difference_dB(float gain) noexcept
	{
		row_col_diff_ = dBm(gain);
	}

	// -------------------------------------------------------------------------

	float get_power_difference() const noexcept
	{
		return row_col_diff_;
	}

	// -------------------------------------------------------------------------

	float get_power_difference_dB() const noexcept
	{
		return to_dBm(get_power_difference());
	}

	// -------------------------------------------------------------------------
	/**
	 * 	Changes sample rate.  Th enew sample rate will be effective after the next
	 *  call to key_press().
	 *
	 * 	@param	[in] freq	New frequency, in Hertz.
	 */

	void set_sample_rate(float freq) noexcept
	{
		sample_rate_ = freq;
	}

	// -------------------------------------------------------------------------

	float get_sample_rate() const noexcept
	{
		return sample_rate_;
	}

	// -------------------------------------------------------------------------

	static std::pair<size_t, size_t> key_to_grid_pos(char key) noexcept
	{
		switch (key)
		{
		case '1': return {0, 0};
		case '2': return {0, 1};
		case '3': return {0, 2};
		case 'A': return {0, 3};

		case '4': return {1, 0};
		case '5': return {1, 1};
		case '6': return {1, 2};
		case 'B': return {1, 3};

		case '7': return {2, 0};
		case '8': return {2, 1};
		case '9': return {2, 2};
		case 'C': return {2, 3};

		case '*': return {3, 0};
		case '0': return {3, 1};
		case '#': return {3, 2};
		case 'D': return {3, 3};
		}
		return {-1, -1};
	}

  private:
	// =========================================================================
	//	signal generation

	sine_generator<float> col_sine_; ///< column sinusoidal generator.
	sine_generator<float> row_sine_; ///< row sinusoidal generator.

	float sample_rate_  = 8000.f;
	float row_col_diff_ = dBm(0.0);
	float gain_         = dBm(-6.0);
	char current_key_   = -1;
	int row_            = -1;
	int col_            = -1;

	// =========================================================================
	//	constants

	/**
	 * 	Frequencies for the 4 telephone keypad rows, in Hertz.
	 */
	static constexpr float row_frequency_[4] = {697, 770, 852, 941};

	/**
	 * 	Frequencies for the 4 telephone keypad columns, in Hertz.
	 */
	static constexpr float col_frequency_[4] = {1209, 1336, 1477, 1633};
};

#endif // ****************************************************************** EOF
