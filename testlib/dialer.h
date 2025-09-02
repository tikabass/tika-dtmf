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
 * @file 	dialer.h
 * @brief 	DTMF dialer for testing
 * @internal
 */

#ifndef TIKA_DIALER_H
#define TIKA_DIALER_H

#include "dtmf_generator.h"
#include <algorithm>
#include <cassert>
#include <string>
#include <string_view>

// *****************************************************************************

class dialer
{
  public:
	// -------------------------------------------------------------------------

	void dial(char c)
	{
		dial_string_ += c;
		dial_next();
	}

	// -------------------------------------------------------------------------

	void dial(const std::string_view& str)
	{
		dial_string_.append(str);
		dial_next();
	}

	// -------------------------------------------------------------------------

	unsigned estimate_dial_duration(std::string_view str) const noexcept
	{
		unsigned result = 0;
		while (!str.empty())
		{
			if (auto c = str.front(); is_valid_key(c))
			{
				result += unsigned(std::floor((keypress_duration_ * get_sample_rate()) + .5f)
				                   + std::floor((spacing_duration_ * get_sample_rate()) + .5f));
			}
			str.remove_prefix(1);
		}
		return result;
	}

	// -------------------------------------------------------------------------

	bool dialing() const noexcept
	{
		return keypress_counter_ != 0 || spacing_counter_ != 0;
	}

	// -------------------------------------------------------------------------

	float operator()() noexcept
	{
		if (keypress_counter_ && --keypress_counter_ == 0)
			release_key();

		if (spacing_counter_ && --spacing_counter_ == 0)
			dial_next();

		return generator_();
	}

	// -------------------------------------------------------------------------

	constexpr static bool is_valid_key(char c) noexcept
	{
		return valid_keys_.find(c) != std::string_view::npos;
	}

	// -------------------------------------------------------------------------

	void set_keypress_duration(float duration) noexcept
	{
		assert(duration > 0);
		keypress_duration_ = duration;
	}

	// -------------------------------------------------------------------------

	float get_keypress_duration() const noexcept
	{
		return keypress_duration_;
	}

	// -------------------------------------------------------------------------

	void set_spacing_duration(float duration) noexcept
	{
		assert(duration > 0);
		spacing_duration_ = duration;
	}

	// -------------------------------------------------------------------------

	float get_spacing_duration() const noexcept
	{
		return spacing_duration_;
	}

	// -------------------------------------------------------------------------

	void set_sample_rate(float freq) noexcept
	{
		generator_.set_sample_rate(freq);
	}

	// -------------------------------------------------------------------------

	float get_sample_rate() const noexcept
	{
		return generator_.get_sample_rate();
	}

	// -------------------------------------------------------------------------

	dtmf_generator& generator() noexcept
	{
		return generator_;
	}

	// -------------------------------------------------------------------------

	const dtmf_generator& generator() const noexcept
	{
		return generator_;
	}

  private:
	// -------------------------------------------------------------------------

	char get_next_key() noexcept
	{
		while (!dial_string_.empty())
		{
			auto c = dial_string_.front();
			dial_string_.erase(dial_string_.begin());

			if (is_valid_key(c))
				return c;
		}
		return 0;
	}

	// -------------------------------------------------------------------------

	void dial_next() noexcept
	{
		if (dialing())
			return;

		if (auto c = get_next_key())
			press_key(c);
	}

	// -------------------------------------------------------------------------

	void press_key(char c) noexcept
	{
		assert(keypress_counter_ == 0 && spacing_counter_ == 0);
		keypress_counter_ = unsigned(std::floor((keypress_duration_ * get_sample_rate()) + .5f));
		generator_.press_key(c);
	}

	// -------------------------------------------------------------------------

	void release_key() noexcept
	{
		assert(keypress_counter_ == 0 && spacing_counter_ == 0);
		spacing_counter_ = unsigned(std::floor((spacing_duration_ * get_sample_rate()) + .5f));
		generator_.release_key();
	}

	// =========================================================================
  private:
	dtmf_generator generator_;
	float keypress_duration_   = .050f;
	float spacing_duration_    = .100f;
	unsigned keypress_counter_ = 0;
	unsigned spacing_counter_  = 0;
	std::string dial_string_;
	static constexpr std::string_view valid_keys_ = ",123A456B789C*0#D";
};

#endif // ****************************************************************** EOF
