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
 * @file dtmf_decoder.h
 * @brief 	Implementation of class tika::dtmf_decoder
 * 
 *	The class and all of its requirements are included in this one file.
 */
/**
 * @mainpage DTMF Decoder
 * 
 * This library implements a DTMF decoder using a modified Goertzel's algorithm 
 * to detect frequencies from keypad input. The performance target is to match 
 * ITU Q.24 requirements. 
 *
 * The library is implemented as a single header file, with no external 
 * dependencies other than the stahndard C++ library.  It requires c++17 or later 
 * to compile. It doesn't allocate memory nor raise eny exceptions.
 * 
 * @section Usage
 * 
 * @code
 * 
 * #define TIKA_SIMD_LEVEL TIKA_SIMD_LEVEL_AVX2	  // this is the default, and can be omitted.
 * #include "dtmf_decoder.h"
 * 
 * struct channel
 * {
 * 		tika::dtmf_decoder dtmf_filter_;
 * 		// ...
 * 
 * 		void on_key_down(char key) {
 * 			// key was pressed
 * 			// process key...
 * 		}
 *
 * 		 void on_key_up(char key) {
 * 			// key was released
 * 			// process key...
 * 		}
 *
 *      void setup() {
 * 			// connect to key events	
 *          dtmf_filter_.on_key_down([this](char key) { this->on_dtmf_key_down(key); });
 * 			dtmf_filter_.on_key_up([this](char key) { this->on_dtmf_key_up(key); });
 * 			// ...
 * 		}
 * 
 * 		// this function is called when audio data is available.  Audio samples can be 
 * 		// either uncompressed PCM 16-bit signed integers or 32 bit floats.
 * 		void on_audio_data(size_t sample_count, const uint16_t* samples) {
 * 			// feed audio data to the filter
 * 			dtmf_filter_(sample_count, samples);
 * 			// do other processing...
 * 		}
 * 
 * 		// ...
 * };
 * 
 * @endcode
 * 
 * @section SIMD CPU acceleration
 * <br>
 * This library uses SIMD CPU instructions to accelerate processing on x86/x64 CPUs.
 * The default is to use AVX2 instructions.  If you want to limit the instruction set
 * to a lower level, define the macro TIKA_SIMD_LEVEL to one of the following values:
 * <br>
 * 	- TIKA_SIMD_LEVEL_AVX2     (default) Use SSE2, SSE4.1, FMA3 and AVX2 instructions. 
 * 	- TIKA_SIMD_LEVEL_FMA3     Use SSE2, SSE4.1 and FMA3 instructions.
 * 	- TIKA_SIMD_LEVEL_SSE41    Use SSE2 and SSE4.1 instructions.
 * 	- TIKA_SIMD_LEVEL_SSE2	 Use SSE2 instructions.
 * 	- TIKA_SIMD_LEVEL_NO_SIMD	 Disable SIMD altogether.
 * 
 * 	Since the library is destined to be used in real-time applications, run-time host 
 *  SIMD support detection is not provided.
 * 
 * @section References
 *	 <br>
 * 	 <A HREF="https://www.itu.int/rec/dologin_pub.asp?lang=e&id=T-REC-Q.24-198811-I!!PDF-E&type=items">
 *   	ITU Q.24 Specification</A><br>
 *	 <A HREF=http://www.ti.com/lit/an/spra066/spra066.pdf>
 *		Modified Goertzel Algorithm in DTMF Detection Using the TMS320C80</A>
 *
 * @section License
 *	 <br>
 * 	This library is licenced under the MIT License, see the LICENSE file for details.
 * 
 */

#ifndef TIKA_DTMF_DECODER_H
#define TIKA_AUDIO_DTMF_DECODER_H
#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>

// *****************************************************************************
// Compilation opetions

#define TIKA_SIMD_LEVEL_NO_SIMD 0
#define TIKA_SIMD_LEVEL_SSE2 1
#define TIKA_SIMD_LEVEL_SSE41 2
#define TIKA_SIMD_LEVEL_FMA3 3
#define TIKA_SIMD_LEVEL_AVX2 4

#ifndef TIKA_SIMD_LEVEL
#define TIKA_SIMD_LEVEL TIKA_SIMD_LEVEL_AVX2
#endif

#if (TIKA_SIMD_LEVEL < TIKA_SIMD_LEVEL_NO_SIMD || TIKA_SIMD_LEVEL_AVX2 < TIKA_SIMD_LEVEL)
#error \
    "TIKA_SIMD_LEVEL must be one of TIKA_SIMD_LEVEL_NO_SIMD, TIKA_SIMD_LEVEL_SSE2, TIKA_SIMD_LEVEL_SSE41, TIKA_SIMD_LEVEL_FMA3, TIKA_SIMD_LEVEL_AVX2"
#endif

#if (TIKA_SIMD_LEVEL > TIKA_SIMD_LEVEL_NO_SIMD)
#include <immintrin.h>
#endif

// *****************************************************************************

namespace tika {

namespace details {

// *****************************************************************************
/**
 * @brief   Transforms 16-bit signed audio PCM data to normalised audio float.
 * 
 * @param [in] src      Pointer to 16-bit signed PCM data to transform
 * @param [in] count    Number of samples to rtransform
 * @param [out] dst     Pointer to space to store the result.  Must be large 
 *                      enough to store @p count 32-bit floats.
 */

inline void pcm_to_float(const int16_t* src, size_t count, float* dst) noexcept
{
	constexpr float SCALING_FACTOR = 1.f / 0x8000;

#if (TIKA_SIMD_LEVEL > TIKA_SIMD_LEVEL_NO_SIMD)

#if (TIKA_SIMD_LEVEL == TIKA_SIMD_LEVEL_AVX2)
	const auto FACTOR = _mm256_set1_ps(SCALING_FACTOR);

	// -------------------------------------------------------------------------

	const auto do_8 = [&FACTOR](const short* s, float* d) noexcept {
		auto in_i  = _mm256_cvtepi16_epi32(_mm_lddqu_si128((const __m128i*)s));
		auto in_ps = _mm256_cvtepi32_ps(in_i);
		auto out   = _mm256_mul_ps(in_ps, FACTOR);
		_mm256_storeu_ps(d, out);
	};

	// -------------------------------------------------------------------------
	// fallback to SSE4.1

	const auto do_4 = [&FACTOR](const int16_t* s, float* d) noexcept {
		auto in_i  = _mm_cvtepi16_epi32(_mm_loadu_si64((const __m128i*)s));
		auto in_ps = _mm_cvtepi32_ps(in_i);
		auto out   = _mm_mul_ps(in_ps, _mm256_castps256_ps128(FACTOR));
		_mm_storeu_ps(d, out);
	};

#elif (TIKA_SIMD_LEVEL >= TIKA_SIMD_LEVEL_SSE41)

	// -------------------------------------------------------------------------

	const auto FACTOR = _mm_set1_ps(SCALING_FACTOR);

	// -------------------------------------------------------------------------

	const auto do_8 = [&FACTOR](const int16_t* s, float* d) noexcept {
		const auto lo_i  = _mm_cvtepi16_epi32(_mm_loadu_si64((const __m128i*)s));
		const auto hi_i  = _mm_cvtepi16_epi32(_mm_loadu_si64((const __m128i*)(s + 4)));
		const auto lo_ps = _mm_cvtepi32_ps(lo_i);
		const auto hi_ps = _mm_cvtepi32_ps(hi_i);
		const auto lo    = _mm_mul_ps(lo_ps, FACTOR);
		const auto hi    = _mm_mul_ps(hi_ps, FACTOR);
		_mm_storeu_ps(d, lo);
		_mm_storeu_ps(d + 4, hi);
	};

	// -------------------------------------------------------------------------

	const auto do_4 = [&FACTOR](const int16_t* s, float* d) noexcept {
		auto in_i  = _mm_cvtepi16_epi32(_mm_loadu_si64((const __m128i*)s));
		auto in_ps = _mm_cvtepi32_ps(in_i);
		auto out   = _mm_mul_ps(in_ps, FACTOR);
		_mm_storeu_ps(d, out);
	};

#elif (TIKA_SIMD_LEVEL >= TIKA_SIMD_LEVEL_SSE2)
	// -------------------------------------------------------------------------

	const auto FACTOR = _mm_set1_ps(SCALING_FACTOR);

	// -------------------------------------------------------------------------

	const auto do_8 = [&FACTOR](const int16_t* s, float* d) noexcept {
		const auto lo_src = _mm_loadu_si64((const __m128i*)s);
		const auto hi_src = _mm_loadu_si64((const __m128i*)(s + 4));
		const auto lo_i   = _mm_srai_epi32(_mm_unpacklo_epi16(lo_src, lo_src), 16);
		const auto hi_i   = _mm_srai_epi32(_mm_unpacklo_epi16(hi_src, hi_src), 16);
		const auto lo_ps  = _mm_cvtepi32_ps(lo_i);
		const auto hi_ps  = _mm_cvtepi32_ps(hi_i);
		const auto lo     = _mm_mul_ps(lo_ps, FACTOR);
		const auto hi     = _mm_mul_ps(hi_ps, FACTOR);
		_mm_storeu_ps(d, lo);
		_mm_storeu_ps(d + 4, hi);
	};

	// -------------------------------------------------------------------------

	const auto do_4 = [&FACTOR](const int16_t* s, float* d) noexcept {
		const auto in_src = _mm_loadu_si64((const __m128i*)s);
		const auto in_i   = _mm_srai_epi32(_mm_unpacklo_epi16(in_src, in_src), 16);
		const auto in_ps  = _mm_cvtepi32_ps(in_i);
		const auto out    = _mm_mul_ps(in_ps, FACTOR);
		_mm_storeu_ps(d, out);
	};
#else
#error "Should not reach"
#endif
	// -------------------------------------------------------------------------

	while (count >= 32)
	{
		do_8(src, dst);
		do_8(src + 8, dst + 8);
		do_8(src + 16, dst + 16);
		do_8(src + 24, dst + 24);

		src += 32;
		dst += 32;
		count -= 32;
	}

	if (count >= 16)
	{
		do_8(src, dst);
		do_8(src + 8, dst + 8);

		src += 16;
		dst += 16;
		count -= 16;
	}

	if (count >= 8)
	{
		do_8(src, dst);

		src += 8;
		dst += 8;
		count -= 8;
	}

	if (count >= 4)
	{
		do_4(src, dst);

		src += 4;
		dst += 4;
		count -= 4;
	}

	if (count)
	{
		const auto N = std::min(count, size_t{4});

		alignas(sizeof(__m128i)) int16_t bufin[4]; // uninitialized OK
		alignas(sizeof(__m128)) float bufout[4];   // uninitialized OK

		std::copy_n(src, N, bufin);
		do_4(bufin, bufout);
		std::copy_n(bufout, N, dst);
	}

#else
	std::transform(src, src + count, dst, [](int16_t s) noexcept -> float { return s * SCALING_FACTOR; });
#endif
}

// *****************************************************************************
/**
 * @brief Get the goertzel coefficient for a specific frequency.
 * 
 * @note This function returns a value for a specific freqyuency.  To get the 
 *        Coefficient for a specific bin#, use the folloowing formula: <br>
 * 		  coef = goertzel_coefficient(float(bin#) / float(DFT_SIZE));
 *
 * @param [in] freq         The normalized frequency of insterest, in the range [0, .5].
 * @return The corresponding coefficient.
 */

inline float goertzel_coefficient(float freq) noexcept
{
	assert(.0f <= freq && freq <= .5f);
	constexpr auto PI = 3.141592653589793238462643383279f;
	return 2.f * std::cos(2 * PI * freq);
}

// *****************************************************************************
/**
 *  Computes power levels for a bank of Goertzel DFT filters.
 *
 *	On x86, this function uses SIMD for better perfomance.
 *
 *	@param [in]  channel_count 		Number of channels.
 *	@param [in]  goertzel_coefs 	Goertzel coefficients, one float per channel.
 *	@param [in]  sample_count 		Number of samples to process, this corresponds to the DFT size.
 *	@param [in]  samples_in 		Input audio samples.
 *	@param [out] power_out 			Pointer to array to store power outputs, must have room for 
 *                                  @p count floats
 *
 */

inline void compute_goertzel_power(size_t channel_count,
                                   const float* goertzel_coefs,
                                   size_t sample_count,
                                   const float* samples_in,
                                   float* power_out) noexcept
{
#if (TIKA_SIMD_LEVEL >= TIKA_SIMD_LEVEL_SSE2)
	// -------------------------------------------------------------------------
	//	processes sample_count samples x 4 channels

	auto do_4 =
	    [](size_t channel_count, const float* goertzel_coefs, size_t sample_count, const float* samples_in, float* power_out) {
		    const auto coefficients = _mm_loadu_ps(goertzel_coefs);
		    auto u0                 = _mm_setzero_ps();
		    auto u1                 = u0;

		    // -----------------------------------------------------------------
		    // processes 1 sample x 4 channels.

		    auto do_1 = [&coefficients, &u0, &u1](__m128 sample) {
#if (TIKA_SIMD_LEVEL >= TIKA_SIMD_LEVEL_FMA3)
			    const auto temp = _mm_fmadd_ps(u0, coefficients, _mm_sub_ps(sample, u1));
#else
			    const auto temp = _mm_add_ps(_mm_mul_ps(u0, coefficients), _mm_sub_ps(sample, u1));
#endif
			    u1 = u0;
			    u0 = temp;
		    };

		    // -----------------------------------------------------------------
		    // loop through all samples.

		    while (sample_count >= 4)
		    {
			    const auto in = _mm_loadu_ps(samples_in);
			    do_1(_mm_shuffle_ps(in, in, _MM_SHUFFLE(0, 0, 0, 0)));
			    do_1(_mm_shuffle_ps(in, in, _MM_SHUFFLE(1, 1, 1, 1)));
			    do_1(_mm_shuffle_ps(in, in, _MM_SHUFFLE(2, 2, 2, 2)));
			    do_1(_mm_shuffle_ps(in, in, _MM_SHUFFLE(3, 3, 3, 3)));
			    samples_in += 4;
			    sample_count -= 4;
		    }

		    while (sample_count > 0)
		    {
			    do_1(_mm_load_ps1(samples_in));
			    ++samples_in;
			    --sample_count;
		    }

		//  power = (u0² + u1²) - (coeff * u0 * u1)
#if (TIKA_SIMD_LEVEL >= TIKA_SIMD_LEVEL_FMA3)
		    const auto t = _mm_fmadd_ps(u0, u0, _mm_mul_ps(u1, u1));
		    u0           = _mm_mul_ps(coefficients, _mm_mul_ps(u0, u1));
#else
		    const auto t = _mm_add_ps(_mm_mul_ps(u0, u0), _mm_mul_ps(u1, u1));
		    u0           = _mm_mul_ps(coefficients, _mm_mul_ps(u0, u1));
#endif
		    _mm_storeu_ps(power_out, _mm_sub_ps(t, u0));
	    };

	// -------------------------------------------------------------------------
	// loop throuh all channels.

	while (channel_count >= 4)
	{
		do_4(channel_count, goertzel_coefs, sample_count, samples_in, power_out);
		goertzel_coefs += 4;
		power_out += 4;
		channel_count -= 4;
	}

	if (channel_count)
	{
		const size_t N = std::min(channel_count, size_t(4));
		alignas(sizeof(__m128)) float c[4]; // uninitialized OK
		alignas(sizeof(__m128)) float r[4]; // uninitialized OK

		std::copy_n(goertzel_coefs, N, c);
		do_4(channel_count, c, sample_count, samples_in, r);
		std::copy_n(r, N, power_out);
	}
#else
	while (channel_count--)
	{
		float u[2] = {};
		auto count = sample_count;
		auto smp   = samples_in;
		auto coef  = *(goertzel_coefs++);
		while (count--)
		{
			auto t = *(smp++) + (coef * u[0]) - u[1];
			u[1]   = u[0];
			u[0]   = t;
		}
		*(power_out++) = (u[0] * u[0]) + (u[1] * u[1]) - (coef * u[0] * u[1]);
	}
#endif
}

} // namespace details

// *****************************************************************************
/**
 *	DTMF decoder, uses a modified Goertzel's algorithm to detect frequencies from keypad input.
 *
 *	The performance target is to match the AT&T detection performance, although the +-3% frequency offset
 *	for non-detection is not met at this time.  This constraint is non-critical for the original application.
 *
 *	Usage:
 *		@li Create a dtmf_decoder object.
 *		@li Connect to asynchronous events on_key_down() and on_key_up()
 *		@li Feed input stream through the process() command.
 *
 *	OR:
 *		@li Create a dtmf_decoder object.
 *		@li get the processing block size by calling get_block_size().
 *		@li Feed input stream through the process() command in blocks of get_block_size() length.
 *		    The function will return the detected, filtered key press.
 *
 *	@see
 *		<A HREF="https://www.itu.int/rec/dologin_pub.asp?lang=e&id=T-REC-Q.24-198811-I!!PDF-E&type=items">ITU Q.24 Specification</A><br>
 *		http://www.ti.com/lit/an/spra066/spra066.pdf
 */

class dtmf_decoder
{
  public:
	// -------------------------------------------------------------------------
	/**
	 * 	Detects DTMF activity in the supplied audio data.
	 *
	 *  @tparam SampleT sample type, can be either 32 bit float or 16 bit signed integer.
	 * 	@param [in] sample_count	Number of audio samples in buffer @p samples_in.
	 *  @param [in] samples_in		Audio data.
	 *
	 * 	@returns ASCII code of the digit being pressed when the end of the 
	 * 			 buffer was reached, or zero if none.
	 */

	template <typename SampleT>
	char operator()(size_t sample_count, const SampleT* samples_in) noexcept
	{
		while (sample_count)
		{
			const auto N = std::min(DFT_SIZE - samples_in_buffer_, sample_count);

			if constexpr (std::is_same_v<SampleT, float>)
			{
				std::copy_n(samples_in, N, samples_buffer_ + samples_in_buffer_);
			}
			else if constexpr (std::is_same_v<SampleT, int16_t> || std::is_same_v<SampleT, short>)
			{
				details::pcm_to_float(samples_in, N, samples_buffer_ + samples_in_buffer_);
			}
			else
			{
				static_assert(false, "SampleT must be either float or signed short");
			}

			samples_in_buffer_ += N;
			samples_in += N;
			sample_count -= N;

			if (samples_in_buffer_ >= DFT_SIZE)
			{
				process_key(process_samples_buffer());
				std::copy_n(samples_buffer_ + BLOCK_SIZE, DFT_SIZE - BLOCK_SIZE, samples_buffer_);
				samples_in_buffer_ -= BLOCK_SIZE;
			}
		}
		return current_key_;
	}

	// -------------------------------------------------------------------------
	/**
	 * @brief 	Registers a callback function for key press events.
	 * 
	 * @param [in] proc User-defined callable.  It is the responsibility of the 
	 * 					caller to ensure this callable does not throw exceptions.
	 */

	void on_key_up(std::function<void(char)> proc)
	{
		on_key_up_ = std::move(proc);
	}

	// -------------------------------------------------------------------------
	/**
	 * @brief 	Registers a callback funct-ion for key up events.
	 * 
	 * @param [in] proc User-defined callable.  It is the responsibility of the 
	 * 					caller to ensure this callable does not throw exceptions.
	 */

	void on_key_down(std::function<void(char)> proc)
	{
		on_key_down_ = std::move(proc);
	}

	// -------------------------------------------------------------------------
	/**
	 * @brief Clears internal buffers and filters.
	 */

	void clear() noexcept
	{
		current_key_       = 0;
		samples_in_buffer_ = 0;
		last_key_          = 0;
		keydown_counter_   = 0;
		keyup_counter_     = 0;
	}

	// -------------------------------------------------------------------------
	/**
	 * @brief Get the proxcessing block size.
	 * 
	 *	The block size is 10 ms of audio.  It is not the same as the DTF size.
	 * 
	 * @return size_t The block size in samples.
	 */

	size_t get_block_size() const noexcept
	{
		return BLOCK_SIZE;
	}

  private:
	// -------------------------------------------------------------------------
	static constexpr char get_key(size_t row, size_t col) noexcept
	{
		if (row >= 4 || col >= 4)
			return 0;
		return KEYPAD[row][col];
	}

	// -------------------------------------------------------------------------
	/**
	 * @brief This function detects key presses using the Goertzel filter bank.
	 * 
	 * @return char The ASCII code for currently detected key, of 0, if none.
	 */

	char process_samples_buffer() noexcept
	{
		// ---------------------------------------------------------------------

		auto find_unique_above_threshold = [](const float* first, const float* last, float threshold) noexcept -> int {
			auto it = std::find_if(first, last, [threshold = threshold](float x) noexcept { return x > threshold; });
			if (it == last)
				return -1;

			if (std::find_if(it + 1, last, [threshold](float x) noexcept { return x > threshold; }) != last)
				return -1;

			return std::distance(first, it);
		};

		// ---------------------------------------------------------------------

		// the minimum threshold is defined in dBs of maximum magnitude
		// this multiplier expresses it as someting closer to mean power
		// as expressed by the DFT.
		constexpr auto power_multiplier = float(BLOCK_SIZE) / (float(DFT_SIZE) * 1.4142135623730950488016887f); // = std::sqrt(2.f);

		alignas(16) float power[4];

		// check for row signal above detection threshold.  This is the only absolute value threshold
		// all other thresholds will be compared relative to that.
		details::compute_goertzel_power(4, &goertzel_coefficients_[0], samples_in_buffer_, samples_buffer_, power);
		const auto max_power_row = *std::max_element(std::begin(power), std::end(power));
		if (max_power_row < MIN_SIGNAL_STRENGTH * power_multiplier)
			return 0;

		const auto row = find_unique_above_threshold(std::begin(power), std::end(power), MAX_SIGNAL_DIFF * max_power_row);
		if (row < 0)
			return 0;

		// checkl for sine signal in column ranges
		details::compute_goertzel_power(4, &goertzel_coefficients_[4], samples_in_buffer_, samples_buffer_, power);
		const auto max_power_col = *std::max_element(std::begin(power), std::end(power));
		if (max_power_row * ROW_COLUM_RATIO > max_power_col || max_power_col * ROW_COLUM_RATIO > max_power_row)
			return 0;

		const auto col = find_unique_above_threshold(std::begin(power), std::end(power), MAX_SIGNAL_DIFF * max_power_col);
		if (col < 0)
			return 0;

		// check for any harmonics above dynamic threshold
		details::compute_goertzel_power(4, &goertzel_coefficients_[8], samples_in_buffer_, samples_buffer_, power);
		if (std::find_if(std::begin(power),
		                 std::end(power),
		                 [t = max_power_row * SECOND_HARMONIC_THRESHOLD](float x) noexcept { return x >= t; })
		    != std::end(power))
		{
			return 0;
		}

		details::compute_goertzel_power(4, &goertzel_coefficients_[12], samples_in_buffer_, samples_buffer_, power);
		if (std::find_if(std::begin(power),
		                 std::end(power),
		                 [t = max_power_col * SECOND_HARMONIC_THRESHOLD](float x) noexcept { return x >= t; })
		    != std::end(power))
		{
			return 0;
		}

		return get_key(row, col);
	}

	// -------------------------------------------------------------------------
	/**
	 *	@brief Filters detected key presses to match Q.23/24 standard.
	 *
	 *	@param [in] key The last key detected by the filter bank 
	 *					(result of process_samples_buffer())
	 *
	 *	@see <A HREF="https://www.itu.int/rec/dologin_pub.asp?lang=e&id=T-REC-Q.24-198811-I!!PDF-E&type=items">
	 *			ITU specification T-REC-Q.24-198811</A>
	 *		 table A1/Q24, ANNEX A, page 3, "signal reception timing".
	 */
	char process_key(char key) noexcept
	{
		if (!key)
		{
			if (++keyup_counter_ >= MAX_INTERRUPTION_DURATION)
			{
				if (current_key_ && on_key_up_)
					on_key_up_(current_key_);

				last_key_        = 0;
				keydown_counter_ = 0;
				keyup_counter_   = MAX_INTERRUPTION_DURATION;
				current_key_     = 0;
			}
		}
		else
		{
			if (key == last_key_ && keyup_counter_ && keyup_counter_ < MAX_INTERRUPTION_DURATION)
				keydown_counter_ += keyup_counter_;

			keyup_counter_ = 0;

			if (key != last_key_)
				keydown_counter_ = 0;

			last_key_ = key;

			if (++keydown_counter_ >= MIN_SIGNAL_DURATION)
			{
				keydown_counter_ = MIN_SIGNAL_DURATION;
				if (key != current_key_)
				{
					if (current_key_ && on_key_up_)
						on_key_up_(current_key_);

					current_key_ = key;
					if (on_key_down_)
						on_key_down_(current_key_);
				}
			}
		}
		return current_key_;
	}

  private:
	// =========================================================================
	//  Keypad characteristics.

	static constexpr size_t NUMROWS          = 4;
	static constexpr size_t NUMCOLS          = 4;
	static constexpr size_t NUMTONES         = NUMROWS * NUMCOLS;        //!< total number of keys
	static constexpr float ROW_FREQ[NUMROWS] = {697, 770, 852, 941};     //!< keypad row frequencies
	static constexpr float COL_FREQ[NUMCOLS] = {1209, 1336, 1477, 1633}; //!< keypad column frequencies

	static constexpr char KEYPAD[NUMROWS][NUMCOLS] = {
	    {'1', '2', '3', 'A'},
	    {'4', '5', '6', 'B'},
	    {'7', '8', '9', 'C'},
	    {'*', '0', '#', 'D'},
	};

	// =========================================================================
	//  DFT filter bank

	static constexpr float SAMPLE_RATE = 8000.f;
	static constexpr size_t DFT_SIZE   = 168;
	static_assert(DFT_SIZE >= 105, "DTMF detection needs at least 105 samples to function at all, much less reliably");

	/**
	 * @brief  Minimal magnitude signal strength threshold.
	 */
	static constexpr float MIN_SIGNAL_STRENGTH = 0.0019952623149688796013524554f; // -27dB;

	/**
	 * @brief  Maximum power difference between frequencies.
	 */
	static constexpr float MAX_SIGNAL_DIFF = 0.15848931924611134852021014f; // -8dB;

	/**
	 * @brief Maximum allowed second harmonic threshold
	 */
	static constexpr float SECOND_HARMONIC_THRESHOLD = 0.15848931924611134852021014f; // -8dB;

	/**
	 * @brief  Maximum power difference between row and column frequencies.
	 */
	static constexpr float ROW_COLUM_RATIO = 0.15848931924611134852021014f; // +-8dB;

	/**
	 * @brief Number of filters in the bank, 4 columns x 4 rows + second harmonics
	 */
	static constexpr size_t NUMCOEFS = 2 * (NUMROWS + NUMCOLS);

	/**
	 * @brief DFT filter coefficients.
	 */
	static const float goertzel_coefficients_[NUMCOEFS];

	// =========================================================================
	//  Keypress filter characteristics

	/**
	 * @brief The block size sets the sample rate of the key detection stream.
	 */
	static constexpr size_t BLOCK_SIZE = size_t((.01 * SAMPLE_RATE) + .5);
	static_assert(BLOCK_SIZE <= DFT_SIZE, "BLOCK_SIZE must be smaller than DFT");

	/**
	 * Rate at which the algorithm produces keypress readings.
	 */
	static constexpr float KEY_SAMPLE_RATE = SAMPLE_RATE / BLOCK_SIZE;

	/**
	 *	Minimum amount of time a key has to be pressed for detection.
	 */
	static constexpr unsigned int MIN_SIGNAL_DURATION = static_cast<unsigned int>((.040f * KEY_SAMPLE_RATE) + .5f);

	/**
	 *	Minimum amount of time with no key pressed to generate a key_up event.
	 */
	static constexpr unsigned int MIN_PAUSE_DURATION = static_cast<unsigned int>((.030f * KEY_SAMPLE_RATE) + .5f);

	/**
	 *	Maximum signal interruption duration
	 */
	static constexpr unsigned int MAX_INTERRUPTION_DURATION = static_cast<unsigned int>((.020f * KEY_SAMPLE_RATE) + .5f);

	// =========================================================================
	//  data

	alignas(64) float samples_buffer_[DFT_SIZE] = {}; //!< DFT imput buffer.
	size_t samples_in_buffer_                   = 0;  //!< Number of samples currently in the buffer
	size_t keydown_counter_                     = 0;  //!< Timer for key down event.
	size_t keyup_counter_                       = 0;  //!< Timer for key up event.
	char current_key_                           = 0;  //!< Current key being detected, after filtering, or zero if there was none.
	char last_key_                              = 0;  //!< Current key being detected, unfiltered, or zero if there was none.
	std::function<void(char)> on_key_down_;           //!< User callback for key down events.
	std::function<void(char)> on_key_up_;             //!< User callback for key up events.
};

// *****************************************************************************

alignas(64) inline const float dtmf_decoder::goertzel_coefficients_[dtmf_decoder::NUMCOEFS] = {
    // base frequencies
    details::goertzel_coefficient(ROW_FREQ[0] / SAMPLE_RATE),
    details::goertzel_coefficient(ROW_FREQ[1] / SAMPLE_RATE),
    details::goertzel_coefficient(ROW_FREQ[2] / SAMPLE_RATE),
    details::goertzel_coefficient(ROW_FREQ[3] / SAMPLE_RATE),

    details::goertzel_coefficient(COL_FREQ[0] / SAMPLE_RATE),
    details::goertzel_coefficient(COL_FREQ[1] / SAMPLE_RATE),
    details::goertzel_coefficient(COL_FREQ[2] / SAMPLE_RATE),
    details::goertzel_coefficient(COL_FREQ[3] / SAMPLE_RATE),

    // harmonics
    details::goertzel_coefficient(2 * ROW_FREQ[0] / SAMPLE_RATE),
    details::goertzel_coefficient(2 * ROW_FREQ[1] / SAMPLE_RATE),
    details::goertzel_coefficient(2 * ROW_FREQ[2] / SAMPLE_RATE),
    details::goertzel_coefficient(2 * ROW_FREQ[3] / SAMPLE_RATE),

    details::goertzel_coefficient(2 * COL_FREQ[0] / SAMPLE_RATE),
    details::goertzel_coefficient(2 * COL_FREQ[1] / SAMPLE_RATE),
    details::goertzel_coefficient(2 * COL_FREQ[2] / SAMPLE_RATE),
    details::goertzel_coefficient(2 * COL_FREQ[3] / SAMPLE_RATE),
};

} // namespace tika

#endif // ****************************************************************** EOF
