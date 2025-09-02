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
 * @file 	dtmf.cpp
 * @brief 	Tests for class tika::dtmf_decoder
 * @internal
 */

#include <dtmf_decoder.h>
#include <gtest/gtest.h>

#include <dialer.h>
#include <dtmf_generator.h>
#include <white_noise_generator.h>

#include <cmath>
#include <fstream>
#include <random>

#ifdef min
#undef min
#endif

#ifdef max
#undef max
#endif

namespace {

constexpr int ITERATIONS = 100;
constexpr auto SMPRATE   = 8000;

// *****************************************************************************

struct ntt_limits
{
	const float pos_freq_tol            = .018f;
	const float neg_freq_tol            = .03f;
	const float min_power_level         = dBm(-24.f);
	const float max_power_level         = dBm(-3.f);
	const float neg_power_level         = dBm(-29.f);
	const float freq_diff_range_min     = dBm(-5.f);
	const float freq_diff_range_max     = dBm(5.f);
	const float signal_duration_min     = .040f;
	const float signal_duration_neg     = .024f;
	const float signal_interruption_max = .010f;
	const float pause_duration_min      = .030f;
	const float signal_velocity_max     = 0.120f;
};

// *****************************************************************************

struct att_limits
{
	const float pos_freq_tol            = .015f;
	const float neg_freq_tol            = .035f;
	const float min_power_level         = dBm(-25.f);
	const float max_power_level         = dBm(0.f);
	const float neg_power_level         = dBm(-55.f);
	const float freq_diff_range_min     = dBm(-8.f);
	const float freq_diff_range_max     = dBm(+4.f);
	const float signal_duration_min     = .040f;
	const float signal_duration_neg     = .023f;
	const float signal_interruption_max = .010f;
	const float pause_duration_min      = .040f;
	const float signal_velocity_max     = 0.093f;
};

// *****************************************************************************

struct danish_limits
{
	const float pos_freq_tol            = .0153f;
	const float neg_freq_tol            = .07f;
	const float min_power_level         = dBm(-30.f);
	const float max_power_level         = dBm(0.f);
	const float neg_power_level         = dBm(-39.f);
	const float freq_diff_range_min     = dBm(-6.f);
	const float freq_diff_range_max     = dBm(+6.f);
	const float signal_duration_min     = .040f;
	const float signal_duration_neg     = .020f;
	const float signal_interruption_max = .020f;
	const float pause_duration_min      = .040f;
	const float signal_velocity_max     = 0.100f;
};

// *****************************************************************************

struct australian_limits
{
	const float pos_freq_tol            = .0156f;
	const float neg_freq_tol            = .07f;
	const float min_power_level         = dBm(-27.f);
	const float max_power_level         = dBm(-5.f);
	const float neg_power_level         = dBm(-30.f);
	const float freq_diff_range_min     = dBm(-10.f);
	const float freq_diff_range_max     = dBm(+10.f);
	const float signal_duration_min     = .040f;
	const float signal_duration_neg     = .025f;
	const float signal_interruption_max = .012f;
	const float pause_duration_min      = .070f;
	const float signal_velocity_max     = 0.125f;
};

// *****************************************************************************

struct brazilian_limits
{
	const float pos_freq_tol            = .018f;
	const float neg_freq_tol            = .03f;
	const float min_power_level         = dBm(-25.f);
	const float max_power_level         = dBm(-3.f);
	const float neg_power_level         = dBm(-50.f);
	const float freq_diff_range_min     = dBm(-9.f);
	const float freq_diff_range_max     = dBm(+9.f);
	const float signal_duration_min     = .040f;
	const float signal_duration_neg     = .020f;
	const float signal_interruption_max = .012f;
	const float pause_duration_min      = .030f;
	const float signal_velocity_max     = 0.120f;
};

// *****************************************************************************

struct test_dialer
{
	// -------------------------------------------------------------------------

	std::vector<float> dial(std::string_view s)
	{
		constexpr auto SMP_RATE = 8000.f;

		dialer_.set_keypress_duration(signal_duration_);
		dialer_.set_spacing_duration(pause_duration_);
		dialer_.generator().set_gain(signal_level_);
		noise_gen_.set_gain(signal_level_ / signal_noise_ratio_);

		const auto prolog_len = unsigned(std::floor((prolog_ * SMP_RATE) + .5f));
		const auto epilog_len = unsigned(std::floor((epilog_ * SMP_RATE) + .5f));

		std::vector<float> result;
		result.reserve(dialer_.estimate_dial_duration(s) + prolog_len + epilog_len);

		for (unsigned i = 0; i < prolog_len; ++i)
			result.push_back(noise_gen_());

		dialer_.dial(s);
		while (dialer_.dialing())
			result.push_back(dialer_() + noise_gen_());

		for (unsigned i = 0; i < epilog_len; ++i)
			result.push_back(noise_gen_());

		return result;
	}

	// -------------------------------------------------------------------------

	void detune(float factor)
	{
		dialer_.set_sample_rate(8000.f * (1.f - factor));
		signal_duration_ *= (1.f + factor);
		pause_duration_ *= (1.f + factor);
	}

	// -------------------------------------------------------------------------

	dialer dialer_;
	white_noise_generator<float> noise_gen_;
	float prolog_             = .1f;
	float epilog_             = .2f;
	float signal_duration_    = .05f; // .4
	float pause_duration_     = .1f;  // .93
	float signal_noise_ratio_ = dBm(12.f);
	float signal_level_       = dBm(-6.f);
};

// *****************************************************************************

std::vector<int16_t> to_pcm(const std::vector<float>& v)
{
	constexpr auto scale = float(0x7FFF);
	auto r               = std::vector<int16_t>{};
	r.reserve(v.size());

	for (const auto& f : v)
		r.push_back(int16_t(floor((f * scale) + .5)));

	return r;
}

// *****************************************************************************

TEST(dtmf_decoder, test_dialer_test)
{
	test_dialer keypad;
	keypad.signal_level_       = dBm(-27.f);
	keypad.prolog_             = .0f;
	keypad.signal_duration_    = .1f;
	keypad.signal_noise_ratio_ = dBm(6.f);

	auto audio = keypad.dial("0123456789*#ABCDEF");

	std::ofstream ofs("test_dialer_output.raw", std::ofstream::binary);
	ofs.write(reinterpret_cast<char*>(audio.data()), audio.size() * sizeof(audio[0]));
	ofs.close();
}

// *****************************************************************************

TEST(dtmf_decoder, small_signal)
{
	test_dialer keypad;
	tika::dtmf_decoder decoder;
	std::string read;
	decoder.on_key_down([&read](char key) noexcept { read += key; });

	keypad.signal_level_       = dBm(-27.f);
	keypad.prolog_             = .1f;
	keypad.signal_duration_    = .05f;
	keypad.signal_noise_ratio_ = dBm(6.f);

	for (int i = 0; i < ITERATIONS; ++i)
	{
		read.clear();
		auto audio = keypad.dial("1");
		decoder(audio.size(), audio.data());
		EXPECT_EQ(read, "1");

		read.clear();
		audio = keypad.dial("2");
		decoder(audio.size(), audio.data());
		EXPECT_EQ(read, "2");

		read.clear();
		audio = keypad.dial("3");
		decoder(audio.size(), audio.data());
		EXPECT_EQ(read, "3");

		read.clear();
		audio = keypad.dial("4");
		decoder(audio.size(), audio.data());
		EXPECT_EQ(read, "4");

		read.clear();
		audio = keypad.dial("5");
		decoder(audio.size(), audio.data());
		EXPECT_EQ(read, "5");

		read.clear();
		audio = keypad.dial("6");
		decoder(audio.size(), audio.data());
		EXPECT_EQ(read, "6");

		read.clear();
		audio = keypad.dial("7");
		decoder(audio.size(), audio.data());
		EXPECT_EQ(read, "7");

		read.clear();
		audio = keypad.dial("8");
		decoder(audio.size(), audio.data());
		EXPECT_EQ(read, "8");

		read.clear();
		audio = keypad.dial("9");
		decoder(audio.size(), audio.data());
		EXPECT_EQ(read, "9");

		read.clear();
		audio = keypad.dial("0");
		decoder(audio.size(), audio.data());
		EXPECT_EQ(read, "0");

		// below are tests using PCM input instead
		read.clear();
		auto audio_pcm = to_pcm(keypad.dial("*"));
		decoder(audio_pcm.size(), audio_pcm.data());
		EXPECT_EQ(read, "*");

		read.clear();
		audio_pcm = to_pcm(keypad.dial("#"));
		decoder(audio_pcm.size(), audio_pcm.data());
		EXPECT_EQ(read, "#");

		read.clear();
		audio_pcm = to_pcm(keypad.dial("A"));
		decoder(audio_pcm.size(), audio_pcm.data());
		EXPECT_EQ(read, "A");

		read.clear();
		audio_pcm = to_pcm(keypad.dial("B"));
		decoder(audio_pcm.size(), audio_pcm.data());
		EXPECT_EQ(read, "B");

		read.clear();
		audio_pcm = to_pcm(keypad.dial("C"));
		decoder(audio_pcm.size(), audio_pcm.data());
		EXPECT_EQ(read, "C");

		read.clear();
		audio_pcm = to_pcm(keypad.dial("D"));
		decoder(audio_pcm.size(), audio_pcm.data());
		EXPECT_EQ(read, "D");
	}
}

// *****************************************************************************

TEST(dtmf_decoder, speed_dial)
{
	test_dialer keypad;
	tika::dtmf_decoder decoder;
	std::string read;
	decoder.on_key_down([&read](char key) noexcept { read += key; });
	keypad.prolog_          = 0.5f;
	keypad.signal_duration_ = .40f;
	keypad.pause_duration_  = .93f;

	std::string dial_string = "123A456B789C*0#D";

	for (int i = 0; i < ITERATIONS; ++i)
	{
		read.clear();
		std::next_permutation(dial_string.begin(), dial_string.end());
		const auto audio = keypad.dial(dial_string);
		decoder(audio.size(), audio.data());
		EXPECT_EQ(read, dial_string);
	}
}

// *****************************************************************************

TEST(dtmf_decoder, speed_dial_synchronous)
{
	test_dialer keypad;

	tika::dtmf_decoder decoder;

	std::string read;

	keypad.prolog_          = 0.5f;
	keypad.signal_duration_ = .40f;
	keypad.pause_duration_  = .93f;

	std::string dial_string = "123A456B789C*0#D";

	for (int i = 0; i < ITERATIONS; ++i)
	{
		read.clear();
		std::next_permutation(dial_string.begin(), dial_string.end());
		auto audio = keypad.dial(dial_string);

		char prev = 0;
		for (size_t j = 0; j < audio.size(); j += decoder.get_block_size())
		{
			auto key = decoder(std::min(decoder.get_block_size(), audio.size() - j), &audio[j]);
			if (key && key != prev)
				read += key;
			prev = key;
		}
		// std::cout << std::endl;
		EXPECT_EQ(read, dial_string);
	}
}

// *****************************************************************************

TEST(dtmf_decoder, detune_plus_1_5_pct)
{
	test_dialer keypad;
	keypad.signal_duration_ = .1f;
	keypad.pause_duration_  = .1f;
	keypad.detune(.0015f);
	tika::dtmf_decoder decoder;
	std::string read;
	decoder.on_key_down([&read](char key) noexcept { read += key; });

	std::string dial_string = "1234567890*#ABCD";

	for (int i = 0; i < ITERATIONS; ++i)
	{
		read.clear();
		std::next_permutation(dial_string.begin(), dial_string.end());
		auto audio = keypad.dial(dial_string);
		decoder(audio.size(), audio.data());
		EXPECT_EQ(read, dial_string);
	}
}

// *****************************************************************************

TEST(dtmf_decoder, detune_minus_1_5_pct)
{
	test_dialer keypad;
	keypad.signal_duration_ = .1f;
	keypad.pause_duration_  = .1f;
	keypad.detune(-.015f);
	tika::dtmf_decoder decoder;
	std::string read;
	decoder.on_key_down([&read](char key) noexcept { read += key; });

	std::string dial_string = "1234567890*#ABCD";

	for (int i = 0; i < ITERATIONS; ++i)
	{
		std::next_permutation(dial_string.begin(), dial_string.end());
		read.clear();
		auto audio = keypad.dial(dial_string);
		decoder(audio.size(), audio.data());
		EXPECT_EQ(read, dial_string);
	}
}

// *****************************************************************************

TEST(dtmf_decoder, detune_plus_3_5_pct)
{
	test_dialer keypad;
	tika::dtmf_decoder decoder;
	std::string read;
	decoder.on_key_down([&read](char key) noexcept { read += key; });

	std::string dial_string = "0123456789*#ABCD";

	keypad.dialer_.generator().set_gain_dB(-3);
	keypad.signal_duration_ = .250f;
	keypad.detune(0.035f);

	for (int i = 0; i < ITERATIONS; ++i)
	{
		read.clear();
		std::next_permutation(dial_string.begin(), dial_string.end());
		auto audio = keypad.dial(dial_string);
		decoder(audio.size(), audio.data());
		EXPECT_TRUE(read.empty());
	}
}

// *****************************************************************************

// TEST(dtmf_decoder, detune_minus_3_5_pct)
// {
// 	test_dialer keypad;
// 	tika::dtmf_decoder decoder;
// 	std::string read;
// 	decoder.on_key_down([&read](char key) noexcept { read += key; });

// 	std::string dial_string = "0123456789*#ABCD";

// 	keypad.dialer_.generator().set_gain_dB(-3);
// 	keypad.detune(-0.035f);
// 	keypad.signal_duration_ = .250f;

// 	for (int i = 0; i < ITERATIONS; ++i)
// 	{
// 		read.clear();
//		std::next_permutation(dial_string.begin(), dial_string.end());
// 		auto audio = keypad.dial(dial_string);
// 		decoder(audio.size(), audio.data());
// 		EXPECT_TRUE(read.empty());
// 	}
// }

} // namespace
// ************************************************************************* EOF
