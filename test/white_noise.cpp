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
 * @file    white_noise.cpp
 * @brief   Unit testing for class white_noise_generator
 * @internal
 */

#include <gtest/gtest.h>

#include <decibel.h>
#include <white_noise_generator.h>

#include <cassert>
#include <cmath>
#include <vector>

namespace {

constexpr size_t TEST_SIZE = 1 << 20;
constexpr float PRECIS     = 1.e-2f; // expected statistical precision

// *****************************************************************************

float get_sum(size_t count, const float* samples) noexcept
{
	auto sum = float{};
	for (size_t i = 0; i < count; ++i)
		sum += samples[i];
	return sum;
}

// *****************************************************************************

float get_sum_of_squares(size_t count, const float* samples) noexcept
{
	auto sum = float{};
	for (size_t i = 0; i < count; ++i)
		sum += (samples[i] * samples[i]);
	return sum;
}

// *****************************************************************************

template <typename Float>
std::vector<Float> make_noise(white_noise_generator<Float>& gen, size_t count)
{
	auto result = std::vector<Float>{};
	result.reserve(count);
	for (size_t i = 0; i < count; ++i)
		result.push_back(gen());
	return result;
}

// *****************************************************************************

TEST(white_noise, statistics)
{
	auto gen = white_noise_generator<float>{};
	gen.set_power_dB(0);

	auto noise = make_noise(gen, TEST_SIZE);

	auto sum            = get_sum(noise.size(), noise.data());
	auto sum_of_squares = get_sum_of_squares(noise.size(), noise.data());

	auto mean               = sum / TEST_SIZE;
	auto mean_power         = sum_of_squares / TEST_SIZE;
	auto standard_deviation = std::sqrt((sum_of_squares - ((sum * sum) / TEST_SIZE)) / (TEST_SIZE - 1));

	EXPECT_NEAR(mean, 0, PRECIS);
	EXPECT_NEAR(mean_power, 1.f, PRECIS);
	EXPECT_NEAR(standard_deviation, 1.f, PRECIS);
}

// *****************************************************************************

TEST(white_noise, gain)
{
	const auto GAIN = dBm(-12.0f);
	const auto PWR  = GAIN * GAIN;

	auto gen = white_noise_generator<float>{};

	gen.set_gain(GAIN);

	EXPECT_EQ(gen.get_gain(), GAIN);
	EXPECT_NEAR(gen.get_power_dB(), to_dBm(PWR), 1.e-6f);

	auto noise = make_noise(gen, TEST_SIZE);

	auto sum            = get_sum(noise.size(), noise.data());
	auto sum_of_squares = get_sum_of_squares(noise.size(), noise.data());

	auto mean       = sum / TEST_SIZE;
	auto mean_power = sum_of_squares / TEST_SIZE;
	auto rms_power  = (sum_of_squares - ((sum * sum) / TEST_SIZE)) / (TEST_SIZE - 1);

	EXPECT_NEAR(mean, 0, PRECIS);
	EXPECT_NEAR(mean_power, PWR, PRECIS);
	EXPECT_NEAR(rms_power, PWR, PRECIS);
}

// *****************************************************************************

TEST(white_noise, power_dB)
{
	constexpr auto PWR_dB = -12.0f;
	const auto PWR        = dBm(PWR_dB);

	auto gen = white_noise_generator<float>{};

	gen.set_power_dB(PWR_dB);

	EXPECT_NEAR(gen.get_power_dB(), PWR_dB, 1.e-6f);
	EXPECT_NEAR(gen.get_gain(), std::sqrt(PWR), 1.e-6f);

	auto noise = make_noise(gen, TEST_SIZE);

	auto sum            = get_sum(noise.size(), noise.data());
	auto sum_of_squares = get_sum_of_squares(noise.size(), noise.data());

	auto mean       = sum / TEST_SIZE;
	auto mean_power = sum_of_squares / TEST_SIZE;
	auto rms_power  = (sum_of_squares - ((sum * sum) / TEST_SIZE)) / (TEST_SIZE - 1);

	EXPECT_NEAR(mean, 0, PRECIS);
	EXPECT_NEAR(mean_power, PWR, PRECIS);
	EXPECT_NEAR(rms_power, PWR, PRECIS);
}

} // namespace
// ************************************************************************* EOF
