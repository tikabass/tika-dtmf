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

#include <benchmark/benchmark.h>

#include <dialer.h>
#include <dtmf_decoder.h>
#include <white_noise_generator.h>

#include <cmath>
#include <vector>

// *****************************************************************************

std::vector<float> make_noise(white_noise_generator<float>& gen, size_t count)
{
	auto result = std::vector<float>{};
	result.reserve(count);
	for (size_t i = 0; i < count; ++i)
		result.push_back(gen());
	return result;
}

// *****************************************************************************

std::vector<int16_t> to_pcm(const std::vector<float>& v)
{
	constexpr auto scale = float(0x7FFF);
	auto r               = std::vector<int16_t>{};
	r.reserve(v.size());

	for (const auto& f : v)
		r.push_back(int16_t(std::floor((f * scale) + .5f)));

	return r;
}

// *****************************************************************************

void throughput_float_no_key(benchmark::State& state)
{
	white_noise_generator<float> noise_gen;
	tika::dtmf_decoder decoder;
	auto data = make_noise(noise_gen, 1 << 20);

	for (auto _ : state)
	{
		benchmark::DoNotOptimize(decoder(data.size(), data.data()));
	}
}

// *****************************************************************************

void throughput_pcm_no_key(benchmark::State& state)
{
	white_noise_generator<float> noise_gen;
	tika::dtmf_decoder decoder;
	auto data = to_pcm(make_noise(noise_gen, 1 << 20));

	for (auto _ : state)
	{
		benchmark::DoNotOptimize(decoder(data.size(), data.data()));
	}
}

// *****************************************************************************

void throughput_float_key_down(benchmark::State& state)
{
	dtmf_generator keypad;
	tika::dtmf_decoder decoder;
	auto data = std::vector<float>{};
	data.resize(1 << 20);

	keypad.press_key('1');
	for (auto& s : data)
		s = keypad();

	for (auto _ : state)
	{
		benchmark::DoNotOptimize(decoder(data.size(), data.data()));
	}
}

// *****************************************************************************

void throughput_pcm_key_down(benchmark::State& state)
{
	dtmf_generator keypad;
	tika::dtmf_decoder decoder;
	auto data = std::vector<float>{};
	data.resize(1 << 20);

	keypad.press_key('1');
	for (auto& s : data)
		s = keypad();

	auto pcm = to_pcm(data);

	for (auto _ : state)
	{
		benchmark::DoNotOptimize(decoder(pcm.size(), pcm.data()));
	}
}

// *****************************************************************************

BENCHMARK(throughput_float_no_key)->Arg(1)->Arg(1)->Arg(1);
BENCHMARK(throughput_float_key_down)->Arg(1)->Arg(1)->Arg(1);
BENCHMARK(throughput_pcm_no_key)->Arg(1)->Arg(1)->Arg(1);
BENCHMARK(throughput_pcm_key_down)->Arg(1)->Arg(1)->Arg(1);

BENCHMARK_MAIN();

// ************************************************************************* EOF
