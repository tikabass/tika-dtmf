# ******************************************************************************
# MIT License
#
# Copyright (c) 2025 tikabass <tika.devel@gmail.com>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
# ******************************************************************************

include(GoogleBenchmark.cmake)

# ******************************************************************************

set(CMAKE_CXX_STANDARD          17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS        OFF)
set(CMAKE_DEBUG_POSTFIX         "d")

# ******************************************************************************

if(${CMAKE_HOST_LINUX})
    set (DTMF_CXX_SSE2_FLAG -msse2)
    set (DTMF_CXX_SSE41_FLAG -msse4.1)
    set (DTMF_CXX_AVX_FLAG -mavx)
    set (DTMF_CXX_AVX2_FLAG -mavx2)
    set (DTMF_CXX_FMA3_FLAG -mfma)
else()
    set (DTMF_CXX_SSE2_FLAG /arch:SSE2)
    set (DTMF_CXX_SSE41_FLAG /arch:SSE4.1)
    set (DTMF_CXX_AVX2_FLAG /arch:AVX2)
    set (DTMF_CXX_FMA3_FLAG /arch:FMA)
endif()

# ******************************************************************************

set(BENCHMARK_INCLUDE_DIRS
    "${googlebenchmark_SOURCE_DIR}"
    "${CMAKE_SOURCE_DIR}/src"
    "${CMAKE_SOURCE_DIR}/testlib"
)

set(BENCHMARK_LINK_DIRS
    "${googlebenchmark_BINARY_DIR}"
)

set(BENCHMARK_LINK_LIBRARIES
    benchmark
    benchmark_main
)

# ******************************************************************************

add_executable(throughput_avx2 
    throughput.cpp
)

target_compile_options(throughput_avx2       
    PRIVATE 
        ${DTMF_CXX_AVX_FLAG} ${DTMF_CXX_AVX2_FLAG} ${DTMF_CXX_FMA3_FLAG}
)

target_include_directories(throughput_avx2   PRIVATE ${BENCHMARK_INCLUDE_DIRS})
target_link_directories(throughput_avx2      PRIVATE ${BENCHMARK_LINK_DIRS})
target_link_libraries(throughput_avx2        PRIVATE ${BENCHMARK_LINK_LIBRARIES})

# ******************************************************************************

add_executable(throughput_fma3 
    throughput.cpp
)

target_compile_options(throughput_fma3       
    PRIVATE 
        ${DTMF_CXX_FMA3_FLAG} ${DTMF_CXX_SSE41_FLAG}
        -DTIKA_SIMD_LEVEL=TIKA_SIMD_LEVEL_FMA3
)

target_include_directories(throughput_fma3   PRIVATE ${BENCHMARK_INCLUDE_DIRS})
target_link_directories(throughput_fma3      PRIVATE ${BENCHMARK_LINK_DIRS})
target_link_libraries(throughput_fma3        PRIVATE ${BENCHMARK_LINK_LIBRARIES})

# ******************************************************************************

add_executable(throughput_sse41
    throughput.cpp
)

target_compile_options(throughput_sse41       
    PRIVATE 
        ${DTMF_CXX_SSE41_FLAG}
        -DTIKA_SIMD_LEVEL=TIKA_SIMD_LEVEL_SSE41
)

target_include_directories(throughput_sse41   PRIVATE ${BENCHMARK_INCLUDE_DIRS})
target_link_directories(throughput_sse41      PRIVATE ${BENCHMARK_LINK_DIRS})
target_link_libraries(throughput_sse41        PRIVATE ${BENCHMARK_LINK_LIBRARIES})

# ******************************************************************************

add_executable(throughput_sse2
    throughput.cpp
)

target_compile_options(throughput_sse2       
    PRIVATE 
        ${DTMF_CXX_SSE2_FLAG}
        -DTIKA_SIMD_LEVEL=TIKA_SIMD_LEVEL_SSE2
)

target_include_directories(throughput_sse2   PRIVATE ${BENCHMARK_INCLUDE_DIRS})
target_link_directories(throughput_sse2      PRIVATE ${BENCHMARK_LINK_DIRS})
target_link_libraries(throughput_sse2        PRIVATE ${BENCHMARK_LINK_LIBRARIES})

# ******************************************************************************

add_executable(throughput_no_simd
    throughput.cpp
)

target_compile_options(throughput_no_simd       
    PRIVATE 
        ${DTMF_CXX_SSE2_FLAG}
        -DTIKA_SIMD_LEVEL=TIKA_SIMD_LEVEL_NO_SIMD
)

target_include_directories(throughput_no_simd   PRIVATE ${BENCHMARK_INCLUDE_DIRS})
target_link_directories(throughput_no_simd      PRIVATE ${BENCHMARK_LINK_DIRS})
target_link_libraries(throughput_no_simd        PRIVATE ${BENCHMARK_LINK_LIBRARIES})

# ************************************************************************** EOF
