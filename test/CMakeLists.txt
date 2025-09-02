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

include(GoogleTest)

# ******************************************************************************

set(CMAKE_CXX_STANDARD          17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS        OFF)
set(CMAKE_DEBUG_POSTFIX         "d")

# ******************************************************************************

if(${CMAKE_HOST_LINUX})
    set(DTMF_CXX_SSE2_FLAG -msse2)
    set(DTMF_CXX_SSE41_FLAG -msse4.1)
    set(DTMF_CXX_AVX_FLAG -mavx)
    set(DTMF_CXX_AVX2_FLAG -mavx2)
    set(DTMF_CXX_FMA3_FLAG -mfma)
else()    
    set(DTMF_CXX_SSE2_FLAG /arch:SSE2)
    set(DTMF_CXX_SSE41_FLAG /arch:SSE4.1)
    set(DTMF_CXX_AVX2_FLAG /arch:AVX2)
    set(DTMF_CXX_FMA3_FLAG /arch:FMA)
endif()

# ******************************************************************************

set(DTMF_TEST_SOURCE_FILES
    dtmf.cpp
)

set(DTMF_TEST_INCLUDE_DIRS
    "${gtest_SOURCE_DIR}" 
    "${CMAKE_SOURCE_DIR}/src"
    "${CMAKE_SOURCE_DIR}/testlib"
)

set(DTMF_TEST_LINK_LIBRARIES
    gtest_main 
    gtest 
)

# ******************************************************************************

add_executable(dtmf_test_avx2 ${DTMF_TEST_SOURCE_FILES})

target_compile_options(dtmf_test_avx2       
    PRIVATE 
        ${DTMF_CXX_AVX_FLAG} ${DTMF_CXX_AVX2_FLAG} ${DTMF_CXX_FMA3_FLAG}
)

target_include_directories(dtmf_test_avx2   PRIVATE ${DTMF_TEST_INCLUDE_DIRS})
target_link_directories(dtmf_test_avx2      PRIVATE ${gtest_BINARY_DIR})
target_link_libraries(dtmf_test_avx2        PRIVATE ${DTMF_TEST_LINK_LIBRARIES})

# ******************************************************************************

add_executable(dtmf_test_fma3 ${DTMF_TEST_SOURCE_FILES})

target_compile_options(dtmf_test_fma3       
    PRIVATE 
        ${DTMF_CXX_FMA3_FLAG}
        ${DTMF_CXX_SSE2_FLAG}
        ${DTMF_CXX_SSE41_FLAG}
        -DTIKA_SIMD_LEVEL=TIKA_SIMD_LEVEL_FMA3
)

target_include_directories(dtmf_test_fma3   PRIVATE ${DTMF_TEST_INCLUDE_DIRS})
target_link_directories(dtmf_test_fma3      PRIVATE ${gtest_BINARY_DIR})
target_link_libraries(dtmf_test_fma3        PRIVATE ${DTMF_TEST_LINK_LIBRARIES})

# ******************************************************************************

add_executable(dtmf_test_sse41 ${DTMF_TEST_SOURCE_FILES})

target_compile_options(dtmf_test_sse41       
    PRIVATE 
        ${DTMF_CXX_SSE2_FLAG}
        ${DTMF_CXX_SSE41_FLAG}
        -DTIKA_SIMD_LEVEL=TIKA_SIMD_LEVEL_SSE41
)

target_include_directories(dtmf_test_sse41   PRIVATE ${DTMF_TEST_INCLUDE_DIRS})
target_link_directories(dtmf_test_sse41      PRIVATE ${gtest_BINARY_DIR})
target_link_libraries(dtmf_test_sse41        PRIVATE ${DTMF_TEST_LINK_LIBRARIES})

# ******************************************************************************

add_executable(dtmf_test_sse2 ${DTMF_TEST_SOURCE_FILES})

target_compile_options(dtmf_test_sse2
    PRIVATE 
        ${DTMF_CXX_SSE2_FLAG}
        -DTIKA_SIMD_LEVEL=TIKA_SIMD_LEVEL_SSE2
)

target_include_directories(dtmf_test_sse2   PRIVATE ${DTMF_TEST_INCLUDE_DIRS})
target_link_directories(dtmf_test_sse2      PRIVATE ${gtest_BINARY_DIR})
target_link_libraries(dtmf_test_sse2        PRIVATE ${DTMF_TEST_LINK_LIBRARIES})

# ******************************************************************************

add_executable(dtmf_test_no_simd ${DTMF_TEST_SOURCE_FILES})

target_compile_options(dtmf_test_no_simd
    PRIVATE 
        -DTIKA_SIMD_LEVEL=TIKA_SIMD_LEVEL_NO_SIMD
)

target_include_directories(dtmf_test_no_simd   PRIVATE ${DTMF_TEST_INCLUDE_DIRS})
target_link_directories(dtmf_test_no_simd      PRIVATE ${gtest_BINARY_DIR})
target_link_libraries(dtmf_test_no_simd        PRIVATE ${DTMF_TEST_LINK_LIBRARIES})

# ******************************************************************************

add_executable(dtmf_test_whitenoise 
    white_noise.cpp
)

target_include_directories(dtmf_test_whitenoise    PRIVATE ${DTMF_TEST_INCLUDE_DIRS})
target_link_directories(dtmf_test_whitenoise       PRIVATE ${gtest_BINARY_DIR})
target_link_libraries(dtmf_test_whitenoise         PRIVATE ${DTMF_TEST_LINK_LIBRARIES})

# ******************************************************************************

enable_testing()

add_test(NAME dtmf_test_avx2        COMMAND dtmf_test_avx2)
add_test(NAME dtmf_test_fma3        COMMAND dtmf_test_fma3)
add_test(NAME dtmf_test_sse41       COMMAND dtmf_test_sse41)
add_test(NAME dtmf_test_sse2        COMMAND dtmf_test_sse2)
add_test(NAME dtmf_test_no_simd     COMMAND dtmf_test_no_simd)  
add_test(NAME dtmf_test_whitenoise  COMMAND dtmf_test_whitenoise)

# ************************************************************************** EOF
