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

set(DOXYGEN_COMMAND  "doxygen")

# ******************************************************************************

add_custom_target(Documentation
    COMMAND 
        ${DOXYGEN_COMMAND} -q Doxyfile
    DEPENDS 
        "${CMAKE_BINARY_DIR}/documentation/Doxyfile"
    WORKING_DIRECTORY 
        "${CMAKE_BINARY_DIR}/documentation"
    COMMENT 
        "Building documentation..."
    SOURCES 
        "${CMAKE_BINARY_DIR}/documentation/Doxyfile"
)

# ******************************************************************************

file(WRITE "${CMAKE_BINARY_DIR}/documentation/Doxyfile.cmake"
[=[
    file(READ "${SOURCE}" TEXT)
    string(REPLACE "{{OUTPUT_DIR}}" "${BINARY_DIR}" TEXT "${TEXT}")
    string(REPLACE "{{SOURCE_DIR}}" "${SOURCE_DIR}" TEXT "${TEXT}")
    string(REPLACE "{{STRIP_FROM_PATH}}" "${STRIP_FROM_PATH}" TEXT "${TEXT}")
    file(WRITE "${TARGET}" "${TEXT}")
]=])

add_custom_command(
    OUTPUT 
        "${CMAKE_BINARY_DIR}/documentation/Doxyfile"
    COMMAND 
        "${CMAKE_COMMAND}"
            "-DSOURCE=${CMAKE_SOURCE_DIR}/documentation/Doxyfile"
            "-DTARGET=${CMAKE_BINARY_DIR}/documentation/Doxyfile"
            "-DSOURCE_DIR=${CMAKE_SOURCE_DIR}/src"
            "-DBINARY_DIR=${CMAKE_BINARY_DIR}/documentation"
            "-DSTRIP_FROM_PATH=${CMAKE_SOURCE_DIR}/src"
            -P "${CMAKE_BINARY_DIR}/documentation/Doxyfile.cmake"
    DEPENDS 
        "${CMAKE_SOURCE_DIR}/documentation/Doxyfile" 
        "${CMAKE_BINARY_DIR}/documentation/Doxyfile.cmake"
)

# ************************************************************************** EOF
