# @copyright (c) 2020-2020 RWTH Aachen. All rights reserved.
#
# ddX software
#
# @file docs/CMakeLists.txt
# @version 1.0.0
# @author Aleksandr Mikhalev
# @date 2020-12-17

# Make lists of sources and headers appropriate for config file
string(REPLACE ";" "\"\\\n                         \"" DOXYGEN_INPUT
    "${SRC}")

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/config.in"
    "${CMAKE_CURRENT_BINARY_DIR}/config")
add_custom_target(docs COMMAND
    "${DOXYGEN_EXECUTABLE}" "${CMAKE_CURRENT_BINARY_DIR}/config")

