################################################################################
# This file is part of pasta::block_tree_experiments
#
# Copyright (C) 2023 Florian Kurpicz <florian@kurpicz.org>
#
# pasta::bit_vector_experiments is free software: you can redistribute
# it and/or modify it under the terms of the GNU General Public
# License as published by the Free Software Foundation, either version
# 3 of the License, or (at your option) any later version.
#
# pasta::bit_vector_experiments is distributed in the hope that it
# will be useful, but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See
# the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with pasta::bit_vector_experiments. If not, see
# <http://www.gnu.org/licenses/>.
#
################################################################################

cmake_minimum_required(VERSION 3.4.1 FATAL_ERROR)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(pasta_bit_vector_experiments)

# Compiler flags for used in experiments, i.e., -03, -march=native, and -DNDEBUG
set(CMAKE_CXX_FLAGS
  "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic -fdiagnostics-color=auto -O3 -march=native -DNDEBUG -fopenmp")


add_subdirectory(extlib/block_tree)
add_subdirectory(extlib/original_block_tree)

add_executable(benchmark benchmarks/block_tree_benchmark.cpp)
target_link_libraries(benchmark
  pasta_block_tree
  blocktree
)

target_include_directories(benchmark
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/extlib/original_block_tree/includes
)

################################################################################
