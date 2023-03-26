cmake_minimum_required(VERSION 3.1)
project (localcondcalc)
set(CMAKE_CXX_STANDARD 11)

#include_directories("${PROJECT_SOURCE_DIR}")
add_executable(testLocalCondCalc testLocalCondCalc.cpp localCondCalc.h)
#add_library(localcondcalc SHARED localCondCalc-impl.h)
#target_compile_features(localcondcalc PRIVATE cxx_variadic_templates)
