# mrock_iEoM_extra_options.cmake
# If the Intel MKL is installed on the system, the iEoM
# library can make use of it to parallelize certain matrix operations
# which tends to make the computations much faster
include_guard()

include(${CMAKE_CURRENT_LIST_DIR}/../../cmake/mrock_message.cmake)

add_library(mrock_iEoM_extra_options INTERFACE)

#
# OpenMP
#

find_package(OpenMP QUIET)
if (OpenMP_FOUND)
    mrock_message("Configuring target mrock_iEoM to use OpenMP")
    target_link_libraries(mrock_iEoM_extra_options INTERFACE OpenMP::OpenMP_CXX)
else()
    mrock_message("OpenMP not found or not enabled; mrock_iEoM will not use OpenMP")
endif()

#
# MKL
#

set(MKL_LINK static)
set(MKL_THREADING gnu_thread)
set(MKL_INTERFACE lp64)

find_package(MKL CONFIG QUIET HINTS $ENV{MKLROOT})
if (MKL_FOUND)
    set(USE_MKL ON CACHE BOOL "Use Intel MKL" FORCE)
else()
    set(USE_MKL OFF CACHE BOOL "Use Intel MKL" FORCE)
endif()

if (USE_MKL)
    mrock_message("Configuring target mrock_iEoM to use Intel MKL")

    target_compile_options(mrock_iEoM_extra_options INTERFACE $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
    target_include_directories(mrock_iEoM_extra_options INTERFACE $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
    target_link_libraries(mrock_iEoM_extra_options INTERFACE $<LINK_ONLY:MKL::MKL>)

	# Let Eigen use MKL
	target_compile_definitions(mrock_iEoM_extra_options INTERFACE EIGEN_USE_MKL_ALL MROCK_IEOM_DO_NOT_PARALLELIZE)
else()
    mrock_message("MKL not found or not enabled; mrock_iEoM will not use Intel MKL")
endif()