cmake_minimum_required(VERSION 3.18)
project(gpu_module LANGUAGES CUDA CXX)
# HELP usage is as follows
# 1) go to root directory where this e.g. root CMakeList.txt is located
# 2) run cmake -B build/ -S. [-D<OPTION_NAME>=<OPTION_VALUE>]
#   * for float enabled compilation and default SIM_TYPE
#       > cmake -B build/ -S. -DUSEFLOAT=ON
#   * for last default precision compilation and EXPANSION SIM_TYPE OPTION
#       > cmake -B build/ -S. -DSIM_TYPE_NAME=CONSTANT
#		DSIM_TYPE_NAME possible values: CONSTANT, EXPANSION, BERENDSEN, MTK, GEAR
# ?) Dont worry about misspelling the SIM_TYPE_NAME you will get error with correct options if you make a mistake
#    alternatively you can use the cmake-gui to configure - Be sure to delete from past changes when it does not work
# to make in paralell, -jN_cores parameter can be used
# 	example for 5 cores: make -j5 argon
# 	
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CUDA_STANDARD 17)

# setting of the constant values used in the mac module for fast and easy check/change
# PREC VALUES
set(FLOAT_VAL 1)
set(DOUBLE_VAL 2)
# SUBSTANCE TYPE VALUES
set(ARGON_VAL 1)
set(NITROGEN_VAL 2)
set(SPCE_VAL 3)
# SIM TYPE VALUES
set(CONSTANT_BOX_VAL 0)
set(CONSTANT_BOX_SEQ_VAL 1) # unused 
set(EXPANSION_BOX_VAL 2)
set(CONSTANT_BOX_BERENDSEN_VAL 3)
set(MTK_VAL 4)
set(GEAR_VAL 5)
# SIM TYPE NAMES
set(CONSTANT_BOX_NAME CONSTANT)
set(CONSTANT_BOX_SEQ_NAME CONSTANT_SEQ)
set(EXPANSION_BOX_NAME EXPANSION)
set(CONSTANT_BOX_BERENDSEN_NAME BERENDSEN)
set(MTK_NAME MTK)
set(GEAR_NAME GEAR)
# LIST OF SIM TYPE NAME
set(SIM_TYPE_OPTIONS "${CONSTANT_BOX_NAME};${EXPANSION_BOX_NAME};${CONSTANT_BOX_BERENDSEN_NAME};${MTK_NAME};${GEAR_NAME}")

macro(set_cuda_properties target substance_num sim_type_name_var)
	# NOTE set the target properties utilized only in target compilation
	# NOTE the PROPERTIES CUDA_ARCHITECTURES is required hint to compile for given architecture
	# GUIDE add your capability after the 50 i.e "50 70" this will compile for both of them

	# set_property(TARGET ${target} PROPERTY CUDA_ARCHITECTURES 50) # celnyd laptop GPU
	set_property(TARGET ${target} PROPERTY CUDA_ARCHITECTURES 61) # GTX 1060 Max-Q (klimam laptop GPU)
	# set_property(TARGET ${target} PROPERTY CUDA_ARCHITECTURES "80;86") # RTX 3090 (VSCHT GPU_PC)

	# NOTE CUDA_RESOLVE_DEVICE_SYMBOLS is mandatory for the correct linking otherwise the cuda specific symbols are not found
	set_target_properties(${target} PROPERTIES CUDA_RESOLVE_DEVICE_SYMBOLS ON)

	# NOTE add compilation options directly to nvcc
	# NOTE relocatable-device-code is mandatory by our implementation
	# NOTE generate-line-info is optional: Generate line-number information for device code
	target_compile_options(${target} PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:
	                       --generate-line-info
	                       --relocatable-device-code=true
						   # --maxrregcount 64 # optimized for RTX 3090
	                       >)
	# NOTE resolve the ouptut name with regard to the substance name and precision
	# NOTE utilizing the case where target names are substance names - easier to make
	if (USEFLOAT)
		target_compile_definitions(${target} PUBLIC PREC=${FLOAT_VAL})
		set(precision "float")		
	else()
		target_compile_definitions(${target} PUBLIC PREC=${DOUBLE_VAL})
		set(precision "double")
	endif()
	
	# NOTE resolve the simulation type from the name given during the configuration
	if (${sim_type_name_var} STREQUAL ${CONSTANT_BOX_NAME} ) 
		set(SIM_TYPE ${CONSTANT_BOX_VAL})
	elseif (${sim_type_name_var} STREQUAL ${EXPANSION_BOX_NAME} )
		set(SIM_TYPE ${EXPANSION_BOX_VAL})
	elseif (${sim_type_name_var} STREQUAL ${CONSTANT_BOX_BERENDSEN_NAME} )
		set(SIM_TYPE ${CONSTANT_BOX_BERENDSEN_VAL})
	elseif (${sim_type_name_var} STREQUAL ${MTK_NAME} )
		set(SIM_TYPE ${MTK_VAL})
	elseif (${sim_type_name_var} STREQUAL ${GEAR_NAME} )
		set(SIM_TYPE ${GEAR_VAL})
	else()
		message( FATAL_ERROR "Invalid simulation type encountered: ${sim_type_name_var} \nChoose one from following instead:\n${SIM_TYPE_OPTIONS}" )
	endif()
	
	# MESSAGE("PRECISION: " ${precision})
	# MESSAGE("SIM_TYPE:  " ${SIM_TYPE})

	target_compile_definitions(${target} PUBLIC SIMULATION_TYPE=${SIM_TYPE})

	# NOTE sets the name of final executable with proper denotion of configuration flags
	set_target_properties(${target} PROPERTIES OUTPUT_NAME "${PROJECT_NAME}_${SIM_TYPE_NAME}_${target}_${precision}")
	
	# NOTE set the compile flag that controls the substance selection
	# add_compile_definitions(SUBSTANCE=${substance_num})
	target_compile_definitions(${target} PUBLIC SUBSTANCE=${substance_num})
endmacro()

# NOTE this does not semms to change situation much for us atm we do not rely on higher standard features
#set(CMAKE_CUDA_STANDARD 14)
# NOTE this is mandatory to separate GPU, CPU code compilation
set(CUDA_SEPARABLE_COMPILATION ON)
# NOTE this is mandatory for the correct linking otherwise the cuda specific symbols are not found

# NOTE se the cmake options for compilation good for FLOAT switching in
# NOTE not so usefull for substances as they are often switched - different targets are better
option(USEFLOAT "Use the float arithmetics")
message(STATUS "USEFLOAT='${USEFLOAT}'")

# NOTE se the cmake options for different type of simulations
# NOTE if none is given than the one specified here is taken instead 
set(SIM_TYPE_NAME ${CONSTANT_BOX_NAME} CACHE STRING "Choose the simulation type during CMake configuration.")

# NOTE this is usefull for cmake-gui settings that one gets the option what to select
set_property(CACHE SIM_TYPE_NAME PROPERTY STRINGS ${SIM_TYPE_OPTIONS})

# NOTE information what configuration SIM_TYPE is selected
message(STATUS "SIM_TYPE='${SIM_TYPE_NAME}'")

# NOTE this set up ${SOURCES_ALL} variable holding all relevant source files
add_subdirectory(src)

# NOTE choose the proper substance for the target name
# NOTE need to choose what to make otherwise all will be done
# LJF 		= 					= 0 - DOES NOT WORK ATM
# argon 	= ${ARGON_VAL} 		= 1 ATM
# nitrogen 	= ${NITROGEN_VAL} 	= 2 ATM
# spce 		= ${SPCE_VAL} 		= 3 ATM

# add_executable("ljf" ${SOURCES_ALL})
# set_cuda_properties("ljf" 0)

add_executable("argon" ${SOURCES_ALL})
set_cuda_properties("argon" ${ARGON_VAL} ${SIM_TYPE_NAME})

add_executable("nitrogen" ${SOURCES_ALL})
set_cuda_properties("nitrogen" ${NITROGEN_VAL} ${SIM_TYPE_NAME})

add_executable("spce" ${SOURCES_ALL})
set_cuda_properties("spce" ${SPCE_VAL} ${SIM_TYPE_NAME})

# NOTE this makes sure usefloat and sim type name does not remain in chache
#      and default values are used instead
unset(USEFLOAT CACHE)
unset(SIM_TYPE_NAME CACHE)
