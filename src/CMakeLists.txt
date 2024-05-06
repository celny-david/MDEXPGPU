# NOTE all files we plan to use in the project have to be located in src and have .h or .cu suffix

FILE(GLOB SOURCES_CU ${CMAKE_SOURCE_DIR}/src/*.cu)
FILE(GLOB SOURCES_H ${CMAKE_SOURCE_DIR}/src/*.h)
SET (SOURCES_ALL ${SOURCES_CU} ${SOURCES_H} PARENT_SCOPE)
