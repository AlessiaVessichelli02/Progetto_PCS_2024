list (APPEND src_sources "${CMAKE_CURRENT_SOURCE_DIR}/Project.cpp")
list (APPEND src_headers "${CMAKE_CURRENT_SOURCE_DIR}/Project.hpp")
list (APPEND src_headers "${CMAKE_CURRENT_SOURCE_DIR}/DFN.hpp")

set (src_sources ${src_sources} PARENT_SCOPE)
set (src_headers ${src_headers} PARENT_SCOPE)
set (src_includes ${CMAKE_CURRENT_SOURCE_DIR} PARENT_SCOPE)

