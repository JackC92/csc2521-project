# This minimal CMake Config Mode file for imporing TensorFlow should be copied to the root of the libtensorflow-(cpu|gpu) folder
get_filename_component(_TensorFlow_PATH "${CMAKE_CURRENT_LIST_FILE}" PATH)

set(_TensorFlow_INCLUDE_DIRS "${_TensorFlow_PATH}/include")
if (MSVC)
    set(_TensorFlow_LIB "${_TensorFlow_PATH}/lib/tensorflow.lib")
    add_library(TensorFlow::tensorflow STATIC IMPORTED)
    set_target_properties(TensorFlow::tensorflow PROPERTIES
        IMPORTED_LOCATION "${_TensorFlow_PATH}/lib/tensorflow.lib")
    add_library(TensorFlow::TensorFlow INTERFACE IMPORTED)
    set_target_properties(TensorFlow::TensorFlow PROPERTIES
            INTERFACE_INCLUDE_DIRECTORIES "${_TensorFlow_INCLUDE_DIRS}"
            INTERFACE_LINK_LIBRARIES TensorFlow::tensorflow)
else()
    set(_TensorFlow_LIB "${_TensorFlow_PATH}/lib/libtensorflow.so" "${_TensorFlow_PATH}/lib/libtensorflow_framework.so")
    add_library(TensorFlow::tensorflow SHARED IMPORTED)
    set_target_properties(TensorFlow::tensorflow PROPERTIES
        IMPORTED_LOCATION "${_TensorFlow_PATH}/lib/libtensorflow.so")
    add_library(TensorFlow::framework SHARED IMPORTED)
    set_target_properties(TensorFlow::framework PROPERTIES
        IMPORTED_LOCATION "${_TensorFlow_PATH}/lib/libtensorflow_framework.so")
    add_library(TensorFlow::TensorFlow INTERFACE IMPORTED)
    set_property(TARGET TensorFlow::TensorFlow PROPERTY
        INTERFACE_INCLUDE_DIRECTORIES "${_TensorFlow_INCLUDE_DIRS}")
    set_property(TARGET TensorFlow::TensorFlow PROPERTY
        INTERFACE_LINK_LIBRARIES TensorFlow::tensorflow TensorFlow::framework)
endif()

if (MSVC)
    add_custom_command(
        TARGET core
        POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${_TensorFlow_PATH}/lib/tensorflow.dll"
        "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/$<CONFIGURATION>")
endif()

set(TensorFlow_INCLUDES _TensorFlow_INCLUDE_DIRS)
set(TensorFlow_LIBRARIES _TensorFlow_LIB)

unset(_TensorFlow_PATH)
