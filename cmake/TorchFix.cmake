if (DEFINED ENV{TORCH_INSTALL_PREFIX})
    set(TORCH_INSTALL_PREFIX $ENV{TORCH_INSTALL_PREFIX})
else()
    # Assume we are in <install-prefix>/share/cmake/Torch/TorchConfig.cmake
    get_filename_component(TORCH_INSTALL_PREFIX "${Torch_DIR}/../../../" ABSOLUTE)
endif()
# caffe2::mkl somehow adds the MKL files by their filename instead of their absolute path, thus causing trouble for linking.
if (TARGET caffe2::mkl)
    # This should be available with Intel MKL enabled.
    set_target_properties(caffe2::mkl PROPERTIES INTERFACE_LINK_LIBRARIES "${mkl_core_file};${mkl_intel_lp64_file};${mkl_intel_thread_file}")
endif()
if (MSVC)
    file(GLOB TORCH_DLLS "${TORCH_INSTALL_PREFIX}/lib/*.dll")
    add_custom_command(
        TARGET core
        POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        ${TORCH_DLLS}
        "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/$<CONFIGURATION>")
endif()
