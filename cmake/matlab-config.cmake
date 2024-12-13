# This minimal CMake Config Mode file for imporing MATLAB can be used anywhere if the environment variable MATLABROOT is set.
add_library(Matlab::Engine UNKNOWN IMPORTED)
set_target_properties(Matlab::Engine PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "$ENV{MATLABROOT}/extern/include"
    IMPORTED_LOCATION "$ENV{MATLABROOT}/extern/lib/win64/microsoft/libMatlabEngine.lib")
    
add_library(Matlab::DataArray UNKNOWN IMPORTED)
set_target_properties(Matlab::DataArray PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "$ENV{MATLABROOT}/extern/include"
    IMPORTED_LOCATION "$ENV{MATLABROOT}/extern/lib/win64/microsoft/libMatlabDataArray.lib")

add_library(Matlab::matlab INTERFACE IMPORTED)
set_property(TARGET Matlab::matlab PROPERTY
    INTERFACE_LINK_LIBRARIES Matlab::Engine Matlab::DataArray)

set(Matlab_INCLUDES $ENV{MATLABROOT}/extern/include)
set(Matlab_LIBRARIES Matlab::matlab)
