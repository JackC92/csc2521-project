set(CMAKE_EXPORT_COMPILE_COMMANDS 1) # Emit a compile flags file to support completion engines, ignored by generators other than Makefile Generators and Ninja.

if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    # Using Clang or GCC
    message(STATUS "Using Clang/gcc compiler flags")
    set(BASE_CXX_FLAGS "-std=c++17 -Wall -Wextra -g3")
    set(DISABLED_WARNINGS " -Wno-unused-parameter -Wno-unused-variable -Wno-unused-function -Wno-deprecated-declarations -Wno-missing-braces")
    set(TRACE_INCLUDES " -H -Wno-error=unused-command-line-argument")
    
    if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
        message(STATUS "Setting Clang specific options")
        set(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} -ferror-limit=5 -fcolor-diagnostics")
        set(CMAKE_CXX_FLAGS_DEBUG "-fsanitize=address -fno-limit-debug-info")
    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        message(STATUS "Setting gcc specific options")
        set(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} -fmax-errors=5")
        set(DISABLED_WARNINGS "${DISABLED_WARNINGS} -Wno-maybe-uninitialized -Wno-format-zero-length -Wno-unused-but-set-parameter -Wno-unused-but-set-variable")
    endif()
    
    set(CMAKE_CXX_FLAGS "${BASE_CXX_FLAGS} ${DISABLED_WARNINGS}")
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TRACE_INCLUDES}") # uncomment if you need to track down where something is getting included from
    set(CMAKE_CXX_FLAGS_DEBUG          "${CMAKE_CXX_FLAGS_DEBUG} -g3")
    set(CMAKE_CXX_FLAGS_MINSIZEREL     "-Os -DNDEBUG")
    set(CMAKE_CXX_FLAGS_RELEASE        "-O3 -DNDEBUG -march=native")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")
    
    add_definitions(-D_USE_MATH_DEFINES)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    # Using Microsoft Visual C++
    message(STATUS "Using Microsoft Visual C++ compiler flags")
    set(BASE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
    set(BASE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP") # Parallel build
    set(DISABLED_WARNINGS "${DISABLED_WARNINGS} /wd\"4267\"")  # ignore conversion to smaller type (fires more aggressively than the gcc version, which is annoying)
    set(DISABLED_WARNINGS "${DISABLED_WARNINGS} /wd\"4244\"")  # ignore conversion to smaller type (fires more aggressively than the gcc version, which is annoying)
    set(DISABLED_WARNINGS "${DISABLED_WARNINGS} /wd\"4305\"")  # ignore truncation on initialization
    set(CMAKE_CXX_FLAGS "${BASE_CXX_FLAGS} ${DISABLED_WARNINGS} /arch:AVX2")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MD")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MDd")
    
    add_definitions(-D_USE_MATH_DEFINES)
else()
    # Unrecognized compiler
    message(FATAL_ERROR "Unrecognized compiler [${CMAKE_CXX_COMPILER_ID}]")
endif()
