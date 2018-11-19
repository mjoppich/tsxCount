#Splitter Application Package - some toolkit to work with gtf/gff files
#Copyright (C) 2015  Markus Joppich
#
#The Splitter Application Package is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#The Splitter Application Package is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.

cmake_minimum_required(VERSION 3.1)
SET (CMAKE_C_COMPILER             "/usr/bin/clang")
SET (CMAKE_C_FLAGS                "-Wall -std=c99")
SET (CMAKE_C_FLAGS_DEBUG          "-g")
SET (CMAKE_C_FLAGS_MINSIZEREL     "-Os -DNDEBUG")
SET (CMAKE_C_FLAGS_RELEASE        "-O4 -DNDEBUG")
SET (CMAKE_C_FLAGS_RELWITHDEBINFO "-O2 -g")

SET (CMAKE_CXX_COMPILER             "/usr/bin/clang++")
SET (CMAKE_CXX_FLAGS                "-Wall")
SET (CMAKE_CXX_FLAGS_DEBUG          "-g")
SET (CMAKE_CXX_FLAGS_MINSIZEREL     "-Os -DNDEBUG")
SET (CMAKE_CXX_FLAGS_RELEASE        "-O4 -DNDEBUG")
SET (CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")
set(CMAKE_CXX_STANDARD 11)


cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)
set(CMAKE_VERBOSE_MAKEFILE ON)
project (tsxCount)




set(CMAKE_CXX_STANDARD 11)

message(STATUS "Checking OpenMP")
find_package(OpenMP REQUIRED)
IF(OPENMP_FOUND)
    message("Found OpenMP!")
    # add flags for OpenMP
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    add_compile_options("-mrtm")
    add_compile_options("-fopenmp")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${OpenMP_SHARED_LINKER_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
ELSE()
    message("Missed OpenMP!")
ENDIF()

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

MESSAGE(WARNING ${CMAKE_SYSTEM_NAME})

if ( CMAKE_SYSTEM_NAME STREQUAL "Windows" OR CYGWIN OR MINGW)

    if (MINGW)

    endif()

else()

endif()

if (MINGW)
else()
endif()

find_package( ZLIB REQUIRED )
if ( ZLIB_FOUND )
    include_directories( ${ZLIB_INCLUDE_DIRS} )
    MESSAGE(WARNING ${ZLIB_LIBRARY})
endif( ZLIB_FOUND )

file(GLOB tsxcountlib_h
        "src/tsxcount/*.h"
)

file(GLOB tsxcountlib_SRC
        "src/tsxcount/*.cpp"
)

MESSAGE(${tsxcountlib_SRC})

add_library(tsxcounter ${tsxcountlib_SRC} ${tsxcountlib_h} )
target_include_directories (tsxcounter PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS tsxcounter DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/lib)
install(FILES ${tsxcountlib_h} DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/include/)

file(GLOB tsxutilslib_h
        "src/tsxutils/*.h"
        "src/utils/*.h"
        "src/fastxutils/*.h"
        )

file(GLOB tsxutilslib_SRC
        "src/tsxutils/*.cpp"
        "src/utils/*.cpp"
        "src/fastxutils/*.cpp"
        )

add_library(tsxutils ${tsxutilslib_SRC} ${tsxutilslib_h})
TARGET_LINK_LIBRARIES(tsxutils z Threads::Threads)

target_include_directories (tsxcounter PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS tsxutils DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/lib)
install(FILES ${tsxutilslib_h} DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/include/)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../splitter/include)
ADD_LIBRARY(splitutils STATIC IMPORTED src/utils/CLParser.cpp src/utils/CLParser.h)
SET_TARGET_PROPERTIES(splitutils PROPERTIES IMPORTED_LOCATION ${CMAKE_CURRENT_SOURCE_DIR}/../splitter/lib/libsplitutils.a)

# Main app
file(GLOB tsxCount_H
        "src/mains/testExecution.h"
        )

file(GLOB tsxCount_SRC
        "src/mains/main.cpp"
        )
ADD_EXECUTABLE(tsxCount ${tsxCount_SRC} ${tsxCount_H})
target_include_directories(tsxCount PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
TARGET_LINK_LIBRARIES(tsxCount tsxcounter tsxutils Threads::Threads)
install(TARGETS tsxCount DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/bin)


# Main app
file(GLOB tsxCountPTHREAD_H
        "src/mains/testExecution.h"
        )

file(GLOB tsxCountPTHREAD_SRC
        "src/mains/main_pthread.cpp"
        )
ADD_EXECUTABLE(tsxCountPThread ${tsxCountPTHREAD_SRC} ${tsxCountPTHREAD_H})
target_include_directories(tsxCountPThread PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
TARGET_LINK_LIBRARIES(tsxCountPThread PRIVATE ${OpenMP_CXX_FLAGS} tsxcounter tsxutils Threads::Threads)
install(TARGETS tsxCountPThread DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/bin)



# Main app
file(GLOB tsxCountOMP_H
        "src/mains/testExecution.h"
        )

file(GLOB tsxCountOMP_SRC
        "src/mains/main_omp.cpp"
        )
ADD_EXECUTABLE(tsxCountOMP ${tsxCountOMP_SRC} ${tsxCountOMP_H})
target_include_directories(tsxCountOMP PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
TARGET_LINK_LIBRARIES(tsxCountOMP PRIVATE ${OpenMP_CXX_FLAGS} tsxcounter tsxutils Threads::Threads)

install(TARGETS tsxCountOMP DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/bin)


# Main app
file(GLOB tsxCountTSX_H
        "src/mains/testExecution.h"
        )

file(GLOB tsxCountTSX_SRC
        "src/mains/main_tsx.cpp"
        )
ADD_EXECUTABLE(tsxCountTSX ${tsxCountTSX_SRC} ${tsxCountTSX_H})
target_include_directories(tsxCountTSX PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)

target_compile_options(tsxCountTSX PUBLIC -Wall ${OpenMP_CXX_FLAGS} -mrtm -fopenmp)

TARGET_LINK_LIBRARIES(tsxCountTSX PUBLIC ${OpenMP_CXX_FLAGS} tsxcounter tsxutils Threads::Threads )
install(TARGETS tsxCountTSX DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/bin)