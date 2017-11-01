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

set(CMAKE_CXX_STANDARD 11)


cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)

project (tsxCount)
set(CMAKE_CXX_STANDARD 11)

find_package(OpenMP)
if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}") #-pg
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} ")
endif()

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
        "src/tsxcount/*.cpp"
)

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
TARGET_LINK_LIBRARIES(tsxutils z)

target_include_directories (tsxcounter PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS tsxutils DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/lib)
install(FILES ${tsxutilslib_h} DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/include/)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../splitter/include)
ADD_LIBRARY(splitutils STATIC IMPORTED src/utils/CLParser.cpp src/utils/CLParser.h)
SET_TARGET_PROPERTIES(splitutils PROPERTIES IMPORTED_LOCATION ${CMAKE_CURRENT_SOURCE_DIR}/../splitter/lib/libsplitutils.a)

# Main app
file(GLOB tsxCount_H
        "src/*.h"
        )

file(GLOB tsxCount_SRC
        "src/*.cpp"
        )
ADD_EXECUTABLE(tsxCount ${tsxCount_SRC} ${tsxCount_H})
target_include_directories(tsxCount PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
TARGET_LINK_LIBRARIES(tsxCount tsxcounter tsxutils)
install(TARGETS tsxCount DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/bin)
