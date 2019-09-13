#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
project(solver)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
# option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
# if(WITH_GEANT4_UIVIS)

#   find_package(Geant4 REQUIRED ui_all vis_all)
# else()
#   find_package(Geant4 REQUIRED)
# endif()

find_package(Threads)

# list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS RIO Net)

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
# include(${Geant4_USE_FILE})
if(ROOT_FOUND)
  include_directories(${ROOT_INCLUDE_DIRS} ${PROJECT_SOURCE_DIR}/include)
else()
  include_directories(${PROJECT_SOURCE_DIR}/include)
endif()


#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cpp)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.h)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#

add_executable( solver ${PROJECT_SOURCE_DIR}/src/solver.cxx ${sources} ${headers})
# target_link_libraries( solver ${Geant4_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} )
target_link_libraries( solver ${CMAKE_THREAD_LIBS_INIT} )

add_executable( drift ${PROJECT_SOURCE_DIR}/src/drift.cxx ${sources} ${headers})
# target_link_libraries( drift ${Geant4_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} )
target_link_libraries( drift ${CMAKE_THREAD_LIBS_INIT} )

add_executable( edep ${PROJECT_SOURCE_DIR}/src/edep.cxx ${sources} ${headers})
# target_link_libraries( edep ${Geant4_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} )
target_link_libraries( edep ${CMAKE_THREAD_LIBS_INIT} )

add_executable( efield ${PROJECT_SOURCE_DIR}/src/efield.cxx ${sources} ${headers})
# target_link_libraries( edep ${Geant4_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} )
target_link_libraries( efield ${CMAKE_THREAD_LIBS_INIT} )

add_executable( induction ${PROJECT_SOURCE_DIR}/src/induction.cxx ${sources} ${headers})
# target_link_libraries( edep ${Geant4_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} )
target_link_libraries( induction ${CMAKE_THREAD_LIBS_INIT} )

if(ROOT_FOUND)
  target_link_libraries(solver ${ROOT_LIBRARIES})
  target_link_libraries(drift ${ROOT_LIBRARIES})
  target_link_libraries(edep ${ROOT_LIBRARIES})
  target_link_libraries(efield ${ROOT_LIBRARIES})
  target_link_libraries(induction ${ROOT_LIBRARIES})
# else()
#   target_link_libraries(solver ${Geant4_LIBRARIES})
#   target_link_libraries(drift ${Geant4_LIBRARIES})
endif()

# install(TARGETS solver DESTINATION bin)
# install(TARGETS drift DESTINATION bin)

file(COPY geometries DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY utils DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
