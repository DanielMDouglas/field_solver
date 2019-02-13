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
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()


# list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
# find_package(ROOT REQUIRED COMPONENTS RIO Net)

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
# if(ROOT_FOUND)
#   include_directories(${ROOT_INCLUDE_DIRS} ${PROJECT_SOURCE_DIR}/include)
# else()
  include_directories(${PROJECT_SOURCE_DIR}/include)
# endif()


#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cpp)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.h)


#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(solver ${sources} ${headers})
# if(ROOT_FOUND)
#   target_link_libraries(G4Sim ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})
# else()
target_link_libraries(solver ${Geant4_LIBRARIES})
# endif()

install(TARGETS solver DESTINATION bin)
