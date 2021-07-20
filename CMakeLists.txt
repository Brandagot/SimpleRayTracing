cmake_minimum_required(VERSION 3.12)

project(ICE4131-202021-Assignment)


# Enable C++ 11################################################################
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)


# JPEG ########################################################################
FIND_PACKAGE(JPEG REQUIRED)
SET(requiredLibs ${requiredLibs} ${JPEG_LIBRARY})


# ASSIMP ######################################################################
OPTION(USE_SYSTEM_ASSIMP "Use system library" OFF)

# Use system library
IF (USE_SYSTEM_ASSIMP)
    FIND_PACKAGE( assimp )
ENDIF (USE_SYSTEM_ASSIMP)

# The system library was not found
# Build ASSIMP instead
IF (NOT assimp_FOUND OR NOT USE_SYSTEM_ASSIMP)
    MESSAGE (STATUS "Build ASSIMP")
    SET (USE_SYSTEM_ASSIMP OFF)
    include(ExternalProject)
    INCLUDE(BuildASSIMP.cmake)
    build_assimp( assimp ${PROJECT_BINARY_DIR} )
ENDIF (NOT assimp_FOUND OR NOT USE_SYSTEM_ASSIMP)


# ZLIB ######################################################################
FIND_PACKAGE( ZLIB REQUIRED )
SET(requiredLibs ${requiredLibs} ${ZLIB_LIBRARY})

# OpenMP ####################################################################
find_package(OpenMP)

# Build RayTracing library ##################################################
add_library(RayTracing
  include/Image.h
  include/Image.inl
  src/Image.cxx
  include/Ray.h
  include/Ray.inl
  src/Ray.cxx
  include/Triangle.h
  include/Triangle.inl
  src/Triangle.cxx
  include/TriangleMesh.h
  include/TriangleMesh.inl
  src/TriangleMesh.cxx
  include/Vec3.h
  include/Vec3.inl
)

TARGET_INCLUDE_DIRECTORIES(RayTracing PRIVATE ${JPEG_INCLUDE_DIR})
TARGET_INCLUDE_DIRECTORIES(RayTracing PUBLIC ${ASSIMP_INCLUDE_DIRS})
TARGET_INCLUDE_DIRECTORIES(RayTracing PRIVATE ${ZLIB_INCLUDE_DIRS})
TARGET_INCLUDE_DIRECTORIES(RayTracing PUBLIC include)

TARGET_LINK_LIBRARIES (RayTracing ${requiredLibs} ${optionalLibs})

IF (NOT USE_SYSTEM_ASSIMP)
    add_dependencies (RayTracing assimp)
ENDIF (NOT USE_SYSTEM_ASSIMP)


# Program ###############################################################
add_executable(main src/main.cxx)
TARGET_LINK_LIBRARIES (main PUBLIC RayTracing ${ASSIMP_LIBRARY})
if(OpenMP_CXX_FOUND)
    TARGET_LINK_LIBRARIES(main PUBLIC OpenMP::OpenMP_CXX)
endif()

#FILE(COPY cloud2.jpg DESTINATION ${CMAKE_BINARY_DIR})
FILE(COPY Bangor_Logo_A1.jpg DESTINATION ${CMAKE_BINARY_DIR})
FILE(COPY dragon.ply DESTINATION ${CMAKE_BINARY_DIR})
SET(CMAKE_THREAD_PREFER_PTHREAD TRUE)
find_package(Threads REQUIRED)
add_executable(main-pthreads src/main-pthreads.cxx)
TARGET_LINK_LIBRARIES(main-pthreads PUBLIC RayTracing ${ASSIMP_LIBRARY} Threads::Threads)

add_executable(main-pthreads-redo src/main-pthreads-redo.cxx)
TARGET_LINK_LIBRARIES(main-pthreads-redo PUBLIC RayTracing ${ASSIMP_LIBRARY} Threads::Threads)

add_executable(main-pthreads-lbuffer src/main-pthreads-lbuffer.cxx)
TARGET_LINK_LIBRARIES(main-pthreads-lbuffer PUBLIC RayTracing ${ASSIMP_LIBRARY} Threads::Threads)

find_package(OpenMP REQUIRED)

add_executable(main-omp src/main-omp.cxx)
TARGET_INCLUDE_DIRECTORIES(main-omp PUBLIC ${ASSIMP_INCLUDE_DIRS})
TARGET_LINK_LIBRARIES(main-omp PUBLIC RayTracing ${ASSIMP_LIBRARY})
if(OpenMP_CXX_FOUND)
    TARGET_LINK_LIBRARIES(main-omp PUBLIC OpenMP::OpenMP_CXX)
endif()

find_package(MPI REQUIRED)
add_executable(main-mpi src/main-mpi.cxx)
TARGET_INCLUDE_DIRECTORIES(main-mpi PUBLIC ${ASSIMP_INCLUDE_DIRS} ${MPI_INCLUDE_PATH})
TARGET_LINK_LIBRARIES(main-mpi PUBLIC RayTracing ${ASSIMP_LIBRARY} ${MPI_CXX_LIBRARIES})

add_executable(main-mpi-omp src/main-mpi-omp.cxx)
TARGET_INCLUDE_DIRECTORIES(main-mpi-omp PUBLIC ${ASSIMP_INCLUDE_DIRS})
TARGET_LINK_LIBRARIES(main-mpi-omp PUBLIC RayTracing ${ASSIMP_LIBRARY})
if(OpenMP_CXX_FOUND)
    TARGET_LINK_LIBRARIES(main-mpi-omp PUBLIC OpenMP::OpenMP_CXX)
endif()
TARGET_INCLUDE_DIRECTORIES(main-mpi-omp PUBLIC ${ASSIMP_INCLUDE_DIRS} ${MPI_INCLUDE_PATH})
TARGET_LINK_LIBRARIES(main-mpi-omp PUBLIC RayTracing ${ASSIMP_LIBRARY} ${MPI_CXX_LIBRARIES})

