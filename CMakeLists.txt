project(HydroCodeSpherical1D)

set(SOURCES
    main_spherical.cpp

    RiemannSolver.hpp
)

add_executable(HydroCodeSpherical1D ${SOURCES})
