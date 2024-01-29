cmake_minimum_required(VERSION 3.5)
project(num-int LANGUAGES C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_COMPILER mpicc)

include_directories(inc)
link_directories(lib)

set(PYTHON_SOURCE
        src/python/atomic_extern_init.c
        src/python/atomicdata.c
        src/python/atomicdata_init.c
        src/python/atomicdata_sub.c
        src/python/python_extern_init.c
        src/python/rdpar.c
        src/python/recipes.c
        src/python/synonyms.c
        src/python/xlog.c
)

add_executable(num-int
        ${PYTHON_SOURCE}
        src/num-int/num_int.c
        src/num-int/alpha_sp.c
        src/num-int/integrate.c
)

add_executable(node-share
        ${PYTHON_SOURCE}
        src/node-share/node_share.c
)

target_link_libraries(num-int m gsl)
target_link_libraries(node-share m gsl)
