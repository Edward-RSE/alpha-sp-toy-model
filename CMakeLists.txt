cmake_minimum_required(VERSION 3.5)
project(num-int LANGUAGES C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_COMPILER gcc)

include_directories(include)
link_directories(lib)

add_executable(num-int
        src/python/atomic_extern_init.c
        src/python/atomicdata.c
        src/python/atomicdata_init.c
        src/python/atomicdata_sub.c
        src/python/python_extern_init.c
        src/python/rdpar.c
        src/python/recipes.c
        src/python/synonyms.c
        src/python/xlog.c
        src/main.c
        src/alpha_sp.c
        include/integrate.h
        src/integrate.c
)

target_link_libraries(num-int m gsl)
