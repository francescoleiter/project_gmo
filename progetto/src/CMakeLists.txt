add_library(${LIB_NAME} SHARED ${SRC})
target_include_directories(${LIB_NAME} PUBLIC ${INCLUDE_DIR})
add_executable(main ${MAIN})
target_include_directories(main PUBLIC ${INCLUDE_DIR})
target_link_libraries(main ${LIB_NAME})

if(${pybind11_FOUND})
    include_directories(SYSTEM ${pybind11_INCLUDE_DIRS})
    pybind11_add_module(${PY_LIB_NAME} ${PY_SRC})
    target_include_directories(${PY_LIB_NAME} PRIVATE include())
    target_link_libraries(${PY_LIB_NAME} PRIVATE ${LIB_NAME})
    set_target_properties(${PY_LIB_NAME} PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${PYLIB_INSTALL_PATH}")
endif()