cmake_minimum_required(VERSION 3.14)
project(MyProject)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_subdirectory(matplotplusplus)
list(APPEND CMAKE_PREFIX_PATH ${CMAKE_SOURCE_DIR}"/matplotplusplus/install")

find_package(pybind11 REQUIRED)

pybind11_add_module(biblioteka main2.cpp biblioteka.cpp)
target_link_libraries(biblioteka PRIVATE matplot)

#Kopiowanie biblioteki do osobnego folderu z programem testowym python

file(COPY ${CMAKE_BINARY_DIR}/Debug/biblioteka.cp313-win_amd64.pyd DESTINATION ${CMAKE_BINARY_DIR}/resources)

# Ścieżki
set(SRC_FILE "${CMAKE_BINARY_DIR}/Debug/biblioteka.cp313-win_amd64.pyd")
set(DST_DIR "D:/OneDrive/szkoła/Projekt pp/python test")

# Komenda kopiująca plik
add_custom_command(
    OUTPUT ${DST_DIR}/biblioteka.cp313-win_amd64.pyd
    COMMAND ${CMAKE_COMMAND} -E make_directory ${DST_DIR}
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_FILE} ${DST_DIR}/biblioteka.cp313-win_amd64.pyd
    DEPENDS ${SRC_FILE}
    COMMENT "Kopiowanie biblioteki do katalogu wyjściowego"
)

# Cel kopiowania
add_custom_target(copy_config ALL
    DEPENDS ${DST_DIR}/biblioteka.cp313-win_amd64.pyd
)
