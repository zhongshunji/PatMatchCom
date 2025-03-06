cmake_minimum_required(VERSION 3.10)

project(PatMatchCom)
set(CMAKE_CXX_STANDARD 11)

if(WIN32)
    find_package(unofficial-brotli CONFIG REQUIRED)
else()
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(BrotliEnc REQUIRED IMPORTED_TARGET libbrotlienc)
    pkg_check_modules(BrotliDec REQUIRED IMPORTED_TARGET libbrotlidec)
    pkg_check_modules(BrotliCommon REQUIRED IMPORTED_TARGET libbrotlicommon)
endif()

set(SOURCES
    main.cpp
    huffman.cpp
    Literal.cpp
)

add_executable(PatMatchCom ${SOURCES})

if(WIN32)
    target_link_libraries(PatMatchCom PRIVATE unofficial::brotli::brotlidec unofficial::brotli::brotlienc unofficial::brotli::brotlicommon)
else()
    target_link_libraries(PatMatchCom PRIVATE PkgConfig::BrotliEnc PkgConfig::BrotliDec PkgConfig::BrotliCommon)
endif()
