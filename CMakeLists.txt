cmake_minimum_required(VERSION 2.6) # compatibilités de CMake
project(C++) # nom du projet

set(CMAKE_BUILD_TYPE debug) # compilation en mode debug
set(CMAKE_CXX_FLAGS "-Wall -Wextra -std=c++11") # options
add_subdirectory(doc)
add_subdirectory(tests)
enable_testing ()
add_test (NAME MyTest  
          COMMAND test_mylib
	  )
file( # liste des fichiers à compiler
	GLOB_RECURSE # recherche récursive
	source_files # fichiers listés dans une variable source_files
	src/main.cpp #src/*.[ch]pp # forme des fichiers à rechercher
)

add_executable( # création de l’exécutable binaire du projet
	create_raster # nom du binaire
	${source_files} # liste des fichiers compilés pour le linking
)


