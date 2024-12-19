cmake_minimum_required(VERSION 3.10)

# Nom du projet
project(CSpp)

cmake_policy(SET CMP0079 NEW)

set(CMAKE_CXX_STANDARD 23)

# Définir un répertoire de sortie pour les fichiers compilés
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Ajouter un répertoire de recherche pour les fichiers d'en-tête
include_directories(${CMAKE_SOURCE_DIR}/inc)
include_directories(${CMAKE_SOURCE_DIR}/inc/neuron)
include_directories(${CMAKE_SOURCE_DIR}/inc/ionChannel)

include_directories(src)

# Ajouter les fichiers sources du projet
set(SOURCE_FILES
    src/ionChannel/IonChannel.cpp
)

# Créer une bibliothèque statique CSpp à partir des fichiers sources
add_library(CSpp STATIC ${SOURCE_FILES})

add_subdirectory(test)

find_package(GTest REQUIRED)

target_link_libraries(run_tests CSpp GTest::GTest GTest::Main pthread)

enable_testing()
add_test(NAME RunUnitTests COMMAND run_tests)
