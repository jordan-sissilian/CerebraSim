set(TEST_FILES
    test_IonChannel.cpp
)

add_executable(run_tests ${TEST_FILES})

find_package(GTest REQUIRED)

target_link_libraries(run_tests CSpp GTest::GTest GTest::Main pthread)

# Activer les tests
enable_testing()
add_test(NAME RunUnitTests COMMAND run_tests)
