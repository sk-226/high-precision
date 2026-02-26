# ---------- テスト・ベンチマークプログラム -----------------------------------

# 基本的なQuadDouble演算テスト
add_executable(test_basic ${CMAKE_SOURCE_DIR}/src/benchmarks/test_basic.cpp)
target_include_directories(test_basic PRIVATE ${COMMON_INCLUDE_DIRS})
target_link_libraries(test_basic PRIVATE ${COMMON_LIBRARIES})
target_compile_features(test_basic PRIVATE cxx_std_17)

# Precision validation test
add_executable(precision_validation_test ${CMAKE_SOURCE_DIR}/src/precision_validation_test.cpp)
target_include_directories(precision_validation_test PRIVATE ${COMMON_INCLUDE_DIRS})
target_link_libraries(precision_validation_test PRIVATE ${COMMON_LIBRARIES})
target_compile_features(precision_validation_test PRIVATE cxx_std_17)

# Check long double properties
add_executable(check_ldbl ${CMAKE_SOURCE_DIR}/src/check_ldbl.cpp)
target_include_directories(check_ldbl PRIVATE ${COMMON_INCLUDE_DIRS})
target_link_libraries(check_ldbl PRIVATE ${COMMON_LIBRARIES})
target_compile_features(check_ldbl PRIVATE cxx_std_17)

# Sanity check program: (sqrt(2))^2 - 2 in DD/DQ/QX
add_executable(prec_sanity ${CMAKE_SOURCE_DIR}/src/prec_sanity.cpp)
target_include_directories(prec_sanity PRIVATE ${COMMON_INCLUDE_DIRS})
target_link_libraries(prec_sanity PRIVATE ${COMMON_LIBRARIES})
target_compile_features(prec_sanity PRIVATE cxx_std_17)

# Scalar precision regression test (DD/DQ/QX)
add_executable(scalar_precision_test scalar_precision_test.cpp)
target_include_directories(scalar_precision_test PRIVATE ${COMMON_INCLUDE_DIRS})
target_link_libraries(scalar_precision_test PRIVATE ${COMMON_LIBRARIES})
target_compile_features(scalar_precision_test PRIVATE cxx_std_17)
add_test(NAME scalar_precision_test COMMAND scalar_precision_test)

# Matrix Market loader smoke test (prints matrix in multiple precisions)
add_executable(matrix_loader_test matrix_loader_test.cpp)
target_include_directories(matrix_loader_test PRIVATE ${COMMON_INCLUDE_DIRS})
target_link_libraries(matrix_loader_test PRIVATE ${COMMON_LIBRARIES})
target_compile_features(matrix_loader_test PRIVATE cxx_std_17)

# Eigen-based arithmetic checks
# add_executable(eigen_arithmetic_test eigen_arithmetic_test.cpp)
# target_include_directories(eigen_arithmetic_test PRIVATE ${COMMON_INCLUDE_DIRS})
# target_link_libraries(eigen_arithmetic_test PRIVATE ${COMMON_LIBRARIES})
# target_compile_features(eigen_arithmetic_test PRIVATE cxx_std_17)

# Eigen-based high-precision arithmetic checks
add_executable(eigen_precision_test eigen_precision_test.cpp)
target_include_directories(eigen_precision_test PRIVATE ${COMMON_INCLUDE_DIRS})
target_link_libraries(eigen_precision_test PRIVATE ${COMMON_LIBRARIES})
target_compile_features(eigen_precision_test PRIVATE cxx_std_17)
add_test(NAME eigen_precision_test COMMAND eigen_precision_test)

# # DQ to QX direct cast test
# add_executable(test_dq_to_qx_cast test_dq_to_qx_cast.cpp)
# target_include_directories(test_dq_to_qx_cast PRIVATE ${COMMON_INCLUDE_DIRS})
# target_link_libraries(test_dq_to_qx_cast PRIVATE ${COMMON_LIBRARIES})
# target_compile_features(test_dq_to_qx_cast PRIVATE cxx_std_17)
# add_test(NAME test_dq_to_qx_cast COMMAND test_dq_to_qx_cast)

# # 1/3 precision comparison (double, DD, DQ)
# add_executable(1div3 1div3.f90)
# target_include_directories(1div3 PRIVATE ${COMMON_INCLUDE_DIRS})
# target_compile_options(1div3 PRIVATE 
#     $<$<COMPILE_LANGUAGE:Fortran>:-I$ENV{DDFUN_DIR} -I$ENV{DQFUN_DIR}>)
# target_link_libraries(1div3 PRIVATE ${COMMON_LIBRARIES})
