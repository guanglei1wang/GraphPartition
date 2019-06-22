include_directories(${PROJECT_INCLUDE_DIR})

file(GLOB MINK_SOURCES ${PROJECT_SOURCE_DIR}/examples/MILP/MinKPart/*.cpp)
#add_executable(gravity101 ${PROJECT_SOURCE_DIR}/examples/Gravity101_main.cpp)
add_executable(min_k_part ${MINK_SOURCES})
#add_executable(fusions ${PROJECT_SOURCE_DIR}/examples/MILP/fusions/Fusion_main.cpp)
target_link_libraries(min_k_part gravity)
