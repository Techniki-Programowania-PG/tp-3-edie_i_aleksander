# Utwórz moduł pybind11 z main.cpp
pybind11_add_module(example
    main.cpp
)

# Dołącz lokalny katalog z nagłówkami
target_include_directories(example PRIVATE ${CMAKE_SOURCE_DIR}/src)

# Dodaj matplot++ do linkowania
target_link_libraries(example PRIVATE matplot)

