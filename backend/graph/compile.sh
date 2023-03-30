
# make shared library
g++ -fPIC -shared -o libelena.so main.cpp include/graph.cpp include/utils.cpp -std=c++17 -stdlib=libc++ -fdeclspec

# compile for run
g++ -o main main.cpp include/graph.cpp include/utils.cpp -std=c++17 -stdlib=libc++ -fdeclspec

# run main.cpp: ./main