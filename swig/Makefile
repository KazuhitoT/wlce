
CXX := g++
SWIG_PATH = /usr/local/bin/
INCLUDE := -I/usr/local/include/eigen3
INCLUDEPYTHON := -I/usr/include/python2.7

SRC_DIR = ../src/
OBJ_DIR = ./build
CXXFLAGS := -std=c++11 -c -fPIC -DEIGEN_NO_DEBUG

UNAME := $(shell uname)

ifeq ($(uname), linux)
	LDCOMMAND = $(CXX) -shared 
else
	LDCOMMAND = ld -bundle -flat_namespace -undefined suppress
endif

ifeq ($(CXX), icpc)
	CXXFLAGS = $(CXXFLAGS) -fast
endif

.PHONY : all
 
all :
	@if [ ! -d $(OBJ_DIR) ]; \
		then echo "mkdir -p $(OBJ_DIR)"; mkdir -p $(OBJ_DIR); \
	fi
	make pyconf_wrap.cpp
	make $(OBJ_DIR)/pyconf.o
	make $(OBJ_DIR)/metroconf.o
	make $(OBJ_DIR)/conf2corr.o
	make $(OBJ_DIR)/input.o
	make $(OBJ_DIR)/parser.o
	make $(OBJ_DIR)/pyconf_wrap.o
	make _pyconf.so

pyconf_wrap.cpp : pyconf.i
	$(SWIG_PATH)/swig -c++ -python -o pyconf_wrap.cpp pyconf.i

$(OBJ_DIR)/pyconf.o : pyconf.cpp
	$(CXX) -c $^ -o $@ $(CXXFLAGS) $(INCLUDE)

$(OBJ_DIR)/metroconf.o : $(SRC_DIR)/metroconf.cpp
	$(CXX) -c $^ -o $@ $(CXXFLAGS) $(INCLUDE)

$(OBJ_DIR)/conf2corr.o : $(SRC_DIR)/conf2corr.cpp
	$(CXX) -c $^ -o $@ $(CXXFLAGS) $(INCLUDE)

$(OBJ_DIR)/input.o : $(SRC_DIR)/input.cpp
	$(CXX) -c $^ -o $@ $(CXXFLAGS) $(INCLUDE)

$(OBJ_DIR)/parser.o : $(SRC_DIR)/parser.cpp
	$(CXX) -c $^ -o $@ $(CXXFLAGS) $(INCLUDE)

$(OBJ_DIR)/pyconf_wrap.o : pyconf_wrap.cpp
	$(CXX) -c $^ -o $@ $(CXXFLAGS) $(INCLUDE) $(INCLUDEPYTHON)

_pyconf.so : $(OBJ_DIR)/*.o
	$(LDCOMMAND) -o _pyconf.so $(OBJ_DIR)/*.o

clean:
	rm -rf $(OBJ_DIR) _pyconf.so
