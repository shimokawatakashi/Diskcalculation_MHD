# 磁気拡散係数 (eta_O, eta_H, eta_A) の nH x T x B テーブル作成
#   make            : eta_table をビルド
#   make table      : ./input/input でテーブルを作成 (INPUT=... で入力ファイル変更)
#   make clean      : obj/, eta_tableを削除

CXX      = g++
# CXX      = clang++
CXXFLAGS = -O3 -std=c++17 -I$(SRCDIR)

SRCDIR  = src
OBJDIR  = obj
TARGET  = eta_table
SOLVER  = calc_resistivity.cpp calculation_tool.cpp
HEADERS = $(addprefix $(SRCDIR)/,calc_resistivity.hpp calculation_tool.hpp typedef.hpp physical_constants.hpp)
SOLVER_OBJECTS = $(patsubst %.cpp,$(OBJDIR)/%.o,$(SOLVER))

INPUT ?= ./input/input

all: $(TARGET)

$(TARGET): $(OBJDIR)/make_eta_table.o $(SOLVER_OBJECTS)
	$(CXX) $^ -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp $(HEADERS)
	@mkdir -p $(OBJDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

table: $(TARGET)
	./$(TARGET) $(INPUT)

clean:
	rm -rf $(OBJDIR) $(TARGET) $(TEST)

.PHONY: all table clean
