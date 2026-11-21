#コンパイラの設定 (MPI 用)
CXX = g++       # MPI/ OPEN MPを使用しないとき
# CXX = mpic++  # MPI を使用するとき
#CXX = g++-14  # OPEN MP のみを使用するとき

# コンパイルフラグ (OpenMP を有効化, LTO 無効化)
#CXXFLAGS = -g -O3 -std=c++17 -I$(SRCDIR)
CXXFLAGS = -g -O3 -std=c++17 -fopenmp -I$(SRCDIR)

# リンクフラグ (MPI & OpenMP & GNU ld)
#LDFLAGS = 
LDFLAGS = -fopenmp # OPEN MP を使用するとき

# ディレクトリ設定
SRCDIR = src
MARCHAND_DIR = Marchand2022/src
OBJDIR = obj
BINDIR = bin
TARGET = $(BINDIR)/run

# ソース & オブジェクトファイル一覧
SOURCES = main.cpp utils.cpp dust_disk.cpp dust_growth.cpp alpha.cpp magnetic_field1D.cpp energy.cpp face_scheme.cpp
MARCHAND_SOURCES = calc_resistivity.cpp calculation_tool.cpp

OBJECTS = $(patsubst %.cpp, $(OBJDIR)/%.o, $(SOURCES))
MARCHAND_OBJECTS = $(patsubst %.cpp, $(OBJDIR)/marchand_%.o, $(MARCHAND_SOURCES))

# 実行ファイル作成
all: $(TARGET)

# MPIを使用するとき
# $(TARGET): $(OBJECTS) $(MARCHAND_OBJECTS)
# 	@mkdir -p $(BINDIR)
# 	$(CXX) $^ -o $@ $(LDFLAGS) -lmpi
# OPEN MPのみを使用するとき
$(TARGET): $(OBJECTS) $(MARCHAND_OBJECTS)
	@mkdir -p $(BINDIR)
	$(CXX) $^ -o $@ $(LDFLAGS)

# コンパイルルール（src/）
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	@mkdir -p $(OBJDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# コンパイルルール（Marchand2022/）
$(OBJDIR)/marchand_%.o: $(MARCHAND_DIR)/%.cpp
	@mkdir -p $(OBJDIR)
	$(CXX) $(CXXFLAGS) -I$(MARCHAND_DIR) -c $< -o $@
	
#実行ルール (スレッド数を指定)
ARGS ?= input
YEAR ?= 320000
# BASE : 計算結果の置き場所 , RUN : 再開するモデル名 (BASE/RUN と BASE/RUN_B を読む)
BASE ?= /mnt/hdd/shimokawa
RUN  ?= test
DUST ?= $(BASE)/$(RUN)/out$(YEAR)
GAS  ?= $(BASE)/$(RUN)/out$(YEAR)_gas
MAG  ?= $(BASE)/$(RUN)_B/out$(YEAR)_B
OPTION ?= -r

run:
	export OMP_PROC_BIND=TRUE
	export OMP_PLACES=cores
	export OMP_NUM_THREADS=4; ./$(TARGET) $(ARGS)

restart:
	export OMP_PROC_BIND=TRUE
	export OMP_PLACES=cores
	export OMP_NUM_THREADS=4; ./$(TARGET) $(OPTION) $(ARGS) $(DUST) $(GAS) $(MAG) $(YEAR)

# 一定時間ごとに休憩を挟みながら自動でチェックポイント再開を繰り返す実行
# (ARGS に input_folder 内のファイルを指定。SEGMENT_HOURS / COOLDOWN_MIN / OMP_NUM_THREADS で上書き可)
supervised: $(TARGET)
	bash scripts/run_supervised.sh $(ARGS)

# 毎日決まった時間帯だけ動かす実行 (内部で supervised と同じ仕組みを使う)
# (ARGS に input_folder 内のファイルを指定。WINDOW_START / WINDOW_END / GRACE_HOURS /
#  OMP_NUM_THREADS(既定2) / CPU_LIST で上書き可)
daytime: $(TARGET)
	bash scripts/run_daytime.sh $(ARGS)


# run: OMP_NUM_THREADS=4 \ OMP_PROC_BIND=TRUE \ OMP_PLACES=cores \ OMP_DYNAMIC=FALSE \ ./$(TARGET)

# run:
# 	./$(TARGET) 

# クリーンアップ (オブジェクトファイルとバイナリ削除)
clean:
	find $(OBJDIR) -type f -name '*.o' -delete
	rm -rf $(BINDIR)/run
