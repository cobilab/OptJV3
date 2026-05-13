CC       := gcc
TARGET   := OptJV3
SRC      := optjv3.c
OBJ      := $(SRC:.c=.o)

STD      := -std=c11
WARN     := -Wall -Wextra
THREADS  := -pthread
OPT      := -O3 -march=native -flto=auto
LIBS     := -lm

CFLAGS   := $(STD) $(WARN) $(THREADS) $(OPT)
LDFLAGS  := $(THREADS) $(OPT)
LDLIBS   := $(LIBS)

DEBUG_CFLAGS := $(STD) $(WARN) $(THREADS) -O0 -g3
PGO_DIR      := profdir

.PHONY: all release debug clean distclean rebuild pgo-gen pgo-use pgo-clean

all: release

release: $(TARGET)

$(TARGET): $(OBJ)
	$(CC) $(LDFLAGS) -o $@ $^ $(LDLIBS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

debug: CFLAGS := $(DEBUG_CFLAGS)
debug: LDFLAGS := $(THREADS)
debug: clean $(TARGET)

rebuild: clean release

clean:
	rm -f $(TARGET) $(OBJ)

distclean: clean pgo-clean

pgo-clean:
	rm -rf $(PGO_DIR)

# Build an instrumented binary for profile generation.
pgo-gen: CFLAGS := $(STD) $(WARN) $(THREADS) -O3 -march=native -flto=auto -fprofile-generate=$(PGO_DIR)
pgo-gen: LDFLAGS := $(THREADS) -O3 -march=native -flto=auto -fprofile-generate=$(PGO_DIR)
pgo-gen: clean $(TARGET)

# After running the instrumented binary on representative workloads,
# rebuild using the collected profiles.
pgo-use: CFLAGS := $(STD) $(WARN) $(THREADS) -O3 -march=native -flto=auto -fprofile-use=$(PGO_DIR) -fprofile-correction
pgo-use: LDFLAGS := $(THREADS) -O3 -march=native -flto=auto -fprofile-use=$(PGO_DIR) -fprofile-correction
pgo-use: clean $(TARGET)
