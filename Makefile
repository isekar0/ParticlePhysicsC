CC = clang
CFlags = -std=c11 # -Wall -Wextra -pedantic-errors 
Debug = 0
Libs  = -lraylib -lGL -lm -lpthread -ldl -lrt -lX11
Src = src
Build = build

.PHONY: all clean

SRCS := $(wildcard $(Src)/*.c)
OBJS := $(patsubst $(Src)/%.c,$(Build)/%.o,$(SRCS))

all: $(Build)/main

$(Build):
	mkdir -p $@

DFlags = 
ifeq ($(Debug), 1)
DFlags += -O0 -ggdb3 -fno-omit-frame-pointer -fsanitize=address # -g enables DWARF symbols for debugging, -O0 is no optimization 
else 
CFlags += -O3 -mavx2 -msse4.2 -march=core-avx2 # Optimization and SIMD instructions
endif

$(Build)/%.o : $(Src)/%.c
	$(CC) $(CFlags) $(DFlags) -c $< -o $@ 


$(Build)/main: $(OBJS)
	$(CC) $^ -o $@ $(Libs)

clean:
	rm $(Build)/* 

