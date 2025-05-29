CC = gcc
CFlags = -Wall -Wextra -std=c99 -pedantic-errors 
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
DFlags += -g -O0 # -g enables DWARF symbols for debugging, -O0 is no optimization 
else 
CFlags += -O2
endif

$(Build)/%.o : $(Src)/%.c
	$(CC) $(CFlags) $(DFlags) -c $< -o $@ 


$(Build)/main: $(OBJS)
	$(CC) $^ -o $@ $(Libs)

clean:
	rm $(Build)/* 

