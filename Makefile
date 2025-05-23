CC = gcc
Flags = -Wall -Wextra -std=C99 -pedantic-errors
Debug = -g
Lib = -lraylib -lGL -lm -lpthread -ldl -lrt -lX11

.PHONY: main clean

main:
	$(CC) main_refactored.c -o main $(CFlags) $(Lib) $(Debug)
	./main

clean:
	rm main


