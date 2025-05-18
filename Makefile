CC = gcc
Flags = -Wall -Wextra -std=C99 -pedantic-errors
Lib = -lraylib -lGL -lm -lpthread -ldl -lrt -lX11

.PHONY: main clean

main:
	$(CC) main.c -o main $(CFlags) $(Lib)	
	./main

clean:
	rm main


