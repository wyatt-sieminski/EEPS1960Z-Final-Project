CC = gcc

CFLAGS = -Wall -Wextra -Werror -pedantic -std=c99 -g

all: main

main: chamber.o
	$(CC) $(CFLAGS) -o chamber chamber.o

clean:
	rm -f chamber chamber.o

.PHONY: all clean