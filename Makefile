CC := gcc
CFLAGS := -std=c99 -W -Wall -Wextra -Wpedantic
CFLAGS += -Wmissing-prototypes -Wstrict-prototypes -Wold-style-definition
CFLAGS += -Wcast-align -Wcast-qual -Wconversion -Winit-self -Wshadow
CFLAGS += -Wsign-conversion -Wswitch-default -Wundef
CFLAGS += -O2
CFLAGS += -g
LDLIBS += -lm

all: convex.o maxrect.o

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

.PHONY: clean
clean:
	$(RM) convex.o maxrect.o
