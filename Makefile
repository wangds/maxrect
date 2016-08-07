CC := gcc
CFLAGS := -std=c99 -W -Wall -Wextra -Wpedantic
CFLAGS += -Wmissing-prototypes -Wstrict-prototypes -Wold-style-definition
CFLAGS += -Wcast-align -Wcast-qual -Wconversion -Winit-self -Wshadow
CFLAGS += -Wsign-conversion -Wswitch-default -Wundef
CFLAGS += -O2
CFLAGS += -g
LDLIBS += -lm

all: demo convex.o maxrect.o

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

demo: demo.o convex.o maxrect.o
	$(CC) $(CFLAGS) -o $@ $^ -lSDL2 $(LDLIBS)

.PHONY: clean
clean:
	$(RM) demo demo.o
	$(RM) convex.o maxrect.o
