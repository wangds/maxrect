#ifndef MAXRECT_H
#define MAXRECT_H

#include <stdbool.h>

enum {
	NUM_VERTICES_IN_QUAD = 4
};

typedef struct {
	double x, y;
} Point;

typedef struct {
	/* Point i is connected to point (i +/- 1) mod 4. */
	Point p[NUM_VERTICES_IN_QUAD];
} Quad;

typedef struct {
	double x1, y1;
	double x2, y2;
} Rect;

extern bool
is_convex_quad(const Quad *quad);

extern Rect
max_aarect_in_quad(const Quad *quad);

#endif
