/* convex.c */

#include "maxrect.h"

#include <assert.h>
#include <math.h>

/*--------------------------------------------------------------*/

static bool
line_line_intersection(
	const Point *a1, const Point *a2,
	const Point *b1, const Point *b2,
	Point *out_p)
{
	const double x1 = a1->x;
	const double y1 = a1->y;
	const double x2 = a2->x;
	const double y2 = a2->y;
	const double x3 = b1->x;
	const double y3 = b1->y;
	const double x4 = b2->x;
	const double y4 = b2->y;
	const double denom = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);

	if (fabs(denom) > 0.0) {
		out_p->x = ((x1 * y2 - y1 * x2) * (x3 - x4) - (x1 - x2) * (x3 * y4 - y3 * x4)) / denom;
		out_p->y = ((x1 * y2 - y1 * x2) * (y3 - y4) - (y1 - y2) * (x3 * y4 - y3 * x4)) / denom;
		return true;
	} else {
		out_p->x = 0.0;
		out_p->y = 0.0;
		return false;
	}
}

/*--------------------------------------------------------------*/

bool
is_convex_quad(const Quad *quad)
{
	bool intersect;
	Point p;
	assert(quad);

	if (!quad)
		return false;

	intersect = line_line_intersection(
			&quad->p[0], &quad->p[2], &quad->p[1], &quad->p[3], &p);
	if (intersect) {
		const bool convex
			= ((quad->p[0].x <= p.x && p.x <= quad->p[2].x) || (quad->p[2].x <= p.x && p.x <= quad->p[0].x))
			&& ((quad->p[0].y <= p.y && p.y <= quad->p[2].y) || (quad->p[2].y <= p.y && p.y <= quad->p[0].y))
			&& ((quad->p[1].x <= p.x && p.x <= quad->p[3].x) || (quad->p[3].x <= p.x && p.x <= quad->p[1].x))
			&& ((quad->p[1].y <= p.y && p.y <= quad->p[3].y) || (quad->p[3].y <= p.y && p.y <= quad->p[1].y));

		return convex;
	} else {
		return false;
	}
}
