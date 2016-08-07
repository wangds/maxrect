/* maxrect.c */

#include "maxrect.h"

#include <assert.h>
#include <math.h>
#include <string.h>

enum {
	NUM_EDGES_IN_QUAD = NUM_VERTICES_IN_QUAD
};

static bool
is_point_inside_quad(const Quad *quad, const Point *p);

static Rect
two_contact_case(const Quad *quad, int v);

static Rect
three_contact_case(const Quad *quad, int i, int j);

static bool
find_same_x_edge_segment(
	const Point *p, const Point *q,
	const Point *a, const Point *b,
	Point *out_a, Point *out_b);

static bool
find_same_y_edge_segment(
	const Point *p, const Point *q,
	const Point *a, const Point *b,
	Point *out_a, Point *out_b);

/*--------------------------------------------------------------*/

static double
clamp(double min, double val, double max)
{
	return (val <= min) ? min : (val >= max) ? max : val;
}

static int
next_vertex_index(int i)
{
	return (i + 1) % NUM_VERTICES_IN_QUAD;
}

static int
opposite_vertex_index(int i)
{
	return (i + 2) % NUM_VERTICES_IN_QUAD;
}

static int
prev_vertex_index(int i)
{
	return (i + 3) % NUM_VERTICES_IN_QUAD;
}

static double
lerp(double a, double b, double s)
{
	return (1.0 - s) * a + s * b;
}

static double
max(double a, double b)
{
	return (a >= b) ? a : b;
}

static double
min(double a, double b)
{
	return (a <= b) ? a : b;
}

/*--------------------------------------------------------------*/

Rect
max_aarect_in_quad(const Quad *quad)
{
	Rect best_rect;
	double best_area = 0.0;
	assert(quad);

	memset(&best_rect, 0, sizeof(best_rect));

	if (!quad)
		return best_rect;

	/* Two contact point cases. */
	for (int i = 0; i < NUM_VERTICES_IN_QUAD; i++) {
		const Rect rect = two_contact_case(quad, i);
		const double area = fabs((rect.x2 - rect.x1) * (rect.y2 - rect.y1));

		if (area > best_area) {
			best_rect = rect;
			best_area = area;
		}
	}

	/* Three contact point cases. */
	for (int i = 0; i < NUM_EDGES_IN_QUAD; i++) {
		for (int j = 0; j < NUM_EDGES_IN_QUAD; j++) {
			if (j == i)
				continue;

			const Rect rect = three_contact_case(quad, i, j);
			const double area = fabs((rect.x2 - rect.x1) * (rect.y2 - rect.y1));

			if (area > best_area) {
				best_rect = rect;
				best_area = area;
			}
		}
	}

	return best_rect;
}

static bool
is_point_inside_quad(const Quad *quad, const Point *p)
{
	const double epsilon = 1.0e-6;
	int num_above = 0;
	int num_below = 0;

	for (int i = 0; i < NUM_VERTICES_IN_QUAD; i++) {
		const Point * const a = &(quad->p[i]);
		const Point * const b = &(quad->p[next_vertex_index(i)]);

		if ((a->x <= p->x && p->x <= b->x)
				|| (b->x <= p->x && p->x <= a->x)) {
			if (a->x == b->x) {
				if ((a->y <= p->y && p->y <= b->y)
						|| (b->y <= p->y && p->y <= a->y)) {
					num_above++;
					num_below++;
				}
			} else {
				const double y = lerp(a->y, b->y, (p->x - a->x) / (b->x - a->x));

				if (y >= p->y - epsilon)
					num_above++;
				if (y <= p->y + epsilon)
					num_below++;
			}
		}
	}

	/* Inside (or on the boundary) if lines cross above and below. */
	return (num_above > 0) && (num_below > 0);
}

static Rect
two_contact_case(const Quad *quad, int contact)
{
	Rect best_rect;
	double best_area = 0.0;

	const int o = opposite_vertex_index(contact);
	const Point * const v = &(quad->p[contact]);
	const Point * const p = &(quad->p[o]);

	memset(&best_rect, 0, sizeof(best_rect));

	for (int i = 0; i < 2; i++) {
		const int qidx = (i == 0) ? prev_vertex_index(o) : next_vertex_index(o);
		const Point * const q = &(quad->p[qidx]);

		const double a = (q->x - p->x) * (q->y - p->y);
		const double b = (q->x - p->x) * (p->y - v->y) + (p->x - v->x) * (q->y - p->y);
		const double c = (p->x - v->x) * (p->y - v->y);
		(void)c;

		double s;

		if (fabs(a) > 0.0) {
			s = -b / (2.0 * a);
		} else if (fabs(q->x - p->x) <= 0.0 && fabs(q->y - p->y) > 0.0) {
			/* Vertical line.  */
			const int j = (i == 0) ? next_vertex_index(o) : prev_vertex_index(o);
			const double yprime = quad->p[j].y;
			s = (yprime - p->y) / (q->y - p->y);
		} else if (fabs(q->y - p->y) <= 0.0 && fabs(q->x - p->x) > 0.0) {
			/* Horizontal line.  */
			const int j = (i == 0) ? next_vertex_index(o) : prev_vertex_index(o);
			const double xprime = quad->p[j].x;
			s = (xprime - p->x) / (q->x - p->x);
		} else {
			continue;
		}

		s = clamp(0.0, s, 1.0);

		const double xprime = lerp(p->x, q->x, s);
		const double yprime = lerp(p->y, q->y, s);
		const double area = fabs((xprime - v->x) * (yprime - v->y));
		const Point test1 = { xprime, v->y };
		const Point test2 = { v->x, yprime };

		if (area > best_area
				&& is_point_inside_quad(quad, &test1)
				&& is_point_inside_quad(quad, &test2)) {
			best_rect.x1 = min(v->x, xprime);
			best_rect.y1 = min(v->y, yprime);
			best_rect.x2 = max(v->x, xprime);
			best_rect.y2 = max(v->y, yprime);
			best_area = area;
		}
	}

	return best_rect;
}

static Rect
three_contact_case(const Quad *quad, int i, int j)
{
	Rect best_rect;
	double best_area = 0.0;

	const Point * const p = &(quad->p[i]);
	const Point * const q = &(quad->p[next_vertex_index(i)]);
	Point xprime1, xprime2;
	double s_lower_x = 0.0;
	double s_upper_x = 1.0;

	memset(&best_rect, 0, sizeof(best_rect));

	/* Find equation for x'(s). */
	if (!find_same_y_edge_segment(
			p, q, &(quad->p[j]), &(quad->p[next_vertex_index(j)]),
			&xprime1, &xprime2)) {
		return best_rect;
	}

	if (fabs(q->y - p->y) > 0.0) {
		const Point * const a = &(quad->p[j]);
		const Point * const b = &(quad->p[next_vertex_index(j)]);
		const double s1 = (a->y - p->y) / (q->y - p->y);
		const double s2 = (b->y - p->y) / (q->y - p->y);

		s_lower_x = max(s_lower_x, min(s1, s2));
		s_upper_x = min(s_upper_x, max(s1, s2));
	}

	/* Take 100 sample points along the edge.
	 * For each x(s) and x'(s), find a third contact point at y'(s).
	 */
	for (int iter = 0; iter <= 100; iter++) {
		const double s = lerp(s_lower_x, s_upper_x, (double)iter / 100.0);

		const double x = lerp(p->x, q->x, s);
		const double y = lerp(p->y, q->y, s);
		const double xprime = lerp(xprime1.x, xprime2.x, s);

		for (int k = 0; k < NUM_EDGES_IN_QUAD; k++) {
			for (int yy = 0; yy < 2; yy++) {
				Point yprime1, yprime2;

				if (yy == 0) {
					if (!find_same_x_edge_segment(
							p, q, &(quad->p[k]), &(quad->p[next_vertex_index(k)]),
							&yprime1, &yprime2)) {
						continue;
					}
				} else {
					if (!find_same_x_edge_segment(
							&xprime1, &xprime2, &(quad->p[k]), &(quad->p[next_vertex_index(k)]),
							&yprime1, &yprime2)) {
						continue;
					}
				}

				const double yprime = lerp(yprime1.y, yprime2.y, s);

				const double area = fabs((xprime - x) * (yprime - y));
				const Point test1 = { x, yprime };
				const Point test2 = { xprime, yprime };

				if (area > best_area
						&& is_point_inside_quad(quad, &test1)
						&& is_point_inside_quad(quad, &test2)) {
					best_rect.x1 = min(x, xprime);
					best_rect.y1 = min(y, yprime);
					best_rect.x2 = max(x, xprime);
					best_rect.y2 = max(y, yprime);
					best_area = area;
				}
			}
		}
	}

	return best_rect;
}

static bool
find_same_x_edge_segment(
	const Point *p, const Point *q,
	const Point *a, const Point *b,
	Point *out_a, Point *out_b)
{
	if (fabs(b->x - a->x) > 0.0) {
		const double c = (p->x - a->x) / (b->x - a->x);
		const double m = (q->x - p->x) / (b->x - a->x);

		out_a->x = p->x;
		out_a->y = a->y + (b->y - a->y) * c;
		out_b->x = q->x;
		out_b->y = (b->y - a->y) * m + out_a->y;

		/* Check: out_a->x == p->x, out_b->x == q->x. */
		/*
		const double epsilon = 0.1;
		out_a->x = a->x + (b->x - a->x) * c;
		out_b->x = (b->x - a->x) * m + out_a->x;
		assert(fabs(p->x - out_a->x) < epsilon);
		assert(fabs(q->x - out_b->x) < epsilon);
		*/
		return true;
	} else {
		out_a->x = 0.0;
		out_a->y = 0.0;
		out_b->x = 0.0;
		out_b->y = 0.0;
		return false;
	}
}

static bool
find_same_y_edge_segment(
	const Point *p, const Point *q,
	const Point *a, const Point *b,
	Point *out_a, Point *out_b)
{
	if (fabs(b->y - a->y) > 0.0) {
		const double c = (p->y - a->y) / (b->y - a->y);
		const double m = (q->y - p->y) / (b->y - a->y);

		out_a->x = a->x + (b->x - a->x) * c;
		out_a->y = p->y;
		out_b->x = (b->x - a->x) * m + out_a->x;
		out_b->y = q->y;

		/* Check: out_a->y == p->y, out_b->y == q->y. */
		/*
		const double epsilon = 0.1;
		out_a->y = a->y + (b->y - a->y) * c;
		out_b->y = (b->y - a->y) * m + out_a->y;
		assert(fabs(p->y - out_a->y) < epsilon);
		assert(fabs(q->y - out_b->y) < epsilon);
		*/
		return true;
	} else {
		out_a->x = 0.0;
		out_a->y = 0.0;
		out_b->x = 0.0;
		out_b->y = 0.0;
		return false;
	}
}
