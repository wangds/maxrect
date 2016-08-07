/* maxrect.c */

#include "maxrect.h"

#include <assert.h>
#include <string.h>

/*--------------------------------------------------------------*/

Rect
max_aarect_in_quad(const Quad *quad)
{
	Rect best_rect;
	assert(quad);

	memset(&best_rect, 0, sizeof(best_rect));

	if (!quad)
		return best_rect;

	return best_rect;
}
