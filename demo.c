/* demo.c */

#include <SDL2/SDL.h>
#include <assert.h>
#include <stdbool.h>
#include "maxrect.h"

enum {
	SCREEN_W = 640,
	SCREEN_H = 480,
};

/*--------------------------------------------------------------*/

static Quad
init_quad(void)
{
	Quad quad;

	/* Parallelogram. */
	quad.p[0] = (Point){ 100, 100 };
	quad.p[1] = (Point){ 300,  90 };
	quad.p[2] = (Point){ 310, 190 };
	quad.p[3] = (Point){ 110, 200 };

	/* Parallelogram, vertical lines. */
	quad.p[0] = (Point){ 200, 158 };
	quad.p[1] = (Point){ 518,  41 };
	quad.p[2] = (Point){ 518, 248 };
	quad.p[3] = (Point){ 200, 343 };

	/* Trapezoid, vertical lines. */
	quad.p[0] = (Point){ 200, 158 };
	quad.p[1] = (Point){ 518,  41 };
	quad.p[2] = (Point){ 518, 408 };
	quad.p[3] = (Point){ 200, 343 };

	return quad;
}

static int
find_closest_vertex(const Quad *quad, int x, int y)
{
	const int threshold = 5;
	int best_dist = threshold * threshold;
	int best_idx = -1;

	for (int i = 0; i < NUM_VERTICES_IN_QUAD; i++) {
		const int dx = (int)quad->p[i].x - x;
		const int dy = (int)quad->p[i].y - y;

		if (dx * dx + dy * dy < best_dist) {
			best_idx = i;
			best_dist = dx * dx + dy * dy;
		}
	}

	return best_idx;
}

int main(void)
{
	bool quit = false;
	SDL_Window *sdl_window;
	SDL_Renderer *sdl_renderer;
	Quad quad;
	int hover = -1;
	int hot = -1;
	int err;

	err = SDL_CreateWindowAndRenderer(SCREEN_W, SCREEN_H, 0,
			&sdl_window, &sdl_renderer);
	if (err != 0)
		return err;

	SDL_SetWindowTitle(sdl_window, "Max Rect");
	quad = init_quad();

	while (!quit) {
		SDL_Event event;
		int timeout = 1000 / 60;
		int read_event = SDL_WaitEventTimeout(&event, timeout);
		if (read_event) {
			switch (event.type) {
				case SDL_QUIT:
					quit = true;
					break;

				case SDL_MOUSEMOTION:
					if (hot != -1) {
						quad.p[hot].x = event.motion.x;
						quad.p[hot].y = event.motion.y;
						hover = hot;
					} else {
						hover = find_closest_vertex(&quad, event.motion.x, event.motion.y);
					}
					break;

				case SDL_MOUSEBUTTONDOWN:
					if (event.button.button == 1)
						hot = find_closest_vertex(&quad, event.button.x, event.button.y);
					break;

				case SDL_MOUSEBUTTONUP:
					if (event.button.button == 1)
						hot = -1;
					break;

				default:
					break;
			}
		}

		{
			const bool convex = is_convex_quad(&quad);

			SDL_SetRenderDrawColor(sdl_renderer, 0x00, 0x00, 0x00, 0xFF);
			SDL_RenderClear(sdl_renderer);

			/* Draw edges. */
			if (convex) {
				SDL_SetRenderDrawColor(sdl_renderer, 0xFF, 0xFF, 0xFF, 0xFF);
			} else {
				SDL_SetRenderDrawColor(sdl_renderer, 0xFF, 0x00, 0x00, 0xFF);
			}

			for (int i = 0; i < NUM_VERTICES_IN_QUAD; i++) {
				const int j = (i + 1) % NUM_VERTICES_IN_QUAD;
				SDL_RenderDrawLine(sdl_renderer,
						(int)quad.p[i].x, (int)quad.p[i].y,
						(int)quad.p[j].x, (int)quad.p[j].y);
			}

			/* Draw hot or hover. */
			if (hot != -1) {
				SDL_Rect rect;
				rect.x = (int)quad.p[hot].x - 5;
				rect.y = (int)quad.p[hot].y - 5;
				rect.w = 2*5 + 1;
				rect.h = 2*5 + 1;

				SDL_SetRenderDrawColor(sdl_renderer, 0xFF, 0x00, 0x00, 0xFF);
				SDL_RenderDrawRect(sdl_renderer, &rect);
			} else if (hover != -1) {
				SDL_Rect rect;
				rect.x = (int)quad.p[hover].x - 5;
				rect.y = (int)quad.p[hover].y - 5;
				rect.w = 2*5 + 1;
				rect.h = 2*5 + 1;

				SDL_SetRenderDrawColor(sdl_renderer, 0x00, 0xFF, 0x00, 0xFF);
				SDL_RenderDrawRect(sdl_renderer, &rect);
			}

			/* Draw max rect. */
			if (convex) {
				const Rect r = max_aarect_in_quad(&quad);
				const int minx = (int)floor((r.x1 <= r.x2) ? r.x1 : r.x2);
				const int maxx = (int) ceil((r.x1 >= r.x2) ? r.x1 : r.x2);
				const int miny = (int)floor((r.y1 <= r.y2) ? r.y1 : r.y2);
				const int maxy = (int) ceil((r.y1 >= r.y2) ? r.y1 : r.y2);
				SDL_Rect rect;

				rect.x = minx;
				rect.y = miny;
				rect.w = (maxx - minx);
				rect.h = (maxy - miny);

				SDL_SetRenderDrawColor(sdl_renderer, 0xA0, 0xA0, 0xFF, 0xFF);
				SDL_RenderFillRect(sdl_renderer, &rect);
			}

			SDL_RenderPresent(sdl_renderer);
		}
	}

	return 0;
}
