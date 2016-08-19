        Max Rect

About
-----

    Determine the maximum area axis-aligned rectangle in a convex quad
    (or a triangle).


Files
-----

convex.c

    Provides a simple test to determine if a quad is convex.

        extern bool
        is_convex_quad(const Quad *quad);

    This is useful if you can not guarantee that the quad is convex.

maxrect.c

    Calculates the maximum axis-aligned rectangle for the given quad.

        extern Rect
        max_aarect_in_quad(const Quad *quad);

    Does not check whether the quad is convex.

maxrect.h

    Exposed interface.


Author
------

David Wang <millimillenary@gmail.com>
