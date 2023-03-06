#if defined(_WIN32)
#ifndef NOMINMAX
#define NOMINMAX
#endif
#endif

#ifdef TINYEXR_USE_MINIZ
#undef TINYEXR_USE_MINIZ
#endif
#define TINYEXR_IMPLEMENTATION
#include "tinyexr.h"
