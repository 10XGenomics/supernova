// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.
//
#ifndef SYSTEM_JEMALLOC_HOOKS_INCLUDED
#define SYSTEM_JEMALLOC_HOOKS_INCLUDED
#include <jemalloc/jemalloc.h>

#ifdef JEMALLOC_HOOKS
void (*__free_hook)(void *ptr) = je_free;
void *(*__malloc_hook)(size_t size) = je_malloc;
void *(*__realloc_hook)(void *ptr, size_t size) = je_realloc;
void *(*__memalign_hook)(size_t alignment, size_t size) = je_memalign;
#endif

#endif
