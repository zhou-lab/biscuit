#ifndef _WQUEUE_H_
#define _WQUEUE_H_

#include <pthread.h>
#include <stdint.h>
#include <stdlib.h>

/* 
   by zhouwanding@gmail.com

   block get when empty
   block put when full

   example usage:

*/

#define DEFINE_WQUEUE(name, _wqueue_elem_t)															\
  typedef struct {																											\
    uint32_t size;																											\
    _wqueue_elem_t *data;																								\
    uint32_t head, tail;																								\
    int full, empty;																										\
    pthread_mutex_t mut;																								\
    pthread_cond_t not_full, not_empty;																	\
  } wqueue_##name##_t;																									\
  static inline wqueue_##name##_t *wqueue_init_##name(uint32_t size) {	\
    wqueue_##name##_t *q = calloc(1, sizeof(wqueue_##name##_t));				\
    q->empty = 1;																												\
    q->size = size;																											\
    q->data = calloc(size, sizeof(_wqueue_elem_t));											\
    pthread_mutex_init(&q->mut, NULL);																	\
    pthread_cond_init(&q->not_full, NULL);															\
    pthread_cond_init(&q->not_empty, NULL);															\
    return q;																														\
  }																																			\
  static inline void wqueue_destroy_##name(wqueue_##name##_t *q) {			\
    free(q->data);																											\
    pthread_mutex_destroy(&q->mut);																			\
    pthread_cond_destroy(&q->not_full);																	\
    pthread_cond_destroy(&q->not_empty);																\
    free(q);																														\
  }																																			\
  static inline void wqueue_get_##name(wqueue_##name##_t *q,						\
																			 _wqueue_elem_t *e) {							\
    pthread_mutex_lock(&q->mut);																				\
    while(q->empty) pthread_cond_wait(&q->not_empty, &q->mut);					\
    *e = q->data[q->head++];																						\
    if (q->head == q->size) q->head = 0;																\
    if (q->head == q->tail) q->empty = 1;																\
    q->full = 0;																												\
    pthread_mutex_unlock(&q->mut);																			\
    pthread_cond_signal(&q->not_full);																	\
  }																																			\
  static inline void wqueue_put_##name(wqueue_##name##_t *q,						\
																			 _wqueue_elem_t *e) {							\
    pthread_mutex_lock(&q->mut);																				\
    while(q->full) pthread_cond_wait(&q->not_full, &q->mut);						\
    q->data[q->tail++] = *e;																						\
    if (q->tail == q->size) q->tail = 0;																\
    if (q->head == q->tail) q->full = 1;																\
    q->empty = 0;																												\
    pthread_mutex_unlock(&q->mut);																			\
    pthread_cond_signal(&q->not_empty);																	\
  }																																			\
  static inline void wqueue_put2_##name(wqueue_##name##_t *q,						\
																				_wqueue_elem_t e) {							\
    pthread_mutex_lock(&q->mut);																				\
    while(q->full) pthread_cond_wait(&q->not_full, &q->mut);						\
    q->data[q->tail++] = e;																							\
    if (q->tail == q->size) q->tail = 0;																\
    if (q->head == q->tail) q->full = 1;																\
    q->empty = 0;																												\
    pthread_mutex_unlock(&q->mut);																			\
    pthread_cond_signal(&q->not_empty);																	\
  }

#define wqueue_t(name) wqueue_##name##_t

#define wqueue_init(name, size) wqueue_init_##name(size)

#define wqueue_destroy(name, q) wqueue_destroy_##name(q)

#define wqueue_get(name, q, e) wqueue_get_##name(q, e)

#define wqueue_put(name, q, e) wqueue_put_##name(q, e)

#define wqueue_put2(name, q, e) wqueue_put2_##name(q, e)

#endif
