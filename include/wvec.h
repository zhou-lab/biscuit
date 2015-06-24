#ifndef _VECTOR_H_
#define _VECTOR_H_

/*
 *  Dynamic Array (or "vector", by C++ convention)
 *  
 *  adapted from Ruan Jue's list.h, Li Heng's kvec.h
 *  
 *  Example usage:
 *  
 *  DEFINE_VECTOR(MyVector, MyElem)
 *  
 *  1. dynamic
 *  MyVector *v = init_MyVector(1023);
 *  push_MyVector(v, elem);
 *  free_MyVector(v)
 *  
 *  2. static
 *  MyVector v;
 *  MyVector_init(&v, 1023);
 *  push_MyVector(&v, elem);
 *  MyVector_free(v);
 *  
 */


#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdio.h>

#define DEFINE_VECTOR_CORE(Vector, Element, Size, inc_size)							\
																																				\
  typedef struct {																											\
    Element* buffer;																										\
    Size size; /* size: logic length*/																	\
    Size cap;  /* capacity: allocated length */													\
  } Vector;																															\
																																				\
  /* MyVector *v = init_MyVector(1023); */															\
  static inline Vector*																									\
  init_##Vector(Size init_size){																				\
    if (init_size == 0)	init_size = 2;																	\
    Vector *vector = (Vector*) malloc(sizeof(Vector));									\
    vector->size = 0;																										\
    vector->cap  = init_size;																						\
    vector->buffer = (Element*) malloc(sizeof(Element) * vector->cap);	\
																																				\
    return vector;																											\
  }																																			\
																																				\
  /* MyVector v; MyVector_init(&v, 1023); */														\
  static inline void																										\
  Vector##_init(Vector *vector,																					\
								Size init_size){																				\
    if (init_size == 0) init_size = 2;																	\
    vector->size = 0;																										\
    vector->cap  = init_size;																						\
    vector->buffer = (Element*) malloc(sizeof(Element) * vector->cap);	\
  }																																			\
																																				\
  /* MyVector v; count_MyVector(&v); */																	\
  static inline Size count_##Vector(Vector *vector){										\
    return vector->size;																								\
  }																																			\
																																				\
  /* MyVector v; clear_MyVector(&v); */																	\
  static inline void clear_##Vector(Vector *vector){										\
    vector->size = 0;																										\
  }																																			\
																																				\
  /* ensure capacity of n more Size */																	\
  /* when capacity < inc_size (1M) increment by doubling */							\
  /* when capacity >= inc_size (1M) increment by inc_size */						\
  static inline void																										\
  encap_##Vector(Vector *vector, Size n){																\
																																				\
    if (vector->size + n <= vector->cap) return;												\
																																				\
    while (vector->size + n > vector->cap) {														\
      if(vector->cap < inc_size){																				\
				vector->cap <<= 1;																							\
      } else {																													\
				vector->cap += inc_size;																				\
      }																																	\
    }																																		\
    vector->buffer = realloc(vector->buffer,														\
														 vector->cap * sizeof(Element));						\
  }																																			\
																																				\
  /* reduce the logic length of the vector by size */										\
  static inline void																										\
  trunc_##Vector(Vector *vector, Size size) {														\
    if(size > count_##Vector(vector)) {																	\
      size = count_##Vector(vector);																		\
    }																																		\
    vector->size -= size;																								\
  }																																			\
																																				\
  static inline void																										\
  set_##Vector##_size(Vector *vector, Size size) {											\
    vector->size = size;																								\
  }																																			\
																																				\
  static inline void																										\
  incre_##Vector(Vector *vector, Size size) {														\
    if(size + vector->size > vector->cap) {															\
      vector->size = vector->cap;																				\
    } else {																														\
      vector->size += size;																							\
    }																																		\
  }																																			\
																																				\
  /* MyVector v; push_MyVector(&v, elem) */															\
  static inline void																										\
  push_##Vector(Vector *vector, Element e) {														\
    encap_##Vector(vector, 1);																					\
    vector->buffer[vector->size++] = e;																	\
  }																																			\
																																				\
  static inline int																											\
  pop_##Vector(Vector *vector, Element*e) {															\
    if(count_##Vector(vector)) {																				\
      vector->size --;																									\
      *e = vector->buffer[vector->size];																\
      return 1;																													\
    } else {																														\
      return 0;																													\
    }																																		\
  }																																			\
																																				\
  static inline Vector*																									\
  dup_##Vector(Vector *vec) {																						\
    Vector *vec2 = calloc(1, sizeof(Vector));														\
    vec2->buffer = malloc(sizeof(Element) * (vec->size+1));							\
    vec2->cap = vec->size+1;																						\
    vec2->size = vec->size;																							\
    memcpy(vec2->buffer, vec->buffer, sizeof(Element)*vec->size);				\
    return vec2;																												\
  }																																			\
																																				\
  /* insert at any location (given by idx) */														\
  /* may be inefficient, be careful */																	\
  static inline void																										\
  insert_##Vector(Vector *vector, Size idx, Element e) {								\
																																				\
    if(idx > vector->size)																							\
      return;																														\
																																				\
    encap_##Vector(vector, 1);																					\
    if(idx == vector->size){																						\
      vector->buffer[vector->size] = e;																	\
    } else {																														\
      memmove(vector->buffer + idx + 1,																	\
							vector->buffer + idx,																			\
							(vector->size - idx) * sizeof(Element));									\
      vector->buffer[idx] = e;																					\
    }																																		\
    vector->size ++;																										\
  }																																			\
																																				\
  static inline Element*																								\
  first_ref_##Vector(Vector *vec) {																			\
    if (vec->size > 0) return vec->buffer;															\
    else return NULL;																										\
  }																																			\
																																				\
  static inline Element*																								\
  last_ref_##Vector(Vector *vec) {																			\
    if (vec->size > 0) return vec->buffer+vec->size-1;									\
    else return NULL;																										\
  }																																			\
																																				\
  /* no check of size, use with caution */															\
  static inline Element																									\
  first_##Vector(Vector *vec) {																					\
    return vec->buffer[0];																							\
  }																																			\
																																				\
  /* no check of size, use with caution */															\
  static inline Element																									\
  last_##Vector(Vector *vec) {																					\
    return vec->buffer[vec->size-1];																		\
  }																																			\
																																				\
  static inline void																										\
  remove_##Vector(Vector *vector, Size idx) {														\
    if(idx >= vector->size) return;																			\
    if(idx + 1 < vector->size){																					\
      memmove(vector->buffer + idx,																			\
							vector->buffer + idx + 1,																	\
							(vector->size - idx - 1) * sizeof(Element));							\
    }																																		\
    vector->size --;																										\
  }																																			\
																																				\
  static inline void																										\
  set_##Vector(Vector *vector, Size idx, Element e) {										\
    vector->buffer[idx] = e;																						\
  }																																			\
																																				\
  static inline Element																									\
  get_##Vector(Vector *vector, Size idx) {															\
    return vector->buffer[idx];																					\
  }																																			\
																																				\
  static inline Element*																								\
  ref_##Vector(Vector *vector, Size idx) {															\
    return vector->buffer + idx;																				\
  }																																			\
																																				\
  /* push an empty element, return its pointer */												\
  static inline Element*																								\
  next_ref_##Vector(Vector *vector) {																		\
    encap_##Vector(vector, 1);																					\
    vector->size++;																											\
    return vector->buffer + vector->size - 1;														\
  }																																			\
																																				\
  static inline Element*																								\
  try_next_##Vector(Vector *vector) {																		\
    encap_##Vector(vector, 1);																					\
    return vector->buffer + vector->size;																\
  }																																			\
																																				\
  static inline void																										\
  commit_next_##Vector(Vector *vector) {																\
    vector->size++;																											\
  }																																			\
																																				\
  /* no encap version of next_ref, memory dangerous! */									\
  static inline Element*																								\
  ref_next_##Vector(Vector *vector) {																		\
    vector->size++;																											\
    return vector->buffer + vector->size - 1;														\
  }																																			\
																																				\
  static inline Element*																								\
  as_array_##Vector(Vector *vector) {																		\
    return vector->buffer;																							\
  }																																			\
																																				\
  static inline void reverse_##Vector(Vector *vector){									\
    Size i, j;																													\
    Element t;																													\
    if (count_##Vector(vector) == 0) return;														\
    i = 0;																															\
    j = count_##Vector(vector) - 1;																			\
    while(i < j){																												\
      t = get_##Vector(vector, i);																			\
      set_##Vector(vector, i, get_##Vector(vector, j));									\
      set_##Vector(vector, j, t);																				\
      i++;																															\
      j--;																															\
    }																																		\
  }																																			\
																																				\
  /* extend a vector by another vector */																\
  static inline void																										\
  extend_##Vector(Vector *vector1, Vector *vector2){										\
    encap_##Vector(vector1, count_##Vector(vector2));										\
		memcpy(vector1->buffer + vector1->size,															\
					 vector2->buffer,																							\
					 sizeof(Element) * vector2->size);														\
		vector1->size += vector2->size;																			\
  }																																			\
																																				\
  /* dump vector to a file */																						\
  static inline Size																										\
  dump_##Vector(Vector *vector, FILE *out){															\
    return fwrite(vector->buffer,																				\
									sizeof(Element),																			\
									count_##Vector(vector),																\
									out);																									\
  }																																			\
																																				\
  /* dynamic free */																										\
  /* free vector buffer also free the vector pointer */									\
  static inline void																										\
  free_##Vector(Vector *vector){																				\
    free(vector->buffer);																								\
    free(vector);																												\
  }																																			\
																																				\
  /* static free */																											\
  /* only free the vector buffer */																			\
  static inline void																										\
  Vector##_free(Vector *vector){																				\
    free(vector->buffer);																								\
    vector->buffer = NULL;																							\
  }																																			\


/* vector extension functions including
 * search, exists, and many others */

#define DEFINE_VECTOR_EXT(Vector, Element, Size, equals)			\
  static inline Size																					\
  delete_##Vector(Vector *vector, Element e){									\
    Size i, ret;																							\
    ret = 0;																									\
    for(i=vector->size;i>0;i--){															\
      if(equals(vector->buffer[i-1], e)){											\
				if(i < vector->size){																	\
					memmove(vector->buffer + i - 1,											\
									vector->buffer + i,													\
									(vector->size - i) * sizeof(Element));			\
				}																											\
				vector->size --;																			\
				ret ++;																								\
      }																												\
    }																													\
    return ret;																								\
  }																														\
																															\
  static inline Size																					\
  occ_##Vector(Vector *vector, Element e){										\
    Size i, n;																								\
    for(i=0,n=0;i<vector->size;i++){													\
      if(equals(vector->buffer[i], e)) n++;										\
    }																													\
    return n;																									\
  }																														\
																															\
  static inline int																						\
  exists_##Vector(Vector *vector, Element e){									\
    Size i;																										\
    for(i=0;i<vector->size;i++){															\
      if(equals(vector->buffer[i], e)) return 1;							\
    }																													\
    return 0;																									\
  }																														\
																															\
  static inline Size																					\
  replace_##Vector(Vector *vector, Element from, Element to){	\
    Size i, ret;																							\
    ret = 0;																									\
    for(i=0;i<vector->size;i++){															\
      if(equals(vector->buffer[i], from)){										\
				vector->buffer[i] = to;																\
				ret ++;																								\
      }																												\
    }																													\
    return ret;																								\
  }																														\
																															\
  static inline Size																					\
  locate_##Vector(Vector *vector, Element e, Size start){			\
    Size i;																										\
    for(i=start;i<vector->size;i++){													\
      if(equals(vector->buffer[i], e)) return i;							\
    }																													\
    return i;																									\
  }

/* the capacity threshold to switch
 * from doubling to linear incrementation
 * is 1M (0xFFFFFU) */
#define DEFINE_VECTOR(Vector, Element)									\
  DEFINE_VECTOR_CORE(Vector, Element, size_t, 0xFFFFFU)

#define NATIVE_NUMBER_EQUALS(e1, e2) ((e1) == (e2))

#define DEFINE_NATIVE_VECTOR(Vector, Element)												\
  DEFINE_VECTOR_CORE(Vector, Element, size_t, 0xFFFFFU);						\
  DEFINE_VECTOR_EXT(Vector, Element, size_t, NATIVE_NUMBER_EQUALS);

DEFINE_NATIVE_VECTOR(u8vector,  uint8_t);
DEFINE_NATIVE_VECTOR(u16vector, uint16_t);
DEFINE_NATIVE_VECTOR(u32vector, uint32_t);
DEFINE_NATIVE_VECTOR(u64vector, uint64_t);

DEFINE_NATIVE_VECTOR(b8vector,  int8_t);
DEFINE_NATIVE_VECTOR(b16vector, int16_t);
DEFINE_NATIVE_VECTOR(b32vector, int32_t);
DEFINE_NATIVE_VECTOR(b64vector, int64_t);

DEFINE_VECTOR(vpvector, void*);

#endif
