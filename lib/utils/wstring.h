/*
 * wstring.h
 * adapted from Ruan Juan's string.h
 * different from the original version
 * the capacity is not necessarily a power of 2
 * removed much irrelevant code
 *
 */
 
#ifndef _WSTRING_H
#define _WSTRING_H

#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <stdlib.h>

#ifndef SWAP_TMP
#define SWAP_TMP
#define swap_tmp(a, b, t) { (t) = (a); (a) = (b); (b) = (t); }
#endif

typedef struct {
  char *string;
  size_t capacity;
} String;

#define uc(ch) (((ch) >= 'a' && (ch) <= 'z')? (ch) + 'A' - 'a' : (ch))
#define lc(ch) (((ch) >= 'A' && (ch) <= 'Z')? (ch) + 'a' - 'A' : (ch))

static inline String* init_string(int cap){
  String *str;
  str = (String*)malloc(sizeof(String));
  str->capacity = cap>1 ? cap : 2;
  str->string = (char*) malloc(sizeof(char) * (str->capacity));
  str->string[0] = 0;
  return str;
}

static inline void encap_string(String *str, int inc){

  if (strlen(str->string) + 1 + inc > str->capacity) {
    while (strlen(str->string) + 1 + inc > str->capacity) {
      if (str->capacity < 0xFFFFF) {
	str->capacity <<= 1;
      } else {
	str->capacity += 0xFFFFF;
      }
    }
    str->string = (char*) realloc(str->string, str->capacity);
  }
}

static inline void uc_string(String *str){
  size_t i;
  for(i=0;i<strlen(str->string);i++){
    if(str->string[i] >= 'a' && str->string[i] <= 'z') str->string[i] = str->string[i] + 'A' - 'a';
  }
}

static inline void lc_string(String *str){
  size_t i;
  for(i=0;i<strlen(str->string);i++){
    if(str->string[i] >= 'A' && str->string[i] <= 'Z') str->string[i] = str->string[i] + 'a' - 'A';
  }
}

static inline char* substr(char *string, size_t start, size_t end, char *dst){
  size_t i, size;
  char *str;
  size = strlen(string);
  if(start > size) start = size;
  if(end > size) end = size;
  if(end < start) size=0;
  else size = end - start;
  if(dst != NULL) str = dst;
  else str = (char*)malloc(sizeof(char) * (size + 1));
  for(i=start;i<end;i++){
    str[i-start] = string[i];
  }
  str[size] = '\0';
  return str;
}

static inline char* catstr(int n_str, ...){
  char *str, *s;
  int i, len;
  va_list params;
	
  len = 0;
  str = NULL;
  va_start(params, n_str);
  for(i=0;i<n_str;i++){
    s = va_arg(params, char*);
    len += strlen(s);
    str = realloc(str, len + 1);
    if(i == 0) str[0] = 0;
    strcat(str, s);
  }
  va_end(params);
  return str;
}

static inline void chomp_string(String *str){
  if(str->string[0] && str->string[strlen(str->string)] == '\n'){
    str->string[strlen(str->string)] = 0;
  }
}


/* static inline void trim_string(String *str){ */
/*   int i, j; */
/*   i = strlen(str->string) - 1; */
/*   while(i >= 0 && (str->string[i] == '\n' || str->string[i] == '\t' || str->string[i] == ' ')) i--;  */
/*   i = 0; */
/*   while(i < strlen(str->string) && (str->string[i] == '\n' || str->string[i] == '\t' || str->string[i] == ' ')) i++; */
/*   if(i){			/\* left shift *\/ */
/*     for(j=i;j<str->size;j++){ str->string[j-i] = str->string[j]; } */
/*     str->size -= i; */
/*   } */
/*   str->string[str->size] = 0; */
/* } */

static inline void
append_string(String *str, const char *src){
  encap_string(str, strlen(src));
  strcat(str->string, src);
}

/* static inline void append_char_string(String *str, char c, int num){ */
/*   encap_string(str, num); */
/*   while(num-- > 0){ str->string[str->size ++] = c; } */
/*   str->string[str->size] = 0; */
/* } */

static inline String* as_string(const char *c_str){
  String *str = init_string(strlen(c_str)+1);
  append_string(str, c_str);
  return str;
}

static inline void putchar_string(String *str, char ch){
  encap_string(str, 1);
  int len = strlen(str->string);
  str->string[len] = ch;
  str->string[len+1] = 0;
}

static inline void
clear_string(String *str){
  str->string[0] = 0;
}

static inline void reverse_string(String *str){
  int i, j;
  char c;
  i = 0;
  j = strlen(str->string) - 1;
  while(i < j){
    swap_tmp(str->string[i], str->string[j], c);
    i++;
    j--;
  }
}

static inline void reverse_str(char *str, int len){
  int i, j;
  char c;
  i = 0;
  j = len - 1;
  while(i < j){
    swap_tmp(str[i], str[j], c);
    i++;
    j--;
  }
}

/* static inline void tidy_string(String *src, String *dst, char ch){ */
/*   int i; */
/*   encap_string(dst, src->size); */
/*   for(i=0;i<src->size;i++){ */
/*     if(src->string[i] != ch){ */
/*       dst->string[dst->size ++] = src->string[i]; */
/*     } */
/*   } */
/*   dst->string[dst->size] = 0; */
/* } */

static inline int occ_str(char *str, int len, char c){
  int i, ret;
  for(i=ret=0;i<len;i++){
    if(str[i] == c) ret ++;
  }
  return ret;
}

static inline void trunc_string(String *str, size_t size){
  if(size >= strlen(str->string)) return;
  str->string[size] = 0;
}

static inline String* clone_string(String *str){
  String *clone;
  clone = init_string(strlen(str->string)+1);
  append_string(clone, str->string);
  return clone;
}

static inline char* clone_native_string(String *str) {
  char *clone = malloc(strlen(str->string)+1);
  strcpy(clone, str->string);
  return clone;
}

static inline void free_string(String *str){ free(str->string); free(str); }

static inline void extend_string(String *str, char *source, int start, int end) {
  encap_string(str, end - start + 1);
  strncat(str->string, source + start, end - start + 1);
}

#endif
