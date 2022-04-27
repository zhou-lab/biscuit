/* visualize bisulfite alignment using ncurses
 * 
 * The MIT License (MIT)
 *
 * Copyright (c) 2016-2020 Wanding.Zhou@vai.org
 *               2021-2022 Jacob.Morrison@vai.org
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:

 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.

 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * 
**/

#include <ncurses.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <assert.h>
#include "wvec.h"
#include "sam.h"
#include "faidx.h"
#include "wzmisc.h"
#include "encode.h"
#include "kstring.h"

static int get_bsstrand(bam1_t *b) {
  uint8_t *s;
  s = bam_aux_get(b, "ZS");     /* bsmap flag */
  if (s) {
    s++;
    if (*s == '+') return 0;
    else if (*s == '-') return 1;
  }

  s = bam_aux_get(b, "YD");     /* bwa-meth flag */
  if (s) {
    s++;
    if (*s == 'f') return 0;
    else if (*s == 'r') return 1;
  }

  s = bam_aux_get(b, "XG");     /* bismark flag */
  if (s) {
    s++;
    if (strcmp((char*)s, "CT")==0) return 0;
    else if (strcmp((char*)s, "GA")) return 1;
  }

  /* otherwise, return failure */
  return -1;
}

/* read information node */
typedef struct __rnode_t {
  bam1_t *b;
  int row;                   /* which row the read should be placed, -1 for skip */
} rnode_t;

DEFINE_VECTOR(rnode_v, rnode_t);

void reset_rnode_v(rnode_v *read_buf) {
  unsigned i;
  for (i=0; i<read_buf->size; ++i) {
    rnode_t *nd = ref_rnode_v(read_buf, i);
    bam_destroy1(nd->b);
  }
  clear_rnode_v(read_buf);
}

/***************************
 * Bisulfite Terminal View *
 ***************************/

typedef struct btview_t {

  /* input files and indexes */
  samFile *fp;
  hts_idx_t *idx;
  bam_hdr_t *header;
  faidx_t *fai;

  /* dimension of the terminal window */
  int mrow, mcol;
  
  /**** plot data buffer ****/
  /* reference sequence and length in the window */
  char *ref; int l_ref;
  /* buffer of reads overlapping current window */
  rnode_v *read_buf;
  /* buffer boundary */
  int buf_flank;
  int buf_tid, buf_left, buf_right;
  /* genomic coordinates of screen */
  int curr_tid, left_pos, row_shift;

  /* display options */
  int show_short_format;      /* show short format of read */
  int show_name;              /* show read name instead of sequence */
  int inverse;                /* inverse color? */
  int color_for;              /* what does color stand for? */
  int base_for;               /* what does base stand for? */
  int is_dot;                 /* use dot instead of showing bases */
  int ins;                    /* whether to display insertion */
  int max_reads_per_pos;      /* max number of reads per position to load */
  char *read_name;            /* target read name to display */

  /* pop-up windows */
  WINDOW *wgoto, *whelp;
  
} btview_t;

#define TV_MIN_ALNROW 2
#define TV_MAX_GOTO  40
#define TV_LOW_MAPQ  10

#define TV_COLOR_MAPQ   0
#define TV_COLOR_BASEQ  1
#define TV_COLOR_NUCL   2
#define TV_COLOR_BSMODE 3

/*********************
 * setup color pairs *
 *********************/
static int btv_init_colors(int inverse) {
  if (inverse) {
    init_pair(1, COLOR_WHITE, COLOR_BLUE);
    init_pair(2, COLOR_BLACK, COLOR_GREEN);
    init_pair(3, COLOR_BLACK, COLOR_YELLOW);
    init_pair(4, COLOR_BLACK, COLOR_WHITE);
    init_pair(5, COLOR_BLACK, COLOR_GREEN);
    init_pair(6, COLOR_BLACK, COLOR_CYAN);
    init_pair(7, COLOR_WHITE, COLOR_MAGENTA);
    init_pair(8, COLOR_WHITE, COLOR_RED);
    init_pair(9, COLOR_WHITE, COLOR_BLUE);
  } else {
    init_pair(1, COLOR_BLUE, COLOR_BLACK);
    init_pair(2, COLOR_GREEN, COLOR_BLACK);
    init_pair(3, COLOR_YELLOW, COLOR_BLACK);
    init_pair(4, COLOR_WHITE, COLOR_BLACK);
    init_pair(5, COLOR_GREEN, COLOR_BLACK);
    init_pair(6, COLOR_CYAN, COLOR_BLACK);
    init_pair(7, COLOR_MAGENTA, COLOR_BLACK);
    init_pair(8, COLOR_RED, COLOR_BLACK);
    init_pair(9, COLOR_BLUE, COLOR_BLACK);
  }

  return 0;
}

/**********************
 * ncurses extensions *
 **********************/
static void vmvprintw(btview_t *tv,int y ,int x,const char* fmt,...) {
  unsigned int size=tv->mcol+2;
  char* str=malloc(size);
  if(str==0) exit(EXIT_FAILURE);
  va_list argptr;
  va_start(argptr, fmt);
  vsnprintf(str,size, fmt, argptr);
  va_end(argptr);
  mvprintw(y,x,str);
  free(str);
}

/* initialize view */
static btview_t *btv_init(const char *fn, const char *ref_fn) {

  /* read input bam and reference fa */
  btview_t *tv = calloc(1, sizeof(btview_t));
 
  if (ref_fn) {
    tv->fai = fai_load(ref_fn);
    if (!tv->fai) wzfatal("Cannot read file: %s.\n", ref_fn);
  }
  tv->fp = sam_open(fn, "r");
  if (!tv->fp) wzfatal("Cannot open sam file: %s.\n", fn);
  tv->header = sam_hdr_read(tv->fp);
  if (!tv->header) wzfatal("Cannot read %s.\n", fn);
  tv->idx = sam_index_load(tv->fp, fn);
  if (!tv->idx) wzfatal("Cannot read index for %s.\n", fn);

  /* read buffer */
  tv->read_buf = init_rnode_v(2);
  tv->buf_left = -1;
  tv->buf_right = -1;
  tv->buf_flank = 0;
  tv->max_reads_per_pos = 50;
  tv->read_name = NULL;

  /* initialize color */
  tv->color_for = TV_COLOR_BSMODE;
  tv->is_dot = 1;
 
  /* initialize screen */
  initscr();
  keypad(stdscr, TRUE);
  clear();
  noecho();
  cbreak();
  getmaxyx(stdscr, tv->mrow, tv->mcol);
  tv->wgoto = newwin(3, TV_MAX_GOTO + 10, 10, 5);
  tv->whelp = newwin(28, 40, 0, 5);
  start_color();
  btv_init_colors(0);
 
  return tv;
}

/* end view */
static void btv_destroy(btview_t *tv) {
  if (tv->fai) fai_destroy(tv->fai);
  if (tv->ref) free(tv->ref);
  reset_rnode_v(tv->read_buf);
  free_rnode_v(tv->read_buf);
  bam_hdr_destroy(tv->header);
  sam_close(tv->fp);
  hts_idx_destroy(tv->idx);
  delwin(tv->wgoto); delwin(tv->whelp);
  endwin();
  free(tv);
}

/* layout which row to put each read
   the following assumes reads are sorted by the start position */
static void set_row_update_endposes(rnode_t *nd, btview_t *tv, int *row_endposes) {
  unsigned i;
  bam1_t *b = nd->b;
  int endpos = bam_endpos(b);
  if (endpos < tv->left_pos || b->core.pos > tv->left_pos + tv->mcol) {
    nd->row = -1;
    return;
  }
  for (i=0; i<=tv->read_buf->size; ++i) {
    if (((b->core.pos > tv->left_pos) ? (b->core.pos - tv->left_pos) : 0) >= row_endposes[i]) {
      nd->row = i+2;
      row_endposes[i] = endpos - tv->left_pos + 5; /* update largest occupied position on the row */
      break;
    }
  }
}

/* assign reads to rows */
static void screen_layout_reads(btview_t *tv) {
  unsigned i;
  int *row_endposes = calloc(tv->read_buf->size, sizeof(int));
  for (i=0; i<tv->read_buf->size; ++i) {
    set_row_update_endposes(ref_rnode_v(tv->read_buf, i), tv, row_endposes);
  }
  free(row_endposes);
}

/* reload data, prepare for drawing */
static void btv_reload_data(btview_t *tv) {

  /* no need to reload buffer if the current window is still within buffer range */
  if (tv->buf_left >= 0 && tv->buf_right >= 0 && tv->curr_tid == tv->buf_tid
      && tv->buf_left+2 <= tv->left_pos && tv->buf_right >= tv->left_pos + tv->mcol + 2)
    return;

  /* reset buffer range */
  tv->buf_tid = tv->curr_tid;
  tv->buf_left = max(0, tv->left_pos-1-tv->buf_flank);
  tv->buf_right = min((int) (tv->header->target_len[tv->curr_tid]), tv->left_pos+tv->mcol+tv->buf_flank);

  /* retrieve reference, set tv->ref
     this is simpler as no concensus calling is drawn, 
     they aren't that useful for bisulfite sequencing */
  if (tv->fai) {
    char *str;
    if (tv->ref) { free(tv->ref); tv->ref = 0; }
    assert(tv->curr_tid>=0);

    str = (char*)calloc(strlen(tv->header->target_name[tv->curr_tid]) + 30, 1);
    assert(str!=NULL);
    sprintf(str, "%s:%d-%d", tv->header->target_name[tv->curr_tid], tv->buf_left+1, tv->buf_right);
    tv->ref = fai_fetch(tv->fai, str, &tv->l_ref);
    free(str);
    if (!tv->ref)
      wzfatal("Could not read the reference sequence.\n");
  }
  
  /* retrieve reads */
  bam1_t *b = bam_init1();

  int n = 1; hts_itr_t *iter;
  iter = sam_itr_queryi(tv->idx, tv->curr_tid, tv->buf_left, tv->buf_right);
  reset_rnode_v(tv->read_buf);
  int prev_pos = -1;
  while (sam_itr_next(tv->fp, iter, b) >= 0) {
    if (b->core.pos != prev_pos) {
      n = 1;
      prev_pos = b->core.pos;
    } else if (tv->read_name==NULL || strcmp(tv->read_name, bam_get_qname(b))!=0) {
      n++; 
      if (n>tv->max_reads_per_pos) continue;
    }
    rnode_t *nd = next_ref_rnode_v(tv->read_buf);
    nd->b = b;
    b = bam_init1();
  }
  bam_destroy1(b);
  hts_itr_destroy(iter);
}

/****************************
 ** main drawing subroutine *
 ****************************/
#define bscall(b, pos) seq_nt16_str[bam_seqi(bam_get_seq(b), pos)]
#define bicall(b, pos) seq_nt16_int[bam_seqi(bam_get_seq(b), pos)]

/* draw the bases of one read */
static void draw_read1(rnode_t *nd, btview_t *tv, int readattr, int bss) {

  bam1_t *b = nd->b; bam1_core_t *c = &b->core;
  uint32_t rpos = c->pos, qpos = 0;
  char qb, rb;
  int attr, i, x, isconv; unsigned j;
  for (i=0; i<c->n_cigar; ++i) {
    uint32_t op = bam_cigar_op(bam_get_cigar(b)[i]);
    uint32_t oplen = bam_cigar_oplen(bam_get_cigar(b)[i]);

    switch(op) {
    case BAM_CMATCH:
      for (j=0; j<oplen; ++j) {
        if (rpos + j < (unsigned) tv->left_pos) continue;
        qb = toupper(bscall(b, qpos + j));
        rb = toupper(tv->ref[rpos+j - tv->buf_left]);

        attr = readattr;
        if (tv->color_for == TV_COLOR_BSMODE) {
          isconv = 0;
          if (rb == 'G' && bss == 1) {
            if (qb == 'G') {    /* RED for retention */
              attr |= COLOR_PAIR(8);
            } else if (qb == 'A') { /* BLUE for conversion */
              isconv = 1;
              attr |= COLOR_PAIR(1);
            }
          } else if (rb == 'C' && bss == 0) {
            if (qb == 'C') {    /* RED for retention */
              attr |= COLOR_PAIR(8);
            } else if (qb == 'T') { /* BLUE for conversion */
              isconv = 1;
              attr |= COLOR_PAIR(1);
            }
          }
          /* when it's not a conversion (error or snp), mark in YELLOW */
          if (!isconv && qb != rb) attr |= COLOR_PAIR(3);
        } else if (tv->color_for == TV_COLOR_NUCL) {
          attr |= COLOR_PAIR(bicall(b, qpos+j) + 5);
        } else if (tv->color_for == TV_COLOR_BASEQ) {
          x = bam_get_qual(b)[qpos + j]/10 + 1;
          if (x > 4) x = 4;
          attr |= COLOR_PAIR(x);
        }

        attron(attr);
        /* exempt the retention under bisulfite mode */
        if (tv->is_dot && qb == rb
            && !(tv->color_for == TV_COLOR_BSMODE &&
                 ((bss == 0 && rb == 'C') || (bss == 1 && rb == 'G')))) {
          mvaddch(nd->row - tv->row_shift, rpos + j - tv->left_pos,
                  bam_is_rev(b) ? ',' : '.');
        } else {
          mvaddch(nd->row - tv->row_shift, rpos + j - tv->left_pos,
                  bam_is_rev(b) ? toupper(qb) : tolower(qb));
        }
        attroff(attr);
      }
      rpos += oplen;
      qpos += oplen;
      break;
    case BAM_CINS:
      qpos += oplen;
      break;
    case BAM_CDEL:
      for (j=0; j<oplen; ++j) {
        mvaddch(nd->row - tv->row_shift, rpos + j - tv->left_pos, '*');
      }
      rpos += oplen;
      break;
    case BAM_CSOFT_CLIP:
      qpos += oplen;
      break;
    case BAM_CHARD_CLIP:
      break;
    default:
      wzfatal("Unknown cigar, %u\n", op);
    }
  }
}

char *sam_short_format1(const bam_hdr_t *h, const bam1_t *b) {

   kstring_t str;
   str.l = str.m = 0; str.s = NULL;

   int i;
   const bam1_core_t *c = &b->core;
   kputw(c->flag, &str); kputc('|', &str); // flag
   if (c->tid >= 0) { // chr
      kputs(h->target_name[c->tid], &str);
      kputc('|', &str);
   } else kputsn("*|", 2, &str);
   kputw(c->pos + 1, &str); kputc('|', &str); // pos
   kputw(c->qual, &str); kputc('|', &str); // qual
   if (c->n_cigar) { // cigar
      uint32_t *cigar = bam_get_cigar(b);
      for (i = 0; i < c->n_cigar; ++i) {
         kputw(bam_cigar_oplen(cigar[i]), &str);
         kputc(bam_cigar_opchr(cigar[i]), &str);
      }
   } else kputc('*', &str);
   kputc('|', &str);
   if (c->mtid < 0) kputsn("*|", 2, &str); // mate chr
   else if (c->mtid == c->tid) kputsn("=|", 2, &str);
   else {
      kputs(h->target_name[c->mtid], &str);
      kputc('|', &str);
   }
   kputw(c->mpos + 1, &str); kputc('|', &str); // mate pos
   kputw(c->isize, &str); kputc('|', &str); // template len
   return str.s;
}

/* draw all alignments */
static void btv_drawaln(btview_t *tv, int re_layout) {

  assert(tv != NULL);
  clear();

  if (re_layout) {
    btv_reload_data(tv);
    screen_layout_reads(tv);
  }
  
  /* draw coordinates and reference */
  int i;
  for (i=1; i<tv->mcol-9; ++i) {
    int pos = tv->left_pos + i;
    if (pos % 20 == 0)
      vmvprintw(tv, 0, i-1, "|%-d", pos);
  }
  if (tv->fai) {
    for (i=0; i<tv->mcol; ++i) {
      int ii = i+tv->left_pos-tv->buf_left;
      unsigned char c = toupper(tv->ref[ii]);
      int attr = 0;
      if (tv->color_for == TV_COLOR_NUCL) {
        attr |= COLOR_PAIR(nt256char_to_nt256int8_table[c]+5);
      } else if (tv->color_for == TV_COLOR_BSMODE) {
        if (c == 'C') {
          if (i+1 < tv->l_ref && toupper(tv->ref[ii+1]) == 'G')
            attr |= COLOR_PAIR(8) | A_UNDERLINE;
          else
            attr |= COLOR_PAIR(1);
        } else if (c == 'G') {
          if (i > 0 && toupper(tv->ref[ii-1]) == 'C')
            attr |= COLOR_PAIR(8) | A_UNDERLINE;
          else
            attr |= COLOR_PAIR(1);
        }
      }
      attron(attr);
      mvaddch(1, i, c);
      attroff(attr);
    }
  } else {
    for (i=0; i<tv->mcol; ++i)
      mvaddch(1, i, 'N');
  }

  /* draw reads */
  unsigned u;
  for (u=0; u<tv->read_buf->size; ++u) {

    rnode_t *nd = ref_rnode_v(tv->read_buf, u);
    bam1_t *b = nd->b; bam1_core_t *c = &b->core;
    int bss = get_bsstrand(b);
    if (nd->row >= 0 && 
        nd->row >= 2+tv->row_shift && 
        nd->row < 2 + tv->row_shift + tv->mrow) {

      /************************
       ** read level coloring *
       ************************/
      int readattr = 0;

      /* color by mapping quality */
      if (tv->color_for == TV_COLOR_MAPQ) {
        int x = c->qual / 10 + 1;
        if (x > 4) x = 4;
        readattr |= COLOR_PAIR(x);
      }

      if (tv->read_name != NULL && strcmp(tv->read_name, bam_get_qname(b))==0) {
        readattr |= A_REVERSE;
      }
      
      /* set underscore for improper pair or secondary mapping */
      if (((c->flag & BAM_FPAIRED) && !(c->flag & BAM_FPROPER_PAIR))
          || (c->flag & BAM_FSECONDARY))
        readattr |= A_UNDERLINE;
      
      if (tv->show_name) {
        attron(readattr);
        mvprintw(nd->row - tv->row_shift, max(c->pos - tv->left_pos, 0), bam_get_qname(b));
        attroff(readattr);
      } else if (tv->show_short_format) {
        attron(readattr);
        mvprintw(nd->row - tv->row_shift, max(c->pos - tv->left_pos, 0), sam_short_format1(tv->header, b));
        attroff(readattr);
      } else {
        draw_read1(nd, tv, readattr, bss);
      }
    }
  }
}

/***************************
 ** the pop-up help window *
 ***************************/
static void btv_win_help(btview_t *tv) {
  int r = 1;
  WINDOW *win = tv->whelp;
  wborder(win, '|', '|', '-', '-', '+', '+', '+', '+');
  mvwprintw(win, r++, 2, "        -=-    Help    -=- ");
  r++;
  mvwprintw(win, r++, 2, "?          This window");
  mvwprintw(win, r++, 2, "Arrows     Small scroll movement");
  /* mvwprintw(win, r++, 2, "h,j,k,l    Small scroll movement"); */
  /* mvwprintw(win, r++, 2, "H,J,K,L    Large scroll movement"); */
  /* mvwprintw(win, r++, 2, "ctrl-H     Scroll 1k left"); */
  /* mvwprintw(win, r++, 2, "ctrl-L     Scroll 1k right"); */
  mvwprintw(win, r++, 2, "space      Scroll one screen");
  mvwprintw(win, r++, 2, "backspace  Scroll back one screen");
  mvwprintw(win, r++, 2, "g          Go to specific location");
  mvwprintw(win, r++, 2, "t          Color for bisulfite mode");
  mvwprintw(win, r++, 2, "m          Color for mapping qual");
  mvwprintw(win, r++, 2, "b          Color for base quality");
  mvwprintw(win, r++, 2, "n          Color for nucleotide");
  mvwprintw(win, r++, 2, ".          Toggle on/off dot view");
  mvwprintw(win, r++, 2, "s          Toggle on/off rd brief");
  mvwprintw(win, r++, 2, "r          Toggle on/off rd name");
  /* mvwprintw(win, r++, 2, "i          Toggle on/off ins"); */
  mvwprintw(win, r++, 2, "v          Inverse video");
  mvwprintw(win, r++, 2, "q          Exit");
  r++;
  mvwprintw(win, r++, 2, "Bisulfite Mode:");
  mvwprintw(win, r++, 2, "Blue:     Conversion;");
  mvwprintw(win, r++, 2, "Red:      Retention;");
  mvwprintw(win, r++, 2, "Yellow:   Other mismatches");
  r++;
  mvwprintw(win, r++, 2, "Underline:      Secondary or orphan");
  mvwprintw(win, r++, 2, "Blue:    0-9    Green: 10-19");
  mvwprintw(win, r++, 2, "Yellow: 20-29   White: >=30");
  wrefresh(win);
  wgetch(win);
}

/****************************
 ** the pop-up goto window **
 ****************************/
static void btv_win_goto(btview_t *tv, int *tid, int *pos) {
  
  char str[256], *p;
  int i, l = 0;

  /* draw window border */
  wborder(tv->wgoto, '|', '|', '-', '-', '+', '+', '+', '+');
  
  mvwprintw(tv->wgoto, 1, 2, "Goto: ");
  for (;;) {
    int invalid = 0;
    int c = wgetch(tv->wgoto);
    wrefresh(tv->wgoto);
    if (c == KEY_BACKSPACE || c == '\010' || c == '\177') {
      if(l > 0) --l;
    } else if (c == KEY_ENTER || c == '\012' || c == '\015') {

      /* parse tid and pos */
      int _tid = -1, _beg, _end;
      if (str[0] == '=') {      /* same tid, only update coordinate */
        _beg = strtol(str+1, &p, 10) - 1;
        if (_beg > 0) {
          *pos = _beg;
          return;
        }
      } else {
        char *name_lim = (char *) hts_parse_reg(str, &_beg, &_end);
        if (name_lim) {
          char name_terminator = *name_lim;
          *name_lim = '\0';
          _tid = bam_name2id(tv->header, str);
          *name_lim = name_terminator;
        } else { /* Unparsable region, but possibly a sequence named "foo:a" */
          _tid = bam_name2id(tv->header, str);
          _beg = 0;
        }

        if (_tid >= 0) {
          *tid = _tid; *pos = _beg;
          return;
        }
      }

      /* If we get here, the region string is invalid */
      invalid = 1;
    } else if (isgraph(c)) {
      if (l < TV_MAX_GOTO) str[l++] = c;
    } else if (c == '\027') {
      l = 0;
    } else if (c == '\033') {
      return;
    }
    str[l] = '\0';
    for (i = 0; i < TV_MAX_GOTO; ++i) mvwaddch(tv->wgoto, 1, 8 + i, ' ');
    if (invalid) mvwprintw(tv->wgoto, 1, TV_MAX_GOTO - 1, "[Invalid]");
    mvwprintw(tv->wgoto, 1, 8, "%s", str);
  }
}

/*************
 * main loop *
 *************/
static int btv_loop(btview_t *tv) {

  int r;                    /* whether to re-layout */
  while (1) {
    r = 0;
    int c = getch();
    switch (c) {
    case '?': btv_win_help(tv); break;
    case '\033':
    case 'q': goto end_loop;
    case '/':
    case 'g': btv_win_goto(tv, &tv->curr_tid, &tv->left_pos); r=1; break;
    case 't': tv->color_for = TV_COLOR_BSMODE; break;
    case 'm': tv->color_for = TV_COLOR_MAPQ; break;
    case 'b': tv->color_for = TV_COLOR_BASEQ; break;
    case 'n': tv->color_for = TV_COLOR_NUCL; break;
    case 'v': btv_init_colors(tv->inverse = !tv->inverse); break;
    case 's': tv->show_short_format = !tv->show_short_format; if (tv->show_short_format) tv->show_name = 0; break;
    case 'r': tv->show_name = !tv->show_name; if (tv->show_name) tv->show_short_format = 0; break;
    case KEY_LEFT:
    case 'h': --tv->left_pos; r=1; break;
    case KEY_RIGHT:
    case 'l': ++tv->left_pos; r=1; break;
    case KEY_SLEFT:
    case 'H': tv->left_pos -= 20; r=1; break;
    case KEY_SRIGHT:
    case 'L': tv->left_pos += 20; r=1; break;
    case '.': tv->is_dot = !tv->is_dot; break;
    case 'i': tv->ins = !tv->ins; break;
    case '\010': tv->left_pos -= 1000; r=1; break;
    case '\014': tv->left_pos += 1000; r=1; break;
    case ' ': tv->left_pos += tv->mcol; r=1; break;
    case KEY_UP:
    case 'j': --tv->row_shift; break;
    case KEY_DOWN:
    case 'k': ++tv->row_shift; break;
    case KEY_PPAGE: tv->row_shift -= 10; break;
    case KEY_NPAGE: tv->row_shift += 10; break;
    case KEY_BACKSPACE:
    case '\177': tv->left_pos -= tv->mcol; r=1; break;
    case KEY_RESIZE: getmaxyx(stdscr, tv->mrow, tv->mcol); r=1; break;
    default: continue;
    }
    if (tv->left_pos < 0) tv->left_pos = 0;
    if (tv->row_shift < 0) tv->row_shift = 0;
    btv_drawaln(tv, r);
  }
 end_loop:
  return 0;
}

static void usage() {
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage: biscuit tview [options] <in.bam> <ref.fa>\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "    -g STR    Go directly to this position\n");
  fprintf(stderr, "    -m INT    Max number of reads to load per position [50]\n");
  fprintf(stderr, "    -n STR    Highlight the read(s) with STR as the read name\n");
  fprintf(stderr, "    -f INT    Flanking sequence length [100]\n");
  fprintf(stderr, "    -h        This help\n");
  fprintf(stderr, "\n");
}

int main_tview(int argc, char *argv[]) {

  char *position = NULL;
  
  int max_reads_per_pos = 50;
  int buf_flank = 100;
  char *read_name = NULL; /* target read name */
  int c;
  if (argc<2) { usage(); return 1; }
  while ((c = getopt(argc, argv, ":g:m:n:f:h")) >= 0) {
      switch (c) {
          case 'g': position = optarg; break;
          case 'm': max_reads_per_pos = atoi(optarg); break;
          case 'n': read_name = strdup(optarg); break;
          case 'f': buf_flank = atoi(optarg); break;
          case 'h': usage(); return 1;
          case ':': usage(); wzfatal("Option needs an argument: -%c\n", optopt); break;
          case '?': usage(); wzfatal("Unrecognized option: -%c\n", optopt); break;
          default: usage(); return 1;
      }
  }

  char *bam_fn; bam_fn = (optind < argc) ? argv[optind++] : NULL;
  char *ref_fn; ref_fn = (optind < argc) ? argv[optind++] : NULL;
  if (!bam_fn) {
    usage();
    wzfatal("No input bam is given.\n");
  }

  if (!ref_fn) {
    usage();
    wzfatal("No reference sequence is given.\n");
  }

  btview_t *tv = btv_init(bam_fn, ref_fn);
  tv->max_reads_per_pos = max_reads_per_pos;
  tv->buf_flank = buf_flank;
  tv->read_name = read_name;
  
  /* if target position is given, parse that */
  if (position) {
    int tid, beg, end;
    char *name_lim = (char *) hts_parse_reg(position, &beg, &end);
    if (name_lim) *name_lim = '\0';
    else beg = 0;               /* region parsing failed */
    tid = bam_name2id(tv->header, position);
    if (tid >= 0) {             /* parsing successful and chromosome name exists */
      tv->curr_tid = tid, tv->left_pos = beg;
    }
  } else if (tv->fai) {
    /* following samtools, find the first sequence present in both BAM and reference */
    int i;
    for (i=0; i<tv->header->n_targets; ++i) {
      if (faidx_has_seq(tv->fai, tv->header->target_name[i])) break;
    }
    if (i == tv->header->n_targets) {
      wzfatal("None of the BAM sequence names present in the fasta file.\n");
    }
    tv->curr_tid = i;
  }

  btv_drawaln(tv, 1);
  btv_loop(tv);
  btv_destroy(tv);
  
  return 0;
}
