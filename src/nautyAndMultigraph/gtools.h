/*****************************************************************************
* This is the main header file for gtools.  nauty version 2.4.               *
* gtools.h.  Generated from gtools-h.in by configure.
*****************************************************************************/

/* The parts between the ==== lines are modified by configure when
creating gtools.h out of gtools-h.in.  If configure is not being
used, it is necessary to check they are correct.
====================================================================*/


#ifndef  _GTOOLS_H_    /* only process this file once */
#define  _GTOOLS_H_

#define HAVE_ERRNO_H  1      /* <errno.h> exists */
#define HAVE_PERROR  1          /* perror() exists */
#define HAVE_PIPE  1          /* pipe() exists */
#define HAVE_WAIT  1          /* wait() exists */
#define HAVE_POPEN  1          /* popen() and pclose() exist */
#define POPEN_DEC  1         /* popen() is declared in stdio.h */
#define FTELL_DEC  1         /* ftell() is declared in stdio.h */
#define FDOPEN_DEC  1        /* fdopen() is declared in stdio.h */
#define SORTPROG  "sort"         /* name of sort program */
#define SORT_NEWKEY 1  /* if -k is supported */
#define HAVE_PID_T 1    /* pid_t is defined */
#define PUTENV_DEC 1   /* putenv() is declared in stdlib.h */
#define SETENV_DEC 1   /* setenv() is declared in stdlib.h */
#define HAVE_PUTENV 1   /* putenv() exists */
#define HAVE_SETENV 1   /* setenv() exists */

/* @edit_msg@ */
/*==================================================================*/

#ifndef MAXN 
#define MAXN  0
#endif

#define SIZELEN(n) ((n)<=SMALLN?1:((n)<=SMALLISHN?4:8))
	/* length of size code in bytes */
#define G6LEN(n)  (((n)*((n)-1)/2+5)/6+SIZELEN(n))
	/* exactly graph6 string length excluding \n\0 */

#include "naututil.h"      /* which includes stdio.h */
#include "nausparse.h"

#if HAVE_ERRNO_H
#include <errno.h>
#else
extern int errno;
#endif

#if HAVE_PERROR
#define ABORT(msg) do {if (errno != 0) perror(msg); exit(1);} while(0)
#else
#define ABORT(msg) do {exit(1);} while(0)
#endif

#if PUTENV_DEC && HAVE_PUTENV
#define SET_C_COLLATION putenv("LC_COLLATE=C")
#elif SETENV_DEC && HAVE_SETENV
#define SET_C_COLLATION setenv("LC_COLLATE","C",1)
#elif HAVE_PUTENV
int putenv(char*);
#define SET_C_COLLATION putenv("LC_COLLATE=C")
#elif HAVE_SETENV
int setenv(const char*,const char*,int);
#define SET_C_COLLATION setenv("LC_COLLATE","C",1)
#else
#define SET_C_COLLATION
#endif

#if HAS_STDIO_UNLOCK && !defined(NAUTY_IN_MAGMA) && !defined(IS_JAVA)
#define FLOCKFILE(f) flockfile(f)
#define FUNLOCKFILE(f) funlockfile(f)
#define GETC(f) getc_unlocked(f)
#undef PUTC
#define PUTC(c,f) putc_unlocked(c,f)
#else
#define FLOCKFILE(f)
#define FUNLOCKFILE(f)
#define GETC(f) getc(f)
#undef PUTC
#define PUTC(c,f) putc(c,f)
#endif

#define BIAS6 63
#define MAXBYTE 126
#define SMALLN 62
#define SMALLISHN 258047
#define TOPBIT6 32
#define C6MASK 63

#define GRAPH6_HEADER ">>graph6<<"
#define SPARSE6_HEADER ">>sparse6<<"
#define PLANARCODE_HEADER ">>planar_code<<"
#define PLANARCODELE_HEADER ">>planar_code le<<"
#define PLANARCODEBE_HEADER ">>planar_code be<<"
#define EDGECODE_HEADER ">>edge_code<<"

#define GRAPH6         1
#define SPARSE6        2
#define PLANARCODE     4
#define PLANARCODELE   8
#define PLANARCODEBE  16
#define EDGECODE      32
#define PLANARCODEANY (PLANARCODE|PLANARCODELE|PLANARCODEBE)
#define UNKNOWN_TYPE 256
#define HAS_HEADER   512

#define ARG_OK 0
#define ARG_MISSING 1
#define ARG_TOOBIG 2
#define ARG_ILLEGAL 3

#define MAXARG 2000000000L
#define NOLIMIT (MAXARG+31L)

#define SWBOOLEAN(c,bool) if (sw==c) bool=TRUE;
#define SWINT(c,bool,val,id) if (sw==c) \
        {bool=TRUE;arg_int(&arg,&val,id);}
#define SWLONG(c,bool,val,id) if (sw==c) \
        {bool=TRUE;arg_long(&arg,&val,id);}
#define SWRANGE(c,sep,bool,val1,val2,id) if (sw==c) \
	{bool=TRUE;arg_range(&arg,sep,&val1,&val2,id);}

#ifdef HELPTEXT2 
#define PUTHELPTEXT printf("\nUsage: %s\n\n%s",USAGE,HELPTEXT1);\
		    printf("%s",HELPTEXT2);
#else
#define PUTHELPTEXT printf("\nUsage: %s\n\n%s",USAGE,HELPTEXT)
#endif

#define HELP if (argc > 1 && (strcmp(argv[1],"-help")==0 \
			   || strcmp(argv[1],"/?")==0 \
			   || strcmp(argv[1],"--help")==0)) \
       { PUTHELPTEXT; return 0;}
#define GETHELP \
fprintf(stderr,"   Use %s -help to see more detailed instructions.\n",argv[0])

#define alloc_error gt_abort

#define CATMSG0(fmt) sprintf(msg+strlen(msg),fmt)
#define CATMSG1(fmt,x1) sprintf(msg+strlen(msg),fmt,x1)
#define CATMSG2(fmt,x1,x2) sprintf(msg+strlen(msg),fmt,x1,x2)
#define CATMSG3(fmt,x1,x2,x3) sprintf(msg+strlen(msg),fmt,x1,x2,x3)
#define CATMSG4(fmt,x1,x2,x3,x4) sprintf(msg+strlen(msg),fmt,x1,x2,x3,x4)
#define CATMSG5(fmt,x1,x2,x3,x4,x5) sprintf(msg+strlen(msg),fmt,x1,x2,x3,x4,x5)
#define CATMSG6(fmt,x1,x2,x3,x4,x5,x6) \
		sprintf(msg+strlen(msg),fmt,x1,x2,x3,x4,x5,x6)

/************************************************************************/


#ifdef __cplusplus
extern "C" {
#endif
 
extern void gtools_check(int,int,int,int);
extern FILE *opengraphfile(char*,int*,boolean,long);
extern void writeline(FILE*,char*);
extern char *gtools_getline(FILE*);
extern int graphsize(char*);
extern void stringcounts(char*,int*,size_t*);
extern void stringtograph(char*,graph*,int);
extern size_t edgecount(char*);
extern graph *readg(FILE*,graph*,int,int*,int*);
extern graph *readg2(FILE*,graph*,int,int*,int*, char*);
extern char *ntog6(graph*,int,int);
extern char *ntos6(graph*,int,int);
extern char *sgtos6(sparsegraph*);
extern void writeg6(FILE*,graph*,int,int);
extern void writes6(FILE*,graph*,int,int);
extern void writes6_sg(FILE*,sparsegraph*);
extern void writepc_sg(FILE*,sparsegraph*);
extern void stringtosparsegraph(char*,sparsegraph*,int*);
extern sparsegraph *read_sg(FILE*,sparsegraph*);
extern sparsegraph *read_sg_loops(FILE*,sparsegraph*,int*);
extern sparsegraph *readpc_sg(FILE*,sparsegraph*);
extern sparsegraph *readpcle_sg(FILE*,sparsegraph*);
extern sparsegraph *read_ec(FILE*,sparsegraph*);
extern void writeec_sg(FILE*,sparsegraph*);
extern void writelast(FILE*);
extern int longval(char**,long*);
extern void arg_int(char**,int*,char*);
extern void arg_long(char**,long*,char*);
extern void arg_range(char**,char*,long*,long*,char*);
extern void writerange(FILE*,int,long,long);
extern void gt_abort(char*);
extern char *stringcopy(char*);
extern boolean strhaschar(char*,int);

extern void fcanonise(graph*,int,int,graph*,char*,boolean);
extern void fcanonise_inv
             (graph*,int,int,graph*,char*,void(*)(graph*,int*,int*,int,
               int,int,permutation*,int,boolean,int,int),int,int,int,boolean);
extern void fcanonise_inv_sg
           (sparsegraph*,int,int,sparsegraph*,char*,void(*)(graph*,int*,int*,
             int,int,int,permutation*,int,boolean,int,int),int,int,int,boolean);
extern void fgroup(graph*,int,int,char*,int*,int*);
extern void fgroup_inv
	     (graph*,int,int,char*,int*,int*,void(*)(graph*,int*,int*,int,
                int,int,permutation*,int,boolean,int,int),int,int,int);
extern int istransitive(graph*,int,int,graph*);
extern void tg_canonise(graph*,graph*,int,int);

extern int readg_code;
extern char *readg_line;
extern size_t ogf_linelen;
extern boolean is_pipe;

#ifdef __cplusplus
}
#endif

#ifdef CPUDEFS
CPUDEFS
#endif

/* @edit_msg@ */

#endif /* _GTOOLS_H_  */
