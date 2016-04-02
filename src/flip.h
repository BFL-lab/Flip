/**************************************************************************** 
*                            flip.h version 2.0.2                           *
*                                                                           *
*                         Copyright (C) 1997 OGMP                           *
*                  (Organelle Genome Megasequencing Project),               *
*                        Departement de Biochimie,                          *
*                         Universite de Montreal,                           *
*                    C.P. 6128, succursale Centre-ville,                    *
*                     Montreal, Quebec, Canada, H3C 2J7                     *
*                                                                           *
*                      Programming: Nicolas Brossard (OGMP)                 *
*               Project management: Gertraud Burger (OGMP)                  *
*               E-Mail information: b.franz.lang@gmail.com                  *
*                                                                           *
*     This software is distributed under the GNU GENERAL PUBLIC LICENSE, as *
* published by the Free Software Foundation. A copy of version 2 of this    *
* license should be included in a file called COPYING. If not, write to the *
* Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   *
****************************************************************************/

/*-----------------------------------------------------------------------------
  
   flip.h - header file for flip.c

  ---------------------------------------------------------------------------*/

/*--------------------------

 CONSTANTS  

  --------------------------*/

#define CODON_LG          3   /* The length (in nucleotide) of a codon       */
#define MAX_LINE_LG    5000   /* The maximum line length of ANY file read by
                                  flip. Longer lines are truncated           */
#define MAX_SEQ_LG  4000000   /* The maximum length a sequence (including all 
                                 ogmp spec. chars and spaces) can have       */
#define NUC_PER_LINE     60   /* Maximum nb of nuc (including spec. ogmp chars
                                  and spaces) to print on a line in all output
                                  files                                      */
#define AA_PER_LINE      60   /* Maximum nb of aa (including *) to print in 
                                  prot.lst and prot.src on a given line      */
#define NB_FRAMES         6   /* Six frames total                            */
#define NB_CHARS        256   /* Number of chars in ascii set                */
#define COMPL             1   /* Cst to indicate that some operation has to be
                                  done on file 'compl'                       */
#define NOCOMPL           0   /* Cst to indicate that some operation has to be
                                  done on file 'nocompl'                     */
#define DIR_LENGTH        3   /* The length of the arrow (==>) in an annot.  */
#define MAX_CTX_LG       12   /* The context is composed of two chunks of nucs
                                 This is the maximum length each of these can
                                 have                                        */
#define VERSION      "2.1.2"  /* Flip's current version number               */


/* Parameters default values */
#define DEF_FASTA_FORMAT       0 /* Default is not to put a '*' after each 
				    contig in prot.lst and prot.src */
#define DEF_GET_DNA            0 /* Default is not to create the file 
				    prot.lst.dna */
#define DEF_NUMBER_PROT        0 /* Default is not to append _<prot_nb> after
				    the contig name in the output files to
				    indicate the protein number */
#define DEF_SILENT             0 /* Default is to be verbose (ie display some
                                      informative messages                   */
#define DEF_CODE_ID            1 /* The default genetic code number (Basic)  */
#define DEF_CIRCULAR           0 /* Whether genome is circular or not (default
                                     linear)                                 */
#define DEF_MIN_CODING_REG_LG 20 /* Default minimum length an orf must have in
                                    order to be reported by flip. The stop 
                                    codon(if any) is ignored during lg calc  */
#define DEF_MIN_PROT_LG        0 /* Default minimum length (calculated from the
                                    first start codon) an orf must have in    
                                    order to be reported by flip. The stop  
                                    codon(if any) is ignored during lg calc.
				    If set to 0, the orf does not even need to
				    have a start to be reported */
#define DEF_TRANS_ALL          0 /* Whether to produce prot.6rf or not 
				    (defaults to no)                         */
#define DEF_FORCE_MET          0 /* The default value of param.force_met (see
				    this field for more description) */

/* Default output file names */
#define PROT_6RF_NAME      "prot.6rf"
#define PROT_LST_DNA_NAME  "prot.lst.dna"
#define PROT_LST_NAME      "prot.lst"
#define PROT_SRC_NAME      "prot.src"
#define NOCOMPL_NAME       "nocompl"
#define COMPL_NAME         "compl"
#define PROT_PRM_NAME      "prot.prm"

/*----------------------
 
  MACROS

  ----------------------*/

#define MIN(a,b)             ( ((a) < (b)) ? (a) : (b) )
#define MAX(a,b)             ( ((a) > (b)) ? (a) : (b) )

/* Is a valid genetic code_id or not */
/* Valid: 1..6, 9..19 */
#define IS_VALID_CODE_ID(a)  ( ((a) > 0) && ((a) <= NB_CODES) && ((a) != 7) \
                                && ((a) != 8) )

/* Returns 1 if two frame ids refer to frames that are on the same strand */
/* Frame ids on forw strand: 0,1,2  */
/* Frame ids on back strand: 3,4,5  */
#define SAME_STRAND(a,b) ( (( (a)<(NB_FRAMES/2) ) && ( (b)<(NB_FRAMES/2) )) ||\
                           (( (a)>=(NB_FRAMES/2) ) && ( (b)>=(NB_FRAMES/2) )) )
                               

/*----------------------

  TYPEDEFS

  ----------------------*/

typedef char CODON[CODON_LG+1];  /* A codon */

typedef struct pos_list{        /* A list of char pos in a string */
  long int         val;         /* The pos (0..str_lg-1) */
  struct pos_list *next;        /* Pointer on next element in list */
} POS_LIST;

typedef struct string_lst{      /* A list of strings */
  char              *str;       /* The string */
  struct string_lst *next;      /* Pointer on next element in list */
} STRING_LST;

typedef struct annot_lst{       /* List of all annots in a contig */
  long int          nb_nuc_bef; /* This includes ALL nucs and OGMP chars that 
                                   lie BEFORE an annot in a contig. SPACES ARE
				   NOT INCLUDED */
  char             *text;       /* The annot's text */
  struct annot_lst *next;       /* Pointer on next element in list */
} ANNOT_LIST;

typedef struct ctg {                 /* A contig */
  char        *name;                 /* The full contig header */
  char         seq[MAX_SEQ_LG+1];    /* The original sequence fetched from the 
					input file associated to the contig. It
					includes all OGMP chars and spaces
					(see function get_next_ctg). Original
					casing in input file is preserved */
  long int     seq_lg;               /* The length of the above string */
  char        *seq_no_blk;           /* The sequence but with all the blanks
					removed */
  long int     seq_no_blk_lg;        /* Length of the above sequence */
  char        *clean_seq;            /* The sequence with all blanks and spec
					OGMP chars removed, uppercased. It 
					should only contain ACGTN */
  long int     clean_seq_lg;         /* Length of the above seq. */
  char        *compl_clean_seq;      /* The cleaned sequence, reversed and
					complemented */
  char        *compl_seq_no_blk;     /* Seq_no_blk reversed and complemented */
  POS_LIST    *seq_pos;              /* Pointer on a list of all the positions
					in the original sequence of spec OGMP 
					chars and spaces */
  POS_LIST    *seq_no_blk_pos;       /* Pointer on a list of all the positions
					in the seq_no_blk of spec OGMP chars */
  POS_LIST    *compl_seq_pos;        /* Pointer on a list of all the positions
					in the reversed compl. sequence of 
					spec OGMP chars and spaces */
  POS_LIST    *compl_seq_no_blk_pos; /* Pointer on a list of all the positions
					in the compl_seq_no_blk of spec OGMP
					chars */
  ANNOT_LIST  *annot_lst;            /* Pointer on the list of (ordered) annots
					for the forw strand */
  ANNOT_LIST  *compl_annot_lst;      /* Pointer on the list of (ordered) annots
					for the backw strand */
} CTG;

typedef struct prot {      /* An orf */
  char         *ctg_name;  /* Full contig header the orf was found in */
  char         *aa;        /* The orf's complete aa sequence */
  long int      aa_lg;     /* Length of the above string */
  char         *dna;       /* The DNA sequence associated to the aa sequence
			      When translating the DNA sequence w/ the corresp.
			      genetic code, you should obtain AA seq, no matter
			      the strand on which it lies */
  long int      dna_lg;    /* Length of the above string */
  long int      length;    /* Length of the above string, but will not include
			      the ending '*' if any */
  struct prot  *next;      /* Pointer on the next orf found */
  long int      begin;     /* Index of the nucleotide (1..ctg_lg) at which the
			      orf begins. If the orf is on the reversed strand,
			      'begin' refers to the nuc at which the orf ENDS*/
  long int      end;       /* Index of the nucleotide (1..ctg_lg) at which the
			      orf ends. If the orf is on the reversed strand, 
			      'end' refers to the nuc at which the orf BEGINS*/
  int           strand;    /* 1 if the orf was found on the reverse strand, 0
			      otherwise */
  char         *comment;    /* All aas upstream from the orf's start codon. If
			       min_prot_lg is 0, then this comment is empty */
  long int      comment_lg; /* The length of the above string */
  char         *context;    /* Some context information (some nucs correspoding
			       to the start and end of the orf, see get_context
			       for exact context definition) */
} PROT;

typedef struct {                  /* The application's parameters */
  int          code_id;           /* The genetic code number used */
  int          circular;          /* Whether the gennome is circular or not */
  long int     min_coding_reg_lg; /* Minimum length an orf must have in order
				     to be reported by flip. The stop codon
				     (if any) is ignored during lg calc   */
  long int     min_prot_lg;       /* Minimum length (calculated from the 
                                     first start codon) an orf must have in   
                                     order to be reported by flip. The stop   
                                     codon(if any) is ignored during lg calc. 
				     If set to 0, the orf does not even need to
				     have a start to be reported */ 
  int          trans_all;         /* 1 if prot.6rf must be produced */
  char        *s_file;            /* Name of the input file */
  int         *starts;            /* Array of 125 ints, one for each codon. 
				     start[x] = 1 iff x is a codon index (see
				     get_codon_idx) of a start codon */
  int         force_met;          /* 1 => if the first codon of a region is a
				     start then it's automatically translated
				     as a metionine ('M'). Default 0. */
  int         silent;             /* Boolean. To display or not some 
				     informative messages */
  int         number_prot;        /* Boolean. Append (or not) "_<protein_nb>"
				     to the contig name in the output files 
				     to indicate the protein number */
  int         get_dna;            /* Boolean. Indicates if prot.lst.dna
				     should be produced */
  int         fasta_format;       /* Boolean. Whether or not to separate the 
				     contigs of compl and nocompl by '*' */
} PARAM;   

typedef struct frame{   /* Data associated to the reading frames */
                        /* 6 frames, indexed 0,1,2,3,4,5. Frame 0 and 5 are
			   always 'in sync'. Same for 1 and 4, 2 and 3 */
  char     *aa;         /* The complete aa sequence obtained when translating
			   the whole frame */
  long int  aa_lg;      /* Length of the above string */
  long int *stop;       /* The positions (0..aa_lg-1) of the aas which are stop
			   codons, according to the current gcode */
  long int  nb_stop;    /* The number of stops in the frame */
  long int *start;      /* The positions (0..aa_lg-1) of the aas which are 
                           start codons, according to the current gcode */
  long int  nb_start;   /* The number of starts in the frame */
  int       offset;     /* The number of leading nucs (in the cleaned (for fw)
			   or compl_clean (for back) sequence) discarded before
			   translation. Frames 0,1,2 always have 0,1,2 (resp.)
			   offset */
} FRAME[NB_FRAMES];

/*----------------------

  PROTOTYPES

  ----------------------*/

int      get_codon_idx( char *s, int idx[] );
int      get_next_ctg( FILE *file, CTG *ctg, int comp[], int is_nuc[],
		       int *read_ctg, PARAM *param, int is_seq_char[] );
int      is_blank_line( char *line );
char    *trim_seq( char string[] );
void     add_contig( CTG **mf, char *name, char *cur_seq, CTG **last_ctg );

void     print_ctg( CTG *ctg );

void     fill_ctg( CTG *ctg, int comp[], int is_nuc[] );

void     mem_err();
void     read_err( char *fname );

void     trans_all( CTG *ctg, char code[], FRAME frame, int idx[],
		    PARAM *param );

void     translate( FRAME f, char *seq, long int start, long int end, 
		    int idx[], PARAM *p, char code[] );
void     init_arrays( int idx[], int comp[], int is_nuc[], int is_aa[], 
		      int is_seq_char[] );

void     open_all_files( FILE **prot_lst_dna, FILE **prot_6rf , FILE **seq_file,
			 FILE **prot_lst, FILE **prot_src, PARAM *param, 
			 FILE **nocompl, FILE **compl, char *seq_file_name );
char    *compl_seq( char seq[], long int lg, int comp[] );
void     update_prot_6rf( FILE *prot_6rf, FRAME frame, CTG *ctg );
void     close_all_files( FILE **prot_lst_dna, FILE **prot_6rf, FILE **seq_file,
			  FILE **prot_lst, FILE **prot_src, FILE **nocompl, 
			  FILE **compl, PARAM *param );
void     find_prot( CTG *ctg, PROT **prot, PARAM *param, FRAME frame );
long int get_p_start( FRAME f, long int r_start, long int r_end, 
		      long int *cur_s_idx );
void     add_prot( PROT **prot, PROT **last_prot, long int r_start, 
		   long int r_end, long int p_start, FRAME f, PARAM *param,
		   CTG *ctg, int f_nb , long int l );
void     update_prot_lst_src( PROT **sorted_prot, int nb_prot, FILE *prot_lst, 
			      FILE *prot_src, int number_prot );
void     free_mem( CTG *ctg, PROT *prot, FRAME frame, PROT **sorted_prot );
void     print_pos( FRAME frame, char *name );
void     add_annot( ANNOT_LIST **lst_start, ANNOT_LIST **lst_end, 
		    ANNOT_LIST **compl_lst_start, ANNOT_LIST **compl_lst_end,
		    char *annot, long int nb_nuc, int *prec_cont );
void     update_seq_file( CTG *ctg, FILE *nocompl, int is_compl, int nb_ctg,
			  int nb_non_empty_ctg, int fasta_format );
void     update_nocompl( CTG *ctg, FILE *compl );
int      get_nb_spec( POS_LIST **pos_ptr, long int start, long int end );
int      ends_with_bksl( char *s );
void     get_context( PROT *prot, CTG *ctg );
void     parse_cmd_line( int argc, char *argv[], PARAM *param, int is_nuc[], 
			 int is_aa[], char all_codes[NB_CODES][NB_CODONS], 
			 int all_starts[NB_CODES][NB_CODONS], char **code,
			 int idx[] );
void     usage();
int      get_starts( char *start_string, PARAM *param, int is_nuc[], int idx[],
		     char *code );
int      get_deviations( char *start_string, int is_nuc[], int is_aa[],
			 char *code );
void     parse_prot_prm( PARAM *param, char **dev_string, 
			 char **starts_string );
int      is_comment( char *string );
void     parse_prot_prm_err( char *line );
long int nuc_lg( PROT *prot, CTG *ctg, long int seq_start, long int seq_end );
void     print_stops( FRAME frame );
void     find_wrap_stops( FRAME frame, CTG *ctg, int idx[], char code[], 
			  PARAM *param );
void     find_wrap_prot( PROT **prot, PROT **last_prot, CTG *ctg, PARAM *param,
		         FRAME frame );
long int get_wrap_rg_lg( CTG *ctg, FRAME frame, int order[], int f );
long int get_wrap_prot_lg( CTG *ctg, FRAME frame, int order[], int start_f, 
			   int end_f, long int r_lg, int *start_o_idx, 
			   long int *start_pos );
void     add_wrap_prot( PROT **prot, PROT **last_prot, CTG *ctg, PARAM *param, 
			FRAME frame, int order[], int f, int last_f, 
			int start_o_idx, long int start_pos, long int r_lg,
			long int prot_lg );
void     get_wrap_aa( PROT *new_prot, PARAM *param, FRAME frame, int f,
		      CTG *ctg, char code[], int idx[], int order[], 
		      long int r, int start_o_idx, long int start_pos );
void     get_wrap_pos( PROT *new_prot, CTG *ctg, FRAME frame, int f, 
		       int last_f, int start_o_idx, long int start_pos,
		       PARAM *param, int order[] );
void     get_wrap_comment( PROT *new_prot, FRAME frame, int f, CTG *ctg,
			   char code[], int idx[], int order[],
			   int start_o_idx, long int start_pos, long int size);
void     get_wrap_context( PROT *prot, CTG *ctg );
char    *find_arrow_ptr( char *annot, char *arrow_ptr, char *compl_arrow_ptr );
int      is_valid_line( char *line, int is_seq_char[] );
void     fix_starts( CTG *ctg, PROT *prot, PARAM *param, int idx[], 
		     int is_nuc[], int comp[] );
int      prot_compare( const void *p1, const void *p2 );
void     sort_prots( PROT *prot, PROT ***sorted_prot, int *nb_prot );
int      is_valid_sfile( char *fname );
void     get_dna_seq( PROT *prot, CTG *ctg );
void     update_prot_lst_dna( PROT **sorted_prot, long int nb_prot, 
			      FILE *prot_lst_dna, int number_prot );
