/**************************************************************************** 
*                            flip.c version 2.0.2                           *
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

/*----------------------------------------------------------------------------

  flip.c

  An orf finder/translator

  ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <ctype.h>
#include "code.h"             /* All the genetic codes & their start codons */
#include "flip.h"             

/*----------------------

  GLOBAL VARIABLES

  ----------------------*/

PARAM  param;                 /* The application's parameters */
CTG    ctg;                   /* Current contig */
FILE  *seq_file=NULL;         /* The input file pointer */
char  *code;                  /* Genetic code used. 125 ints (one foreach codon
				 idx). */
FRAME  frame;                 /* The 6 reading frames */
int    idx[NB_CHARS];         /* Used to get the index of each nuc (see 
				 get_codon_idx)  */
int    comp[NB_CHARS];        /* To find the complement of a nucleotide */
int    is_nuc[NB_CHARS];      /* Whether a char is a nuc or not */
int    is_aa[NB_CHARS];       /* Whether a char is a valid amino acid or not */
int    is_seq_char[NB_CHARS]; /* If a char is a sequencing char or not (ie
				 aAcCgGtTnNxX#@!+- or any whitespace*/
FILE  *prot_6rf=NULL;         /* Pointer on the file of all 6 reading frames */
FILE  *prot_lst_dna=NULL;     /* Pointer on the file prot.lst.dna */
FILE  *prot_lst=NULL;         /* Pointer on the file prot.lst */
FILE  *prot_src=NULL;         /* Pointer on the file prot.src */
PROT  *prot;                  /* Pointer on the list of proteins found */
FILE  *nocompl=NULL;          /* Pointer on the file nocompl */
FILE  *compl=NULL;            /* Pointer on the file compl */
int    read_ctg = 0;          /* If a ctg has been read or not */
int    nb_non_empty_ctg = 0;  /* Nb of non-empty (ie at least one nuc or OGMP
				 spec char) contig written to compl and nocompl
				 so far */
int    nb_ctg = 0;            /* Nb of contigs (empty or non-emptyz) written to
				 compl and nocompl so far */
PROT **sorted_prot;           /* array of sorted proteins */
int    nb_prot;               /* the number of found proteins/orf (elements in
			         sorted_prot) */
 


/*----------------------

  MAIN PROGRAM

  ----------------------*/

int main( int argc, char *argv[] ){
  printf( "Turbo flip %s\n", VERSION );
  init_arrays( idx, comp, is_nuc, is_aa, is_seq_char );

  parse_cmd_line( argc, argv, &param, is_nuc, is_aa, all_codes, all_starts,
		  &code, idx );
  open_all_files( &prot_lst_dna, &prot_6rf, &seq_file, &prot_lst, &prot_src,
		  &param, &nocompl, &compl, param.s_file );
  while( get_next_ctg( seq_file, &ctg, comp, is_nuc, &read_ctg, &param, 
		       is_seq_char ) )
  {
      #ifdef DEBUG
      print_ctg( &ctg ); 
      #endif

      prot        = NULL;
      sorted_prot = NULL;
      if( ctg.seq_lg )
      {
	trans_all( &ctg, code, frame, idx, &param );
        
        #ifdef DEBUG
        print_stops( frame ); 
        #endif

	if( param.trans_all ) update_prot_6rf( prot_6rf, frame, &ctg );
	find_prot( &ctg, &prot, &param, frame );
	sort_prots( prot, &sorted_prot, &nb_prot );
	if( param.force_met ) 
	  fix_starts( &ctg, prot, &param, idx, is_nuc, comp );
	update_prot_lst_src( sorted_prot, nb_prot, prot_lst, prot_src, 
			     param.number_prot );
	if( param.get_dna ) 
	  update_prot_lst_dna( sorted_prot, nb_prot, prot_lst_dna, 
			       param.number_prot );
      }
      update_seq_file( &ctg, nocompl, NOCOMPL, nb_ctg, nb_non_empty_ctg,
		       param.fasta_format );
      update_seq_file( &ctg, compl, COMPL, nb_ctg, nb_non_empty_ctg,
		       param.fasta_format );
      //free_mem( &ctg, prot, frame, sorted_prot );
      nb_ctg++;
      if( ctg.seq_lg ) nb_non_empty_ctg++;
  }
  
  close_all_files( &prot_lst_dna, &prot_6rf, &seq_file, &prot_lst, &prot_src,
		   &nocompl, &compl, &param);
  exit(0);
}

/*----------------------

  SUBROUTINES

  ----------------------*/

/*-----------------------------------------------------------------------------

  Name: print_stops

  Description: Prints all the positions (0..frame_lg-1) of the stops and starts
               for all the frames.
  
  Returns: Nothing
  
  ---------------------------------------------------------------------------*/


void print_stops( 
		 FRAME frame  /* All the six frames */
 )
{
  int i,j;  /* Loop counters */

  /* Print all stops */
  for( i=0; i < NB_FRAMES ; i++ ){
    printf( "Stop pos for frame %d (offset %d):\n\t", i, frame[i].offset );
    for( j=0; j < frame[i].nb_stop ; j++ )
      printf( "%ld ", frame[i].stop[j] );
    printf( "\n" );
  }

  /* Print all starts */
  printf( "\n" );
  for( i=0; i < NB_FRAMES ; i++ ){
    printf( "Start pos for frame %d (offset %d):\n\t", i, frame[i].offset );
    for( j=0; j < frame[i].nb_start ; j++ )
      printf( "%ld ", frame[i].start[j] );
    printf( "\n" );
  }
}

/*-----------------------------------------------------------------------------

  Name: mem_err

  Description: Prints an error message that a memory error occured and exits
               with an exit status of 1
  
  Returns: Nothing
  
  ---------------------------------------------------------------------------*/

void mem_err()
{
  fprintf( stderr, "No more memory. Aborting.\n" );
  exit(1);
}

/*-----------------------------------------------------------------------------

  Name: get_codon_idx

  Description: Returns the index associated to a given codon. All 125 codons
               (yes, codons with 'N' are included and only uppercased codons
	       have a defined index) are numbered from 0 to 124. Here are the 
	       first few:
                           AAA=0, AAC=1, AAG=2, AAT=3, AAN=4,
			   ATA=5, ATC=6, ATG=7, ATT=8, ATN=9, etc...
               It is assumed that s points to a string that contains at least 3
	       chars

  Returns: The index
  
  ---------------------------------------------------------------------------*/

int get_codon_idx( 
		  char *s,       /* The codon beginning                     */
		  int idx[]      /* The index of all nucs (see global vars) */
)
{
  return( ( 25 * idx[*s] ) + ( 5 * idx[ *(s+1) ] ) + idx[ *(s+2) ] );
}


/*-----------------------------------------------------------------------------

  Name: get_next_ctg

  Description: Gets the next contig in the input file and fills the current
               contig structure. The input file is taken to be in masterfile or
	       fasta or extended fasta format. There's a limit on the length of
	       the lines in the file (see constants). Contig headers can be 
	       anything (even duplicates). Verifies that there is a unique
	       contig if genome is circular.

  Returns: 1 if there are no more contigs left to read, 0 otherwise
  
  ---------------------------------------------------------------------------*/


int get_next_ctg( 
	FILE *file,        /* The input file's pointer */
	CTG *ctg,          /* Adress of the current contig structure */
	int compl[],       /* Array of complement foreach nuc */
	int is_nuc[],      /* Array telling if char is a nuc or not */
	int *read_ctg,     /* 1 if a contig was ALREADY read, else 0 */
	PARAM *param,      /* Adress of the param structure */
	int is_seq_char[]  /* Array telling if a char is one of 
		              aAcCgGtTnNxX#@!+- or any whitespace*/
)
{
  char        line[MAX_LINE_LG+1]; /* Buffer. A line */
  char       *seq_start;           /* Will point the start of the ctg's nuc 
				      sequence */
  int         c;                   /* A nucleotide */
  int         line_lg;             /* Length of the line read */
  long int    nb_nuc;              /* nb of nucs read so far */
  ANNOT_LIST *last_annot;          /* Pointer on the last annot read */
  ANNOT_LIST *compl_last_annot;    /* pointer on the last annot of the list of
				      annots on the compl strand */
  int         first_line;          /* 1 if current line is the 1st one read */
  int         prec_cont;           /* 1 if the last line read was a line ending
				      with a "\" (continuation line) */

  /* init section */
  ctg->name = NULL;
  ctg->annot_lst = last_annot = ctg->compl_annot_lst = compl_last_annot = NULL;
  ctg->seq[0] = '\0';
  nb_nuc = 0;
  first_line = 1;
  prec_cont = 0;

  /* While not eof */
  while( ( c = getc( file ) ) != EOF ){

    /* A contig header */
    if( c == '>' ) {
      ungetc( c, file );

      /* if we already read one, stop */
      if( ctg->name != NULL ) break;
      else

        /* else, store it */
	if( first_line )
	  if( fgets( line, MAX_LINE_LG, file ) != NULL ){
	    line_lg = strlen( line );
	    if( !(ctg->name = (char *) calloc( line_lg+1, sizeof(char) )) ) 
	      mem_err();
	    strcpy( ctg->name, line );
            ctg->name[line_lg-1] = '\0';  /* Remove '\n' */
	  }
	  else read_err( param->s_file );
        else break;
    }     
    /* Then it's a nucleotide line (with optional leading digits) or annot or
       blank line */
    else{
      ungetc( c, file );
      if( fgets( line, MAX_LINE_LG, file ) != NULL ){
	if( !is_blank_line( line ) )
	  if( line[0] != ';' ) {
            if( prec_cont ){
	      fprintf( stderr, "Invalid continuation line:\t%s\n", line );
	      exit(1);
	    }
	    prec_cont = 0;
	    seq_start = trim_seq( line );
            if( !is_valid_line( seq_start, is_seq_char ) ){
	      fprintf( stderr, "Invalid line in %s:\n%s\n", param->s_file,
		       line );
	      exit(1);
	    }
	    nb_nuc += strlen( seq_start );
	    strcat( ctg->seq, seq_start );
	  }
          else add_annot( &(ctg->annot_lst), &last_annot, 
			  &(ctg->compl_annot_lst), &compl_last_annot, line,
			  nb_nuc, &prec_cont );
      }
      else read_err( param->s_file );
    }   /* else */
    
    first_line = 0;
  } /* while */  

  /* Check that contig is unique if genome is circular */
  if( param->circular && *read_ctg && ctg->name ){
    fprintf( stderr, "Multiple contigs for circular genome. Aborting.\n" );
    exit(1);
  }

  /* If contig is not empty and no header was found */
  if( strlen( ctg->seq ) && !ctg->name ){
    ctg->name = (char *) calloc( strlen( param->s_file ) + 2, sizeof(char) );
    ctg->name[0] = '>';
    strcpy( ctg->name + 1, param->s_file );
  }

  /* Fill missing info / adjust some info in ctg */
  fill_ctg( ctg, comp, is_nuc );

  if( strlen( ctg->seq ) ) *read_ctg = 1;
  return( !( feof(file) && !ctg->name ) );

} /* read_seq_file */


/*-----------------------------------------------------------------------------

  Name: is_blank_line

  Description: Verifies if a line is composed entirely of whitespace chars
               The null string is considered to be a blank line.
  
  Returns: 1 if it is so, 0 otherwise
  
  ---------------------------------------------------------------------------*/


int is_blank_line( 
		  char line[] /* Pointer on an array of characters */
)
{
  int i;

  for( i=0 ; i < strlen( line ) ; i++ ){
    if( !isspace( line[i] ) ) return(0);
  }
  return(1);
}

  
/*-----------------------------------------------------------------------------

  Name: read_err

  Description: Prints an error message that a read error occured while reading
               a given file and exit with an exit status of 1
  
  Returns: Nothing
  
  ---------------------------------------------------------------------------*/

void read_err( 
	      char *fname  /* The name of the file in which a reading error
			      occured */
)
{
  fprintf( stderr, "Unexpected read error on file %s. Aborting\n" , fname );
  exit(1);
}


/*-----------------------------------------------------------------------------

  Name: trim_seq

  Description: Eliminates all the leading and trailing spaces in a string
  
  Returns: A pointer on the beginning of the trimed sequence
  
  ---------------------------------------------------------------------------*/

char * trim_seq( 
		char string[]  /* A pointer on the string to trim */
)
{
  int i;  /* loop counter */

  /* this should trim the leading spaces and/or tabs */
  while( (*string) != '\0' && ( isspace( *string ) || isdigit( *string ) ) ){
    string++;
  }

  /* position just before the '\0' and "trim" the trailing spaces */
  i = strlen( string ) - 1;
  while( i >= 0 && isspace( string[i] ) ){
    i--;
  }
  string[i+1] = '\0'; /* readjust the end */

  return( string );
}


/*-----------------------------------------------------------------------------

  Name: print_ctg

  Description: Debug routine to print all the information stored in memory 
               regarding the current contig
  
  Returns: Nothing
  
  ---------------------------------------------------------------------------*/

void print_ctg( 
	       CTG *ctg  /* pointer on the current contig */
)
{
  char       *seq_ptr;               /* Tmp sequence ptr */
  char        chunk[NUC_PER_LINE+1]; /* Some nucleotides */
  ANNOT_LIST *ptr;                   /* pointer on an annot list */
  POS_LIST   *pos_ptr;               /* pointer on a list of positions */


  /* Print original sequence, NUC_PER_LINE nucs per line */
  if( ctg->name ) printf( "%s\n", ctg->name );
  printf( "Seq_lg = %ld,clean_seq_lg = %ld,", ctg->seq_lg, ctg->clean_seq_lg );
  printf( "seq_w_no_blanks_lg = %ld\n", ctg->seq_no_blk_lg );
  if( ctg->seq )
    for( seq_ptr = ctg->seq; strlen( seq_ptr ) > 0; ){
      strncpy( chunk, seq_ptr, NUC_PER_LINE );
      if( NUC_PER_LINE <= strlen( seq_ptr ) ) chunk[NUC_PER_LINE] = '\0';
      printf( "%s\n", chunk );
      seq_ptr += 
	strlen( seq_ptr ) > NUC_PER_LINE ? NUC_PER_LINE : strlen( seq_ptr );
    }

  /* Print cleaned sequence, NUC_PER_LINE nucs per line */  
  if( ctg->name ) printf( "%s (cleaned)\n", ctg->name );
  if( ctg->clean_seq )
    for( seq_ptr = ctg->clean_seq; strlen( seq_ptr ) > 0; ){
      strncpy( chunk, seq_ptr, NUC_PER_LINE );
      if( NUC_PER_LINE <= strlen( seq_ptr ) ) chunk[NUC_PER_LINE] = '\0';
      printf( "%s\n", chunk );
      seq_ptr += 
	strlen( seq_ptr ) > NUC_PER_LINE ? NUC_PER_LINE : strlen( seq_ptr );
    }
  
  printf( "\n" );
  
  /* Print all the positions of the non-nucs (ie spaces and spec. ogmp chars)
     in original sequence */
  printf( "The non nuc pos (forw):\n" );
  for( pos_ptr = ctg->seq_pos; pos_ptr; pos_ptr = pos_ptr->next )
    printf( "%ld  ", pos_ptr->val );
  printf( "\n" );

  /* Print all the non-nucs chars per se (ie spaces and spec. ogmp chars)
     in original sequence */
  printf( "The non nuc chars (forw):\n" );
  for( pos_ptr = ctg->seq_pos; pos_ptr; pos_ptr = pos_ptr->next )
    printf( "%c,", ctg->seq[ pos_ptr->val ] );
  printf( "\n" );
 
  /* Print all the positions of the non-nucs (ie spaces and spec. ogmp chars)
     in reversed complemented sequence */
  printf( "The non nuc pos (back):\n" );
  for( pos_ptr = ctg->compl_seq_pos; pos_ptr; pos_ptr = pos_ptr->next )
    printf( "%ld  ", pos_ptr->val );
  printf( "\n" );
  
  /* Print all the non-nucs chars per se (ie spaces and spec. ogmp chars)
     in reversed complemented sequence */
  printf( "The non nuc chars (back):\n" );
  for( pos_ptr = ctg->compl_seq_pos; pos_ptr; pos_ptr = pos_ptr->next )
    printf( "%c,", ctg->seq[ ctg->seq_lg - pos_ptr->val - 1 ] );
  printf( "\n" );

  /* Print all the positions of the non-nucs (ie spec. ogmp chars)
     in original sequence whose blanks have already been removed */
  printf( "The seq_no_blk non nuc pos (forw):\n" );
  for( pos_ptr = ctg->seq_no_blk_pos; pos_ptr; pos_ptr = pos_ptr->next )
    printf( "%ld  ", pos_ptr->val );
  printf( "\n" );

  /* Print all the positions of the non-nucs (ie spec. ogmp chars)
     in reversed complemented  sequence whose blanks have already been 
     removed */
  printf( "The compl_seq_no_blk non nuc pos (back):\n" );
  for( pos_ptr = ctg->compl_seq_no_blk_pos; pos_ptr; pos_ptr = pos_ptr->next )
    printf( "%ld  ", pos_ptr->val );
  printf( "\n" );

  /* print list of forw and back annotations */
  printf( "Here is the list of all annotations (forw):\n" );
  for( ptr = ctg->annot_lst ; ptr ; ptr = ptr->next )
    printf( "%ld --- Text: %s", ptr->nb_nuc_bef, ptr->text );
  printf( "Here is the list of all annotations (back):\n" );
  for( ptr = ctg->compl_annot_lst ; ptr ; ptr = ptr->next )
    printf( "%ld --- Text: %s", ptr->nb_nuc_bef, ptr->text );
}


/*-----------------------------------------------------------------------------

  Name: fill_ctg

  Description: Fills the remaining missing info in the current contig structure
               (ie calculate seq_lg, complemented and reversed sequence, get
	       the position of all the non-nucs, etc...)
  
  Returns: Nothing
  
  ---------------------------------------------------------------------------*/

void fill_ctg( 
	      CTG *ctg,     /* Pointer on current contig structure */
	      int comp [],  /* Array of complemented nucs */
	      int is_nuc[]  /* Array telling if a char is a nuc or not */
)
{
  POS_LIST   *ptr_int_list  = NULL;      /* beginning of list of non-nuc pos in
					    orig sequence (forw) */
  POS_LIST   *last_int      = NULL;      /* end of list of non-nuc pos in
					    original sequence (forw) */
  POS_LIST   *ptr_int_list2 = NULL;      /* beginning of list of OGMP chars
					    in orig sequence w/no blanks */
  POS_LIST   *last_int2     = NULL;      /* end of list of OGMP chars in 
					    orig sequence w/no blanks */
  char       *seq_ptr       = ctg->seq;  /* beginning of original sequence */
  POS_LIST   *new_int, *new_int2;        /* tmp new positions */
  POS_LIST   *cur_pos;                   /* pointer on a pos element */
  char       *clean_seq;                 /* pointer on the beginning of the 
					    cleaned sequence */
  char       *seq_no_blk;                /* pointer on the beginning of the 
					    sequence w/no blanks */
  char       *cur_ptr1, *cur_ptr2;       /* tmp seq ptrs */
  long int    i;                         /* loop counter */
  long int    nb_spaces;                 /* nb of spaces seen so far in orig 
					    sequence */
  char        c;                         /* a nuc */
  ANNOT_LIST *ptr;                       /* tmp pointer on an annot */

  /* Allocate max space for cleaned seq. */
  ctg->seq_lg          = strlen( seq_ptr );
  clean_seq  = (char *) calloc( ctg->seq_lg + 1, sizeof( char ) );
  seq_no_blk = (char *) calloc( ctg->seq_lg + 1, sizeof( char ) );

  /* Compute the uppercased 'cleaned' (ie no other nuc. than ACGTN) sequence */
  cur_ptr1 = clean_seq;
  cur_ptr2 = seq_no_blk;
  for( i=0 ; i < ctg->seq_lg ; i++ ){
    c = seq_ptr[i];

    /* if it's a nucleotide, add it in uppercase in clean sequence and as is
       in seq w/no blanks */
    if( is_nuc[c] ){
      *cur_ptr1++ = toupper(c);
      *cur_ptr2++ = c;
    }
    /* if not a nucleotide */
    else{
      /* add position of non-nuc in list of non-nuc associated to original
	 sequence*/
      new_int = (POS_LIST *) malloc( sizeof( POS_LIST ) );
      new_int->val = i, new_int->next = NULL;

      if( last_int == NULL ) ptr_int_list = last_int = new_int;
      else last_int->next = new_int, last_int = new_int;

      /* If not a space add char to seq_no_blk. Add pos to list of non-nuc pos
         for sequence w/no blanks */
      if( !isspace(c) ){
	*cur_ptr2++ = c; 
        new_int2 = (POS_LIST *) malloc( sizeof( POS_LIST ) );
	new_int2->val = cur_ptr2 - seq_no_blk - 1, new_int2->next = NULL;
	
	if( last_int2 == NULL ) ptr_int_list2 = last_int2 = new_int2;
	else last_int2->next = new_int2, last_int2 = new_int2;
      }
    }
  }
  /* terminate cleaned sequence string and seq without blanks string */
  *cur_ptr1 = '\0';    
  *cur_ptr2 = '\0';    

  /* update the info pointed to by ctg */
  ctg->clean_seq       = clean_seq;
  ctg->clean_seq_lg    = strlen( clean_seq );

  ctg->seq_no_blk      = seq_no_blk;
  ctg->seq_no_blk_lg   = strlen( seq_no_blk );

  ctg->compl_seq_no_blk = 
    compl_seq( ctg->seq_no_blk, ctg->seq_no_blk_lg, comp );
  ctg->compl_clean_seq  = compl_seq( ctg->clean_seq, ctg->clean_seq_lg, comp );

  ctg->seq_pos          = ptr_int_list;
  ctg->seq_no_blk_pos   = ptr_int_list2;
  
  /* Update positions in list of non-nuc positions for complementary seq */
  ctg->compl_seq_pos = NULL;
  for( cur_pos = ctg->seq_pos; cur_pos ; cur_pos = cur_pos->next ){
    new_int = (POS_LIST *) malloc( sizeof( POS_LIST ) );
    new_int->val = ctg->seq_lg - cur_pos->val - 1;
    new_int->next = ctg->compl_seq_pos;
    ctg->compl_seq_pos = new_int;
  }

  /* Update positions in list of non-nuc positions for complementary seq w/no
     blanks */
  ctg->compl_seq_no_blk_pos = NULL;
  for( cur_pos = ctg->seq_no_blk_pos; cur_pos ; cur_pos = cur_pos->next ){
    new_int = (POS_LIST *) malloc( sizeof( POS_LIST ) );
    new_int->val = ctg->seq_no_blk_lg - cur_pos->val - 1;
    new_int->next = ctg->compl_seq_no_blk_pos;
    ctg->compl_seq_no_blk_pos = new_int;
  }

  /* Modify the nb_nuc_bef field of all the annots for ctg 'cause it includes
     the spaces */
  for( ptr = ctg->annot_lst ; ptr ; ptr = ptr->next ){
    nb_spaces = 0;
    for( cur_pos = ctg->seq_pos ; cur_pos ; cur_pos = cur_pos->next ){
      if( cur_pos->val > ptr->nb_nuc_bef - 1 ) break;
      if( isspace( ctg->seq[ cur_pos->val ] ) ) nb_spaces++;
    }
    ptr->nb_nuc_bef -= nb_spaces;
  }

  /* Adjust the field nb_nuc_bef of annots on compl strand (right now, they're 
     forw AND they take spaces into account !) */
  for( ptr = ctg->compl_annot_lst ; ptr ; ptr = ptr->next ){
    nb_spaces = 0;
    for( cur_pos = ctg->compl_seq_pos ; cur_pos ; cur_pos = cur_pos->next ){
      if( cur_pos->val > ctg->seq_lg - ptr->nb_nuc_bef - 1 ) break;
      if( isspace( ctg->seq[ ctg->seq_lg - cur_pos->val - 1 ] ) ) nb_spaces++;
    }
    ptr->nb_nuc_bef = ctg->seq_lg - ptr->nb_nuc_bef - nb_spaces;
  }
}


/*-----------------------------------------------------------------------------

  Name: trans_all

  Description: Translates the current contig in all six reading frames 
               Gets the stops and starts
  
  Returns: Nothing
  
  ---------------------------------------------------------------------------*/

void trans_all(
	       CTG *ctg,     /* pointer on current contig struct. */
	       char code[],  /* genetic code used */
	       FRAME frame,  /* the six reading frames */
	       int idx[],    /* idx of all nucleotides */
	       PARAM *param  /* application's parameters */
)
{
  int       i;       /* loop counter */
  int       offset;  /* the offset associated to a given frame (see FRAME 
			typedef) */
  char     *seq;     /* a tmp sequence pointer */
  long int  aa_lg;   /* length of the aa seq associated to a frame */
  for( i=0; i < NB_FRAMES ; i++ )
  {
    offset = i < (NB_FRAMES/2) ? 
      i % CODON_LG : ( ctg->clean_seq_lg - (NB_FRAMES - i - 1) ) % CODON_LG;
    frame[i].offset = offset;
    seq = ( i >= NB_FRAMES/2 ) ? ctg->compl_clean_seq : ctg->clean_seq;
    aa_lg = 
      ctg->clean_seq_lg < offset ? 0 : (ctg->clean_seq_lg - offset) / CODON_LG;
    frame[i].aa_lg = aa_lg;
    if( !(frame[i].aa = (char *) calloc( aa_lg + 1, sizeof(char) ) ) ) 
      mem_err();     
    if( !(frame[i].stop = (long int*) calloc( aa_lg + 1, sizeof(long int) ) ) )
      mem_err();     
    if( !(frame[i].start = (long int*) calloc( aa_lg + 1, sizeof(long int) ) ))
      mem_err();     
    translate( &(frame[i]), seq, offset, ctg->clean_seq_lg - 1, idx, param, 
	       code );
  }

  /* This will find all the wrap codons that yield a stop in a frame */
  if( param->circular ) find_wrap_stops( frame, ctg, idx, code, param );
}


/*-----------------------------------------------------------------------------

  Name: find_wrap_stops

  Description: Examines all the wrap-around codons to see if there is a stop or
               a start that crosses the split point. If so, the pos is put in 
	       the frame where the stop/start begins. Note that this pos refers
	       to an aa that DOESN'T EXIST in the frame. It's always gonna be
	       aa_lg, where aa_lg is the current number of aa in the frame
  
  Returns: Nothing
  
  ---------------------------------------------------------------------------*/

void find_wrap_stops( 
		     FRAME frame,  /* the six reading frame */
		     CTG *ctg,     /* pointer on the current contig struct */
		     int idx[],    /* index of all nucs */
		     char code[],  /* genetic code used for translation */
		     PARAM *param  /* application's parameters */
)
{
  char  *seq;        /* tmp sequence pointer */
  CODON  codon;      /* a codon */
  int    f;          /* a frame id (0..5) */
  int    nb_ending;  /* the number of trailing nuc in the cleaned sequence for
			a given reading frame (ex: cleaned = ACGTACGT, 
			nb_ending = 2 for fr. 0, 1 for fr. 1 and 0 for fr. 2 */
  int    code_idx;   /* a given codon index */

  /* there cannot be a wrap-around prot for this case, so why bother finding
     wrap-around codons */
  if( ctg->clean_seq_lg < CODON_LG ) return;

  /* Scan every frame */
  for( f=0; f<NB_FRAMES; f++ ){
    seq = ( f >= (NB_FRAMES/2) ) ? ctg->compl_clean_seq : ctg->clean_seq;
    nb_ending = ( ctg->clean_seq_lg - frame[f].offset ) % CODON_LG;

    /* If there are trailing nucs, there's a wrap-around codon */
    if( nb_ending ){
      strcpy( codon, seq + ctg->clean_seq_lg - nb_ending );
      strncat( codon, seq, CODON_LG - nb_ending ) ;
      codon[CODON_LG] = '\0';
      code_idx = get_codon_idx( codon, idx );

      /* Don't forget here that we store the pos of the AA in the frame that
         would yield a stop */
      if( code[ code_idx ] == '*' ){
	frame[f].stop[ frame[f].nb_stop++ ] = 
	  ( ctg->clean_seq_lg - frame[f].offset ) / CODON_LG;
      }

      /* Don't forget here that we store the pos of the AA in the frame that
         would yield a start */
      if( param->starts[ code_idx ] ){
	frame[f].start[ frame[f].nb_start++ ] = 
	  ( ctg->clean_seq_lg - frame[f].offset ) / CODON_LG;
      }
    } /* if */

  } /* for */
  
}


/*-----------------------------------------------------------------------------

  Name: translate

  Description: Translate a given nucleotide sequence according to a given 
               genetic code. The routine translates each codon. Sequence is
	       assumed to be uppercased and to contain only ACGTN. If 
	       seq_lg % 3 != 0, the trailing nucs are ignored.
  
  Returns: Nothing
  
  ---------------------------------------------------------------------------*/

void translate( 
	       FRAME f,        /* a given frame */
	       char *seq,      /* pointer on a given nuc sequence */
	       long int start, /* the first nuc of the region to translate 
				  (0.seq_lg-1) */
	       long int end,   /* the last nuc of the region to translate 
				  (0.seq_lg-1) */
	       int idx[],      /* index of all nucs */
	       PARAM *param,   /* application's parameters */
	       char code[]     /* the genetic code used */
)
{
 char     *last_ptr;   /* ptr on the last nuc of the last codon to translate */
 char     *cur_ptr;    /* ptr on first nuc of first codon to translate */
 char     *aa_ptr;     /* pointer on beginning of resulting aa sequence  */
 char      aa;         /* a tmp amino acid */
 int       code_idx;   /* the codon idx of a codon */
 long int  nb_stop;    /* the nb of stops found so far */
 long int  nb_start;   /* the nb of starts found so far */
 long int  pos;        /* a given aa position */

 /* init */
 aa_ptr = f->aa;
 nb_stop = 0;
 nb_start = 0;
 pos = 0;

 /* if there are at least CODON_LG nucs to translate */
 if( start + CODON_LG - 1 <= end ) {
   last_ptr = 
     seq + end - ( CODON_LG - 1 ) - ( ( end - start + 1 ) % CODON_LG );
   cur_ptr = seq + start;

   /* translate all */
   while( cur_ptr <= last_ptr ){ 
     code_idx = get_codon_idx( cur_ptr, idx );
     aa = code[ code_idx  ];    
     *aa_ptr++ = aa;

     /* store stop/start pos if current aa is stoip/start */
     if( aa == '*' ) f->stop[nb_stop++] = pos;
     if( param->starts[code_idx] ) f->start[nb_start++] = pos;

     cur_ptr += CODON_LG;
     pos++;
   }
 }

 /* update */
 *aa_ptr     = '\0';
 f->nb_start = nb_start;
 f->nb_stop  = nb_stop;
   
}

/*-----------------------------------------------------------------------------

  Name: init_arrays

  Description: performs some initialization of the basic arrays needed by
              flip
  
  Returns: Nothing
  
  ---------------------------------------------------------------------------*/

void init_arrays( 
		 int idx[],    /* Each nuc has an index, used to compute the
				  codon index of all codons. This is the array
				  of all nuc indexes */
		 int comp[],   /* Each nuc has a complement (A->T T->A, etc...)
				  This is the array of all complements
				  (including the complement of all OGMP spec.
				  chars) */
		 int is_nuc[], /* Array telling, foreach char, if it's a 
				  nucleotide (valid nucs are ACGTNXacgtnx) */
		 int is_aa[],  /* Array telling, foreach char, if it's a valid
				  amino acid (valid aas are:
				     - every letter (upper or lower) except:
				         bBjJoOqQsSuUxXzZ
                                     - *                                     */
                 int is_seq_char[] /* will tell if a char is a sequencing char
				    (ie aAcCgGtTnNxX#@!+- or any whitespace) */
)
{
  int i;   /* loop counter */

  /* Init is_seq_char array */
  for( i=0; i < NB_CHARS ; i++ )
    is_seq_char[i] = isspace(i);
  is_seq_char['a'] = is_seq_char['A'] = 1;
  is_seq_char['c'] = is_seq_char['C'] = 1;
  is_seq_char['g'] = is_seq_char['G'] = 1;
  is_seq_char['t'] = is_seq_char['T'] = 1;
  //Modified by David To 12th April 2005
  //Adding in "x" and "X"
  is_seq_char['n'] = is_seq_char['N'] = is_seq_char['x'] = is_seq_char['X'] = 1;
  is_seq_char['#'] = is_seq_char['@'] = is_seq_char['!'] = 1;
  is_seq_char['+'] = is_seq_char['-'] = 1;
  

 /* Init idx array */
 idx['A'] = 0;
 idx['C'] = 1;
 idx['G'] = 2;
 idx['T'] = 3;
 idx['N'] = 4;
 //Modified by David To 12th April 2005
 //Adding in "X"
 idx['X'] = 4;

 /* Init is_aa array */
 for( i=0; i < NB_CHARS ; i++ )
   is_aa[i] = isalpha(i);
 is_aa['b'] = is_aa['B'] = 0;
 is_aa['j'] = is_aa['J'] = 0;
 is_aa['o'] = is_aa['O'] = 0;
 is_aa['q'] = is_aa['Q'] = 0;
 is_aa['s'] = is_aa['S'] = 0;
 is_aa['u'] = is_aa['U'] = 0;
 is_aa['x'] = is_aa['X'] = 0;
 is_aa['z'] = is_aa['Z'] = 0;
 is_aa['*'] = 1;

 /* Init comp array */
  for( i = 0; i < NB_CHARS; i++ ){
    comp[i] = i;
  }
  
  comp['a'] = 't', comp['A'] = 'T';
  comp['t'] = 'a', comp['T'] = 'A';
  comp['c'] = 'g', comp['C'] = 'G';
  comp['g'] = 'c', comp['G'] = 'C';
  comp['+'] = '-';
  comp['-'] = '+'; 
  comp['@'] = '#'; 
  comp['#'] = '@';

  /* Init is_nuc array */ 
  for( i=0; i < NB_CHARS ; i++ ) is_nuc[i] = 0 ;
  is_nuc['a'] = is_nuc['A'] = 1;
  is_nuc['c'] = is_nuc['C'] = 1;
  is_nuc['g'] = is_nuc['G'] = 1;
  is_nuc['t'] = is_nuc['T'] = 1;
  is_nuc['n'] = is_nuc['N'] = 1;
  //Modified by David To 12th April 2005
  //Adding in "X"
  is_nuc['x'] = is_nuc['X'] = 1;
}


/*-----------------------------------------------------------------------------

  Name: compl_seq

  Description: Complements and reverses a given sequence. Sequence is assumed
               to composed of aAcCgGtTnNxX@#+- or any whitespace only
  
  Returns: The complemented reversed sequence
  
  ---------------------------------------------------------------------------*/

char *compl_seq( 
		char seq[],  /* a tmp nuc pointer */
		long int lg, /* the length of the sequence to complement */
		int comp[]   /* array of all the nuc  complements */
)
{
  char *compl_seq;  /* will contain the complemented reversed sequence */
  long int i;       /* loop counter */

  compl_seq = (char *) calloc( lg+1, sizeof( char ) );
  seq += lg-1;

  /* reverse and complement */
  for( i=0; i < lg ; i++, seq-- ){
    compl_seq[i] = comp[ *seq ];
  }
  compl_seq[i] = '\0';

  return( compl_seq );
}

/*-----------------------------------------------------------------------------

  Name: open_all_files

  Description: opens all the files needed by flip. Exits flip with an exit 
               status of 1 if anything goes wrong 
  
  Returns: nothing
  
  ---------------------------------------------------------------------------*/

void open_all_files( 
		    FILE **prot_lst_dna, /* prot.lst.dna file pointer */
		    FILE **prot_6rf,     /* prot.6rf file ptr */
		    FILE **seq_file,     /* input file ptr */
		    FILE **prot_lst,     /* prot.lst file ptr */
		    FILE **prot_src,     /* prot.src file ptr */
		    PARAM *param,        /* param struct pointer  */
		    FILE **nocompl,      /* nocompl file ptr */
		    FILE **compl,        /* compl file ptr */
		    char *seq_file_name  /* input file name */
)
{
 /* Open prot.lst.dna only if user wants the DNA seq. of each orf */
 if( param->get_dna ) 
   if( !( *prot_lst_dna = fopen( PROT_LST_DNA_NAME, "w+" ) ) ){
     fprintf(stderr, "File %s could not be opened. Aborting.\n", 
	     PROT_LST_DNA_NAME);
     exit(1);
   }

 /* Open prot.6rf only if user wants all six reading frames fully translated */
 if( param->trans_all ) 
   if( !( *prot_6rf = fopen( PROT_6RF_NAME, "w+" ) ) ){
     fprintf(stderr, "File %s could not be opened. Aborting.\n", 
	     PROT_6RF_NAME);
     exit(1);
   }

 /* open prot.lst. Exit if error */
 if( !( *prot_lst = fopen( PROT_LST_NAME, "w+" ) ) ){
   fprintf(stderr, "File %s could not be opened. Aborting.\n", PROT_LST_NAME);
   exit(1);
 }

 /* open prot.src. Exit if error. Print the leading header (filename) now */
 if( !( *prot_src = fopen( PROT_SRC_NAME, "w+" ) ) ){
   fprintf(stderr, "File %s could not be opened. Aborting.\n", PROT_SRC_NAME);
   exit(1);
 }
 else fprintf( *prot_src, ">%s\n", seq_file_name );

 /* open input file. Exit if error*/
 if( !( *seq_file = fopen( param->s_file, "r" ) ) ){
   fprintf(stderr, "File %s could not be opened. Aborting.\n", param->s_file );
   exit(1);
 }

 /* open nocompl. Exit if error*/
 if( !( *nocompl = fopen( NOCOMPL_NAME, "w+" ) ) ){
   fprintf(stderr, "File %s could not be opened. Aborting.\n", NOCOMPL_NAME );
   exit(1);
 }

 /* open compl. Exit if error*/
 if( !( *compl = fopen( COMPL_NAME, "w+" ) ) ){
   fprintf(stderr, "File %s could not be opened. Aborting.\n", COMPL_NAME );
   exit(1);
 }
}

/*-----------------------------------------------------------------------------

  Name: update_prot_6rf

  Description: updates the file prot.6rf by writing to it the fully translated
               six reading frames associated to the current contig 
  
  Returns: Nothing
  
  ---------------------------------------------------------------------------*/


void update_prot_6rf( 
		     FILE *prot_6rf, /* pointer on prot.6rf */
		     FRAME frame,    /* the six frames of the current ctg */
		     CTG *ctg        /* ptr on the current contig */
)
{
  char      chunk[NUC_PER_LINE+1]; /* a chunk (max 60 nucs) of nucleotides */
  long int  i;                     /* the number of nucs printed so far */
  long int  seq_lg;                /* length of the clenaed sequence */
  char     *seq;                   /* pointer on start of cleaned seq. */
  int       f;                     /* tmp frame id (0..5) */
  int       chunk_lg;              /* length of the current printed chunk */
  int       nb_spaces;             /* tmp for formatting */
  int       f_lg;                  /* tmp */
  char     *cur_aa[NB_FRAMES];     

  /* tmp for pretty formatting of the frames */
  int       spaces[3][3] = { { 1, 0, 2 }, { 2, 1, 0 }, { 0, 2, 1 } };

  fprintf( prot_6rf, "%s\n\n", ctg->name );

  i=0;
  seq = ctg->clean_seq;
  seq_lg = ctg->clean_seq_lg;
  for( f=0; f < NB_FRAMES; f++ ){
    cur_aa[f] =
      f >= (NB_FRAMES/2) ? frame[f].aa + frame[f].aa_lg - 1 : frame[f].aa;
  }

  /* Print all frames */
  while( *seq != '\0' ){
    chunk_lg = seq_lg > NUC_PER_LINE ? NUC_PER_LINE : seq_lg;

    /* forward frames */
    for( f=0 ; f < (NB_FRAMES/2) ; f++ ){
      fprintf( prot_6rf, "%s    F:%d  ", ctg->name, f+1 );
      nb_spaces = (f == 2) && (i==0) ? 3 : spaces[f][ i % CODON_LG ];
      fprintf( prot_6rf, "%*s", nb_spaces, "" );
      f_lg = ( chunk_lg - nb_spaces - 1 ) / CODON_LG + 1;
      while( f_lg-- && *cur_aa[f] ) 
	fprintf( prot_6rf, "%c  ", *cur_aa[f]++ );
      fprintf( prot_6rf, "\n" );
    }

    /* print nuc sequence */
    fprintf( prot_6rf, "%*ld  ", (int) strlen( ctg->name ) + 7, i+1 );
    strncpy( chunk, seq, chunk_lg );
    chunk[chunk_lg] = '\0'; 
    fprintf( prot_6rf, "%s\n", chunk );

    /* backward frames */
    for( f=NB_FRAMES-1 ; f >= (NB_FRAMES/2)  ; f-- ){
      fprintf( prot_6rf, "%s    F:%d  ", ctg->name, f+1 );
      nb_spaces = !i && ( f == NB_FRAMES/2 ) ? 3 : spaces[5-f][ i % CODON_LG ];
      fprintf( prot_6rf, "%*s", nb_spaces, "" );
      f_lg = ( chunk_lg - nb_spaces - 1 ) / CODON_LG + 1;
      while( f_lg-- && ( cur_aa[f] >= frame[f].aa ) )
	fprintf( prot_6rf, "%c  ", *cur_aa[f]-- );
      fprintf( prot_6rf, "\n" );
    }

    fprintf( prot_6rf, "\n\n\n" );

    seq    += chunk_lg;
    seq_lg -= chunk_lg;
    i      += chunk_lg;
  }
}

/*-----------------------------------------------------------------------------

  Name: close_all_files

  Description: closes all the files that flip used 
  
  Returns: Nothing
  
  ---------------------------------------------------------------------------*/

void close_all_files( 
		     FILE **prot_lst_dna, /* prot.lst.dna file pointer */
		     FILE **prot_6rf,     /* prot.6rf file ptr         */
		     FILE **seq_file,     /* input file ptr            */
		     FILE **prot_lst,     /* prot.lst file ptr         */
		     FILE **prot_src,     /* prot.src file ptr         */
		     FILE **nocompl,      /* nocompl file ptr          */
		     FILE **compl,        /* compl file ptr            */
		     PARAM *param         /* param struct pointer      */
)
{
  if( *prot_6rf!=NULL )      fclose( *prot_6rf );
  if( *prot_lst_dna!=NULL )  fclose( *prot_lst_dna );
  if( *seq_file!=NULL )  fclose( *seq_file );
  if( *prot_lst!=NULL )  fclose( *prot_lst );
  if( *prot_src!=NULL )  fclose( *prot_src );
  if( *nocompl!=NULL )   fclose( *nocompl );
  if( *compl!=NULL )     fclose( *compl );

  if( !param->silent ){
    printf( "Wrote %s\n", NOCOMPL_NAME);
    printf( "Wrote %s\n", COMPL_NAME);
    printf( "Wrote %s\n", PROT_LST_NAME);
    printf( "Wrote %s\n", PROT_SRC_NAME);
    if( param->trans_all ) printf( "Wrote %s\n", PROT_6RF_NAME );
    if( param->get_dna ) printf( "Wrote %s\n", PROT_LST_DNA_NAME );
  }
}

/*-----------------------------------------------------------------------------

  Name: find_prot

  Description: finds all the proteins in all the frames, including the 
               wrap-around proteins and saves them .
  
  Returns: Nothing
  
  ---------------------------------------------------------------------------*/

void find_prot( 
	       CTG *ctg,     /* ptr on current contig */
	       PROT **prot,  /* pointer on the protein list */
	       PARAM *param, /* ptr on the application's params */
	       FRAME frame   /* all six reading frames */
)
{
  PROT     *last_prot;       /* pointer on the last found prot in the list */
  int       f;               /* a frame_id (0..5) */
  int       s;               /* index of a stop in the list of stops */
  long int  r_start, r_end;  /* two consecutive stop positions (r_end might 
				point to something that's not a stop for the 
				last region) */
  long int  r_lg;            /* length of a region (between two consecutive
				stops, or aa after last stop) */
				
  long int  prot_start;
  long int  p_lg;
  long int  cur_s_idx;

  last_prot = *prot = NULL;  

  /* scan all frames  */
  for( f=0; f<NB_FRAMES; f++ ){
    cur_s_idx = 0;
    /* scan all regions */
    for( s=0; s <= frame[f].nb_stop; s++ ){
      
      /* If genome is circular, first and last regions are skipped */
      if( ( !s || (s == frame[f].nb_stop) ) && param->circular ) continue;
      
      r_start = s ? frame[f].stop[s-1] + 1 : 0;
      r_end = s == frame[f].nb_stop ? frame[f].aa_lg - 1 : frame[f].stop[s];
      r_lg = s == frame[f].nb_stop ? r_end - r_start + 1 : r_end - r_start;

      /* does the orf respect the min_codong_reg_lg  ? */
      if( r_lg >= param->min_coding_reg_lg ){

        /* if there is a minimum length (from start codon) the orf must have */
	if( param->min_prot_lg ){
	  prot_start = get_p_start( frame+f, r_start, r_end, &cur_s_idx );
	  if( prot_start != -1 ){
	    p_lg = s == frame[f].nb_stop ?
	      r_end - prot_start + 1 : r_end - prot_start; 
	    if( p_lg >= param->min_prot_lg )
	      add_prot( prot, &last_prot, r_start, r_end, prot_start,
			frame + f, param, ctg, f, p_lg );
	  }
	}
	else add_prot( prot, &last_prot, r_start, r_end, -1, frame + f, 
		       param, ctg, f, r_lg );
      }
    } /* for */
  }

  if( param->circular ) find_wrap_prot( prot, &last_prot, ctg, param, frame );
}

/*-----------------------------------------------------------------------------

  Name: find_wrap_prot

  Description: Assuming the genome is circular, finds the proteins that span 
               across the split point (these don't include proteins for which
	       only the last stop codon crosses the split poitn, which are 
	       reported by find_prot)
             
  Returns: Nothing
  
  ---------------------------------------------------------------------------*/

void find_wrap_prot( 
		    PROT **prot,      /* pointer on the ptr on the first prot
					 in the list of found prot */
		    PROT **last_prot, /* pointer on the ptr on the last prot
					 in the list of found prot */
		    CTG *ctg,         /* pointer on the current ctg struct */
		    PARAM *param,     /* pointer on the application's params
					 structure */
		    FRAME frame       /* the siz reading frames for the ctg */
)
{
 int       f;                   /* a frame_id (0..NB_FRAMES-1) */
 int       cur_f;               /* the current frame index (0..NB_FRAMES-1) */
 int       nb_ending;           /* the number of trailing nucs that are ignored
				   during linear translation of a seq. for a
				   given frame (ex, if seq=ACTG, then nb_ending
				   is 1 for frame 0) */
 int       f_idx;               /* a given frame index (0..NB_FRAMES-1) */
 int       next_f;              /* frame index of the frame following a given 
				   one during circular translation */
 int       order[NB_FRAMES/2];  /* this array gives the order in which the 
				   frame indexes are to be read, starting from
				   a given frame index. Ex: seq=TAAC and the 
				   current frame is 0. order=(2,1,0,-1) (the -1
				   is to signal end of reading) */
 int       nb_pass;             /* nb of times we entered a loop */
 long int  r_lg;                /* length of the wrap-around orf */
 long int  prot_lg;             /* length of the wrap-around orf, counted from 
				   the first start codon */
 int       start_o_idx;         /* the index, in the order array, of the frame
				   where the first start codon is found */
 long int  start_pos;           /* the position of the start aa in the frame 
				   refered to by start_o_idx */

  /* Clearly, if clean_seq_lg is < 3, there cannot be ANY wrap-around prots */
  if( ctg->clean_seq_lg < CODON_LG ) return;
  
  /* scan all frames */
  for( f=0; f < NB_FRAMES ; f++ ){

    /* if no stop, then a wrap-around prot cannot start in the current frame */
    if( !frame[f].nb_stop ) continue;

    /* find next frame (wrap-around-wise) that has a stop */  
    cur_f  = f;  
    for( f_idx=0; f_idx < NB_FRAMES ; f_idx++ ) order[f_idx] = -1;
    nb_pass = 0;
    while( 1 ){
      nb_ending = ( ctg->clean_seq_lg - frame[cur_f].offset ) % CODON_LG;

      /* get frame following current one (wrap-around-wise) */
      if( !nb_ending ) nb_ending = CODON_LG;
      next_f = -1;
      for( f_idx=0; f_idx < NB_FRAMES ; f_idx++ )
	if( SAME_STRAND( f_idx, cur_f ) )  
	  if( frame[f_idx].offset == ( CODON_LG - nb_ending ) ){
	    next_f = f_idx;
	    break;
	  }

      if( next_f == -1 )
	fprintf( stderr, "Impossible, no next frame !\n" ), exit(1);

      order[nb_pass++] = next_f;

      /* we might have a wrap-around prot here */
      if( frame[ next_f].nb_stop ){
/*        printf( "Found wrap-around region starting in frame %d and ", f );
	printf( "ending in frame %d\n", next_f );  
	printf( "Order = " );  
        for( f_idx=0; f_idx < (NB_FRAMES/2) ; f_idx++ )    
	  printf( "%d ", order[f_idx] );   
	printf( "\n" ); */

        r_lg = get_wrap_rg_lg( ctg, frame, order, f );
/*	printf( "The region length is %ld\n", r_lg ); */
	if( r_lg >= param->min_coding_reg_lg )
	  if( param->min_prot_lg ){
/*	    printf( "Seeking a wrap-around prot\n" );*/
	    prot_lg = get_wrap_prot_lg( ctg, frame, order, f, next_f, r_lg, 
					&start_o_idx, &start_pos );
	    if( prot_lg >= param->min_prot_lg ){
/*	      printf( "A valid prot was found\n" );*/
/*	      printf( "The calculated protein length is %ld\n", prot_lg );*/
	      add_wrap_prot( prot, last_prot, ctg, param, frame, order, f, 
			     next_f, start_o_idx, start_pos, r_lg, prot_lg );
	    }
	  }
	  else 
	    add_wrap_prot( prot, last_prot, ctg, param, frame, order, f, 
			   next_f, start_o_idx, start_pos, r_lg, prot_lg );  
 	break;
      }
      else cur_f = next_f;
    }

  }

}

/*-----------------------------------------------------------------------------

  Name: add_wrap_prot

  Description: Custom routine to add a wrap-around protein in the list of all
               proteins found.
  
  Returns: Nothing
  
  ---------------------------------------------------------------------------*/

void add_wrap_prot( 
		   PROT **prot,        /* ptr on the ptr on the first prot */
		   PROT **last_prot,   /* ptr on the ptr on the last prot */
		   CTG *ctg,           /* ptr on the current ctg struct */
		   PARAM *param,       /* ptr on the param struct */
		   FRAME frame,        /* the six reading frames */
		   int order[],        /* the order in which the frames are to
					  be read (circular-wise, see 
					  find_wrap_prot) */
		   int f,              /* the frame in which the orf starts 
					  (ie leading stop) */
		   int last_f,         /* the frame in which the orf ends (ie
					  last stop) */
		   int start_o_idx,    /* index, in order array, of the frame 
					  where the first start was found */
		   long int start_pos, /* pos, in the frame referred to by 
					  start_o_idx, where the first stop aa
					  was found */
		   long int r_lg,      /* length of the wrap-around orf  */
		   long int prot_lg    /* length of the wrap-aournd orf, 
					  counted from the first start codon */
)
{
  PROT     *new_prot; /* a new protein structure */
  int       on_rev;   /* 1 if protein is on reverse strand, 0 otherwise */
  long int  length;   /* the length of the prot (might be from the first start
			 codon if param->min_prot_lg > 0). Always includes last
			 stop. */

  if( !( new_prot = (PROT *) malloc( sizeof(PROT) ) ) ) mem_err();

  /* get prot length (including last stop) */
  length = param->min_prot_lg ? prot_lg + 1 : r_lg + 1;

  on_rev = f >= (NB_FRAMES/2);
  new_prot->ctg_name = ctg->name;
  new_prot->aa_lg    = length;
  new_prot->length   = length-1;
  new_prot->strand   = on_rev;


  /* Get the aa composing the protein */
  if( !( new_prot->aa = (char *) calloc( length + 1, sizeof( char ) ) ) ) 
    mem_err();

  /* get info associated to wrap-around prot */
  get_wrap_aa( new_prot, param, frame, f, ctg, code, idx, order, r_lg, 
	       start_o_idx, start_pos );
  get_wrap_pos( new_prot, ctg, frame, f, last_f, start_o_idx, start_pos, 
		param, order );
  if( param->min_prot_lg ) 
    get_wrap_comment( new_prot, frame, f, ctg, code, idx, order, start_o_idx, 
		      start_pos, r_lg - prot_lg );
  get_wrap_context( new_prot, ctg );

  /* Get orf's DNA sequence if needed */
  if( param->get_dna ) get_dna_seq( new_prot, ctg );
  else new_prot->dna = NULL;

  /* Add in list. Update pointers to start and end of list */
  if( *prot ){
    (*last_prot)->next = new_prot;
    *last_prot = new_prot;
  }
  else *last_prot = *prot = new_prot;
  new_prot->next = NULL;
}

/*-----------------------------------------------------------------------------

  Name: get_wrap_context

  Description: custom routine to get the context associated to a wrap-around 
               protein (see get_context for a complete description of the 
	       context). Note that for wrap-around prots, there are ALWAYS 
	       "..." in the context  

  Returns: Nothing
  
  ---------------------------------------------------------------------------*/

void get_wrap_context( 
		      PROT *prot, /* ptr on wrap-around prot */
		      CTG *ctg    /* ptr on current contig */
)
{
 long int  clean_seq_start = prot->begin - 1; /* first nuc in cleaned seq. 
						 assoc to orf */
 long int  clean_seq_end   = prot->end - 1;   /* last nuc in cleaned seq. assoc
						 to orf */
 long int  i;                                 /* loop counter */
 long int  seq_start;                         /* the position (in the original
						 seq) of the nuc whose pos is
						 clean_seq_start in the cleaned
						 seq. */
 long int  seq_end;                           /* the position (in the original
						 seq) of the nuc whose pos is
						 clean_seq_end in the cleaned
						 seq. */
 long int  chunk_lg;                          /* a chunk of context */
 long int  old_lg;                            /* tmp */
 POS_LIST *pos_ptr;                           /* ptr on the list of non-nuc pos
						 either in orig or rev compl. 
						 seq */

 /* get start in orig sequence */
 i = 0;
 for( pos_ptr = ctg->seq_pos; pos_ptr ; pos_ptr = pos_ptr->next, i++ ){
   if( pos_ptr->val > i + clean_seq_start ) break;
 }
 seq_start = clean_seq_start + i;

 /* get end in original sequence */
 i = 0;
 for( pos_ptr = ctg->seq_pos; pos_ptr ; pos_ptr = pos_ptr->next, i++ ){
   if( pos_ptr->val > i + clean_seq_end ) break;
 }
 seq_end = clean_seq_end + i;
 
 /* Get C-term and N-term */
 prot->context = (char *) calloc( 2*MAX_CTX_LG + 4, sizeof( char ) );
 
 /* Get C-term (or N-term for backw prot) */
 if( !prot->strand ){
   strncpy( prot->context, ctg->seq + seq_start, MAX_CTX_LG );
   prot->context[MAX_CTX_LG] = '\0';
 }
 else{
   chunk_lg = MIN( MAX_CTX_LG, seq_start + 1 );
   strncpy( prot->context, ctg->seq + seq_start - chunk_lg + 1, chunk_lg );
   prot->context[chunk_lg] = '\0';
 }

 /* Put separator between the two chunks */
 strcat( prot->context, "..." );
 old_lg = strlen( prot->context );

 /* N-term (or C-term for back) */
 if( !prot->strand ){
   chunk_lg = MIN( MAX_CTX_LG, seq_end + 1 );
   strncat( prot->context, ctg->seq + seq_end - chunk_lg + 1, chunk_lg );
   prot->context[old_lg + chunk_lg] = '\0';
 }
 else{ 
   strncat( prot->context, ctg->seq + seq_end, MAX_CTX_LG );
   prot->context[ old_lg + MAX_CTX_LG ] = '\0';
 }
 
}

/*-----------------------------------------------------------------------------

  Name: get_wrap_comment

  Description: Gets the aa upstream from the first start codon in a wrap-around
               orf
  
  Returns: Nothing
  
  ---------------------------------------------------------------------------*/

void get_wrap_comment(
	    PROT *new_prot,     /* ptr on the new wrap-around prot */
            FRAME frame,        /* six reading frames */
	    int f,              /* the frame in which the orf starts (ie 
				   leading stop) */
	    CTG *ctg,           /* ptr on the current ctg struct */
	    char code[],        /* genetic code to use during translation */
	    int idx[],          /* index of each nuc */
	    int order[],        /* the order in which the frames are to be read
				   (circular-wise, see find_wrap_prot) */
	    int start_o_idx,    /* index, in order array, of the frame where 
				   the first start was found */
	    long int start_pos, /* pos, in the frame referred to by 
				   start_o_idx, where the first stop aa was
				   found */
	    long int size       /* length of the wrap-around orf, counted from 
				   the first start codon */
)
{
  CODON     codon;      /* a codon that lies across the split point */
  int       i;          /* loop counter */
  long int  lg;         /* tmp */
  long int  old_lg;     /* tmp */
  long int  last_stop;  /* the pos of the lst stop in f */
  int       cur_f;      /* a frame index (0.NB_FRAMES-1) */
  int       nb_ending;  /* the number of trailing nucs that are ignored during
			   linear translation of a seq. for a given frame (ex,
			   if seq=ACTG, then nb_ending is 1 for frame 0) */
  char     *seq_ptr;    /* a tmp sequence ptr */

  /* init */
  new_prot->comment = (char *) calloc( size + 1, sizeof(char) );
  new_prot->comment[0] = '\0';

  /* get all comment chunks in all the frames */
  for( i = -1 ;; i++ ){

    /* the first frame, f, is not in order, so it's a special case */
    cur_f = i == -1 ? f : order[i];

    /* that's the end of the context */
    if( i == start_o_idx ){

      /* if the start fr. is also the frame where orf starts, get portion from
	 (excluding) last stop of frame to (excluding) first start after this 
	 stop */
      if( cur_f == f ){ 
	last_stop = frame[f].stop[ frame[f].nb_stop -1 ];
        old_lg = strlen( new_prot->comment );
	lg = start_pos - last_stop - 1;
	strncat( new_prot->comment, frame[f].aa + last_stop + 1, lg );
	new_prot->comment[old_lg + lg] = '\0';
      }
      /* get everything up to (but excluding) the first start in the current
	 frame */
      else{
	old_lg = strlen( new_prot->comment );
	strncat( new_prot->comment, frame[cur_f].aa, start_pos );
	new_prot->comment[old_lg + start_pos] = '\0';
      }

      break;
    }
    /* part of the chunk, but not the last part */
    else{
   
      /* if we are processing the frame where orf starts, get portion from
	 (excluding) last stop of frame to end of aa in frame */
       if( cur_f == f ){ 
	last_stop = frame[f].stop[ frame[f].nb_stop -1 ];
	if( last_stop < frame[f].aa_lg - 1 )
	  strcpy( new_prot->comment, frame[f].aa + last_stop + 1 );
       }
       /* if it's not the first chunk, get the whole translated frame */
       else strcat( new_prot->comment, frame[cur_f].aa );
 
       /* this will add in the comment the wrap-around codon (if any) */
       nb_ending = ctg->clean_seq_lg - 
	 ( CODON_LG * frame[cur_f].aa_lg + frame[cur_f].offset );
       /* watch out ! If there's a wrap-around, make sure it's not the stop
	  that's initiating the orf */
       if( nb_ending && !( (cur_f == f) && (last_stop >= frame[f].aa_lg - 1) ))
	 {
	   seq_ptr = new_prot->strand ? ctg->compl_clean_seq : ctg->clean_seq;
	   strcpy( codon, seq_ptr + ctg->clean_seq_lg - nb_ending );
	   strncat( codon, seq_ptr, CODON_LG - nb_ending );
	   old_lg = strlen( new_prot->comment );
	   new_prot->comment[old_lg] = code[ get_codon_idx( codon, idx ) ];
	   new_prot->comment[old_lg+1] = '\0';
	 }

    } /* else */
  } /* for */
  
  new_prot->comment_lg = strlen( new_prot->comment );
}

/*-----------------------------------------------------------------------------

  Name: get_wrap_aa

  Description: gets the amino acid sequence associated to a wrap-around protein
  
  Returns: Nothing
  
  ---------------------------------------------------------------------------*/

void get_wrap_aa( 
       PROT *new_prot,     /* ptr on the prot for which we want aa sequence */
       PARAM *param,       /* ptr on the application's parameter structure */
       FRAME frame,        /* the six reading frames */
       int f,              /* the frame idx (0..5) where the orf's initiating
			      stop was found */
       CTG *ctg,           /* ptr on current contig struct */
       char code[],        /* genetic code used during translation */
       int idx[],          /* index of each nuc */
       int order[],        /* this array gives the order in which the frame
			      indexes are to be read, starting from a given
			      frame index. Ex: seq=TAAC and the current frame
			      is 0. order=(2,1,0,-1) (the -1 is to signal end
			      of reading) */
       long int r_lg,      /* length of the orf (from first to last stop 
			      (excluding the stops themselves) */
       int start_o_idx,    /* the index, in the order array, of the frame where
			      the first start codon is found */
       long int start_pos  /* the position of the start aa in the frame 
			      referred to by start_o_idx */
)
{
  long int  last_stop;    /* position of the last stop of frame f */
  int       nb_ending;    /* the number of trailing nucs that are ignored
			     during linear translation of a seq. for a given
			     frame (ex, if seq=ACTG, then nb_ending is 1 for
			     frame 0) */
  char     *seq_ptr;      /* a tmp seq ptr */
  long int  old_lg;       /* a tmp length */
  int       i;            /* loop counter */
  int       cur_f;        /* a frame index (0..5) */
  int       on_rev;       /* 1 if prot is on backw strand, 0 otherwise */
  CODON     codon;        /* a tmp wrap-around codon */
  char     *aa_ptr;       /* a tmp aa seq ptr */
  long int  lg;           /* position of the last stop in the frame in which 
			     orf ends */

  new_prot->aa[0] = '\0';
  on_rev = f >= ( NB_FRAMES/2 );  

  /* if min_prot_lg was specified */
  if( param->min_prot_lg )
    for( i = start_o_idx ; i < (NB_FRAMES/2); i++ ){
      cur_f = i == -1 ? f : order[i];
      if( cur_f == -1 ) break;

      /* last chunk...*/
      if( order[i+1] == -1 ){
	old_lg = strlen( new_prot->aa );
        aa_ptr = i == start_o_idx ? 
	  frame[cur_f].aa + start_pos : frame[cur_f].aa;

	/* if first aa is not a stop, get aa before it */
	if( frame[cur_f].stop[0] ){
	  lg = frame[cur_f].stop[0];
	  if( i == start_o_idx ) lg -= start_pos;
	  strncat( new_prot->aa, aa_ptr, lg );
	}
	else lg = 0;

        new_prot->aa[old_lg + lg ] = '*';
        new_prot->aa[old_lg + lg + 1 ] = '\0';
      }
      /* it's not the last chunk */
      else
	/* first chunk */
	if( i == start_o_idx ){
          /* if the start is not the last wrap-around codon, get it and all the
	     aas following it in the same frame */
          if( start_pos < frame[cur_f].aa_lg ){
	    aa_ptr = frame[cur_f].aa +  start_pos;
	    strcpy( new_prot->aa, aa_ptr );
	  }

	  /* Add the wrap-around codon (if any) */
	  nb_ending = ctg->clean_seq_lg - 
	    ( CODON_LG * frame[cur_f].aa_lg + frame[cur_f].offset );
	  if( nb_ending ){
	    seq_ptr = on_rev ? ctg->compl_clean_seq : ctg->clean_seq;
	    strcpy( codon, seq_ptr + ctg->clean_seq_lg - nb_ending );
	    strncat( codon, seq_ptr, CODON_LG - nb_ending );
	    old_lg = strlen( new_prot->aa );
	    new_prot->aa[old_lg] = code[ get_codon_idx( codon, idx ) ];
	    new_prot->aa[old_lg+1] = '\0';
	  }
	} 

        /* other chunks (ie not first and not last): we get the whole frame */
	else{
	  strcat( new_prot->aa, frame[cur_f].aa );

 	  /* Add wrap-around codon (if any) */
	  nb_ending = ctg->clean_seq_lg - 
	    ( CODON_LG * frame[cur_f].aa_lg + frame[cur_f].offset );	  
	  if( nb_ending ){
	    seq_ptr = on_rev ? ctg->compl_clean_seq : ctg->clean_seq;
	    strcpy( codon, seq_ptr + ctg->clean_seq_lg - nb_ending );
	    strncat( codon, seq_ptr, CODON_LG - nb_ending );
	    old_lg = strlen( new_prot->aa );
	    new_prot->aa[old_lg] = code[ get_codon_idx( codon, idx ) ];
	    new_prot->aa[old_lg+1] = '\0';
	  }

	} /* else */
    } /* for */

  /* If min_prot_lg was not specified */
  else{
    /* if the last stop in f is not a wrap-around codon, get all aas following
       it (if any) */
    if( frame[f].stop[ frame[f].nb_stop - 1 ] < frame[f].aa_lg ){ 
      last_stop = frame[f].nb_stop - 1;
      strcat( new_prot->aa, frame[f].aa + frame[f].stop[ last_stop ] + 1 );

      /* Add wrap-around codon (if any) */
      nb_ending = 
	ctg->clean_seq_lg - ( CODON_LG * frame[f].aa_lg + frame[f].offset );
      if( nb_ending ){
	seq_ptr = on_rev ? ctg->compl_clean_seq : ctg->clean_seq;
	strcpy( codon, seq_ptr + ctg->clean_seq_lg - nb_ending );
	strncat( codon, seq_ptr, CODON_LG - nb_ending );
        old_lg = strlen( new_prot->aa );
        new_prot->aa[old_lg] = code[ get_codon_idx( codon, idx ) ];
	new_prot->aa[old_lg+1] = '\0';
      }

    }

    /* get remaining aa chunks */
    for( i=0; i < (NB_FRAMES/2); i++ ){
      cur_f = order[i];
      if( cur_f == -1 ) break;
      /* if there's a stop in the current frame, then it's the last chunk. Get
	 everyhting up to this stop (excluding)*/
      if( frame[cur_f].nb_stop ){
	strncat( new_prot->aa, frame[cur_f].aa, frame[cur_f].stop[0] );
           
        /* add OB stop */
	new_prot->aa[r_lg] = '*';
	new_prot->aa[r_lg+1] = '\0';
      }
      /* no stop in current frame, then we take all the aas in the frame */
      else{
	strcat( new_prot->aa, frame[cur_f].aa );

	/* Add wrap-around codon (if any) */
	nb_ending = ctg->clean_seq_lg - 
	  ( CODON_LG * frame[cur_f].aa_lg + frame[cur_f].offset );
	if( nb_ending ){
	  seq_ptr = on_rev ? ctg->compl_clean_seq : ctg->clean_seq;
	  strcpy( codon, seq_ptr + ctg->clean_seq_lg - nb_ending );
	  strncat( codon, seq_ptr, CODON_LG - nb_ending );
	  old_lg = strlen( new_prot->aa );
	  new_prot->aa[old_lg] = code[ get_codon_idx( codon, idx  ) ];
	  new_prot->aa[old_lg+1] = '\0';
	}

      } /* else */

    } /* for */  
  } /* else */
 
}  /* get_wrap_aa */

/*-----------------------------------------------------------------------------

  Name: get wrap_region_lg

  Description: Gets the length (in ass) of a wrap-around orf. This length 
               doesn't include the stops
  
  Returns: the computed length
  
  ---------------------------------------------------------------------------*/

long int get_wrap_rg_lg(
      CTG *ctg,     /* ptr on the current contig */
      FRAME frame,  /* the six reading frames */
      int order[],  /* this array gives the order in which the frame indexes 
		       are to be read, starting from a given frame index. Ex:
		       seq=TAAC and the current frame is 0. order=(2,1,0,-1)
		       (the -1 is to signal end of reading) */
      int f        /* the frame in which the orf starts (ie leading stop) */
)
{
  long int r_lg;   /* will be the sought after region length */
  int      o_idx;  /* loop counter */
  int      cur_f;  /* frme index (0..5) */

  /* Get length of chunk of region on first frame */
  r_lg = frame[f].aa_lg - frame[f].stop[ frame[f].nb_stop - 1 ] - 1;

  /* Add the wrap-around codon, if any */
  if( ( ctg->clean_seq_lg - frame[f].offset ) % 3 ) r_lg++;

  /* This is to take into account the wrap-around stops whose pos is >= aa_lg,
     and so might yield negative region length */
  if( r_lg < 0 ) r_lg = 0;

  /* Add the length of the other chunks on the other frames */
  for( o_idx=0 ; o_idx < (NB_FRAMES/2) ; o_idx++ ){
    cur_f = order[o_idx];
    if( cur_f == -1 ) break;

     /* this indicates the last chunk (ie from first aa to (excluding) the 1st
	stop found */
    if( frame[cur_f].nb_stop ){
      r_lg += frame[cur_f].stop[0];
    }
    /* get everything */
    else{
      r_lg += frame[cur_f].aa_lg;
      
      /* Add the wrap-around codon if there is one */
      if( ( ctg->clean_seq_lg - frame[cur_f].offset ) % CODON_LG ) r_lg++;
    }
  }

  return( r_lg );
}

/*-----------------------------------------------------------------------------

  Name: get_wrap_prot_lg

  Description: Gets the length (counted from the first start codon) of an orf
               that is assumed to be wrapped around the split point.
               The stops are never part of this length computation, but the 
	       start are. Also, sets start_o_idx to be the index, in the order
               array of the frame where the first start was found. Also sets 
               start_pos to be the position of the start (in the frame referred
               to by start_o_idx).
  
  Returns: The computed length if the orf has at least one start codon, -1
           otherwise
  
  ---------------------------------------------------------------------------*/

long int get_wrap_prot_lg( 
      CTG *ctg,            /* ptr on the contig struct */
      FRAME frame,         /* the six reading frames */
      int order[],         /* this array gives the order in which the frame
			      indexes are to be read, starting from a given
			      frame index. Ex: seq=TAAC and the current frame
			      is 0. order=(2,1,0,-1) (the -1 is to signal end
			      of reading) */
      int begin_f,         /* the frame in which the orf starts (ie leading
			      stop) */
      int end_f,           /* the frame in which the orf ends (ie last stop)*/
      long int r_lg,       /* whole length of the orf (between the two stops)*/
      int *start_o_idx,    /* index, in order array, of the frame where the
			      first start will be found. If the 1st start lies
			      after the last stop of f and before the end of f,
			      then start_o_idx will be set to -1 */
      long int *start_pos  /* pos, in the frame referred to by start_o_idx, 
			      where the first stop aa was found */
)
{
  long int last_stop_pos;  /* position of the last stop in frame f */
  long int last_start_idx; /* position of the last start in f */
  long int pos;            /* a tmp start pos */
  int      i;              /* tmp loop counter */
  int      cur_f;          /* a frame index (0..5) */
  int      o_idx;          /* loop counter */
  int      tmp_idx;        /* tmp order index */
  long int prot_lg;        /* sought length */
  long int stop_idx;       /* an index of a stop in the stop array */

  last_stop_pos = frame[begin_f].stop[ frame[begin_f].nb_stop - 1 ];

  /* Check to see if there's a start that's after this stop in the same frm */
  if( frame[begin_f].nb_start ){
    last_start_idx = frame[begin_f].nb_start - 1;
    
    *start_pos = -1;
    for( i = last_start_idx ; i >= 0 ; i-- ){
      pos = frame[begin_f].start[i];
      if( pos < last_stop_pos ) break;
      *start_pos = pos;
    }

    /* return the whole orf length minus the portion between the first stop
       and the found start */
    if( *start_pos != -1 ){
      *start_o_idx = -1;
      return( r_lg - ( *start_pos - last_stop_pos - 1 ) );
    }
  }

  /* seek a start */
  *start_pos = -1;
  for( o_idx = 0 ; o_idx < (NB_FRAMES/2) ; o_idx++ ){
    cur_f = order[o_idx];
    if( cur_f == -1 ) break; 

    /* we found a start */
    if( frame[cur_f].nb_start )

      /* if it's in the same frame as the one in ehich the orf ends, make sure
	 that this start is before the first stop */
      if( cur_f == end_f ){
	if( frame[cur_f].start[0] < frame[end_f].stop[0] ){
	  *start_o_idx = o_idx;
	  *start_pos = frame[cur_f].start[0];
	  break;
	}
      }
    /* we got what we wanted */
      else{
	*start_o_idx = o_idx;
	*start_pos = frame[cur_f].start[0];
	break;
      }
  }  
  
  /* if a stop was found */
  if( *start_pos != -1 ){

    /* Get the first chunk (might also be the last) */
    /* if it's in the frame where the orf ends */
    if( cur_f == end_f ){ 
      if( frame[cur_f].nb_stop ){
        stop_idx = -1;
	for( i = 0 ; i < frame[cur_f].nb_stop; i++ ){
	  if( frame[cur_f].stop[i] > *start_pos ){
	    stop_idx = i; break;
	  } 
	}
         /* if there's a stop after this start, return the length of the chunk
	    between (including) this start and the stop (excluding) */
        if( stop_idx == -1 ) 
	  prot_lg = frame[cur_f].aa_lg - *start_pos;
        /* else first chunk = from start to end of frame */
	else
	  prot_lg = frame[cur_f].stop[stop_idx] - *start_pos;
      }
      /* else first chunk = from start to end of frame */
      else
	prot_lg = frame[cur_f].aa_lg - *start_pos;
    }
    /* if the start is not in the frame where the orf ends */
    else{
      prot_lg = frame[cur_f].aa_lg - *start_pos;
      
      /* Add the wrap-around aa, if any */
      if( ( ctg->clean_seq_lg - frame[cur_f].offset ) % 3 ) prot_lg++;
    }

    /* get length of remaining chunks */
    for( tmp_idx = o_idx+1; tmp_idx < (NB_FRAMES/2) ; tmp_idx++ ){
      cur_f = order[tmp_idx];
      if( cur_f == -1 ) break;

      /* that's the last chunk */
      if( frame[cur_f].nb_stop ){
	prot_lg += frame[cur_f].stop[0];
      }
      /* not last and not first, so get every aa */
      else{
	prot_lg += frame[cur_f].aa_lg;
	
	/* Add the wrap-around codon if there is one */
	if( ( ctg->clean_seq_lg - frame[cur_f].offset ) % CODON_LG ) prot_lg++;
      }
    }

    return( prot_lg );
  }
  /* No start codon found... */
  else return( -1 );
}  

/*-----------------------------------------------------------------------------

  Name: get_p_start

  Description: gets the position of the first aa which is a start aa in the
               orf
  
  Returns: The position of the first start aa if there's a start aa in the orf
           -1 otherwise
  
  ---------------------------------------------------------------------------*/

long int get_p_start(
        FRAME f,             /* a given frame struct */
	long int r_start,    /* the position of the first aa of the orf (not
				the leading stop!) */
	long int r_end,      /* the position of the first aa of the orf (if orf
				ends w/stop, then it's the stop pos) */
	long int *cur_s_idx  /* index of the last start, in the starts array, 
				that was examined to find the previous 
				protein's start */
)
{
  long int  pos;        /* a tmp start pos */
  long int *start_ptr;  /* ptr on the starts array assoc. to a frame */

  /* find the first start in [r_start, r_end[  */
  start_ptr = f->start; 
  for( ; (*cur_s_idx) < f->nb_start ; (*cur_s_idx)++ ){
    pos = start_ptr[*cur_s_idx];
    if( pos > r_end ) return( -1 );
    if( pos >= r_start ) return( pos );
  }
  /* if no start, return -1 */
  return( -1 );
}

/*-----------------------------------------------------------------------------

  Name: add_prot

  Description: Adds an orf/protein in list of current detected orfs/proteins
  
  Returns: Nothing
  
  ---------------------------------------------------------------------------*/

void add_prot(
    PROT **prot,       /* ptr on the ptr on the beginning of the list */
    PROT **last_prot,  /* ptr on the ptr on the end of the list */
    long int r_start,  /* the position of the first aa of the orf (not the
			  leading stop!) */
    long int r_end,    /* the position of the first aa of the orf (not the
			  leading stop!) */
    long int p_start,  /* position of the first start in the orf (if any) */
    FRAME f,           /* frame where the orf was found */
    PARAM *param,      /* ptr on the application's parameters structure */
    CTG *ctg,          /* ptr on the current contig structure */
    int f_nb,          /* frame index of the frame where the orf was found */
    long int l         /* length of the orf (counted from the first start codon
			  if param->min_prot_lg ) */
)
{
  PROT     *new_prot;    /* a new prot struct */
  long int  prot_begin;  /* will be r_start if param->min_prot_lg and p_start
			    otherwise (ie actual orf's start) */
  long int  p_lg;        /* actual prot_lg (why this ? we laready have l...) */
  int       on_rev;      /* 1 if prot is on the reverse strand, 0 otherwise */
  long int  c_lg;        /* length of comment to put in prot.lst (if any) */

  prot_begin = param->min_prot_lg ? p_start : r_start;
  p_lg = r_end - prot_begin + 1;

  if( !( new_prot = (PROT *) malloc( sizeof(PROT) ) ) ) mem_err();

  if( !( new_prot->aa = (char *) calloc( p_lg + 2, sizeof( char ) ) ) ) 
    mem_err();
  on_rev = f_nb >= (NB_FRAMES/2);

  new_prot->ctg_name = ctg->name;
  new_prot->aa_lg    = p_lg;
  new_prot->length   = l;
  new_prot->strand   = on_rev;
  strncpy( new_prot->aa, f->aa + prot_begin, p_lg );

  /* For prots that are linear but the last stop crosses the split-point */
  if( prot_begin + p_lg >= f->aa_lg )
    strcat( new_prot->aa, "*" );

  new_prot->aa[p_lg] = '\0'; 
  new_prot->begin    = CODON_LG * prot_begin + f->offset + 1;
  new_prot->end      = CODON_LG * r_end      + CODON_LG + f->offset;

  /* For prots that are linear but the last stop crosses the split-point */ 
  if( new_prot->end > ctg->clean_seq_lg ) 
    new_prot->end = new_prot->end - ctg->clean_seq_lg;

  if( on_rev ){
    new_prot->begin = ctg->clean_seq_lg - new_prot->begin + 1;
    new_prot->end   = ctg->clean_seq_lg - new_prot->end + 1;
  }

  /* Get context info */
  get_context( new_prot, ctg );
  
  /* Get DNA seq of AA seq, if needed. */
  if( param->get_dna ) get_dna_seq( new_prot, ctg );
  else new_prot->dna = NULL; 

  /* Get comment: all aa upstream from the start codon. If no minimum length
     specified for proteins, then no comment */
  if( param->min_prot_lg ){
    c_lg = p_start - r_start;
    if( !( new_prot->comment = (char *) calloc( c_lg + 1, sizeof( char ) ) ) ) 
      mem_err();
    strncpy( new_prot->comment, f->aa + r_start, c_lg );
    new_prot->comment[c_lg] = '\0';
    new_prot->comment_lg = c_lg;
  }
  else{
    new_prot->comment = NULL;
    new_prot->comment_lg = 0;
  }

  /* Add in list. Update pointers to start and end of list */
  if( *prot ){
    (*last_prot)->next = new_prot;
    *last_prot = new_prot;
  }
  else *last_prot = *prot = new_prot;
  new_prot->next = NULL;
}

/*-----------------------------------------------------------------------------

  Name: get_wrap_pos

  Description: gets the positions, in the cleaned sequence, of the first and
                last nucleotides associated to a wrap-around protein. Put them
		in the prot structure 
  
  Returns: nothing
  
  ---------------------------------------------------------------------------*/

void get_wrap_pos(
      PROT *new_prot,     /* pointer on the new prot */
      CTG *ctg,           /* pointer on the current contig struct */
      FRAME frame,        /* the six reading frames */
      int f,              /* frame index of frame in which orf begins (ie 
			     leading stop)  */
      int last_f,         /* frame index of frame in which orf begins (ie 
			     trailing stop)  */
      int start_o_idx,    /* index, in the order array, of the frame in which 
			     the first start was found (if any) (might be -1
			     if first start is in f) */
      long int start_pos, /* pos of the first start (in frame pointed to by 
			     start_o_idx) */
      PARAM *param,       /* ptr on the application's param structure */
      int order[]         /* this array gives the order in which the frame
			      indexes are to be read, starting from a given
			      frame index. Ex: seq=TAAC and the current frame
			      is 0. order=(2,1,0,-1) (the -1 is to signal end
			      of reading) */
)
{
  int       start_f;       /* index of the frame where first start was found */
  long int  last_stop_idx; /* index of the last stop (in stop array) of 
			      last_f */
  int       nb_ending;


  /* Get orf begin pos */
  /* if param->min_prot_lg, pos is pos of found start */
  if( param->min_prot_lg ){
    start_f = start_o_idx == -1 ? f : order[start_o_idx];
    new_prot->begin = CODON_LG * start_pos + frame[start_f].offset + 1;
  }
  /* else pos is pos of first nucleotide after leading stop */
  else{
    last_stop_idx = frame[f].nb_stop - 1;
    new_prot->begin = CODON_LG * ( frame[f].stop[last_stop_idx] + 1 ) + 
      frame[f].offset + 1;

    /* the "%" if for when the leading stop is a wrap-around stop */
    if( new_prot->begin > ctg->clean_seq_lg )
      new_prot->begin = new_prot->begin % ctg->clean_seq_lg;
  }

  /* Get the orf's end pos */
  /* if the first stop is a wrap-around stop */
  if( frame[last_f].stop[0] >= frame[last_f].aa_lg ){
    nb_ending = ctg->clean_seq_lg - 
      ( frame[last_f].offset + CODON_LG * frame[last_f].aa_lg ); 
    new_prot->end = CODON_LG - nb_ending;
  }
  /* else if the nuc in clean seq assoc. w/ last nuc of stop */
  else 
    new_prot->end = 
      CODON_LG * ( frame[last_f].stop[0] + 1 ) + frame[last_f].offset;

  /* if on compl strand, reverse and complement positions */
  if( new_prot->strand ){
    new_prot->begin = ctg->clean_seq_lg - new_prot->begin + 1;
    new_prot->end    = ctg->clean_seq_lg - new_prot->end   + 1;
  }
}

/*-----------------------------------------------------------------------------

  Name: update_prot_lst_src

  Description: puts all orfs found in the current contig in either prot.lst
               and prot.src
  
  Returns: Nothing
  
  ---------------------------------------------------------------------------*/

void update_prot_lst_src(
     PROT **sorted_prot,      /* array of sorted proteins */
     int    nb_prot,          /* number of elements in the sorted array */
     FILE  *prot_lst,         /* the prot.lst file pointer */
     FILE  *prot_src,         /* the prot.src file pointer */
     int   number_prot        /* boolean. Indicates if we should append a
				 "_<protein_nb>" to indicate the protein number
				 in the output files */
)
{
  char      chunk[AA_PER_LINE+1]; /* a chunk of aa sequence */
  long int  i;                    /* tmp loop counter */
  int       chunk_lg;             /* the length of a chunk of aa sequence */
  PROT     *prot;                 /* pointer on a protein */
  int       prot_idx;

  /* if there were no proteins found, leave */
  if( !nb_prot ) return;

  for( prot_idx = 0; prot_idx < nb_prot; prot_idx++ ){
    prot = sorted_prot[prot_idx];
 
    /* Print prot headers in both files */
    fprintf( prot_src, ";" );
    
    if( !number_prot ){
      fprintf( prot_lst, "%s; ", prot->ctg_name );
      fprintf( prot_src, "%s; ", prot->ctg_name );
    }
    else{
      fprintf( prot_lst, "%s_%d; ", prot->ctg_name, prot_idx+1 );
      fprintf( prot_src, "%s_%d; ", prot->ctg_name, prot_idx+1 );      
    }

    fprintf( prot_lst, "%s ", prot->strand ? "compl." : "orig." );
    fprintf( prot_src, "%s ", prot->strand ? "compl." : "orig." );

    fprintf( prot_lst, "%6ld to %6ld ;", prot->begin, prot->end );
    fprintf( prot_src, "%6ld to %6ld ;", prot->begin, prot->end );

    fprintf( prot_lst, " %s\n", prot->context );
    fprintf( prot_src, " %s\n", prot->context );
    
    /* Print comment line if any */
    if( prot->comment_lg )
      for( i = 0; i < prot->comment_lg ; i += AA_PER_LINE ){
	strncpy( chunk, prot->comment + i, AA_PER_LINE );
	chunk[AA_PER_LINE] = '\0';
	fprintf( prot_lst, ";%s\n", chunk );
	fprintf( prot_src, ";%s\n", chunk );
      }
    
    /* Print aa sequence formatted to AA_PER_LINE amino acids per line */
    /* In prot.lst, the ending stop If any) is printed.  */
    for( i = 1; i <= prot->aa_lg ; i += AA_PER_LINE ){
      strncpy( chunk, prot->aa + i - 1, AA_PER_LINE );
      chunk[AA_PER_LINE] = '\0';
      fprintf( prot_lst, "%6ld  %s\n", i, chunk );
    }
    
    /* Print aa sequence formatted to AA_PER_LINE amino acids per line */
    /* In prot.src, the ending stop (if any) is NOT printed */
    for( i = 1; i <= prot->length; i += AA_PER_LINE ){
      chunk_lg = MIN( AA_PER_LINE, prot->length - i + 1 );
      strncpy( chunk, prot->aa + i - 1, chunk_lg );
      chunk[ chunk_lg ] = '\0';
      fprintf( prot_src, "%s\n", chunk );
    }

    /* Print aa sequence length (only in prot.lst) */
    fprintf( prot_lst, "%6ld\n", prot->length );

    fprintf( prot_lst, "\n" );
    fprintf( prot_src, "\n" );
  }
}

static void  FreePOSList( POS_LIST *pos )
{
  POS_LIST   *p, *next;

  p = pos;
  while( p != NULL )
  {
    next = p->next;
    free( p );
    p = next;
  }

  return;
}

static void  FreeAnnotList( ANNOT_LIST *annot )
{
  ANNOT_LIST   *p, *next;

  p = annot;
  while( p != NULL )
  {
    next = p->next;
    free( p );
    p = next;
  }

  return;
}

/*-----------------------------------------------------------------------------

  Name: free_mem

  Description: frees the memory occupied by the prot list, the current contig
               structure. Frees the aa sequences associated to each frame along
	       with their stop and start list 
  
  Returns: Nothing
  
  ---------------------------------------------------------------------------*/


void free_mem( 
    CTG   *ctg,        /* ptr on the current ctg structure */
    PROT  *prot,       /* ptr on the current protein list */
    FRAME  frame,      /* all six reading frames */
    PROT **sorted_prot /* array of sorted proteins */
)
{
  PROT       *p, *next; /* tmp pos list pointers */
  int         f;        /* a frame index (0...5) */

  /* Free space occupied by current contig */
 
  if( ctg->name != NULL )             free( ctg->name );
  if( ctg->seq_no_blk != NULL )       free( ctg->seq_no_blk );
  if( ctg->clean_seq != NULL )        free( ctg->clean_seq );
  if( ctg->compl_clean_seq != NULL )  free( ctg->compl_clean_seq );
  if( ctg->compl_seq_no_blk != NULL ) free( ctg->compl_seq_no_blk );
  FreePOSList( ctg->seq_pos );
  FreePOSList( ctg->seq_no_blk_pos );
  FreePOSList(  ctg->compl_seq_pos );
  FreePOSList( ctg->compl_seq_no_blk_pos );
  FreeAnnotList( ctg->annot_lst );
  FreeAnnotList( ctg->compl_annot_lst );
  memset( ctg, 0, sizeof(CTG) );

  /* Protein link list */

  p = prot;
  while( p != NULL )
  {
    if( p->aa != NULL )  free( p->aa );
    
    /* If we had the orfs DNA seq. free the string */
    
    if( p->dna != NULL )      free( p->dna );
    if( p->comment != NULL )  free( p->comment );
    if( p->context != NULL )  free( p->context );

    next = p->next;
    free( p );
    p = next;
  }

  /* Free space occupied by all frame's aa sequences, stop list and start 
     list */
  
  for( f=0 ; f < NB_FRAMES ; f++ )
  {
    if( frame[f].aa != NULL )    free( frame[f].aa );
    if( frame[f].stop != NULL )  free( frame[f].stop );
    if( frame[f].start != NULL ) free( frame[f].start );
  }

  /* Free space occupied by current array of sorted proteins */
  
  if( sorted_prot != NULL ) free( sorted_prot );
}

/*-----------------------------------------------------------------------------

  Name: print_pos

  Description: debug routine to print various informations on all frames 
  
  Returns: Nothing
  
  ---------------------------------------------------------------------------*/


void print_pos( 
   FRAME frame, /* all the frames */
   char *name   /* name of the current contig */
)
{
  int f, i;  /* tmps */

  for( f=0 ; f < NB_FRAMES ; f++ ){
    printf( "Contig %s frame %d:\n", name, f+1 );
    printf( "\t - AA length is %ld\n", strlen( frame[f].aa ) );
    printf( "\t - Offset is %d\n", frame[f].offset );
    printf( "\t - Number of stops is %ld", frame[f].nb_stop );
    printf( "\t - Stop pos: " ); 
    for( i=0; i< frame[f].nb_stop ; i++ ) printf( "%ld,", frame[f].stop[i] );
    printf( "\n" );
    printf( "\t - Number of starts is %ld", frame[f].nb_start );
    printf( "\t - Start pos: " );
    for( i=0; i< frame[f].nb_start ; i++ ) printf( "%ld,", frame[f].start[i] );
    printf( "\n\n" );
  }
}

/*-----------------------------------------------------------------------------

  Name: add_annot

  Description: adds an annotation in the list of annotations on the forw and
               backw strand. While adding the annot in the list of annots on 
	       the backw strand, the arrow is reversed. The routine supports
	       annots that end with a backslash (indicating an annot that 
	       continues on the next line) and will insert them in both lists
	       taking this into account. 

  Returns: nothing
  
  ---------------------------------------------------------------------------*/

void add_annot( 
    ANNOT_LIST **lst_start,       /* ptr on ptr of beginning of annot list on
				     forw strand */		
    ANNOT_LIST **lst_end,         /* ptr on ptr of end of annot list on forw
				     strand */
    ANNOT_LIST **compl_lst_start, /* ptr on ptr of beginning of annot list on
				     backw strand */
    ANNOT_LIST **compl_lst_end,   /* ptr on ptr of end of annot list on backw
				     strand */
    char *annot,                  /* the annotation's text */
    long int nb_nuc,              /* the nb of nucleotides (does not include
				     spec. OGMP chars or spaces) that lie
				     before the annot on the forw. strand */
    int  *prec_cont               /* ptr on a var which = 1 if the last line */
)
{
  ANNOT_LIST *new_annot;        /* a new annot on the forw strand */
  ANNOT_LIST *compl_new_annot;  /* a new annot on the backw strand */
  char       *dir_start;        /* tmp */
  char       *tmp;              /* tmp */
  int         nb_stars;         /* the number of stars before/after the arrow 
				   in an annotation */

  /* Get new annotation info: text, number of nucleotide preceding the annot */
  new_annot             = (ANNOT_LIST *) malloc( sizeof(ANNOT_LIST) ); 
  new_annot->nb_nuc_bef = nb_nuc;
  new_annot->text = (char *) calloc( strlen(annot) + 1, sizeof(char) ); 
  strcpy( new_annot->text, annot );
  new_annot->next = NULL;
  
  /* Add new annot in list of annot on forw strand */
  if( *lst_start ) (*lst_end)->next = new_annot;
  else *lst_start = new_annot;  
  *lst_end = new_annot;

  /* Get new annotation info: text, nb of nucleotide preceeding the annot */
  compl_new_annot = (ANNOT_LIST *) malloc( sizeof(ANNOT_LIST) ); 

  /* The field nb_nuc_bef will be set as if the annot was forward (to be 
     changed later when we know the length of seq_no_blk */
  compl_new_annot->nb_nuc_bef = nb_nuc;

  /* Get compl annotation text */
  compl_new_annot->text = (char *) calloc( strlen(annot) + 1, sizeof(char) );
  strcpy( compl_new_annot->text, annot );

  /* Reverse direction of arrow if not a comment */
  if( strlen(annot) > 2 && annot[1] != ';' )
/*    if( (dir_start = strstr( compl_new_annot->text, "==>" )) ){*/
    if( (dir_start = find_arrow_ptr( compl_new_annot->text, "==>", "<==" )) ){
      nb_stars = 0;
      while( (dir_start > compl_new_annot->text) && *(dir_start-1) == '*' )
	nb_stars++, dir_start--;
      strncpy( dir_start, "<==", DIR_LENGTH );
      tmp = dir_start + DIR_LENGTH - 1;
      while( nb_stars ) *(tmp + nb_stars) = '*', nb_stars--;
    }
    else
      if( (dir_start = find_arrow_ptr( compl_new_annot->text, "<==", "==>" )) )
	{
	  nb_stars = 0;
	  tmp = dir_start + DIR_LENGTH;
	  while( *(tmp + nb_stars) == '*' ) nb_stars++;
	  tmp = dir_start + nb_stars;
	  while( nb_stars ) *(dir_start + nb_stars - 1) = '*', nb_stars--;
	  strncpy( tmp, "==>", DIR_LENGTH );
	}
  
  /* Add new annot in list of annots on back strand */
  if( *compl_lst_end )
    if( ends_with_bksl( (*compl_lst_end)->text ) ){
      compl_new_annot->next  = (*compl_lst_end)->next;
      (*compl_lst_end)->next = compl_new_annot;
    }
    else{
      compl_new_annot->next = *compl_lst_start;
      *compl_lst_start = compl_new_annot;
    }
  else {
    *compl_lst_start = compl_new_annot;
    compl_new_annot->next = NULL;
  }

  *prec_cont = ends_with_bksl( annot );

  *compl_lst_end = compl_new_annot;
}

/*-----------------------------------------------------------------------------

  Name: update_seq_file

  Description: Updates either nocompl or compl w/the information we have on the
               current contig
  
  Returns: Nothing
  
  ---------------------------------------------------------------------------*/

void update_seq_file(  
    CTG  *ctg,              /* ptr on current contig struct */
    FILE *nocompl,          /* pointer on the file to update (not necess.
			       nocompl) */
    int   is_compl,         /* 1 if the file to update is compl, 0 otherwise */
    int   nb_ctg,           /* Nb of contigs (empty or not) written to compl
			       and nocompl so far */
    int   nb_non_mepty_ctg, /* Nb of non-empty (ie at least one nuc or spec 
			       OGMP char) contig written to compl and nocompl
			       so far */
    int   fasta_format      /* Boolean. Whether to separate the contigs by '*'
			       or not */
)
{
  char       *seq_ptr;                           /* tmp seq ptr */
  ANNOT_LIST *cur_annot;                         /* tmp annot list ptr */
  POS_LIST   *cur_non_nuc = ctg->seq_pos;        /* tmp ptr on beginning of 
						    list of non-nuc pos on orig
						    strand */
  long int    nb_char     = 0;                   /* nb of chars (including spec
						    OGMP ch.) printed so far*/
  long int    nb_nuc      = 0;                   /* nb of nucs printed */
  long int    seq_lg      = ctg->seq_no_blk_lg;  /* lg of orig seq. w/blanks
						    stripped */
  char        chunk[NUC_PER_LINE+1];             /* a chunk of nuc sequence */
  int         chunk_lg;                          /* length of this chunk */
  int         cur_nb_nuc;                        /* nb of nuc in cur chunk*/
  long int    tmp;                               /* tmp */

  seq_ptr     = is_compl ? ctg->compl_seq_no_blk     : ctg->seq_no_blk;
  cur_annot   = is_compl ? ctg->compl_annot_lst      : ctg->annot_lst;
  cur_non_nuc = is_compl ? ctg->compl_seq_no_blk_pos : ctg->seq_no_blk_pos;

  if( nb_ctg ) 
    fprintf( nocompl, "%s\n\n", 
	     (fasta_format && nb_non_empty_ctg) ? "*": "" );
      

  /* Print ctg name, which might be empty (for header, which is considered as
     a contig containing only comments */
  if( ctg->name ) fprintf( nocompl, "%s\n", ctg->name );

  /* Print all annots that are before the first nuc printed. Set cur_annot to
     point after the last annot printed */
  for( ; cur_annot && !(cur_annot->nb_nuc_bef) ; cur_annot = cur_annot->next )
    fprintf( nocompl, "%s", cur_annot->text );

  /* If there are no seq. char to print, leave */
  if( !seq_ptr ) return ;

  /* Print all nuc w/intermixed annots */
  while( *seq_ptr ){
    if( cur_annot ){
      tmp = cur_annot->nb_nuc_bef - nb_char;
      chunk_lg = NUC_PER_LINE > tmp ? tmp : NUC_PER_LINE;
    }
    else{
      tmp = seq_lg - nb_char;
      chunk_lg = NUC_PER_LINE > tmp ? tmp : NUC_PER_LINE;
    }
    strncpy( chunk, seq_ptr, chunk_lg );
    chunk[chunk_lg] = '\0';

    /* Print seq chars. with nuc count in front of them */
    cur_nb_nuc =  
      chunk_lg - get_nb_spec( &cur_non_nuc, nb_char, nb_char + chunk_lg - 1 );
    fprintf( nocompl, "%6ld  %s\n", nb_nuc+1, chunk );

    /* Print all annots that lie just after the last seq. char. printed */
    while( cur_annot && ( cur_annot->nb_nuc_bef == nb_char + chunk_lg ) ){
      fprintf( nocompl, "%s", cur_annot->text );
      cur_annot = cur_annot->next;
    }

    /* Update counters */
    seq_ptr += chunk_lg;
    nb_nuc  += cur_nb_nuc;
    nb_char += chunk_lg;
  }

  /* Print annots after last seq. char */
  for( ; cur_annot && !(cur_annot->nb_nuc_bef) ; cur_annot = cur_annot->next )
    fprintf( nocompl, "%s", cur_annot->text );

  /* Print total contig length (in nb of nucleotide) */
  if( nb_nuc ) { fprintf( nocompl, ";;%6ld\n", nb_nuc); }
}

/*-----------------------------------------------------------------------------

  Name: get_nb_spec

  Description: gets the number of special (OGMP) characters in a given sequence
  
  Returns: the number of spec chars found
  
  ---------------------------------------------------------------------------*/

int get_nb_spec( POS_LIST **pos_ptr, long int start, long int end )
{
  int nb = 0;

  while( *pos_ptr && ((*pos_ptr)->val >= start) && ((*pos_ptr)->val <= end) ){
    nb++;
    *pos_ptr = (*pos_ptr)->next;
  }
  
  return( nb );
}

/*-----------------------------------------------------------------------------

  Name: ends_with_backslash

  Description: Verifies if a given string's last non-whitespace character is
               a backslash
  
  Returns: 1 if the last non-whitespace character of the string is '\', 0 
           otherwise
  
  ---------------------------------------------------------------------------*/

int ends_with_bksl( 
	   char *s   /* a string */
)
{
  char *cur_char;

  if( !strlen(s) ) return 0;
  for( cur_char = s + strlen(s) - 1 ; cur_char >= s ; cur_char-- ){
    if( *cur_char == '\\' ) return (1);
    if( !isspace( *cur_char ) ) break;
  }
  return(0);
}

/*-----------------------------------------------------------------------------

  Name: get_context

  Description: Gets the context of a linear orf. The context returned always
               refers to chunks of chars on the forw original strand.
	       Suppose that orf starts at nuc x and ends at nuc y (if orf on 
	       forw strand y > x, otherwise x > y (except for orfs for which
	       the last stop crosses the split point).
          
               -if the orf doesn't end with a stop that crosses the split pnt:
                 - If abs(x-y) + 1 > 24:
	            - if the prot is on forw strand:
		         return: (x,x+11)...(y-11,y)
			      (where (x,x+11) means chars (including spec. OGMP
			       chars and spaces) x thru x+11)
                    - if the prot is on back strand:
                         return: (x-11,x)..(y,y+11)

                 - If abs(x-y) + 1 <= 24
                    - if orf on forw. strand return (x,y)
		    - if orf on backw strand return (y,x)

               -if orf ends with a stop that crosses the split point:
	           the same proc is applied but:

	           -if the orf is on the forw. strand, it's taken to end at the
		    last char of the sequence.
	             
	           -if the orf is on the back. strand, y is set to be the first
		    char of the sequence.


  Returns: Nothing
  
  ---------------------------------------------------------------------------*/

void get_context( 
	PROT *prot,  /* ptr on the current protein */
	CTG *ctg     /* ptr on the current contig */
)
{
 long int  clean_seq_start = prot->begin - 1; /* start of prot in cleaned seq 
						 string */
 long int  clean_seq_end;                     /* end of prot in cleaned seq 
						 string */
 long int  seq_start;                         /* the position (in the original
						 seq) of the nuc whose pos is
						 clean_seq_start in the cleaned
						 seq. */
 long int  seq_end;                           /* the position (in the original
						 seq) of the nuc whose pos is
						 clean_seq_end in the cleaned
						 seq. */
 long int  i;                                 /* loop counter */
 long int  chunk_lg, old_lg;                  /* tmp */
 POS_LIST *pos_ptr;                           /* tmp pos list ptr */
 long int  prot_lg;   /* length of the nuc sequence (including spec. OGMP chars
                         and spaces) associated to prot in orig seq. */

 /* This is for prots for which the last stop crosses the split point */
 /* The last nuc of this stop will not be given in the context */
 if( prot->strand )
   if( prot->begin < prot->end ) clean_seq_end = 0;
   else clean_seq_end = prot->end - 1;
 else
   if( prot->begin > prot->end ) clean_seq_end = ctg->clean_seq_lg;
   else clean_seq_end = prot->end - 1;

 /* get start of sequence in orig forw sequence */
 i = 0;
 for( pos_ptr = ctg->seq_pos; pos_ptr ; pos_ptr = pos_ptr->next, i++ ){
   if( pos_ptr->val > i + clean_seq_start ) break;
 }
 seq_start = clean_seq_start + i;

 /* get start of sequence in orig back sequence */
 i = 0;
 for( pos_ptr = ctg->seq_pos; pos_ptr ; pos_ptr = pos_ptr->next, i++ ){
   if( pos_ptr->val > i + clean_seq_end ) break;
 }
 seq_end = clean_seq_end + i;

 prot_lg = nuc_lg( prot, ctg, seq_start, seq_end );

 
 if( prot_lg >= 2*MAX_CTX_LG ){

   /* Get C-term and N-term */
   prot->context = (char *) calloc( 2*MAX_CTX_LG + 4, sizeof( char ) );
   
   /* C-term (or N-term for backw) */
   if( !prot->strand ){
     strncpy( prot->context, ctg->seq + seq_start, MAX_CTX_LG );
     prot->context[MAX_CTX_LG] = '\0';
   }
   else{
     chunk_lg = MIN( MAX_CTX_LG, seq_start + 1 );
     strncpy( prot->context, ctg->seq + seq_start - chunk_lg + 1, chunk_lg );
     prot->context[chunk_lg] = '\0';
   }
   
   /* Separator */
   strcat( prot->context, "..." );
   old_lg = strlen( prot->context );
   
   /* N-term (or C-term for back) */
   if( !prot->strand ){
     chunk_lg = MIN( MAX_CTX_LG, seq_end + 1 );
     strncat( prot->context, ctg->seq + seq_end - chunk_lg + 1, chunk_lg );
     prot->context[old_lg + chunk_lg] = '\0';
   }
   else{ 
     strncat( prot->context, ctg->seq + seq_end, MAX_CTX_LG );
     prot->context[ old_lg + MAX_CTX_LG ] = '\0';
   }
 }
 else{
   /* Get whole nuc. seq. (including spec. chars & spaces ) assoc. to prot */
   prot->context = (char *) calloc( prot_lg + 1, sizeof( char ) );
   strncpy( prot->context, ctg->seq + MIN( seq_start, seq_end ), prot_lg );
   prot->context[ prot_lg ] = '\0';
 }
}


/*-----------------------------------------------------------------------------

  Name: parse_cmd_line

  Description: parses the command line and fills the param structure. If there
               were no switches on the cmd line, tries to read prot.prm. 
	       All unset params are set to the defaults. If min_coding_reg_lg
	       is unset and min_prot_lg is set, then we set min_coding_reg_lg
	       to min_prot_lg
  
  Returns: Nothing
  
  ---------------------------------------------------------------------------*/

void parse_cmd_line(
      int    argc,        /* nb of cmd line args */
      char  *argv[],      /* all cmd line args */
      PARAM *param,       /* param structure */
      int    is_nuc[],    /* array telling if a char is a nuc */
      int    is_aa[],     /* array telling if a char is a valid aa */
      char   all_codes[NB_CODES][NB_CODONS], /* all NCBI's genetic codes */
      int    all_starts[NB_CODES][NB_CODONS], /* all starts for every genetic 
					      code */
      char **code,                         /* will contain the current code to 
					      use */
      int idx[]                            /* index of each nuc */
)
{
  int   i;                    /* loop counter */
  char *end;                  /* tmp */
  char *arg_ptr;              /* ptr on an arg */
  char *starts_string = NULL; /* the starts string passed on the command line*/
  char *dev_string    = NULL; /* the dev string passed on the command line*/
  int   has_switches = 0;     /* if there were any switches or not */

  /* Init param structure to it's default values */
  param->code_id           = DEF_CODE_ID;
  param->circular          = DEF_CIRCULAR;
  param->min_coding_reg_lg = -1;
  param->min_prot_lg       = -1;
  param->trans_all         = DEF_TRANS_ALL;
  param->s_file            = NULL;
  param->starts            = NULL;
  param->force_met         = DEF_FORCE_MET;
  param->silent            = DEF_SILENT;
  param->number_prot       = DEF_NUMBER_PROT;
  param->get_dna           = DEF_GET_DNA;
  param->fasta_format      = DEF_FASTA_FORMAT;

  for( i = 1 ; i < argc-1 ; i++ ){
    if( argv[i][0] == '-' ){
      if( strlen( argv[i] ) <= 1 ) usage();
      switch( argv[i][1] )
        {
	  
        /* Tell that if first codon is a start, translate it as 'M' */
	case 'm':
	  if( strlen( argv[i] ) > 2 ) usage();
	  else param->force_met = 1;
	  has_switches = 1;
	  break;

	/* Set the genome to be circular */
	case 'c': 
	  if( strlen( argv[i] ) > 2 ) usage();
	  else param->circular = 1;
	  has_switches = 1;
	  break;

	/* Set the amount of messages to display */
	case 'S': 
	  if( strlen( argv[i] ) > 2 ) usage();
	  else param->silent = 1;
	  break;

	/* n = we want proteins to be "numbered" in the output files */
	case 'n':  
	  if( strlen( argv[i] ) > 2 ) usage();
	  else param->number_prot = 1;
	  break;

	/* -D = we want the file prot.lst.dna to be produced */
	case 'D':  
	  if( strlen( argv[i] ) > 2 ) usage();
	  else param->get_dna = 1; 
	  break;

	/* -F = separate each contig of compl and nocompl by '*' */
	case 'F':  
	  if( strlen( argv[i] ) > 2 ) usage();
	  else param->fasta_format = 1; 
	  break;

        /* User wants to produce prot.6rf (ie translate all frames) */
	case 'f': 
	  if( strlen( argv[i] ) > 2 ) usage();
	  else param->trans_all = 1;
	  has_switches = 1;
	  break;

        /* User wants to set the minimum coding region length */
	case 'l':
	  if( strlen( argv[i] ) <= 2 )
	    if( i < argc-1 ) arg_ptr = argv[++i];
	    else usage();
	  else arg_ptr = argv[i] + 2;

	  param->min_coding_reg_lg = strtol( arg_ptr, &end, 10 );
	  if( *end != '\0' ) usage();
	  if( param->min_coding_reg_lg <= 0 ) usage();
	  else{
	    has_switches = 1;
	    break;
	  }

        /* User wants to set the minimum protein length */
	case 'L': 
	  if( strlen( argv[i] ) <= 2 )
	    if( i < argc-1 ) arg_ptr = argv[++i];
	    else usage();
	  else arg_ptr = argv[i] + 2;

	  param->min_prot_lg = strtol( arg_ptr, &end, 10 );
	  if( *end != '\0' ) usage();
	  if( param->min_prot_lg <= 0 ) usage();
	  else{
	    has_switches = 1;
	    break;
	  }

        /* Set the genetic code to be used during translation */
        case 'g':
	  if( strlen( argv[i] ) <= 2 )  
	    if( i < argc-1 ) arg_ptr = argv[++i];
	    else usage();
	  else arg_ptr = argv[i] + 2;

	  param->code_id = strtol( arg_ptr, &end, 10 );
	  if( *end != '\0' ) usage();

	  if( param->code_id <= 0 ) usage();
	  if( (param->code_id > NB_CODES) || (param->code_id == 7) || 
	      (param->code_id == 8) ) usage();
	  else{
	    break;
	  }

        /* Redefine the start codons */
	case 's':
	  if( strlen( argv[i] ) <= 2 )  
	    if( i < argc-1 ) arg_ptr = argv[++i];
	    else usage();
	  else arg_ptr = argv[i] + 2;

	  starts_string = (char *) calloc( strlen(arg_ptr) + 1, sizeof(char) );
          strcpy(starts_string, arg_ptr);
	  has_switches = 1;
	  break;

        /* Impose deviations to the genetic code used */
	case 'd':
	  if( strlen( argv[i] ) <= 2 )  
	    if( i < argc-1 ) arg_ptr = argv[++i];
	    else usage();
	  else arg_ptr = argv[i] + 2;

          dev_string = (char *) calloc( strlen(arg_ptr) + 1, sizeof(char) );
	  strcpy( dev_string, arg_ptr );
	  has_switches = 1;
	  break;
 
	default: usage();
	}
    }
    else break;
  }
 
  if( i != argc-1 ) usage();

  /* if no switches were found, read prot.prm (if it exists ) */
  if( !has_switches ) parse_prot_prm( param, &dev_string, &starts_string );
  param->s_file = argv[argc-1];

  *code = all_codes[ param->code_id - 1 ];
  param->starts = all_starts[ param->code_id - 1 ];

  if( dev_string )
    if( get_deviations( dev_string, is_nuc, is_aa, *code ) ){
      fprintf( stderr, "Invalid deviation string %s. Aborting.\n", 
	       dev_string );
      usage();
    }

  if( starts_string )
    if( get_starts( starts_string, param, is_nuc, idx, *code ) ){
      fprintf( stderr, "Invalid starts string %s. Aborting.\n", 
	       starts_string );
      usage();
    }

  /* If min_prot_lg and min_coding_reg_lg both set: error */
  /* If min_prot_lg set but not min_coding_reg_lg, set the later to 
     min_prot_lg */
  /* If min_coding_reg_lg set but not min_prot_lg, set the later to 0 */
  /* If neither is set, set both to defaults */
  if( ( param->min_coding_reg_lg != -1 ) && ( param->min_prot_lg != -1 ) )
    usage();
  if( param->min_prot_lg != -1 ) param->min_coding_reg_lg = param->min_prot_lg;
  else
    if( param->min_coding_reg_lg == -1 ) 
      param->min_coding_reg_lg = DEF_MIN_CODING_REG_LG;

  if( param->min_prot_lg == -1 ) param->min_prot_lg = DEF_MIN_PROT_LG;
 

  if( !is_valid_sfile(param->s_file) ){
    fprintf(stderr, "Invalid sequence file name %s. Aborting\n", 
	    param->s_file);
    usage();
  }

  /* print params if not silent */
  if( !param->silent ){
    printf("\nGenetic code used           : %d\n", param->code_id );
    if( dev_string ){
      printf("Deviation(s)                : %s\n", dev_string); 
    }
    if( starts_string ){
      printf("Start codon(s) to use       : %s\n", starts_string); 
    }
    printf("Circular genome ?           : %s\n", 
	   param->circular ? "Yes" : "No");
    printf("Minimum orf length          : %ld\n", param->min_coding_reg_lg );
    printf("Minimum protein length      : %ld\n", param->min_prot_lg );
    printf("Force start codons as 'M' ? : %s\n", 
	   param->force_met ? "Yes" : "No");
    printf("\n-----------------\n\n");
  }
}

/*-----------------------------------------------------------------------------

  Name: usage 

  Description: prints a usage error message and exits w/ an exit status of 1
  
  Returns: Nothing
  
  ---------------------------------------------------------------------------*/

void usage()
{
  fprintf(stderr, "Usage: flip [-S] [-c] [-f] [-l length | -L length] [-m]");
  fprintf(stderr, " [-d deviation]\n         [-s start_codons] [-g gen_code]");
  fprintf(stderr, " [-n] [-D] [-F] file\n");
  exit(1);
}

/*-----------------------------------------------------------------------------

  Name: get_starts

  Description: parse the start string read on command line or prot.prm and
               updates the starts array. The new starts will replace the old
	       ones
  
  Returns: 1 if the starts string is invalid, 0 otherwise
  
  ---------------------------------------------------------------------------*/


int get_starts( 
     char *start_string,  /* the starts string */ 
     PARAM *param,        /* ptr on the param struct */
     int is_nuc[],        /* array telling if a char is a nuc or not */
     int idx[],           /* index of all nucs */
     char *code           /* current genetic code */
)
{
  char     codon[CODON_LG+1];
  long int str_lg;
  int      start_idx;
  int      i;

  if( !start_string ) return(1);

  /* Reinitialize all codons as not being start codons */
  for( i=0; i < NB_CODONS ; i++ ){
    param->starts[i] = 0;
  }

/*  while( start_string[i] != '\0' ) printf( "%c", start_string[i++] );*/

  str_lg = strlen(start_string);
  while( str_lg >= CODON_LG ){
    strncpy(codon, start_string, CODON_LG);
    codon[CODON_LG] = '\0';
    
    /* get a (possible) codon */
    for( i=0 ; i< CODON_LG ; i++ ) 
      if( is_nuc[ codon[i] ] ) codon[i] = toupper( codon[i] );
      else return(1);
    
    if( str_lg > CODON_LG )
      if( str_lg != (CODON_LG + 1) && ( *(start_string + CODON_LG) == ',' ) ){
	str_lg -= CODON_LG+1;
	start_string += CODON_LG +1;
      }
      else return(1);
    else{
      str_lg -= CODON_LG;
      start_string += CODON_LG;
    }

    /* Verify that codon is not already a stop */
    start_idx = get_codon_idx( codon, idx );
    if( code[ start_idx ] == '*' ){
      fprintf( stderr, "%s defined as both start AND stop codon\n", codon );
      return(1);
    }
    else param->starts[ start_idx ] = 1;
  }
  
  if( str_lg ) return(1);
  return(0);
}

/*-----------------------------------------------------------------------------

  Name: get_deviation

  Description: sets deviations to the genetic code according to a deviation
               string passed on the command line or read in prot.prm. Casing of
	       aas are maintained, while codons are all uppercased.
  
  Returns: 1 if the dev string is invalid, 0 otherwise
  
  ---------------------------------------------------------------------------*/

int get_deviations( 
	char *dev_string,   /* the dev string */
	int is_nuc[],       /* array telling if a char is a nuc or not */
	int is_aa[],        /* array telling if a char is a aa or not */
	char code[]         /* current genetic code */
)
{
  long int str_lg;            /* length of the dev string */
  char     codon[CODON_LG+1]; /* a codon */
  char     aa;                /* an aa */
  int      i;                 /* loop counter */

  if( !dev_string ) return(1);

  /* parse dev string */
  str_lg = strlen( dev_string );
  while( str_lg >= CODON_LG + 2 ){
    strncpy( codon, dev_string, CODON_LG );
    codon[CODON_LG] = '\0';

    /* Verify that codon is valid */    
    for( i=0 ; i< CODON_LG ; i++ ) 
      if( is_nuc[ codon[i] ] ) codon[i] = toupper( codon[i] );
      else return(1);
   
    if( dev_string[CODON_LG] != '=' ) return(1);
    aa = dev_string[ CODON_LG + 1 ];
    if( !is_aa[aa] ) return(1);

    if( str_lg > CODON_LG + 2 )
      if( str_lg != (CODON_LG + 3) && ( dev_string[ CODON_LG + 2 ] == ',' ) ){
	str_lg -= CODON_LG + 3;
	dev_string += CODON_LG + 3;
      }
      else return(1);
    else{
      str_lg -= CODON_LG + 2;
      dev_string += CODON_LG + 2;
    }

    /* Note that casing of the aa is maintained, while codon is uppercased */
    code[ get_codon_idx( codon, idx ) ] = aa;
  }

  /* If str_lg is not 0, dev_string is invalid */
  return( str_lg );
}

/*-----------------------------------------------------------------------------

  Name: parse_prot_prm

  Description: Parses prot.prm to set application's parameters
  
  Returns: Nothing
  
  ---------------------------------------------------------------------------*/


void parse_prot_prm( 
     PARAM *param,          /* ptr on application's params */
     char **dev_string,     /* will contain the dev string read */
     char **starts_string   /* will contain the starts string read */
)
{
  FILE *prm;                       /* ptr on prot.prm */
  char  buffer[ MAX_LINE_LG + 1 ]; /* a line in prot.prm */
  char  f_str[ MAX_LINE_LG + 1 ];  /* tmp */
  char  field[ MAX_LINE_LG + 1 ];  /* a field id */
  char  value[ MAX_LINE_LG + 1 ];  /* the value of a field */
  char  rest[ MAX_LINE_LG + 1 ];   /* the things after the value */
  char *end;                       /* tmp */
  int   n;                         /* tmp  */
  int   i;                         /* loop counter */

  if( (prm = fopen( PROT_PRM_NAME, "r" )) ){
    if( !param->silent ) printf( "Reading %s\n", PROT_PRM_NAME );
    while( fgets( buffer, MAX_LINE_LG, prm ) != NULL ){
      if( !is_blank_line( buffer ) && !is_comment( buffer ) ){

        /* make sure that line has the format FIELD = VALUE */
	if( sscanf( buffer, "%[ \ta-zA-Z0-9]=%s%s", f_str, value, rest ) == 2 )
	  {
	    /* Get rid of leading and trailing blanks */
	    if( (n = sscanf( f_str, "%s%s", field, rest )) != 1 )
	      parse_prot_prm_err( buffer );
	    
	    /* Uppercase field id */
	    for(i=0 ; i < strlen(field) ; i++ ) field[i] = toupper( field[i] );
	    
	    /* to set the genetic code number */
	    if( !strcmp( field, "GENCODE" ) ){
	      param->code_id = strtol( value, &end, 10 );
	      if( *end != '\0' ) parse_prot_prm_err( buffer );
	      if( !IS_VALID_CODE_ID( param->code_id ) ) 
		parse_prot_prm_err( buffer );
	      continue;
	    }
	    
	    /* to set the circular genome indicator */
	    if( !strcmp( field, "CIRCULAR" ) ){
	      if( strcmp( value, "0" ) && strcmp( value, "1" ) )
		parse_prot_prm_err( buffer );
	      param->circular = atoi( value );
	      continue;
	    }
	    
	    /* to indicate that prot.6rf is to be produced */
	    if( !strcmp( field, "6RF" ) ){
	      if( strcmp( value, "0" ) && strcmp( value, "1" ) )
		parse_prot_prm_err( buffer );
	      param->trans_all = atoi( value );
	      continue;
	    }
	    
	    /* to set the min orf length */
	    if( !strcmp( field, "MINORFLENGTH" ) ){
	      param->min_coding_reg_lg = strtol( value, &end, 10 );
	      if( *end ) parse_prot_prm_err( buffer );
	      continue;
	    }
	    
	    /* to set the min prot length (from first start codon) */
	    if( !strcmp( field, "MINPROTLENGTH" ) ){
	      param->min_prot_lg = strtol( value, &end, 10 );
	      if( *end ) parse_prot_prm_err( buffer );
	      continue;
	    }
	    
	    /* to set new starts */
	    if( !strcmp( field, "STARTS" ) ){
	      *starts_string = 
		(char *) calloc( strlen(value) + 1, sizeof(char) );
	      strcpy( *starts_string, value );
	      continue;
	    }
	    
	    /* to set deviations */
	    if( !strcmp( field, "DEVIATION" ) ){
	      *dev_string = (char *) calloc( strlen(value) + 1, sizeof(char) );
	      strcpy( *dev_string, value );
	      continue;
	    }
	    
	    /* If we're here, the line read is invalid */
	    parse_prot_prm_err( buffer );
	  }
	else parse_prot_prm_err( buffer );
      }
    }
    fclose( prm );
  }
  else printf( "Using default parameters\n" );
}

/*-----------------------------------------------------------------------------

  Name: is_comment

  Description: detects if a given line in prot.prm is a comment (first non-
               whitespace is '#')
  
  Returns:1 if line is a comment, 0 otherwise
  
  ---------------------------------------------------------------------------*/


int is_comment(
	char *string  /* a line in prot.prm */
)
{
  for( ; *string; string++ ){
    if( *string == '#' ) return(1);
    if( !isspace( *string ) ) return(0);
  }
  return(0);
}

/*-----------------------------------------------------------------------------

  Name: parse_prot_prm_err

  Description: Signals to the user thatprot.prm is in invalid format and
               prints an error message including the faulty line, then exits
	       with an eit status of 1
  
  Returns: Nothing
  
  ---------------------------------------------------------------------------*/


void parse_prot_prm_err( char *line )
{
 
  /* Notice that *line is modified... */
  if( line[ strlen(line) - 1 ] == '\n' ) line[ strlen(line) - 1 ] = '\0';
  fprintf(stderr, "Invalid line \"%s\" in %s. Aborting\n", line,
	  PROT_PRM_NAME);
  exit(1);
}


/*-----------------------------------------------------------------------------

  Name: nuc_lg

  Description: Gets the length of the nuc. seq.( including spec. OGMP chars 
               and spaces ) associated to prot in orig seq.
  
  Returns: the computed length 
  
  ---------------------------------------------------------------------------*/

long int nuc_lg( 
    PROT *prot,         /* pointer on the protein */
    CTG *ctg,           /* pointer on current contig */
    long int seq_start, /* start of prot in original sequence string */ 
    long int seq_end    /* end of prot in original sequence string */
)
{
  long int lg;   /* the length */

  /* if length <= total length of the cleaned sequence */
  if( seq_start >= seq_end )
    if( prot->strand ) lg = seq_start - seq_end + 1;
    else lg = ctg->seq_lg - seq_start + 1 + seq_end + 1;
  else
    if( prot->strand ) lg = ctg->seq_lg - seq_end + 1 + seq_start + 1;
    else lg = seq_end - seq_start + 1;

  if( CODON_LG * prot->aa_lg > 2 * ctg->clean_seq_lg ){
    lg += 2 * ctg->seq_lg;
    if( seq_start == seq_end ) lg--;
  }
  else
    if( CODON_LG * prot->aa_lg > ctg->clean_seq_lg ){
      lg += ctg->seq_lg;
      if( seq_start == seq_end ) lg--;
    }

  return( lg );
}

/*-----------------------------------------------------------------------------

  Name: find_arrow_ptr

  Description: Finds a pointer on an arrow string (either '==>' or '<==') in 
               an annotation. It makes sure that the arrow is actually part
	       of the annotation and not part of a possible comment at the
	       end of the annot (by checking that if an arrow is found and
	       that if the is a semi-colon after the first one of the annot, 
	       then the arrow occurs BEFORE the semi-colon. If there are no
	       semi-colon aside from the leading one, this check is not done)
	       It also checks that the arrow is not part of a qualifier text
	       string.
  
  Returns: The ptr on the valid arrow found, or NULL if no valid arrow was 
           found
  
  ---------------------------------------------------------------------------*/

char *find_arrow_ptr(
    char *annot,            /* the annot to search in */
    char *arrow_ptr,        /* ptr on the arrow string to search in 'annot' */
    char *compl_arrow_ptr   /* ptr on the compl. arrow (for '==>' it's '<==', 
			       and vice-versa)  */
)
{
  char *ptr;
  char *semi_colon_ptr;
  char *compl_ptr;

  ptr = strstr( annot, arrow_ptr );
  if( !ptr ) return(NULL);

  semi_colon_ptr = strstr( annot + 1, ";" );
  if( semi_colon_ptr && (semi_colon_ptr < ptr) ) return( NULL );

  /* try to find the complement arrow (for '==>' it's '<==', and vice-versa)
     If there's one before ptr, then return NULL. This will handle the case
     ";  G-nad5 <== /substitution="G==>A"  ", for which, without this loop, a
     ptr on the forw arrow would be returned if we sought '==>' */
  else{ 
    compl_ptr = strstr( annot, compl_arrow_ptr );
    if( compl_ptr && (compl_ptr < ptr) ) return( NULL );
    else return( ptr );
  }
} 

/*-----------------------------------------------------------------------------

  Name: is_valid_line

  Description: Determines if a line in the masterfile is valid or not
  
  Returns: 1 for valid lines, 0 otherwise
  
  ---------------------------------------------------------------------------*/

int is_valid_line( 
      char *line,       /* pointer on the beginning of the line */
      int is_seq_char[] /* array telling if a char is one of aAcCgGtTnNxX#@!+- or
			   any whitespace */
)
{
  while( isdigit( *line ) ) line++;
  while( isspace( *line ) ) line++;

  /* Only sequencing chars are allowed now */
  for( ; *line; line++ )
    if( !is_seq_char[*line] ) return(0);
  return(1);
}

/*-----------------------------------------------------------------------------

  Name: fix_starts

  Description: Looks at every protein. If one starts with a valid start codon
               then the first amino acid of this protein is automatically
               changed to 'M'
  
  Returns: Nothing
  
  ---------------------------------------------------------------------------*/

void fix_starts( 
     CTG    *ctg, 
     PROT   *prot, 
     PARAM  *param,
     int     idx[],
     int     is_nuc[],
     int    comp[]
)
{
  CODON     codon;
  PROT     *cur_prot;
  int       i;
  int       inc;
  long int  pos;
  char      nuc;

  for( cur_prot = prot; cur_prot ; cur_prot = cur_prot->next ){
    i = 0;
    inc = cur_prot->strand ? -1 : 1;
    pos = cur_prot->begin - 1;
    while( i != CODON_LG ){
      if( is_nuc[ ctg->clean_seq[pos] ] ){
	nuc = 
	  cur_prot->strand ? comp[ ctg->clean_seq[pos] ] : ctg->clean_seq[pos];
	codon[i++] = toupper( nuc );
      }
      pos += inc;
      if( pos < 0 ) pos = ctg->clean_seq_lg - 1;
      if( pos >= ctg->clean_seq_lg ) pos = 0;
    }
    codon[i] = '\n';

    if( param->starts[ get_codon_idx( codon, idx ) ] )
      *( cur_prot->aa ) = 'M';
  }
}

/*-----------------------------------------------------------------------------

  Name: sort_prots

  Description: Sorts the proteins according the the function prot_compare (see
               below). Takes the list of proteins whose start is pointed to by 
	       prot and builds an array of pointers to prots. After that, the
	       pointers are sorted according to prot_compare.
  
  Returns: Nothing
  
  ---------------------------------------------------------------------------*/

void sort_prots( 
   PROT   *prot,         /* The list of all proteins in a contig */
   PROT ***sorted_prot,  /* Will contain the list of sorted prot pointers */
   int    *nb_prot       /* Will contain the nb of proteins in the contig */
)
{ 
  int   i;          /* counter */
  PROT *cur_prot;   /* prot pointer */

  /* Calculate number of proteins */
  *sorted_prot = NULL; 
  *nb_prot = 0;
  for( cur_prot = prot; cur_prot ; cur_prot = cur_prot->next ) (*nb_prot)++;

  /* If there are any */
  if( *nb_prot ){
    /* stock each protein pointer in an array */
    *sorted_prot = (PROT **) calloc( sizeof(PROT *), *nb_prot );
    for( i = 0, cur_prot = prot; cur_prot ; cur_prot = cur_prot->next, i++ ){
      (*sorted_prot)[i] = cur_prot;
    }
    
    /* sort the pointers, according to sort function */
    qsort( *sorted_prot, *nb_prot, sizeof(PROT *), prot_compare );
  }
}

/*-----------------------------------------------------------------------------

  Name: prot_compare

  Description: Used to sort proteins: proteins on original strand before 
               before proteins on compl. strand. For proteins on the original
               strand, sort them in increasing order of start position. For
               proteins on compl. str., sort them in decreasing order of start
               position
  
  Returns: 1, 0 or -1 depending on the comparison between the two protein
           pointers
  
  ---------------------------------------------------------------------------*/

int prot_compare( 
       const void *p1,
       const void *p2 
)
{
  PROT *prot1     = * ((PROT **)p1);
  PROT *prot2     = * ((PROT **)p2);
  int strand1     = prot1->strand;
  int strand2     = prot2->strand;
  long int begin1 = prot1->begin;
  long int begin2 = prot2->begin;

  if( strand1 != strand2 ){
    return( strand1 - strand2 );
  }
  else
    if( strand1 ){ return( begin2 - begin1 ); }
    else{ return( begin1 - begin2 ); }
}

/*-----------------------------------------------------------------------------

  Name: is_valid_sfile

  Description: Verifies that the sequence file name is valid, ie that it is
               not the name of one of the file that flip produces.
  
  Returns: 0 if the file name is invalid, 1 otherwise 
  
  ---------------------------------------------------------------------------*/

int is_valid_sfile( 
		   char *fname   /* pointer on the filename */
)
{
  char *tmp;  /* temporary char pointer */

  tmp = PROT_6RF_NAME;
  if( !strcmp(fname, tmp) ) return 0;

  tmp = PROT_LST_DNA_NAME;
  if( !strcmp(fname, tmp) ) return 0;

  tmp = PROT_LST_NAME;
  if( !strcmp(fname, tmp) ) return 0;

  tmp = PROT_SRC_NAME;
  if( !strcmp(fname, tmp) ) return 0;

  tmp = NOCOMPL_NAME;
  if( !strcmp(fname, tmp) ) return 0;

  tmp = COMPL_NAME;
  if( !strcmp(fname, tmp) ) return 0;

  tmp = PROT_PRM_NAME;
  if( !strcmp(fname, tmp) ) return 0;

  return 1;
}

/*-----------------------------------------------------------------------------

  Name: get_dna_seq

  Description: Gets the DNA sequence associated to a given orf
               Fills the 'dna' field of the prot structure
  
  Returns: Nothing
  
  ----------------------------------------------------------------------------*/

void get_dna_seq(
		 PROT *prot, /* Pointer on an orf */
		 CTG  *ctg   /* Pointer on the contig the orf is in */
)
{
  long int  dna_start; /* Orf's start (1..ctg_lg) in orig (or reverse) strand. 
			  We always have dna_start <= dna_end for ANY orf, 
			  except the circular ones for which the relationship 
			  might fail */
  long int  dna_end;   /* Orf's end (1..ctg_lg) in orig (or reverse) strand. 
			  We always have dna_start <= dna_end for ANY orf, 
			  except the circular ones for which the relationship 
			  might fail */   
  char     *whole_dna; /* Seq from which the DNA is fetched (orig or rev
			  strand) */
  long int  cur_lg;    /* Cur lg of DNA chunk fetched */
  long int  cur_start; /* DNA chunk seq start  */
  long int  cur_end;   /* DNA chunk seq end  */
  
  
  prot->dna_lg = CODON_LG * prot->aa_lg;
  prot->dna = (char *) calloc( prot->dna_lg + 1, sizeof( char ) );

  /* Start/end of the dna sequence, but on the ORIGINAL strand */
  dna_start = prot->begin;
  dna_end   = prot->end;
  
  /* If on reverse strand, we want pos. on rev compl. strand so rev. and compl.
     start/end pos */
  if( prot->strand ){
    dna_start = ctg->clean_seq_lg - dna_start + 1;
    dna_end   = ctg->clean_seq_lg - dna_end   + 1;
  }

  whole_dna = prot->strand ? ctg->compl_clean_seq : ctg->clean_seq;

  /* Get the DNA seq of the orf. Will work even for circular orfs */
  cur_lg       = 0;
  cur_start    = dna_start;
  cur_end      = dna_end;
  prot->dna[0] = '\0';
  while( cur_lg + cur_end - cur_start + 1 <= prot->dna_lg ){
    if( cur_lg + cur_end - cur_start + 1 == prot->dna_lg ){
      strncat( prot->dna, whole_dna + cur_start - 1, cur_end - cur_start + 1 );
      break;
    }
    else{
      strcat( prot->dna, whole_dna + cur_start - 1 );
      cur_lg += ctg->clean_seq_lg - cur_start + 1;
      cur_start = 1;
    }
  }
  
  /* Terminate string properly */
  prot->dna[prot->dna_lg] = '\0';
}

/*-----------------------------------------------------------------------------

  Name: update_prot_lst_dna

  Description: Adds all the orf's DNA (for current contig) in prot.lst.dna
  
  Returns: Nothing
  
  ----------------------------------------------------------------------------*/

void update_prot_lst_dna(
			 PROT     **sorted_prot,  /* The orfs, in appropriate
						     order */
			 long int   nb_prot,      /* Nb of prots found */
			 FILE      *prot_lst_dna, /* Pointer on prot.lst.dna 
						     file */
			 int        number_prot   /* boolean. Indicates if we
						     should append a 
						     "_<protein_nb>" to indicate
						     the protein number in the
						     output files */
)
{
  long int prot_idx;
  long int i;
  char     chunk[NUC_PER_LINE+1]; /* a chunk of aa sequence */
  
  /* if there were no proteins found, leave */
  if( !nb_prot ) return;

  for( prot_idx = 0; prot_idx < nb_prot; prot_idx++ ){
    prot = sorted_prot[prot_idx];

    /* Print prot headers */
    if( !number_prot )
      fprintf( prot_lst_dna, "%s; ", prot->ctg_name );
    else
      fprintf( prot_lst_dna, "%s_%d; ", prot->ctg_name, prot_idx+1 );      
 
    fprintf( prot_lst_dna, "%s ", prot->strand ? "compl." : "orig." );
    fprintf( prot_lst_dna, "%6ld to %6ld\n", prot->begin, prot->end );
    
    /* Print aa sequence formatted to NUC_PER_LINE nucleotides per line */
    /* In prot.lst, the ending stop If any) is printed.  */
    for( i = 1; i <= prot->dna_lg ; i += NUC_PER_LINE ){
      strncpy( chunk, prot->dna + i - 1, NUC_PER_LINE );
      chunk[NUC_PER_LINE] = '\0';
      fprintf( prot_lst_dna, "%6ld  %s\n", i, chunk );
    }

    fprintf( prot_lst_dna, "\n" );
  }
}
