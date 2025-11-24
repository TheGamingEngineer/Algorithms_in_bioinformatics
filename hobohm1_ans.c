/* M.Nielsen May 2008, mniel@cbs.dtu.dk */

/* 
Copyright (C) 2008-2015 Danish Technical University

This suite of programs and library routine is free software. You can 
redistribute it and/or modify it under the terms of the GNU General Public 
License as published by the Free Software Foundation; either version 2 of the
License, or (at your option) any later version.

In other words, you are free to modify, copy, or redistribute this
source code and its documentation in any way you like, but you must
distribute all derivative versions as free software under the same
terms that I've provided my code to you (i.e. the GNU General Public
License). This precludes any use of the code in proprietary or
commercial software unless your source code is made freely available.

If you wish to use the code under a different Open Source license
that's not compatible with the GPL (like the Artistic License, BSD
license, or the Netscape Public License), please contact me
(Morten Nielsen, mniel@cbs.dtu.dk) for permission.

Incorporation into commercial software under non-GPL terms is possible
by obtaining a specially licensed version from The Danish Technical University.
Contact Morten Nielsen (mniel@cbs.dtu.dk) to arrange licensing terms.

This software is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.
*/


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"

float   p_gapf; 
float   p_gapn;
int     p_minall;
FILENAME	p_blmat;
int	p_verbose;

PARAM   param[] =
{
	"-gf", VFLOAT p_gapf, "penalty for first gap", "11",
	"-gn", VFLOAT p_gapn, "penalty for next gap", "1",
	"-mal", VINT p_minall, "Minimum alignment lenght", "1",
	"-m", VFNAME	p_blmat, "Alignment matrix", "$ALGOHOME/data/BLOSUM50",
	"-v", VSWITCH	p_verbose, "Verbose mode", "0",
	0
};

float	sw_alignment( float **m,		/* Scoring matrix */
	int l1,			/* Length of query sequence */
	int l2,			/* Length of database sequence */
	float fg,		/* Penalty for first gap */
	float ng,		/* Penalty for each of subsequent gaps */
	float **S,		/* match scores, called D in lecture */
	int *firsti,		/* Offset in query sequence */
	int *firstj,		/* Offset in database sequence */
	char *qal,		/* Query alignment */
	char *dal, 		/* DB algnemnt */
	int *alen,		/* alignemnt length */
	char *qseq,		/* Query sequence */
	char *dseq		/* Database sequence */
	)

{
	float **P;		/* P matrix */	
	float **Q;		/* Q matrix */
	int     i, j;
	float   temp1, temp2, temp;
	float   sij, pij, qij;
	float  *Si, *Sp;
	float  *Pi;
	float  *Qi, *Qp;
	float  *Mi;
	float  *Oi;
	int	**eij;
	int	e;
	int	*Ei;
	float	score;
	int	keep_going;
	int	best;
	int	k;

	score = 0;
	(*firsti) = -1;
	(*firstj) = -1;

	eij = imatrix(0, l1, 0, l2);
	P = fmatrix(0, l1, 0, l2);
	Q = fmatrix(0, l1, 0, l2);
	S[l1][l2] = 0.0;
	P[l1][l2] = 0.0;
	Q[l1][l2] = 0.0;

	/* Remember

                              S[m+1][n+1] + m[m][n]
             S[m][n] = max{   P[m][n]
                              Q[m][n]
                              0

                              S[m+1][n] + fg
             P[m][n] = max{   P[m+1][n] + ng


                              S[m][n+1] + fg
             Q[m][n] = max{   Q[m][m+1] + ng

        */

	/* Initialize edges of S, Q and P matrices. */
	/* Here they are all set to 0! */

	for (j = l2 - 1; j >= 0; j--) {
		sij = S[l1][j + 1];
		S[l1][j] = sij;
		P[l1][j] = sij;  /* Here one can penalize ends and set P[l1][j] = sij - fg */
		Q[l1][j] = sij;  /* Here one can penalize ends and set Q[l1][j] = sij - fg */
	}

	for (j = l1 - 1; j >= 0; j--) {
		sij = S[j+1][l2];
		S[j][l2] = sij;
		P[j][l2] = sij;  /* Here one can penalize ends and set P[j][l2] = sij - fg */
		Q[j][l2] = sij;	 /* Here one can penalize ends and set Q[j][l2] = sij - fg */
	}

	/* Loop over Query sequence */
	for (i = l1 - 1; i >= 0; i--) {

		/* indirect array access to speed up code */
		/* Mi, Si, Pi etc points to row i in the corresponding matrix */

		Mi = m[i];
		Si = S[i];
		Pi = P[i];
		Qi = Q[i];

		Sp = S[i + 1];
		Qp = Q[i + 1];

		Ei = eij[i];

		/* eij is the backtrack direction matrix
			eij = 0 stop back tracking
			eij = 1 match
			eij = 2 gap-opening database
			eij = 3 gap-extension database
			eij = 4 gap-opening query
			eij = 5 gap-extension query
		*/

		for (j = l2 - 1; j >= 0; j--) {

			/* Try match state */

			sij = Sp[j + 1] + Mi[j]; /* sij = S[i+1][j+1] + m[i][j] */

			/* Try Gap in Query sequence (insertion in database sequence) */

			temp1 = Qp[j] - ng; /* Gap extension temp1 = Q[i+1][j] - ng */
                        temp2 = Sp[j] - fg; /* Gap opening   temp2 = S[i+1][j] - fg */

			if ( temp1 > temp2 ) { /* extension best */
				qij = temp1;
				e = 3;
			}
			else { /* gap opening best */
				qij = temp2; 
				e = 2;
			}

			/* Select if match or gap scores best */

			if ( qij > sij ) 
				sij = qij; /* Gap is best */
			else 
				e = 1; /* Match best */

			/* Try Gap in database sequence (insertion in Query sequence) */

			temp1 = Pi[j + 1] - ng; /* Gap extension temp1 = P[i][j+1] - ng */
			temp2 = Si[j + 1] - fg; /* Gap opening   temp2 = S[i][j+1] - fg */

			if ( temp1 > temp2 ) { /* extension best */

				pij = temp1;

				if ( temp1 > sij ) {
                                        sij = temp1;
                                        e = 5;
                                }  

			}
			else { /* gap opening best */

				pij = temp2;

				if ( temp2 > sij ) {
                                        sij = temp2;
                                        e = 4;
                                } 
			}

			if (sij > score) {
				score = sij;
				(*firsti) = i;
				(*firstj) = j;
			}

			if ( sij <= 0 ) {
				sij = 0.0;
				e = 0;
			}

			Si[j] = sij; /* S[i][j] = sij */
			Qi[j] = qij; /* Q[i][j] = qij */
			Pi[j] = pij; /* P[i][j] = pij */
			Ei[j] = e;   /* eij[i][j] = e */
		}
	}

	fmatrix_free(P, 0, l1, 0, l2);
	fmatrix_free(Q, 0, l1, 0, l2);

	/* Do back tracking */

	if (*firsti < 0 || *firstj < 0 ) 
		printf( "No alignment found. Exit\n" );

	(*alen) = 0;

	i = *firsti;
	j = *firstj;
	qal[(*alen)] = qseq[i];
	dal[(*alen)] = dseq[j];
	i++;
	j++;

	(*alen)++;

	keep_going = 1;

#define SMALL 1e-10

	while ((i < l1 ) && (j < l2 ) && keep_going ) {

		if ( eij[i][j] == 0 ) {
			keep_going = 0;
		} 
		else if ( eij[i][j] == 1 ) { /* Match */

			qal[(*alen)] = qseq[i];
			dal[(*alen)] = dseq[j];
			i++;
			j++;
			(*alen)++;

		}
		else if ( eij[i][j] == 4 ) { /* gap opening in Query */

			qal[(*alen)] = '-';
			dal[(*alen)] = dseq[j];
			j++;
			(*alen)++;

		}
		else if (  eij[i][j] == 5 ) { /* gap extension in Query */

			best = j+2;
			Si = S[i];
			temp = Si[best] - fg - ( best-j-1) * ng - Si[j];

			while ( temp*temp > SMALL ) {
				best++;
				temp = Si[best] - fg - ( best-j-1) * ng - Si[j];
			}

			for ( k=j; k<best; k++ ) {
				qal[(*alen)] = '-';
				dal[(*alen)] = dseq[j];
				j++;
				(*alen)++;
			}

		}
		else if ( eij[i][j] == 2 ) { /* gap opening in Database */

			qal[(*alen)] = qseq[i];
			dal[(*alen)] = '-';
			i++;
			(*alen)++;

		}
		else if (  eij[i][j] == 3 ) { /* gap extension in Database */

			best = i+2;
			temp = S[best][j] - fg - ( best-i-1) * ng - S[i][j];

			while ( temp*temp > SMALL ) {
				best++;
				temp = S[best][j] - fg - ( best-i-1) * ng - S[i][j];
			}

			for ( k=i; k<best; k++ ) {
				qal[(*alen)] = qseq[i];
				dal[(*alen)] = '-';
				i++;
				(*alen)++;
			}
		}
	}

	qal[(*alen)] = 0;
	dal[(*alen)] = 0;

	imatrix_free( eij, 0, l1, 0, l2);

	return( score );
}

ALN    *align(float **m, FSALIST *q, FSALIST *d, float gapf, float gapn)

{
	float **sco;
	float   score;
	int     firsti, firstj;
	int     alen;
	ALN    *new = NULL;
	int     i;
	char	qpal[15000], dpal[15000];

	sco = fmatrix(0, q->len, 0, d->len);

	score = sw_alignment(m, q->len, d->len, gapf, gapn, sco, &firsti, &firstj, 
		qpal, dpal, &alen, q->seq, d->seq );

	fmatrix_free(sco, 0, q->len, 0, d->len);

	if (alen < p_minall)
		return( NULL );

	if ((new = aln_alloc()) == NULL) {
		printf("Cannot alloc new ALN\n");
		exit(1);
	}

	new->alen = alen;

	new->score = score;

	new->mlen = 0;
	new->ngap = 0;
	new->nid = 0;

	for (i = 0; i < new->alen; i++) {
		if (qpal[i] != '-' && dpal[i] != '-') {
			new->mlen++;
		} else {
			new->ngap++;
		}
		if (qpal[i] == dpal[i])
			new->nid++;
	}

	new->qof = firsti;
	new->qal = cvector(0, strlen(qpal));
	strcpy(new->qal, qpal);

	new->dof = firstj;
	new->dal = cvector(0, strlen(dpal));
	strcpy(new->dal, dpal);

	strcpy(new->qname, q->name);
	new->qlen = q->len;

	strcpy(new->dname, d->name);
	new->dlen = d->len;

	strcpy(new->type, "SW_ALN");

	new->rscore = -new->score;

	return (new);

}

float   **score_mat( FSALIST *q, FSALIST *d, float **blmat )

{
        float   **scomat;
        int     i,j;
        int     ix, jx;

        scomat = fmatrix( 0, q->len-1, 0, d->len-1 );

        for ( i=0; i<q->len; i++ ) {

                ix = q->i[i];

                if ( ix < 0 ) {
                        printf( "Error. Unconventional amino acid i query sequence %c %s\n", q->seq[i], q->name );
                        exit( 1 );
                }

		if ( ix > 19 )
			continue;
                 
                for ( j=0; j<d->len; j++ ) {

                        jx = d->i[j];

                        if ( jx < 0 ) {
                                printf( "Error. Unconventional amino acid i query sequence %c %s\n", q->seq[i], q->name );
                                exit( 1 );                                                                             
                        }
               
			if ( jx > 19 )
				continue;
               
                        scomat[i][j] = blmat[ix][jx];
                }
        }      
               
        return( scomat );
          
} 

FSALIST     *homolog( FSALIST *keep, FSALIST *fsa, float *score, float **blmat )

{
        ALN     *aln;
        int     alen, nid;
        float   s;
        FSALIST     *k;
	float	**scomat;
      
	/* Loop over all elements in kepp list */
 
        for ( k=keep; k; k=k->next ) {

		/* Make alignment scoring matrix */
		scomat = score_mat( k, fsa, blmat );	
       
		/* Align the two sequences */
                aln = align( scomat, k, fsa, p_gapf, p_gapn );

		/* Free dynamically allocated memory */
		fmatrix_free( scomat, 0, k->len-1, 0, fsa->len-1 );
       
		/* %id > 290/sqrt(alen) => fid * alen > 2.9*sqrt(alen) => 

			nid > 2.9*sqrt(alen) 

		*/

                if ( aln ) {

                        alen = aln->alen;
                        nid = aln->nid;

                        aln_free( aln );

			/* XXXXXXXX Here you must write the missing code. s is the score that defines if nid
                        is too large to include fsa in the list of unique sequences */

			/*                        
                        s = XXXX;
			*/
       
                        s = 2.9*sqrt(alen);
       
                        if ( s < nid ) {
                                *score = s;
                                return( k );
                        }
                }
        }
       
        return( NULL );
}      

main(int argc, char *argv[])

{
	FSALIST		*db_fsa, *d, *keep, *last, *d_next, *hom;
	float   	gapf, gapn;
	float           **blmat, score;
	char            *alphabet;
	int		n, nc;

	pparse(&argc, &argv, param, 1, "db");

	/* Read Blosum substutution scoring matrix from file */

        blmat = read_blosummat( p_blmat, &alphabet );

	if ( ! blmat ) {
		printf( "Error. Cannot read BLOSUM matrix from file %s. Exit\n", p_blmat );
		exit( 1 );
	}

	/* Read the fasta file */

	if ( ( db_fsa = fsalist_read( argv[1] ) ) == NULL ) {
                printf("Cannot read fasta file %s\n", argv[1] );
                exit(1);
        }

	db_fsa = fsalist_check_names( db_fsa );

	/* Assign Blosum alphabet order to variable ->i in all faste entries in db_fsa */

	fsalist_iassign_profile_order( db_fsa );

	gapf = p_gapf;
	gapn = p_gapn;

	printf( "# Gap penalties. fgap: %f. ngap: %f\n", gapf, gapn);

	/* Initiate keep list to be empty */

	keep = NULL;
	n = 0; /* Initiate number of sequenecs */
	nc = 0; /* Initiate number of clusters */
         
	for ( d = db_fsa; d; d=d_next ) {

		/* Remove d from db_fsa list */
         
                d_next = d->next;
                d->next = NULL;
               
                if ( ( hom = homolog( keep, d, &score, blmat ) ) == NULL ) {

			/* Insert d in keep list */

                        if ( keep == NULL ) 	/* Insert d as first element in keep list */
                                keep = d;
                        else 			/* Insert d in end of keep list */
                                last->next = d;
               
                        last = d;
               
                        printf( "# Unique. %i %i %s\n", n, nc, d->name );

			nc++;
                }
                else {
                        printf( "# Not unique. %i %s is homolog to %s %f\n", n, d->name, hom->name, score );
                }

		n++;
        }

	for ( d=keep; d; d=d->next )
		fsa_print( d );

	printf( "# Number of unique clusters %i\n", nc );

	exit(0);
}
