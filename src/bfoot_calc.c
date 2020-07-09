
/* BaGFoot_calc.c :

 Programmed by Songjoon Baek 2017
 
 Last Update : 01.05.2017
*/


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>

#define max(a,b) ((a)>(b) ? (a):(b))
#define min(a,b) ((a)<(b) ? (a):(b))



#define MAXNUMCHROM 100
#define BUFFERSIZE 2000
#define MAXCHRLENGTH 50


typedef char fixed_string[MAXCHRLENGTH];



SEXP intersect_intervals(SEXP X1,SEXP X2,SEXP Y1,SEXP Y2, SEXP Isorted1, SEXP Isorted2)
{

  int i, j,k, s, starter;
  
  SEXP I1;
  
  int len1 = LENGTH(Isorted1);
  int len2 = LENGTH(Isorted2);   
  
  PROTECT(I1=allocVector(INTSXP, len1+len2));
  int result;
  
  double *pY1 = REAL(Y1);
  double *pY2 = REAL(Y2);
  double *pX1 = REAL(X1);
  double *pX2 = REAL(X2);
  int *pI1 = INTEGER(I1);
  
  //int *pI2 = INTEGER(I2);
  int *pIsorted1 = INTEGER(Isorted1);
  int *pIsorted2 = INTEGER(Isorted2);
  
   k=0;starter=0;
   for (int ii=0; ii<len1;ii++) {
     i=pIsorted1[ii]-1;
     s=starter;
     for (int jj=s; jj<len2;jj++) {
        j=pIsorted2[jj]-1;
        result = (pY1[j]> pX2[i]) ?  -1 : ((pY2[j]< pX1[i]) ? 1: 0);
        if (result==0)
        {
          pI1[k] = i+1;
          k++;
        }
        if (result < 0) break;
      
        if (result > 0) starter = jj;
     }
    }
    SET_LENGTH(I1, k);
    UNPROTECT(1);
    return(I1);
}


SEXP intersect_intervals_index(SEXP X1,SEXP X2,SEXP Y1,SEXP Y2, SEXP Isorted1, SEXP Isorted2)
{
  int i, j,k, s, starter;
  
  SEXP I1,I2;
  
  int len1 = LENGTH(Isorted1);
  int len2 = LENGTH(Isorted2);   
  
  PROTECT(I1=allocVector(INTSXP, len1+len2));
  PROTECT(I2=allocVector(INTSXP, len1+len2));
  
  int result;
  
  double *pY1 = REAL(Y1);
  double *pY2 = REAL(Y2);
  double *pX1 = REAL(X1);
  double *pX2 = REAL(X2);
  int *pI1 = INTEGER(I1);
  int *pI2 = INTEGER(I2);
  
  const char *names[] = {"I1", "I2",""}; 
  
  SEXP res = PROTECT(mkNamed(VECSXP, names));
  
  int *pIsorted1 = INTEGER(Isorted1);
  int *pIsorted2 = INTEGER(Isorted2);
  
   k=0;starter=0;
   for (int ii=0; ii<len1;ii++) {
     i=pIsorted1[ii]-1;
     s=starter;
     for (int jj=s; jj<len2;jj++) {
        j=pIsorted2[jj]-1;
        result = (pY1[j]> pX2[i]) ?  -1 : ((pY2[j]< pX1[i]) ? 1: 0);
        if (result==0)
        {
          pI1[k] = i+1;
          pI2[k] = j+1;
          k++;
        }
        if (result < 0) break;
      
        if (result > 0) starter = jj;
     }
    }
   SET_LENGTH(I1, k);
   SET_LENGTH(I2, k);
   
   SET_VECTOR_ELT(res, 0, I1);       /* numeric(1) */ 
   SET_VECTOR_ELT(res, 1, I2);   /* numeric(<some length>) */

   UNPROTECT(3);
   return(res);
}

/* readcutcount */

int max_intArray(int * a, int length) {
    int M;
    M = a[1];
    for (int i = 1; i <= length; i++) {
        if (M < a[i])
            M = a[i];
    }
    return M;
}

int min_intArray(int * a, int length) {
    int m;
    m = a[1];
    for (int i = 1; i <= length; i++) {
        if (m > a[i])
            m = a[i];
    }
    return m;
}


int count_line_num(char *filename) {
    FILE *f;
    char c;
    
    int lines = 0;
    f = fopen(filename, "r");

    if (f == NULL) {
        printf("Error: can't open file:%s\n", filename);
        exit(-1);
    }

    int line = 0;
    while ((c = fgetc(f)) != EOF) {
        if (c != ' ' && c != '\t' && c != '\n') {
            line = 1;
        }
        if (c == '\n') {
            lines = lines + line;
            line = 0;
        }
    };
    fclose(f);
    return lines;
}

int read_data(char* filename, int *pos, int *posend, char *str, int nline) {
    // returns number of reads

    char chrNum[10];
    
    printf("reading file %s...\n", filename);

    FILE *f;
    f = fopen(filename, "r");

    for (int ii = 1; ii <= nline; ii++) {
        fscanf(f, "%s %d %d %c", chrNum, &(pos[ii]), &(posend[ii]), &(str[ii]));
        //       printf("%s %d %d %c\n", chrNum, pos[ii], posend[ii], str[ii]);
    }

    fclose(f);

    return nline;
}

void readcutcount(char **pdatafilepath, char **pchr, int* span, int* count) {
   
    char *datafilepath = *pdatafilepath;
    char *chr = *pchr;

    int* pos;
    int* posend;
    char* str;


    char filename[500];
    
    int MM;
    int nline,nreads;
        
    sprintf(filename, "%s_%s.txt", datafilepath, chr);
     
    nline = count_line_num(filename);       
    
    if (nline < 1)
     return;   
     
    printf("line # is %d\n", nline);
    
    pos = (int *) malloc( (nline+1) * sizeof(int));
    posend = (int *) malloc( (nline+1) * sizeof(int));
    str = (char *) malloc( (nline+1) * sizeof(char));


    nreads = read_data(filename, pos, posend, str, nline);

    if (nreads < 0) {
            printf("File not found: %s\n", filename);
            return;
    }

    
    MM = max_intArray(posend, nreads);
    
    *span = MM; // required length of memory
       
    int p,q;
    
    for (int i = 0; i <= *span; i++)
            count[i] = 0;

    for (int i = 1; i <= nreads; i++) {
            p = pos[i];
            q = posend[i];

            if (str[i] == 'F') 
                count[p]++;
            else
                count[q+1]++;
    };
        
    free(pos);
    free(posend);
    free(str);
}

/* fastaggregation_fixedwidth */

SEXP CalcAggregationPerChromFixedWidthInC(SEXP x1,SEXP x2,SEXP y1,SEXP y2, SEXP dir, SEXP dvalue)
{
		
	int *px1 = INTEGER(x1);//  x1 = density$st;
	int *px2 = INTEGER(x2);//  x2 = density$ed;
	int *py1 = INTEGER(y1);//  y1 = region$st;
	int *py2 = INTEGER(y2);//  y2 = region$ed;
	int *pdir = INTEGER(dir);
	double *pdvalue = REAL(dvalue);

	// numregion = nrow(region);
	SEXP hmap, indstart, indend; 

	
	int numx = LENGTH(x1);//68028
	int numy = LENGTH(y1); //16
	int im,ii,jj,kk,istart,iend,b1,start,end,aa,bb;
	
	int hmapwidth = py2[0]-py1[0]+1;
		
	PROTECT(indstart=allocVector(INTSXP, numy));// indstart = rep(0, numy);
	PROTECT(indend=allocVector(INTSXP, numy)); // indend = rep(0, numy);
	PROTECT(hmap = allocMatrix(REALSXP, numy, hmapwidth));
	
	int *pindstart = INTEGER(indstart);
	int *pindend = INTEGER(indend);
	double *phmap = REAL(hmap);
	double x0,tmp;
	
	for (int i0=0; i0<numy; i0++) {
		for (int j0=0; j0<hmapwidth; j0++) {
			phmap[i0 + numy*j0] = 0.0;
		}
	}
	
	im=0;

	for (int i0=0; i0<numy; i0++) { // for (i0 in 1:numy) {	
//		printf("i0=%d : \n", i0);
		x0 = py1[i0];
		istart = 0; iend = numx-1;
		while (istart < iend-1) {
		   im = (istart+iend)/2;  // binary search		#   cat(sprintf('step %d: [%d, %d]  middle=%d \n', step, istart, iend, im));
		   if (x0 < px1[im]) {
			iend =im;
		   } else {
			istart =im;
		   }
		}
		
		pindstart[i0] = max(0, im-1);
		ii=pindstart[i0];
		while (ii< numx) {
		    if (px1[ii] > py2[i0])
				break;
		    ii++;
		}
		pindend[i0] = min(ii,numx-1);
	}

	for (int i0=0; i0<numy; i0++) {
		b1=py1[i0]; //b2=py2[i0];

		for (jj=pindstart[i0]; jj<=pindend[i0]; jj++) {
			 start=max(0,px1[jj]-b1);
			 end = min(px2[jj]-b1, hmapwidth-1);

				
			 if (start<=end) {
				 for (kk=start; kk<=end; kk++) {
					phmap[i0+ numy * kk] = pdvalue[jj];
				}
			 }
		}
		if (pdir[i0]!=1) {
			for (jj=0; jj< hmapwidth/2; jj++) {
				aa= i0+jj*numy;
				bb= i0+ (hmapwidth-jj-1)*numy;
				
				tmp=phmap[aa];
				phmap[aa] = phmap[bb];
				phmap[bb] = tmp;
			}
		}
	}
	
	UNPROTECT(3);
	return(hmap);
}


//extern "C" {
		
	SEXP scanCompact(SEXP fname) {
			
	//    PROTECT(filename=AS_CHARACTER(fname));
		
	//    SEXP chroms, SEXP offsets,  SEXP numChrom
		SEXP chroms;
		SEXP linestart;
		SEXP result;
		SEXP resnames;
		int nprotect=0;
		
		PROTECT(fname = AS_CHARACTER(fname));nprotect++;
		const char* filename;
		filename = STRING_VALUE(fname);
		PROTECT(result = NEW_LIST(3)); nprotect++;
				
				
		
		FILE *fp;

		int chromchars=20;
		
		int lstart[MAXNUMCHROM];
		fixed_string lchrom[MAXNUMCHROM];
		
		char* linebuffer;
		char* chr;
		char  current_chr[chromchars];
		char buffer1[BUFFERSIZE];
		char buffer2[BUFFERSIZE];
		
		size_t buffersize=(size_t) BUFFERSIZE;
		int numChrom;
		strcpy(current_chr,"");
		
		linebuffer = (char*) malloc(buffersize+1);
	   
		numChrom=0;
		
		if ((fp=fopen(filename,"r"))==NULL) {
			return 0;  /* file not found */
		}
		int nch;
		int coord;
		int nline;
		
		nline=0;	
		while ((nch= getline(&linebuffer, &buffersize, fp))>0) {
			coord = ftell(fp);
			chr= strtok(linebuffer," \t");
			nline++;
				if (strcmp(current_chr,chr)!=0) {
					strcpy(current_chr, chr);
					printf("%s %ld\n", chr, nline);
					
					strcpy(lchrom[numChrom], chr);
					lstart[numChrom] = nline;
					
					numChrom++;
					if (numChrom>MAXNUMCHROM) {
						printf("Error: # Chromosomes > %d\n", MAXNUMCHROM);
						break;
					}
				}
				//sscanf(linebuffer, "%s %d %d %c", chrNum, &(pos), &(posend), &(str));
		}
		
		PROTECT(linestart = NEW_INTEGER(numChrom)); nprotect++;    
		PROTECT(chroms = NEW_CHARACTER(numChrom)); nprotect++;
	

		
		for (int ii=0; ii<numChrom; ii++) {
			SET_STRING_ELT(chroms, ii, mkChar(lchrom[ii]));
			INTEGER(linestart)[ii]=lstart[ii];
		}

		SET_ELEMENT(result,0, chroms);
		SET_ELEMENT(result,1, linestart);
		SET_ELEMENT(result,2, fname);			
			
					
		PROTECT(resnames= NEW_CHARACTER(3)); nprotect++;
		SET_STRING_ELT(resnames,0,mkChar("chr"));
		SET_STRING_ELT(resnames,1,mkChar("linenum"));
		SET_STRING_ELT(resnames,2,mkChar("filename"));
		SET_NAMES(result, resnames);
			
		  //  cout << "Total # reads =" << numtotal << endl;
		   // return 0;
		UNPROTECT(nprotect);
		
		
		free(linebuffer);
		fclose(fp); 
		return result;
		
	}
		
	
	void writeBGRperChromosomeInt(char** filename, char** chrom, int* count, int *span, int * binsize, int * thr) {

	    FILE *fp;
	    char* fname= *filename;
	    char* chr = *chrom;
	    //printf("Filename=[%s]   Chr=[%s]\n",fname, chr);
	    fp = fopen(fname, "w");
	    
	    int ii = 0;
	    int jj = 0;
	    
	    int c, prevc;
	    int  prevj; 
	    
//	    while (ii < *span) {
	    while (ii < *span) {
			c = count[ii];
			if (c > 0) {
				jj = ii + 1;
				while (count[jj] == c && jj < *span)
				    jj++;
				if (jj < *span) {
	                fprintf(fp, "%s %d %d %d\n", chr, ii+1, jj+1, c);
				} else 
					fprintf(fp, "%s %d %d %d\n", chr, ii+1, span+1,  c);
				ii = jj;
			} else {
				ii++;
			}
	    }
	    fclose(fp);
	}
//}
