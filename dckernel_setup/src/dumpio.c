/******************************************************************//**
 *functions
 *<utilities>
 *void           dumpio_ednchg                   : endian changer
 *void           dumpio_mk_fname                 : filename generator
 *void           dumpio_set_str, dumpio_trim_str    : string preprocess
 *int32_t        dumpio_syscheck                 : check endian
 *<file r/w>
 *int32_t        dumpio_fopen                    : open file IO stream
 *int32_t        dumpio_fclose                   : close file IO stream
 *int32_t        dumpio_write_data               : write data array
 *int32_t        dumpio_read_data                : read data array
 **********************************************************************/
#define _LARGEFILE_SOURCE
#define _FILE_OFFSET_BITS 64
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define H_SHORT  16
#define H_MID    64
#define H_LONG  256

/* action type */
#define H_FREAD   0
#define H_FWRITE  1
#define H_FAPPEND 2

/* return type */
#define ERROR_CODE   -1
#define SUCCESS_CODE  1

typedef float (real32_t);
typedef double(real64_t);

/* file structure */
typedef struct{
  char fname[H_LONG];
  int32_t opened;
  FILE *fp;
} fileinfo_t;

extern fileinfo_t *dumpio_finfo;

/******************************************************************************/
/* C functions                                                                */
/******************************************************************************/

/** endian change *****************************************************/
extern void dumpio_ednchg( void* avp_pointer,
                           const int32_t ai_size,
                           const int32_t ai_num   );

/** filename generator ************************************************/
extern void dumpio_mk_fname( char *fname,
                             const char *base,
                             const char *ext,
                             int32_t i,
                             int32_t y );

/** string preprocess *************************************************/
extern void dumpio_set_str( char *_str,
                            const char *str,
                            const int str_len );

/** string postprocess *************************************************/
extern void dumpio_trim_str( char *_str,
                             const char *str,
                             const int str_len );

/** check system & initialze ******************************************/
extern int32_t dumpio_syscheck( void );

/** open file IO stream ***********************************************/
extern int32_t dumpio_fopen( const char *vname_in, int32_t mode );

/** close file IO stream **********************************************/
extern int32_t dumpio_fclose( int32_t fid );

/** write data array **************************************************/
extern int32_t dumpio_write_data( int32_t fid,
                                  int32_t idxsize,
                                  void *data );

/** read data array (full size) ***************************************/
extern int32_t dumpio_read_data( int32_t fid,
                                 int32_t idxsize,
                                 void *data );

/* file ID counter */
int32_t dumpio_num_of_file = 0;

/* package+data+status container */
fileinfo_t *dumpio_finfo = NULL;

/* system information */
int32_t dumpio_system_ednchg = 0;

/** endian change *****************************************************/
void dumpio_ednchg( void* avp_pointer,
                    const int32_t ai_size,
                    const int32_t ai_num   )
{
  int ai_csize, ai_cnum;
  char ac_buf[16];
  char *acp_tmp;
  char *acp_local;

  memset(ac_buf,'\0',sizeof(ac_buf));

  acp_tmp   = avp_pointer;
  acp_local = avp_pointer;
  for( ai_cnum=0; ai_cnum<ai_num; ai_cnum++ ) {
    memcpy(ac_buf, acp_local, ai_size);
    for( ai_csize=0; ai_csize<ai_size; ai_csize++ ) {
      *acp_local = ac_buf[ai_size-ai_csize-1];
      acp_local++;
    }
    acp_tmp += ai_size;
  }
}

/** filename generator ************************************************/
void dumpio_mk_fname( char *fname,
                      const char *base,
                      const char *ext,
                      int32_t i,
                      int32_t y )
{
  char _fname[H_LONG];

  switch (y) {
  case 4 :
    sprintf(_fname,"%s.%s%04d",base,ext,i);
    break;
  case 5 :
    sprintf(_fname,"%s.%s%05d",base,ext,i);
    break;
  case 6 :
    sprintf(_fname,"%s.%s%06d",base,ext,i);
    break;
  default :
    break;
  }

  dumpio_trim_str( fname, _fname, H_LONG-1 );
}

/** string preprocess *************************************************/
void dumpio_set_str( char *_str,
                     const char *str,
                     const int str_len )
{
  int i;

  strncpy(_str, str, str_len);

  _str[str_len] = '\0'; /* [fix] H.Yashiro 20120621 */
  for( i=str_len-1; i>=0; i-- ) {
    if( _str[i] == ' ' ) {
      _str[i] = '\0';
    } else {
      break;
    }
  }
}

/** string postprocess *************************************************/
void dumpio_trim_str( char *_str,
                      const char *str,
                      const int str_len )
{
  int i;

  strncpy(_str, str, str_len);

  _str[str_len] = ' ';
  for( i=str_len-1; i>=0; i-- ) {
    if( _str[i] == '\0' ) {
      _str[i] = ' ';
    } else {
      break;
    }
  }
}

/** check system & initialze ******************************************/
int32_t dumpio_syscheck( void )
{
  int32_t i=1;

  if ( (sizeof(real32_t)!=4) || (sizeof(real64_t)!=8) ) {
    printf("Data type (real) is inconsistent!\n");
    exit(1);
  }

  if ( *(char*)&i ) {
    dumpio_system_ednchg = 1;
  } else {
    dumpio_system_ednchg = 0;
  }

  return(SUCCESS_CODE);
}

/** open file IO stream ***********************************************/
int32_t dumpio_fopen( const char *fname, int32_t mode )
{
  int32_t fid;

  /* get file ID */
  fid = dumpio_num_of_file++;

  /* memory re-allocation (expand by new dumpio_num_of_file) */
  dumpio_finfo=(fileinfo_t *)realloc(dumpio_finfo,sizeof(fileinfo_t)*(dumpio_num_of_file));

  /* intitialize */
  dumpio_set_str(dumpio_finfo[fid].fname,fname,H_LONG-1);
  dumpio_finfo[fid].opened = 0;
  dumpio_finfo[fid].fp     = NULL;

  if ( mode==H_FWRITE ) {
    if ( (dumpio_finfo[fid].fp=fopen(dumpio_finfo[fid].fname,"wb"))==NULL ) {
      fprintf(stderr,"Can not open file : %s!\n",dumpio_finfo[fid].fname);
      exit(1);
    }
  } else if( mode==H_FREAD ) { /* [mod] H.Yashiro 20110907 avoid overwrite action */
    if ( (dumpio_finfo[fid].fp=fopen(dumpio_finfo[fid].fname,"rb"))==NULL ) {
      fprintf(stderr,"Can not open file : %s!\n",dumpio_finfo[fid].fname);
      exit(1);
    }
  } else if( mode==H_FAPPEND ) { /* [add] H.Yashiro 20110907 overwrite mode */
    if ( (dumpio_finfo[fid].fp=fopen(dumpio_finfo[fid].fname,"r+b"))==NULL ) {
      fprintf(stderr,"Can not open file : %s!\n",dumpio_finfo[fid].fname);
      exit(1);
    }
  }
  dumpio_finfo[fid].opened = 1;

  return(fid);
}

/** close file IO stream **********************************************/
int32_t dumpio_fclose( int32_t fid )
{
  fclose(dumpio_finfo[fid].fp);
  dumpio_finfo[fid].opened = 0;

  return(SUCCESS_CODE);
}

/** write data array **************************************************/
int32_t dumpio_write_data( int32_t fid,
                           int32_t idxsize,
                           void *data )
{
  int64_t datasize;
  void *_data;

  datasize = idxsize * 8;

  /* data */
  if(dumpio_system_ednchg){
    _data = malloc(datasize);
    memcpy( _data, data, datasize);
    dumpio_ednchg(_data,8,idxsize);
  }else{
    _data = data;
  }

  fwrite(_data,datasize,1,dumpio_finfo[fid].fp);

  if(dumpio_system_ednchg) { free(_data); }

  return(SUCCESS_CODE);
}

/** read data array (full size) ***************************************/
int32_t dumpio_read_data( int32_t fid,
                          int32_t idxsize,
                          void *data )
{
  int64_t datasize;

  datasize = idxsize * 8;

  fread(data,datasize,1,dumpio_finfo[fid].fp);
  if(dumpio_system_ednchg){
    dumpio_ednchg(data,8,idxsize);
  }

  return(SUCCESS_CODE);
}

/******************************************************************************/
/* Fortran interface                                                          */
/******************************************************************************/

extern void dumpio_mk_fname_( char *fname,
                              const char *base,
                              const char *ext,
                              int32_t *i,
                              int32_t *y,
                              int32_t fname_len,
                              int32_t base_len,
                              int32_t ext_len    );

extern void dumpio_syscheck_( void );

extern void dumpio_fopen_( int32_t *fid, const char *fname, int32_t *mode, int32_t fname_len );

extern void dumpio_fclose_( int32_t *fid );

extern void dumpio_write_data_( int32_t *fid, int32_t *idxsize, void *data );

extern void dumpio_read_data_( int32_t *fid, int32_t *idxsize, void *data );

/** filename generator ************************************************/
void dumpio_mk_fname_( char *fname,
                    const char *base,
                    const char *ext,
                    int32_t *i,
                    int32_t *y,
                    int32_t fname_len,
                    int32_t base_len,
                    int32_t ext_len    )
{
  char _base[H_LONG];
  char _ext[H_SHORT];
  int32_t _base_len = base_len;
  int32_t _ext_len  = ext_len;
  if( _base_len > H_LONG-1  ){ _base_len = H_LONG-1;  }
  if( _ext_len  > H_SHORT-1 ){ _ext_len  = H_SHORT-1; }

  dumpio_set_str( _base, base, _base_len );
  dumpio_set_str( _ext,  ext,  _ext_len  );

  dumpio_mk_fname( fname, _base, _ext, *i, *y );
}

/** check system & initialze ******************************************/
void dumpio_syscheck_( void )
{
  int32_t ierr;

  ierr = dumpio_syscheck();
}

/** register new file *************************************************/
void dumpio_fopen_( int32_t *fid, const char *fname, int32_t *mode, int32_t fname_len )
{
  int32_t ierr;

  char _fname[H_LONG];
  int32_t _fname_len = fname_len;
  if( _fname_len > H_LONG-1 ){ _fname_len = H_LONG-1; }

  dumpio_set_str( _fname, fname, _fname_len );

  *fid = dumpio_fopen( _fname, *mode );
}

/** close file IO stream **********************************************/
void  dumpio_fclose_( int32_t *fid )
{
  int32_t ierr;

  ierr = dumpio_fclose( *fid );
}

/** write data array **************************************************/
void dumpio_write_data_( int32_t *fid, int32_t *idxsize, void *data )
{
  int32_t ierr;

  ierr = dumpio_write_data( *fid, *idxsize, data );
}

/** read data array (full size) ***************************************/
void dumpio_read_data_( int32_t *fid, int32_t *idxsize, void *data )
{
  int32_t ierr;

  ierr = dumpio_read_data( *fid, *idxsize, data );
}
