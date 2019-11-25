/**
 * Writer: Xiaolin.liu
 * xiaolin.liu@mail.bnu.edu.cn
 *
 * This head file define basic datatype in the math lib.
 * Referenced from LAL
 * Based on glib.h and comples.h
 *
 * 2019.05.21, UWA
 **/


#include "myUtils.h"


/** < FUNCTIONS > **/

/* Used to print error standard error message and output message */
static INT myPrintVstderr(const CHAR *fmt, va_list ap)
{
    return vfprintf(stderr, fmt, ap);
}


INT print_err(const CHAR *fmt, ...)
{
    INT n = 0;
    va_list ap;
    va_start(ap, fmt);
    n = myPrintVstderr(fmt, ap);
    va_end(ap);
    return n;
}


INT print_log(CHAR *fmt, ...)
{
    INT n = 0;
    INT len = strlen(fmt);
    CHAR paste[len+25];
    memset(paste, '\0', sizeof(paste));
    CHAR log[] = "\033[1;34;47mLOG\033[0m:";
    strcat(paste, log);
    strcat(paste, fmt);

    va_list ap;
    va_start(ap, fmt);
    n = myPrintVstderr(paste, ap);
    va_end(ap);
    return n;
}

INT print_debug(CHAR *fmt, ...)
{
    INT n = 0;
    INT len = strlen(fmt);
    CHAR paste[len+25];
    memset(paste, '\0', sizeof(paste));
    CHAR debug[] = "\033[1;32;47mDEBUG\033[0m:";
    strcat(paste, debug);
    strcat(paste, fmt);
    
    va_list ap;
    va_start(ap, fmt);
    n = myPrintVstderr(paste, ap);
    va_end(ap);
    return n;
}


INT print_warning(CHAR *fmt, ...)
{
    INT n = 0;
    INT len = strlen(fmt);
    CHAR paste[len+25];
    memset(paste, '\0', sizeof(paste));
    CHAR war[] = "\033[4;31mWARNING\033[0m:";
    strcat(paste, war);
    strcat(paste, fmt);
    
    va_list ap;
    va_start(ap, fmt);
    n = myPrintVstderr(paste, ap);
    va_end(ap);
    return n;
}


static INT myPrintVstdout(const CHAR *fmt, va_list ap)
{
    return vfprintf(stdout, fmt, ap);
}


INT print_out(const CHAR *fmt, ...)
{
    INT n = 0;
    va_list ap;
    va_start(ap, fmt);
    n = myPrintVstdout(fmt, ap);
    va_end(ap);
    return n;
}
/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */

/** < create & destroy struct > **/

/* REAL8Vector */
REAL8Vector* CreateREAL8Vector(UINT length)
{
    REAL8Vector* vector;
    vector = (REAL8Vector*)malloc(sizeof(*vector));
    vector->length = length;
    if ( ! length ) /* zero length: set data pointer to be NULL */
    {
        vector->data = NULL;
    }
    else /* non-zero length: allocate memory for data */
    {
        vector->data = (REAL8 *)malloc(length * sizeof( *vector->data));
    }
    
    return vector;
}

/* Destroy */
void DestroyREAL8Vector(REAL8Vector* vector)
{
    if(NULL == vector)
    {
        return;
    }
    if(vector->data)
    {
        vector->length = 0;
        free(vector->data);
    }
    vector->data = NULL;
    free(vector);
    vector = NULL;
    return;
}

REAL8VectorSequence * CreateREAL8VectorSequence ( UINT length, UINT veclen )
{
	
	REAL8VectorSequence *seq;
		
	seq = (REAL8VectorSequence *)malloc( sizeof( *seq ) );
	if ( ! seq )
	  return NULL;
	
	seq->length = length;
	seq->vectorLength = veclen;
	
	if ( ! length || ! veclen )
	  seq->data = NULL;
	else
	{
	  seq->data = (REAL8 *)malloc( length * veclen * sizeof( *seq->data ) );
	  if ( ! seq )
	  {
		free( seq );
		return NULL;
	  }
	}
	
	return seq;
}

void DestroyREAL8VectorSequence ( REAL8VectorSequence * vseq )
{
	if ( ! vseq )
    return;
  if ( ( ! vseq->length || ! vseq->vectorLength ) && vseq->data )
	  return;
  if ( ! vseq->data && ( vseq->length || vseq->vectorLength ) )
	  return;
  if ( vseq->data )
    free( vseq->data );
  vseq->data = NULL; /* leave lengths as they are to indicate freed vector */
  free( vseq );
  return;
}


/* CHARVector */
CHARVector* CreateCHARVector(UINT length, UINT STR_LEN)
{
    CHARVector* vector;
    vector = (CHARVector*)malloc(sizeof(*vector));
    vector->length = length;
    if (STR_LEN > MAX_STR_LEN || !STR_LEN)
    {
        print_err("Warning -%s: STR_LEN = %d is invalid reset it to MAX_STR_LEN = %d\n", __func__, STR_LEN, MAX_STR_LEN);
        STR_LEN = MAX_STR_LEN;
    }
    vector->STR_LEN = STR_LEN;
    if ( ! length ) /* zero length: set data pointer to be NULL */
    {
        vector->data = NULL;
    }
    else /* non-zero length: allocate memory for data */
    {
        vector->data = (CHAR **)malloc(length * sizeof( *vector->data));
        UINT i;
        for(i=0; i < length; ++i)
        {
            vector->data[i] = (CHAR *)malloc(STR_LEN * sizeof(CHAR));
        }
    }
    return vector;
}

/* Destroy */
void DestroyCHARVector(CHARVector* vector)
{
    if(NULL == vector)
    {
        return;
    }
    if(vector->data)
    {
        UINT i;
        for (i=0; i < vector->length; ++i)
        {
            free(vector->data[i]);
            vector->data[i] = NULL;
        }
        vector->length = 0;
        free(vector->data);
    }
    vector->data = NULL;
    free(vector);
    vector = NULL;
    return;
}


/* COMPLEX16Vector */
COMPLEX16Vector* CreateCOMPLEX16Vector(UINT length)
{
    COMPLEX16Vector* vector;
    vector = (COMPLEX16Vector*)malloc(sizeof(*vector));
    vector->length = length;
    if ( ! length ) /* zero length: set data pointer to be NULL */
    {
        vector->data = NULL;
    }
    else /* non-zero length: allocate memory for data */
    {
        vector->data = (COMPLEX16 *)malloc(length * sizeof( *vector->data));
    }
    
    return vector;
}

/* Destroy */
void DestroyCOMPLEX16Vector(COMPLEX16Vector* vector)
{
    if(NULL == vector)
    {
        return;
    }
    if(vector->data)
    {
        vector->length = 0;
        free(vector->data);
    }
    vector->data = NULL;
    free(vector);
    vector = NULL;
    return;
}

/* UINTVector */
UINTVector* CreateUINTVector(UINT length)
{
    UINTVector* vector;
    vector = (UINTVector*)malloc(sizeof(*vector));
    vector->length = length;
    if ( ! length ) /* zero length: set data pointer to be NULL */
    {
        vector->data = NULL;
    }
    else /* non-zero length: allocate memory for data */
    {
        vector->data = (UINT *)malloc(length * sizeof( *vector->data));
    }
    
    return vector;
}
/* Destroy */
void DestroyUINTVector(UINTVector* vector)
{
    if(NULL == vector)
    {
        return;
    }
    if(vector->data)
    {
        vector->length = 0;
        free(vector->data);
    }
    vector->data = NULL;
    free(vector);
    vector = NULL;
    return;
}



/* REAL8Array */
REAL8Array *CreateREAL8Array_comm(UINT *dims, UINT ndim)
{
    UINTVector dimLength;
    UINT size = 1;
    REAL8Array *arr;
    dimLength.length = ndim;
    dimLength.data = dims;
    
    UINT i;
    for(i = 0; i < ndim; ++i)
        size *= dimLength.data[i];
    arr = (REAL8Array *)malloc(sizeof(*arr));
    arr->size = size;
    arr->dimLength = CreateUINTVector(ndim);
    if (! arr->dimLength)
    {
        free(arr);
        arr = NULL;
        return NULL;
    }
    
    memcpy(arr->dimLength->data, dimLength.data, ndim*sizeof(*arr->dimLength->data));
    
    arr->data = (REAL8 *)malloc( size * sizeof(*arr->data));
    if(!arr->data)
    {
        DestroyUINTVector(arr->dimLength);
        free(arr);
        arr = NULL;
        return NULL;
    }
    return arr;

}

REAL8Array *CreateREAL8Array(UINT ndim, ...)
{
    enum { maxdim = 16 };
    if(ndim > maxdim)
    {
        print_err("Error %s: dimension input %d is greater than max dimension %d", __func__, ndim, maxdim);
        return NULL;
    }
    va_list ap;
    UINT dims[ndim];
    UINT dim;
    
    va_start(ap, ndim);
    for(dim = 0; dim<ndim; ++dim)
        dims[dim] = va_arg(ap, UINT);
    va_end(ap);
    return CreateREAL8Array_comm(dims, ndim);
}

/* Destroy */
void DestroyREAL8Array(REAL8Array* array)
{
    DestroyUINTVector( array->dimLength );
    free( array->data );
    array->data = NULL;
    free( array );
    array = NULL;
}


/* COMPLEX16Array */
COMPLEX16Array *CreateCOMPLEX16Array_comm(UINT *dims, UINT ndim)
{
    UINTVector dimLength;
    UINT size = 1;
    COMPLEX16Array *arr;
    dimLength.length = ndim;
    dimLength.data = dims;
    
    UINT i;
    for(i = 0; i < ndim; ++i)
        size *= dimLength.data[i];
    arr = (COMPLEX16Array *)malloc(sizeof(*arr));
    arr->size = size;
    arr->dimLength = CreateUINTVector(ndim);
    if (! arr->dimLength)
    {
        free(arr);
        arr = NULL;
        return NULL;
    }
    
    memcpy(arr->dimLength->data, dimLength.data, ndim*sizeof(*arr->dimLength->data));
    
    arr->data = (COMPLEX16 *)malloc( size * sizeof(*arr->data));
    if(!arr->data)
    {
        DestroyUINTVector(arr->dimLength);
        free(arr);
        arr = NULL;
        return NULL;
    }
    return arr;
    
}

COMPLEX16Array *CreateCOMPLEX16Array(UINT ndim, ...)
{
    enum { maxdim = 16 };
    if(ndim > maxdim)
    {
        print_err("Error %s: dimension input %d is greater than max dimension %d", __func__, ndim, maxdim);
        return NULL;
    }
    va_list ap;
    UINT dims[ndim];
    UINT dim;
    
    va_start(ap, ndim);
    for(dim = 0; dim<ndim; ++dim)
        dims[dim] = va_arg(ap, UINT);
    va_end(ap);
    return CreateCOMPLEX16Array_comm(dims, ndim);
}

/* Destroy */
void DestroyCOMPLEX16Array(COMPLEX16Array* array)
{
    DestroyUINTVector( array->dimLength );
    free( array->data );
    array->data = NULL;
    free( array );
    array = NULL;
}


/* REAL8TimeSeries */
REAL8TimeSeries *CreateREAL8TimeSeries (const REAL8 epoch, REAL8 deltaT, UINT length)
{
    REAL8TimeSeries *r8ts;
    REAL8Vector *r8data;
    
    r8ts = (REAL8TimeSeries *)malloc(sizeof(*r8ts));
    r8data = CreateREAL8Vector(length);
    if(!r8ts || !r8data) {
        free(r8ts);
        r8ts = NULL;
        DestroyREAL8Vector (r8data);
        return NULL;
    }
    
    r8ts->epoch = epoch;
    r8ts->deltaT = deltaT;
    r8ts->data = r8data;
    
    return r8ts;
}

/* Destroy */
void DestroyREAL8TimeSeries( REAL8TimeSeries * series )
{
    if(series)
        DestroyREAL8Vector(series->data);
    free(series);
    series = NULL;
}

/* COMPLEX16TimeSeries */
COMPLEX16TimeSeries *CreateCOMPLEX16TimeSeries (const REAL8 epoch, REAL8 deltaT, UINT length)
{
    COMPLEX16TimeSeries *c16ts;
    COMPLEX16Vector *c16data;
    
    c16ts = (COMPLEX16TimeSeries *)malloc(sizeof(*c16ts));
    c16data = CreateCOMPLEX16Vector(length);
    if(!c16ts || !c16data) {
        free(c16ts);
        c16ts = NULL;
        DestroyCOMPLEX16Vector (c16data);
        return NULL;
    }
    
    c16ts->epoch = epoch;
    c16ts->deltaT = deltaT;
    c16ts->data = c16data;
    
    return c16ts;
}

/* Destroy */
void DestroyCOMPLEX16TimeSeries( COMPLEX16TimeSeries * series )
{
    if(series)
        DestroyCOMPLEX16Vector(series->data);
    free(series);
    series = NULL;
}

/* REAL8FrequencySeries */
REAL8FrequencySeries *CreateREAL8FrequencySeries(const REAL8 epoch, REAL8 f0, REAL8 deltaF, UINT length)
{
    REAL8FrequencySeries *r8fs;
    REAL8Vector *r8data;
    
    r8fs = (REAL8FrequencySeries *)malloc(sizeof(*r8fs));
    r8data = CreateREAL8Vector(length);
    if(!r8fs || !r8data)
    {
        free(r8fs);
        r8fs = NULL;
        DestroyREAL8Vector(r8data);
        return NULL;
    }
    
    r8fs->epoch = epoch;
    r8fs->deltaF = deltaF;
    r8fs->f0 = f0;
    r8fs->data = r8data;
    
    return r8fs;
}

/* Destroy */
void DestroyREAL8FrequencySeries( REAL8FrequencySeries * series )
{
    if (series)
    {
        DestroyREAL8Vector(series->data);
    }
    free(series);
    series = NULL;
}

/* COMPLEX16FrequencySeries */
COMPLEX16FrequencySeries *CreateCOMPLEX16FrequencySeries(const REAL8 epoch, REAL8 f0, REAL8 deltaF, UINT length)
{
    COMPLEX16FrequencySeries *c16fs;
    COMPLEX16Vector *c16data;
    
    c16fs = (COMPLEX16FrequencySeries *)malloc(sizeof(*c16fs));
    c16data = CreateCOMPLEX16Vector(length);
    if(!c16fs || !c16data)
    {
        free(c16fs);
        c16fs = NULL;
        DestroyCOMPLEX16Vector(c16data);
        return NULL;
    }
    
    c16fs->epoch = epoch;
    c16fs->deltaF = deltaF;
    c16fs->f0 = f0;
    c16fs->data = c16data;
    
    return c16fs;
}

/* Destroy */
void DestroyCOMPLEX16FrequencySeries( COMPLEX16FrequencySeries * series )
{
    if (series)
    {
        DestroyCOMPLEX16Vector(series->data);
    }
    free(series);
    series = NULL;
}

/* gsl Vector Array */
/* valid when ndim >= 2 */
VectorArray *CreateVectorArray_comm(UINT *dims, UINT ndim)
{
    if (ndim <= 1)
    {
        print_err("Error %s: The dimensions of VectorArray is at least 2, not %d\n", __func__, ndim);
        return NULL;
    }
    VectorArray *VA;
    VA = (VectorArray *)malloc(sizeof(*VA));
    UINTVector *dimvec;
    dimvec = CreateUINTVector(ndim - 1);
    memcpy(dimvec->data, dims, (ndim - 1) * sizeof(*dimvec->data));
    UINT vsize = 1;
    UINT i;
    for(i=0; i < ndim-1; ++i)
    {
        vsize *= dims[i];
    }
    gsl_vector **vector;
    vector = (gsl_vector **)malloc(vsize * sizeof(gsl_vector *));
    VA->size = vsize;
    VA->dimLength = dimvec;
    UINT nv = dims[ndim - 1];
    for (i=0 ; i<vsize; ++i)
    {
        vector[i] = gsl_vector_alloc(nv);
    }
    VA->vlength = nv;
    VA->data = vector;
    return VA;
}

VectorArray *CreateVectorArray(UINT ndim, ...)
{
    UINT dims[ndim];
    UINT i;
    
    va_list ap;
    va_start(ap, ndim);
    for(i=0; i < ndim; ++i)
    {
        dims[i] = va_arg(ap, UINT);
    }
    va_end(ap);
    return CreateVectorArray_comm(dims, ndim);
}

void DestroyVectorArray(VectorArray *varr)
{
    UINT i;
    for (i=0; i < varr->size; ++i)
    {
        gsl_vector_free(varr->data[i]);
    }
    DestroyUINTVector( varr->dimLength );
    free( varr->data );
    varr->data = NULL;
    free( varr );
    varr = NULL;
}

/* gsl Matrix Array */
/* valid on ndim >= 3 */
MatrixArray *CreateMatrixArray_comm(UINT *dims, UINT ndim)
{
    if (ndim <= 2)
    {
        print_err("Error %s: The dimensions of MatrixArray is at least 3, not %d\n", __func__, ndim);
        return NULL;
    }
    MatrixArray *MA;
    MA = (MatrixArray *)malloc(sizeof(*MA));
    UINTVector *dimvec;
    dimvec = CreateUINTVector(ndim - 2);
    memcpy(dimvec->data, dims, (ndim - 2) * sizeof(*dimvec->data));
    UINT msize = 1;
    UINT i;
    for(i=0; i<ndim - 2; i++)
    {
        msize *= dims[i];
    }
    gsl_matrix **matrix;
    matrix = (gsl_matrix **)malloc(msize * sizeof(gsl_matrix *));
    MA->size = msize;
    MA->dimLength = dimvec;
    UINT n1 = dims[ndim-2], n2 = dims[ndim-1];
    UINT *mdim;
    mdim = (UINT *)malloc(2 * sizeof(UINT));
    mdim[0] = n1;
    mdim[1] = n2;
    for (i=0; i < msize; ++i)
    {
        matrix[i] = gsl_matrix_alloc(n1, n2);
    }
    MA->mdim = mdim;
    MA->data = matrix;
    return MA;
}

MatrixArray *CreateMatrixArray(UINT ndim, ...)
{
    UINT dims[ndim];
    UINT i;
    
    va_list ap;
    va_start(ap, ndim);
    for(i=0; i < ndim; ++i)
    {
        dims[i] = va_arg(ap, UINT);
    }
    va_end(ap);
    return CreateMatrixArray_comm(dims, ndim);
}

void DestroyMatrixArray(MatrixArray *marr)
{
    UINT i;
    for (i=0; i < marr->size; ++i)
    {
        gsl_matrix_free(marr->data[i]);
    }
    DestroyUINTVector( marr->dimLength );
    free( marr->mdim );
    marr->mdim = NULL;
    free( marr->data );
    marr->data = NULL;
    free( marr );
    marr = NULL;
}


/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */
/** Vector Cut **/
REAL8Vector *CutREAL8Vector(REAL8Vector *vec,
                            size_t first,
                            size_t length)
{
    REAL8Vector *new = CreateREAL8Vector(length);
    if(!new)
    {
        return NULL;
    }
    memcpy(new->data, vec->data + first, length * sizeof(*new->data));
    
    return new;
}

COMPLEX16Vector *CutCOMPLEX16Vector(COMPLEX16Vector *vec,
                                    size_t first,
                                    size_t length)
{
    COMPLEX16Vector *new = CreateCOMPLEX16Vector(length);
    if(!new)
    {
        return NULL;
    }
    memcpy(new->data, vec->data + first, length * sizeof(*new->data));
    
    return new;
}


/** TimeSeries Cut **/
REAL8TimeSeries *CutREAL8TimeSeries(const REAL8TimeSeries *series,
                                    size_t first,
                                    size_t length)
{
    REAL8TimeSeries *new;
    REAL8Vector *vec;
    
    new = malloc(sizeof(*new));
    vec = CutREAL8Vector(series->data, first, length);
    if(!new || !vec)
    {
        free(new);
        DestroyREAL8Vector(vec);
        return NULL;
    }
    
    *new = *series;
    new->data = vec;
    new->epoch = new->epoch + first * new->deltaT;
    return new;
}

COMPLEX16TimeSeries *CutCOMPLEX16TimeSeries(const COMPLEX16TimeSeries *series,
                                            size_t first,
                                            size_t length)
{
    COMPLEX16TimeSeries *new;
    COMPLEX16Vector *vec;
    
    new = malloc(sizeof(*new));
    vec = CutCOMPLEX16Vector(series->data, first, length);
    if(!new || !vec)
    {
        free(new);
        DestroyCOMPLEX16Vector(vec);
        return NULL;
    }
    
    *new = *series;
    new->data = vec;
    new->epoch = new->epoch + first * new->deltaT;
    return new;
}


/** < REAL8Array > **/
/* Array element: for [i,j,k],
 * get corresponding dim from Array->dimLength [0],[1],[2],
 * data = Array->data[i + j*di + k*di*dj]
 */

static inline UINT Cood2Index(UINT *cood, UINT *shape, UINT ndim)
{
    UINT index = 0, i;
    UINT cum = 1;
    for(i=0; i < ndim; ++i)
    {
        cum *= i == 0 ? 1 : shape[ndim - i];
        index += cood[ndim - i - 1] * cum;
    }
    return index;
}

static inline INT Index2Cood(UINT **cood ,UINT index, UINT *shape, UINT ndim, UINT size)
{
    // The length of cood should be ndim, be careful!!
    UINT cum = size, i, tmp;
    for(i=0; i < ndim - 1; ++i)
    {
        cum /= shape[i];
        tmp = index / cum;
        (*cood)[i] = tmp;
        index -= cum * tmp;
    }
    (*cood)[ndim - 1] = index;
    return CEV_SUCCESS;
}

/* get value for REAL8Array */
INT getArrayValue_comm(REAL8 *val, REAL8Array *arr, UINT *cood)
{
    UINT index;
    index = Cood2Index(cood, arr->dimLength->data, arr->dimLength->length);
    *val = arr->data[index];
    return CEV_SUCCESS;
}

INT getArrayValue(REAL8 *val, REAL8Array *arr, UINT ndim,...)
{
    enum { maxdim = 16 };
    if(ndim > maxdim)
    {
        print_err("Error %s: dimension input %d is greater than max dimension %d", __func__, ndim, maxdim);
        return CEV_FAILURE;
    }

    UINT dims[ndim];
    UINT i;
    
    va_list ap;
    va_start(ap, ndim);
    for(i=0; i < ndim; ++i)
    {
        dims[i] = va_arg(ap, UINT);
    }
    va_end(ap);

    if (ndim != arr->dimLength->length)
    {
        print_err("Error %s: The shape of input cood is not equal to the shape of the input array.\n", __func__);
        return CEV_FAILURE;
    }
    return getArrayValue_comm(val, arr, dims);
}

/* set value for REAL8Array */
void setArrayValue_comm(REAL8 val, REAL8Array *arr, UINT *cood)
{
    UINT index;
    index = Cood2Index(cood, arr->dimLength->data, arr->dimLength->length);
    arr->data[index] = val;
}

INT setArrayValue(REAL8 val, REAL8Array *arr, UINT ndim,...)
{
    enum { maxdim = 16 };
    if(ndim > maxdim)
    {
        print_err("Error %s: dimension input %d is greater than max dimension %d", __func__, ndim, maxdim);
        return CEV_FAILURE;
    }
    UINT dims[ndim];
    UINT i;
    
    va_list ap;
    va_start(ap, ndim);
    for(i=0; i < ndim; ++i)
    {
        dims[i] = va_arg(ap, UINT);
    }
    va_end(ap);

    if (ndim != arr->dimLength->length)
    {
        print_err("Error %s: The shape of input cood is not equal to the shape of the input array.\n", __func__);
        return CEV_FAILURE;
    }

    setArrayValue_comm(val, arr, dims);
    return CEV_SUCCESS;
}

/* set value for VectorArray */
void setVectorArrayValue_comm(REAL8 val, VectorArray *varr, UINT *cood)
{
    UINT index;
    index = Cood2Index(cood, varr->dimLength->data, varr->dimLength->length);
    gsl_vector_set(varr->data[index], cood[varr->dimLength->length], val);
}

INT setVectorArrayValue(REAL8 val, VectorArray *varr, UINT ndim, ...)
{
    UINT dims[ndim];
    UINT i;
    
    va_list ap;
    va_start(ap, ndim);
    for(i=0; i < ndim; ++i)
    {
        dims[i] = va_arg(ap, UINT);
    }
    va_end(ap);
    
    if (ndim != varr->dimLength->length + 1 || dims[ndim-1] >= varr->vlength)
    {
        print_err("Error %s: Incorrect index for VectorArray.\n", __func__);
        return CEV_FAILURE;
    }
    setVectorArrayValue_comm(val, varr, dims);
    return CEV_SUCCESS;
}

/* get value from VectorArray */
INT getVectorArrayValue_comm(REAL8 *val, VectorArray *varr, UINT *cood)
{
    UINT index;
    index = Cood2Index(cood, varr->dimLength->data, varr->dimLength->length);
    *val = gsl_vector_get(varr->data[index], cood[varr->dimLength->length]);
    return CEV_SUCCESS;
}

INT getVectorArrayValue(REAL8 *val, VectorArray *varr, UINT ndim, ...)
{
    UINT dims[ndim];
    UINT i;
    
    va_list ap;
    va_start(ap, ndim);
    for(i=0; i < ndim; ++i)
    {
        dims[i] = va_arg(ap, UINT);
    }
    va_end(ap);
    
    if (ndim != varr->dimLength->length + 1 || dims[varr->dimLength->length] >= varr->vlength)
    {
        print_err("Error %s: Incorrect index for VectorArray.\n", __func__);
        return CEV_FAILURE;
    }
    return getVectorArrayValue_comm(val, varr, dims);
}


/* set value for MatrixArray */
void setMatrixArrayValue_comm(REAL8 val, MatrixArray *marr, UINT *cood)
{
    UINT index;
    index = Cood2Index(cood, marr->dimLength->data, marr->dimLength->length);
    //print_err("IDX : %d, cood = [%d, %d]\n", index, cood[marr->dimLength->length], cood[marr->dimLength->length + 1]);
    gsl_matrix_set(marr->data[index], cood[marr->dimLength->length], cood[marr->dimLength->length + 1], val);
}

INT setMatrixArrayValue(REAL8 val, MatrixArray *marr, UINT ndim, ...)
{
    UINT dims[ndim];
    UINT i;
    
    va_list ap;
    va_start(ap, ndim);
    for(i=0; i < ndim; ++i)
    {
        dims[i] = va_arg(ap, UINT);
    }
    va_end(ap);
    
    if (ndim != marr->dimLength->length + 2 || dims[ndim-2] >= marr->mdim[0] || dims[ndim-1] >= marr->mdim[1])
    {
        print_err("Error %s: Incorrect index for MatrixArray.\n", __func__);
        return CEV_FAILURE;
    }
    setMatrixArrayValue_comm(val, marr, dims);
    return CEV_SUCCESS;
}

/* get value from MatrixArray */
INT getMatrixArrayValue_comm(REAL8 *val, MatrixArray *marr, UINT *cood)
{
    UINT index;
    index = Cood2Index(cood, marr->dimLength->data, marr->dimLength->length);
    *val = gsl_matrix_get(marr->data[index], cood[marr->dimLength->length], cood[marr->dimLength->length + 1]);
    return CEV_SUCCESS;
}

INT getMatrixArrayValue(REAL8 *val, MatrixArray *marr, UINT ndim, ...)
{
    UINT dims[ndim];
    UINT i;
    
    va_list ap;
    va_start(ap, ndim);
    for(i=0; i < ndim; ++i)
    {
        dims[i] = va_arg(ap, UINT);
    }
    va_end(ap);
    
    if (ndim != marr->dimLength->length + 2 || dims[marr->dimLength->length] >= marr->mdim[0] || dims[marr->dimLength->length + 1] >= marr->mdim[1])
    {
        print_err("Error %s: Incorrect index for MatrixArray.\n", __func__);
        return CEV_FAILURE;
    }
    return getMatrixArrayValue_comm(val, marr, dims);
}

/* MatrixArray Transpose */
INT MatrixArrayTranspose(MatrixArray *marr)
{
    if (!marr)
    {
        print_err("Error %s: The input MatrixArray is empty!\n", __func__);
        return CEV_FAILURE;
    }
    UINT i;
    for(i=0; i< marr->size; ++i)
    {
        gsl_matrix_transpose(marr->data[i]);
    }
    SWAP(marr->mdim[0], marr->mdim[1]);
    return CEV_SUCCESS;
}

/* For array with dimLength >= 3, we can convert it to a gsl matrix array.
 * e.g. Array(2,6,3,4) would be converted to a MatrixArray(2,6),
 *      where each element is a (3,4) gsl matrix.
 */
INT Array2gslMatrixArray(REAL8Array *arr, MatrixArray **marr)
{
    MatrixArray *MA;
    MA = CreateMatrixArray_comm(arr->dimLength->data, arr->dimLength->length);
    if(!MA)
    {
        print_err("Error %s: Cannot create MatrixArray.\n", __func__);
        return CEV_FAILURE;
    }
    // print_err("mdim : [%d %d]\n", MA->mdim[0], MA->mdim[1]);
    UINT idx, midx;
    UINT *cood;
    //UINT *mcood; // debug
    cood = (UINT *)malloc(arr->dimLength->length * sizeof(UINT));
    //mcood = (UINT *)malloc(MA->dimLength->length * sizeof(UINT));
    for(idx=0; idx < arr->size; ++idx)
    {
        // Be Careful
        Index2Cood( &cood, idx, arr->dimLength->data, arr->dimLength->length, arr->size);
        //PRINT_LEN(cood, arr->dimLength->length, " %d ");
        midx = Cood2Index(cood, MA->dimLength->data, MA->dimLength->length);
        setMatrixArrayValue_comm(arr->data[idx], MA, cood);
    }
    *marr = MA;
    free(cood);
    cood = NULL;
    return CEV_SUCCESS;
}

INT gslMatrixArray2Array(MatrixArray *marr, REAL8Array **arr)
{
    REAL8Array *array;
    UINT *dims, mlength;
    mlength = marr->dimLength->length;
    dims = (UINT *)malloc((mlength + 2) * sizeof(UINT));
    memcpy(marr->dimLength->data, dims, mlength * sizeof(UINT));
    dims[mlength] = marr->mdim[0];
    dims[mlength + 1] = marr->mdim[1];
    array = CreateREAL8Array_comm(dims, mlength + 2);
    if (!array)
    {
        print_err("Error %s: Cannot create REAL8Array.\n", __func__);
        return CEV_FAILURE;
    }
    UINT idx;
    UINT *cood;
    cood = (UINT *)malloc(array->dimLength->length * sizeof(UINT));
    for(idx=0; idx < array->size; ++idx)
    {
        Index2Cood( &cood, idx, array->dimLength->data, array->dimLength->length, array->size);
        getMatrixArrayValue_comm(array->data + idx, marr, cood);
    }
    (*arr) = array;
    free(cood);
    cood = NULL;
    return CEV_SUCCESS;
}

/* Math */
COMPLEX16 CX16polar(REAL8 r, REAL8 phi)
{
    REAL8 re,im;
    re = r*GET_COS(phi);
    im = r*GET_SIN(phi);
    
    
    COMPLEX16 z=re+I*im;
    
    return z;
}

/* String */
int StringNCaseCompare(const char *s1, const char *s2, size_t n)
{
    
    /* ctype replacements w/o locale */
    const char upper_chars[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    const char lower_chars[] = "abcdefghijklmnopqrstuvwxyz";
    
    int c1 = 0, c2 = 0;
    
    /* based on implementation of strncmp() in glibc */
    while (n-- > 0) {
        
        c1 = (s1 == NULL) ? 0 : *s1++;
        c2 = (s2 == NULL) ? 0 : *s2++;
        
        /* convert c1 to lower case */
        if (c1) {
            char *p = strchr(upper_chars, c1);
            if (p) {
                int offset = p - upper_chars;
                c1 = lower_chars[offset];
            }
        }
        
        /* convert c2 to lower case */
        if (c2) {
            char *p = strchr(upper_chars, c2);
            if (p) {
                int offset = p - upper_chars;
                c2 = lower_chars[offset];
            }
        }
        
        /* compare characters */
        if (c1 == '\0' || c1 != c2) {
            return c1 - c2;
        }
        
    }
    
    return c1 - c2;
    
}


char * StringCaseSubstring(const char *haystack, const char *needle)
{
    size_t haystack_length;
    size_t needle_length = strlen(needle);
    
    /* return haystack if needle is empty */
    if (needle_length == 0)
        return (char *)(intptr_t)(haystack);
    
    haystack_length = strlen(haystack);
    while (needle_length <= haystack_length) {
        if (StringNCaseCompare(haystack, needle, needle_length) == 0)
            return (char *)(intptr_t)(haystack);
        --haystack_length;
        ++haystack;
    }
    
    /* needle not found in haystack */
    return NULL;
}

int SphericalToCartesian(REAL8 qCart[],
                         REAL8 pCart[],
                         const REAL8 qSph[],
                         const REAL8 pSph[])
{
    
    REAL8 r;
    REAL8 pr, pTheta, pPhi;
    
    REAL8 x, y, z;
    REAL8 pX, pY, pZ;
    
    r      = qSph[0];
    pr     = pSph[0];
    pTheta = pSph[1];
    pPhi   = pSph[2];
    
    x = r;
    y = 0.0;
    z = 0.0;
    
    pX = pr;
    pY = pPhi / r;
    pZ = - pTheta / r;
    
    /* Copy these values to the output vectors */
    qCart[0] = x;
    qCart[1] = y;
    qCart[2] = z;
    pCart[0] = pX;
    pCart[1] = pY;
    pCart[2] = pZ;
    
    return CEV_SUCCESS;
}

int CartesianToSpherical(REAL8 qSph[],
                         REAL8 pSph[],
                         const REAL8 qCart[],
                         const REAL8 pCart[])
{
    
    REAL8 r;
    REAL8 pr, pTheta, pPhi;
    
    REAL8 x; //, y, z;
    REAL8 pX, pY, pZ;
    
    x  = qCart[0];
    //y  = qCart[1];
    //z  = qCart[2];
    pX = pCart[0];
    pY = pCart[1];
    pZ = pCart[2];
    
    r  = x;
    pr = pX;
    pTheta = - r * pZ;
    pPhi   =   r * pY;
    
    /* Fill the output vectors */
    qSph[0] = r;
    qSph[1] = CST_PI_2;
    qSph[2] = 0.0;
    pSph[0] = pr;
    pSph[1] = pTheta;
    pSph[2] = pPhi;
    
    return CEV_SUCCESS;
}

