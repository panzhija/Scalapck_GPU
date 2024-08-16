#include <stdio.h>
#include <stdlib.h>
//#include <hipblas/hipblas.h>
#include <hipblas.h>
#include <hip/hip_runtime.h>
#include "magma_v2.h"

#ifndef		__INCLUDE_PB_HIPBLAS_H__
#define		__INCLUDE_PB_HIPBLAS_H__

#define		NG	0
#define		OK	1

#define 	checkhipErrors(err)	__checkhipErrors( err, #err, __FILE__, __LINE__ )
#define		checkhipblasErrors(err)	__checkhipblasErrors( err, #err, __FILE__, __LINE__ )
#define     checkMagmaErrors(err) __checkMagmaErrors( err, #err, __FILE__, __LINE__ ) 

typedef void  (*FUNCPTR)( char* transa, char* transb,
		                  int m, int n, int k, char* alpha, char* beta,
						  char* A, int lda, char* B, int ldb,
						  char* C, int idc );

static hipblasOperation_t hipblas_trans_const(char* key)
{
	switch(*key)
	{
		case 'N':	return HIPBLAS_OP_N;
		case 'T':	return HIPBLAS_OP_T;
		case 'C':	return HIPBLAS_OP_C;
	}
	return HIPBLAS_OP_N;
}

static magma_trans_t magma_trans_const(char* key)
{
	switch(*key)
	{
		case 'N':	return MagmaNoTrans;
		case 'T':	return MagmaTrans;
		case 'C':	return MagmaConjTrans;
	}
	return MagmaNoTrans;
}

static char* hipblasGetErrorString(hipblasStatus_t status)
{
	switch(status)
	{
		case HIPBLAS_STATUS_NOT_INITIALIZED:	return "HIPBLAS_STATUS_NOT_INITIALIZED";
		case HIPBLAS_STATUS_ALLOC_FAILED:	return "HIPBLAS_STATUS_ALLOC_FAILED";
		case HIPBLAS_STATUS_INVALID_VALUE:	return "HIPBLAS_STATUS_INVALID_VALUE";
		case HIPBLAS_STATUS_ARCH_MISMATCH:	return "HIPBLAS_STATUS_ARCH_MISMATCH";
		case HIPBLAS_STATUS_MAPPING_ERROR:	return "HIPBLAS_STATUS_MAPPING_ERROR";
		case HIPBLAS_STATUS_EXECUTION_FAILED:return "HIPBLAS_STATUS_EXECUTION_FAILED";
		case HIPBLAS_STATUS_INTERNAL_ERROR:	return "HIPBLAS_STATUS_INTERNAL_ERROR";
		case HIPBLAS_STATUS_NOT_SUPPORTED:	return "HIPBLAS_STATUS_NOT_SUPPORTED";
	}
	return "unknown error";
}

static inline int __checkhipErrors(int err, const char* func, const char* file, const int line)
{
	if( hipSuccess != err )
	{
		printf("hip error at %s:%d code=%d(%s) \"%s\" \n",
				file, line, (int)err, hipGetErrorString(err), func);
		return NG;
	}
	return OK;
}

static inline int __checkhipblasErrors(int err, const char* func, const char* file, const int line)
{
	if( HIPBLAS_STATUS_SUCCESS != err )
	{
		char* hipblasStatus = hipblasGetErrorString(err);
		printf("HIPBLAS error at %s:%d code=%d(%s) \"%s\" \n",
				file, line, (int)err, hipblasStatus, func);
		return NG;
	}
	return OK;
}

static inline int __checkMagmaErrors(int err, const char* func, const char* file, const int line)
{
	if (err != MAGMA_SUCCESS)
    {
		printf("MAGMA error at %s:%d code=%s \"%s\" \n",
				file, line, magma_strerror(err), func);
		return NG;
	}
	return OK;
}

#endif		// __INCLUDE_PB_HIPBLAS_H__