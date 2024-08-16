#include "../pblas.h"
#include "../PBpblas.h"
#include "../PBtools.h"
#include "../PBblacs.h"
#include "../PBblas.h"
#include "../PB_hipBLAS.h"
//#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef __STDC__
void hipPB_CpgemmAB( PBTYP_T * TYPE, char * DIRECA, char * DIRECB,
                  char * TRANSA, char * TRANSB, int M, int N, int K,
                  char * ALPHA, char * A, int IA, int JA, int * DESCA,
                  char * B, int IB, int JB, int * DESCB, char * BETA,
                  char * C, int IC, int JC, int * DESCC )
#else
void hipPB_CpgemmAB( TYPE, DIRECA, DIRECB, TRANSA, TRANSB, M, N, K, ALPHA,
                  A, IA, JA, DESCA, B, IB, JB, DESCB, BETA, C, IC, JC,
                  DESCC )
/*
*  .. Scalar Arguments ..
*/
   char           * DIRECA, * DIRECB, * TRANSA, * TRANSB;
   int            IA, IB, IC, JA, JB, JC, K, M, N;
   char           * ALPHA, * BETA;
   PBTYP_T        * TYPE;
/*
*  .. Array Arguments ..
*/
   int            * DESCA, * DESCB, * DESCC;
   char           * A, * B, * C;
#endif
{
/*
*  .. Local Scalars ..
*/
   char           Aroc, Broc, TrA, TrB, * one, * tbeta, * zero;
   int            ABrocs, Abufld, AcurrocR, Afr, Afwd, AiD, AiR, AiiD, AiiR,
                  AinbD, AinbR, Ainb1D, Ainb1R, AisR, AkkR, Ald, AmyprocD,
                  AmyprocR, AnbD, AnbR, AnpD, AnpR, AnprocsD, AnprocsR, Aoff,
                  ArocD, ArocR, AsrcR, Bbufld, BcurrocR, Bfr, Bfwd, BiD, BiR,
                  BiiD, BiiR, BinbD, BinbR, Binb1D, Binb1R, BisR, BkkR, Bld,
                  BmyprocD, BmyprocR, BnbD, BnbR, BnpD, BnpR, BnprocsD,
                  BnprocsR, Boff, BrocD, BrocR, BsrcR, Ccol, Cii, Cimb1, Cinb1,
                  Cjj, Cld, Cmb, Cmp, Cnb, Cnq, Crow, WAfr, WAsum, WBfr, WBsum,
                  Wkbb=0, ctxt, k, kb, kbb, lcmb, maxp, maxpm1, maxq, mycol,
                  myrow, ncpq, nota, notb, npcol, npq=0, nprow, nrpq, p=0, q=0,
                  size, tmp;
   GEMM_T         gemm;
/*
*  .. Local Arrays ..
*/
   PB_VM_T        VM;
   int            Cd0[DLEN_], DBUFA[DLEN_], DBUFB[DLEN_], WAd0[DLEN_],
                  WBd0[DLEN_];
   char           * Abuf = NULL, * Bbuf = NULL, * Cptr = NULL, * WA = NULL,
                  * WB   = NULL;
/* ..
*  .. Executable Statements ..
*
*/
   // struct timeval start, end;
   // double timeuse = 0.0;
   if(!checkMagmaErrors(magma_init()))  return;
   char           * WA_BUF = NULL, * WB_BUF = NULL;
   int buf_tag = 1;
   int W_ir = 4;
   int kbb_old;
   int WA_size, WB_size;
   int d_A_size, d_B_size, d_C_size;
   magma_queue_t myqueue;
   //magma_queue_t myqueue_buf;
   // hipStream_t mystream;
   // hipblasHandle_t handle;
   double* d_A = NULL;
   double* d_B = NULL;
   double* d_C = NULL;
   magma_device_t device_id;

   // if(!checkhipblasErrors(hipblasCreate(&handle)))  return;
   // hipStreamCreate(&mystream);
   // hipblasSetStream(handle, mystream);

   magma_getdevice(&device_id);
   magma_queue_create(device_id, &myqueue);
   //magma_queue_create(device_id, &myqueue_buf);

/*
*  Retrieve process grid information
*/
   Cblacs_gridinfo( ( ctxt = DESCC[CTXT_] ), &nprow, &npcol, &myrow, &mycol );

   nota = ( ( TrA = Mupcase( TRANSA[0] ) ) == CNOTRAN );
   notb = ( ( TrB = Mupcase( TRANSB[0] ) ) == CNOTRAN );
   TrA  = ( ( TrA == CCOTRAN ) ? CCONJG : CNOCONJG );
   TrB  = ( ( TrB == CCOTRAN ) ? CCONJG : CNOCONJG );

   size = TYPE->size;
/*
*  Retrieve local information for sub( A ), sub( B ) and sub( C )
*/
   if( nota )
   {
      AiR   = JA;          Aroc = CCOLUMN;     AnprocsR = npcol;
      AinbR = DESCA[INB_]; AnbR = DESCA[NB_ ]; AsrcR    = DESCA[CSRC_];
   }
   else
   {
      AiR   = IA;          Aroc = CROW;        AnprocsR = nprow;
      AinbR = DESCA[IMB_]; AnbR = DESCA[MB_ ]; AsrcR    = DESCA[RSRC_];
   }

   if( notb )
   {
      BiR   = IB;          Broc = CROW;        BnprocsR = nprow;
      BinbR = DESCB[IMB_]; BnbR = DESCB[MB_ ]; BsrcR    = DESCB[RSRC_];
   }
   else
   {
      BiR   = JB;          Broc = CCOLUMN;     BnprocsR = npcol;
      BinbR = DESCB[INB_]; BnbR = DESCB[NB_ ]; BsrcR    = DESCB[CSRC_];
   }
/*
*  Retrieve sub( C )'s local information: Aii, Ajj, Arow, Acol ...
*/
//PB_Cdescribe返回子矩阵的全局描述符。该例程还计算从IA, JA指向的入口全局开始的子矩阵的起始局部索引II, JJ。
//这个例程返回拥有全局索引I、J的矩阵项的进程的网格坐标，即PROW和PCOL。真正的全局块大小IMB, INB, MB和NB也会返回。
   PB_Cdescribe( M, N, IC, JC, DESCC, nprow, npcol, myrow, mycol, &Cii, &Cjj,
                 &Cld, &Cimb1, &Cinb1, &Cmb, &Cnb, &Crow, &Ccol, Cd0 );
//PB_Cnumroc返回矩阵处理的本地行/列数,如果我们从全局索引I开始给出N行/列，PROC将会得到。
   Cmp = PB_Cnumroc( M, 0, Cimb1, Cmb, myrow, Crow, nprow );
   Cnq = PB_Cnumroc( N, 0, Cinb1, Cnb, mycol, Ccol, npcol );
/*
*  When sub( A ) and sub( B ) do not span more than one process row or column,
*  there is no need to pack the data.
*/
//PB_Cspan返回1，如果行(resp。I:I+N-1跨出多个进程行(响应。列)，否则为0。
   if( !( PB_Cspan( K, AiR, AinbR, AnbR, AsrcR, AnprocsR ) ) &&
       !( PB_Cspan( K, BiR, BinbR, BnbR, BsrcR, BnprocsR ) ) )
   {
      //PB_CInV返回一个指向数组的指针，该数组包含一个仅输入的一维子向量，该子向量在DESCA描述的子矩阵的行或列上复制。
      //子向量是在这个例程的输入上指定的，只要可能就重用它。返回时，子向量由指向某些数据的指针、描述其布局的描述符数组和指示该局部数据是否已由该函数动态分配的逻辑值指定。
      //这个例程专门为传统的2级设计，如PBLAS操作，使用仅输入的向量，如PxGER, PxSYR…
      hipPB_CInV( TYPE, &TrA, COLUMN, M, N, Cd0, K, A, IA, JA, DESCA, &Aroc, &WA,
               WAd0, &WAfr );
      hipPB_CInV( TYPE, &TrB, ROW,    M, N, Cd0, K, B, IB, JB, DESCB, &Broc, &WB,
               WBd0, &WBfr );
/*
*  Perform the local update if I own some of sub( C )
*/
        // Allocate
      //   if(!checkhipErrors(hipMalloc((void**)&d_A, sizeof(double)*Cmp*K  )))  return;
      //   if(!checkhipErrors(hipMalloc((void**)&d_B, sizeof(double)*K*Cnq  )))  return;
      //   if(!checkhipErrors(hipMalloc((void**)&d_C, sizeof(double)*Cmp*Cnq)))  return;

      // //   // HtoD
      //   if(!checkhipblasErrors(hipblasSetMatrix(Cmp, K  , sizeof(double), WA, WAd0[LLD_], (double*)d_A, Cmp))) return;
      //   if(!checkhipblasErrors(hipblasSetMatrix(K  , Cnq, sizeof(double), WB, WBd0[LLD_], (double*)d_B, K  ))) return;
      //   if(!checkhipblasErrors(hipblasSetMatrix(Cmp, Cnq, sizeof(double), Mptr( C, Cii, Cjj, Cld, size ), Cld, (double*)d_C, Cmp))) return;

      // //   // Run
      //   if(!checkhipblasErrors(hipblasDgemm(handle, hipblas_trans_const(NOTRAN), hipblas_trans_const(NOTRAN),
      //                                       Cmp, Cnq, K,
      //                                       (double*)ALPHA, (double*)d_A, Cmp,
      //                                       (double*)d_B, K,
      //                                       (double*)BETA,  (double*)d_C, Cmp)))    return;

      // //   // DtoH
      //   if(!checkhipblasErrors(hipblasGetMatrix(Cmp, Cnq, sizeof(double), (double*)d_C, Cmp, Mptr( C, Cii, Cjj, Cld, size ), Cld))) return;

      // //   // Deallocate
      //   hipFree(d_A);
      //   hipFree(d_B);
      //   hipFree(d_C);
      //   hipblasDestroy(handle);
      //   hipStreamDestroy(mystream);
          TYPE->Fgemm( C2F_CHAR( NOTRAN ), C2F_CHAR( NOTRAN ), &Cmp, &Cnq, &K,
                              ALPHA, WA, &WAd0[LLD_], WB, &WBd0[LLD_], BETA, Mptr( C,
                              Cii, Cjj, Cld, size ), &Cld );
      if( WAfr ) free( WA );
      if( WBfr ) free( WB );
      return;
   }
/*
*  sub( A ) and sub( B ) span more than one process row or column.
*/
   Afwd  = ( Mupcase( DIRECA[0] ) == CFORWARD );
   Bfwd  = ( Mupcase( DIRECB[0] ) == CFORWARD );

   one   = TYPE->one; zero = TYPE->zero; tbeta = BETA; gemm = TYPE->Fgemm;
   //reture逻辑快的大小
   kb    = pilaenv_( &ctxt, C2F_CHAR( &TYPE->type ) );
/*
*  Compute local information for sub( A ) and sub( B )
*/
   if( nota )
   {
      AiD      = IA;          AinbD    = DESCA[IMB_]; AnbD     = DESCA[MB_];
      Ald      = DESCA[LLD_]; AmyprocD = myrow;       AmyprocR = mycol;
      AnprocsD = nprow;
      //PB_Cinfog2l计算在I, J指向的条目处全局开始的子矩阵的起始局部索引II, JJ。该例程返回拥有全局索引I, J的矩阵条目的进程在网格中的坐标，即PROW和PCOL。
      PB_Cinfog2l( IA, JA, DESCA, AnprocsD, AnprocsR, AmyprocD, AmyprocR,
                   &AiiD, &AiiR, &ArocD, &ArocR );
   }
   else
   {
      AiD      = JA;          AinbD    = DESCA[INB_]; AnbD     = DESCA[NB_];
      Ald      = DESCA[LLD_]; AmyprocD = mycol;       AmyprocR = myrow;
      AnprocsD = npcol;
      PB_Cinfog2l( IA, JA, DESCA, AnprocsR, AnprocsD, AmyprocR, AmyprocD,
                   &AiiR, &AiiD, &ArocR, &ArocD );
   }
   //pb_confirstnb返回第一个块的矩阵行或列的全局数量，如果从全局索引i开始给出N行或列。注意，如果N等于0，这个例程返回0。
   Ainb1D = PB_Cfirstnb( M, AiD, AinbD, AnbD );
   AnpD   = PB_Cnumroc( M, 0, Ainb1D, AnbD, AmyprocD, ArocD, AnprocsD );
   Ainb1R = PB_Cfirstnb( K, AiR, AinbR, AnbR );
   AisR   = ( ( AsrcR < 0 ) || ( AnprocsR == 1 ) );

   if( notb )
   {
      BiD      = JB;          BinbD    = DESCB[INB_]; BnbD     = DESCB[NB_];
      Bld      = DESCB[LLD_]; BmyprocD = mycol;       BmyprocR = myrow;
      BnprocsD = npcol;
      PB_Cinfog2l( IB, JB, DESCB, BnprocsR, BnprocsD, BmyprocR, BmyprocD,
                   &BiiR, &BiiD, &BrocR, &BrocD );
   }
   else
   {
      BiD      = IB;          BinbD    = DESCB[IMB_]; BnbD     = DESCB[MB_];
      Bld      = DESCB[LLD_]; BmyprocD = myrow;       BmyprocR = mycol;
      BnprocsD = nprow;
      PB_Cinfog2l( IB, JB, DESCB, BnprocsD, BnprocsR, BmyprocD, BmyprocR,
                   &BiiD, &BiiR, &BrocD, &BrocR );
   }
   Binb1D = PB_Cfirstnb( N, BiD, BinbD, BnbD );
   BnpD   = PB_Cnumroc( N, 0, Binb1D, BnbD, BmyprocD, BrocD, BnprocsD );
   Binb1R = PB_Cfirstnb( K, BiR, BinbR, BnbR );
   BisR   = ( ( BsrcR < 0 ) || ( BnprocsR == 1 ) );
/*
*  When sub( A ) is not replicated and backward pass on sub( A ), find the
*  virtual process q owning the last row or column of sub( A ).
*/
//当子(A)未被复制并向后传递子(A)时，找到拥有子(A)最后一行或最后一列的虚拟进程q。
   if( !( AisR ) && !( Afwd ) )
   {
      //PB_Cindxg2p计算具有全局索引IG指定的矩阵项的流程坐标。
      tmp = PB_Cindxg2p( K - 1, Ainb1R, AnbR, ArocR, ArocR, AnprocsR );
      q   = MModSub( tmp, ArocR, AnprocsR );
   }
/*
*  When sub( B ) is not replicated and backward pass on sub( B ), find the
*  virtual process p owning the last row or column of sub( B ).
*/
//当子(B)未被复制并向后传递子(B)时，找到拥有子(B)最后一行或最后一列的虚拟进程p。
   if( !( BisR ) && !( Bfwd ) )
   {
      tmp = PB_Cindxg2p( K - 1, Binb1R, BnbR, BrocR, BrocR, BnprocsR );
      p   = MModSub( tmp, BrocR, BnprocsR );
   }

   if( Cmp > 0 && Cnq > 0 ) Cptr = Mptr( C, Cii, Cjj, Cld, size );
/*
*  Allocate work space in process rows and columns spanned by sub( C )
*/
//PB_COutV返回一个指向数组的指针，该数组包含一个一维输出零子向量，该子向量在DESCA描述的子矩阵的行或列上复制。
//在返回时，子向量由指向某些数据的指针、描述其布局的描述符数组、指示该局部数据是否已由该函数动态分配的逻辑值、指定是否应进行和减的逻辑值指定。
//这个例程专门为传统的2级和3级PBLAS操作设计，使用仅输出向量，如PxTRMV或PxTRMM。
   hipPB_COutV( TYPE, COLUMN, NOINIT, M, N, Cd0, kb, &WA, WAd0, &WAfr, &WAsum, &WA_size);
   hipPB_COutV( TYPE, ROW,    NOINIT, M, N, Cd0, kb, &WB, WBd0, &WBfr, &WBsum, &WB_size);
/*
*  Loop over the virtual process grid induced by the sub( A ) and sub( B )
*/
//由子(A)和子(B)诱导的虚拟进程网格上的循环
//PB_Clcm计算并返回两个正整数M和N的最小公倍数(LCM)，实际上，该例程计算最大公约数(GCD)，并使用M*N = GCD*LCM的属性。
   lcmb   = PB_Clcm( ( maxp = ( BisR ? 1 : BnprocsR ) ) * BnbR,
                     ( maxq = ( AisR ? 1 : AnprocsR ) ) * AnbR );
   maxpm1 = maxp - 1;
/*
*  Find out process coordinates corresponding to first virtual process (p,q)
*/
//求出第一个虚拟进程(p,q)对应的进程坐标
   AcurrocR = ( AisR ? -1 : MModAdd( ArocR, q, AnprocsR ) );
//PB_Cg2lrem计算全局索引IG所指向的矩阵项的局部索引。请注意，当MYPROC不是拥有该条目的进程时，该例程返回与IG对应的最近的较大本地索引，就像例程PB_Cinfog2l一样。
   AkkR     = PB_Cg2lrem( AiR, AinbR, AnbR, AcurrocR, AsrcR, AnprocsR );
   AnpR     = PB_Cnumroc( K, 0, Ainb1R, AnbR, AcurrocR, ArocR, AnprocsR );

   BcurrocR = ( BisR ? -1 : MModAdd( BrocR, p, BnprocsR ) );
   BkkR     = PB_Cg2lrem( BiR, BinbR, BnbR, BcurrocR, BsrcR, BnprocsR );
   BnpR     = PB_Cnumroc( K, 0, Binb1R, BnbR, BcurrocR, BrocR, BnprocsR );
/*
*  Find out how many diagonals this virtual process (p,q) has
*/
//PB_CVMinit使用相对坐标进程(MRROW, MRCOL)拥有的m × n本地数组的信息初始化一个虚拟矩阵。
   PB_CVMinit( &VM, 0, BnpR, AnpR, Binb1R, Ainb1R, BnbR, AnbR, p, q,
               maxp, maxq, lcmb );
//PB_CVMnpq计算VM指定的虚拟mma中的对角线条目数。
   npq = PB_CVMnpq( &VM );

   // for hipBLAS (Initialize)
   if( Cmp > 0 && Cnq > 0 && kb > 0)
   {
      // Allocate
      // d_A_size = Cmp * kb;
      // d_B_size = kb * Cnq;
      // d_C_size = Cmp * Cnq;
      // if(!checkhipErrors(hipMalloc((void**)&d_A, sizeof(double) * d_A_size * 2)))  return;
      // if(!checkhipErrors(hipMalloc((void**)&d_B, sizeof(double) * d_B_size * 2)))  return;
      // if(!checkhipErrors(hipMalloc((void**)&d_C, sizeof(double) * d_C_size * 2)))  return;
      if(!checkMagmaErrors(magma_dmalloc(&d_A, Cmp * kb)))  return;
      if(!checkMagmaErrors(magma_dmalloc(&d_B, kb * Cnq)))  return;
      if(!checkMagmaErrors(magma_dmalloc(&d_C, Cmp * Cnq)))  return;
      WA_BUF = hipPB_Cmalloc( WA_size * size);
      WB_BUF = hipPB_Cmalloc( WB_size * size);
      // memset(WA_BUF,0,WA_size * size);
      // memset(WB_BUF,0,WB_size * size);
      if(!checkhipErrors(hipHostRegister(WA, sizeof(double)*WA_size, hipHostRegisterDefault)))  return;
      if(!checkhipErrors(hipHostRegister(WB, sizeof(double)*WB_size, hipHostRegisterDefault)))  return;
      if(!checkhipErrors(hipHostRegister(WA_BUF, sizeof(double)*WA_size, hipHostRegisterDefault)))  return;
      if(!checkhipErrors(hipHostRegister(WB_BUF, sizeof(double)*WB_size, hipHostRegisterDefault)))  return;
      // if(!checkhipErrors(hipMalloc((void**)&d_A, sizeof(double)*Cmp*kb)))  return;
      // if(!checkhipErrors(hipMalloc((void**)&d_B, sizeof(double)*kb*Cnq)))  return;
      // if(!checkhipErrors(hipMalloc((void**)&d_C, sizeof(double)*Cmp*Cnq)))    return;
      // HtoD
      //if(!checkhipblasErrors(hipblasSetMatrix(Cmp, Cnq, sizeof(double), Cptr, Cld, (double*)d_C, Cmp))) return;
      //magma_setvector(batch_size, sizeof(double *), h_C_array, 1, d_C_array, 1, myqueue);
      //magma_setmatrix(Cmp, Cnq, sizeof(double), Cptr, Cld, d_C, Cmp, myqueue);

   }

   for( k = 0; k < K; k += kb )
   {
      kbb = K - k; kbb = MIN( kbb, kb );
      while( Wkbb != kbb )
      {
/*
*  Ensure that the current virtual process (p,q) has something to contribute
*  to the replicated buffers WA and WB.
*/
//确保当前虚拟进程(p,q)对复制缓冲区WA和WB有贡献。
         while( npq == 0 )
         {
            if( ( Bfwd      && ( p == maxpm1 ) ) ||
                ( !( Bfwd ) && ( p == 0      ) ) )
               q = ( Afwd ? MModAdd1( q, maxq ) : MModSub1( q, maxq ) );
            p = ( Bfwd ? MModAdd1( p, maxp ) : MModSub1( p, maxp ) );

            AcurrocR = ( AisR ? -1 : MModAdd( ArocR, q, AnprocsR ) );
            AkkR     = PB_Cg2lrem(  AiR, AinbR,  AnbR, AcurrocR, AsrcR,
                                   AnprocsR );
            AnpR     = PB_Cnumroc( K, 0, Ainb1R, AnbR, AcurrocR, ArocR,
                                   AnprocsR );

            BcurrocR = ( BisR ? -1 : MModAdd( BrocR, p, BnprocsR ) );
            BkkR     = PB_Cg2lrem(  BiR, BinbR,  BnbR, BcurrocR, BsrcR,
                                   BnprocsR );
            BnpR     = PB_Cnumroc( K, 0, Binb1R, BnbR, BcurrocR, BrocR,
                                   BnprocsR );

            PB_CVMinit( &VM, 0, BnpR, AnpR, Binb1R, Ainb1R, BnbR, AnbR,
                        p, q, maxp, maxq, lcmb );
            npq = PB_CVMnpq( &VM );
         }
/*
*  Current virtual process (p,q) has something, find out how many rows or
*  columns could be used: ABrocs.
*/
         if( Wkbb == 0 ) { ABrocs = ( npq < kbb ? npq : kbb ); }
         else            { ABrocs = kbb - Wkbb; ABrocs = MIN( ABrocs, npq ); }
/*
*  Find out how many rows or columns of sub( A ) and sub( B ) are contiguous
*/
         PB_CVMcontig( &VM, &nrpq, &ncpq, &Boff, &Aoff );

         if( nota )
         {
/*
*  Compute the descriptor DBUFA for the buffer that will contained the packed
*  columns of sub( A ).
*/
            if( ( Afr = ( ncpq < ABrocs ) ) != 0 )
            {
/*
*  If columns of sub( A ) are not contiguous, then allocate the buffer and
*  pack the ABrocs columns of sub( A ).
*/
               Abufld = MAX( 1, AnpD );
               if( AisR || ( AmyprocR == AcurrocR ) )
               {
                  Abuf = hipPB_Cmalloc( AnpD * ABrocs * size );
                  PB_CVMpack( TYPE, &VM, COLUMN, &Aroc, PACKING, NOTRAN,
                              ABrocs, AnpD, one, Mptr( A, AiiD, AkkR, Ald,
                              size ), Ald, zero, Abuf, Abufld );
               }
            }
            else
            {
/*
*  Otherwise, re-use sub( A ) directly.
*/
               Abufld = Ald;
               if( AisR || ( AmyprocR == AcurrocR ) )
                     Abuf = Mptr( A, AiiD, AkkR + Aoff, Ald, size );
            }
            PB_Cdescset( DBUFA, M, ABrocs, Ainb1D, ABrocs, AnbD, ABrocs,
                         ArocD, AcurrocR, ctxt, Abufld );
         }
         else
         {
/*
*  Compute the descriptor DBUFA for the buffer that will contained the packed
*  rows of sub( A ).
*/
            if( ( Afr = ( ncpq < ABrocs ) ) != 0 )
            {
/*
*  If rows of sub( A ) are not contiguous, then allocate the buffer and
*  pack the ABrocs rows of sub( A ).
*/
               Abufld = ABrocs;
               if( AisR || ( AmyprocR == AcurrocR ) )
               {
                  Abuf = hipPB_Cmalloc( AnpD * ABrocs * size );
                  PB_CVMpack( TYPE, &VM, COLUMN, &Aroc, PACKING, NOTRAN,
                              ABrocs, AnpD, one, Mptr( A, AkkR, AiiD, Ald,
                              size ), Ald, zero, Abuf, Abufld );
               }
            }
            else
            {
/*
*  Otherwise, re-use sub( A ) directly.
*/
               Abufld = Ald;
               if( AisR || ( AmyprocR == AcurrocR ) )
                  Abuf = Mptr( A, AkkR + Aoff, AiiD, Ald, size );
            }
            PB_Cdescset( DBUFA, ABrocs, M, ABrocs, Ainb1D, ABrocs, AnbD,
                         AcurrocR, ArocD, ctxt, Abufld );
         }

         if( notb )
         {
/*
*  Compute the descriptor DBUFB for the buffer that will contained the packed
*  rows of sub( B ).
*/
            if( ( Bfr = ( nrpq < ABrocs ) ) != 0 )
            {
/*
*  If rows of sub( B ) are not contiguous, then allocate the buffer and
*  pack the ABrocs rows of sub( B ).
*/
               Bbufld = ABrocs;
               if( BisR || ( BmyprocR == BcurrocR ) )
               {
                  Bbuf = hipPB_Cmalloc( BnpD * ABrocs * size );
                  PB_CVMpack( TYPE, &VM, ROW,    &Broc, PACKING, NOTRAN,
                              ABrocs, BnpD, one, Mptr( B, BkkR, BiiD, Bld,
                              size ), Bld, zero, Bbuf, Bbufld );
               }
            }
            else
            {
/*
*  Otherwise, re-use sub( B ) directly.
*/
               Bbufld = Bld;
               if( BisR || ( BmyprocR == BcurrocR ) )
                  Bbuf = Mptr( B, BkkR + Boff, BiiD, Bld, size );
            }
            PB_Cdescset( DBUFB, ABrocs, N, ABrocs, Binb1D, ABrocs, BnbD,
                         BcurrocR, BrocD, ctxt, Bbufld );
         }
         else
         {
/*
*  Compute the descriptor DBUFB for the buffer that will contained the packed
*  columns of sub( B ).
*/
            if( ( Bfr = ( nrpq < ABrocs ) ) != 0 )
            {
/*
*  If columns of sub( B ) are not contiguous, then allocate the buffer and
*  pack the ABrocs columns of sub( B ).
*/
               Bbufld = MAX( 1, BnpD );
               if( BisR || ( BmyprocR == BcurrocR ) )
               {
                  Bbuf = hipPB_Cmalloc( BnpD * ABrocs * size );
                  PB_CVMpack( TYPE, &VM, ROW,    &Broc, PACKING, NOTRAN,
                              ABrocs, BnpD, one, Mptr( B, BiiD, BkkR, Bld,
                              size ), Bld, zero, Bbuf, Bbufld );
               }
            }
            else
            {
/*
*  Otherwise, re-use sub( B ) directly.
*/
               Bbufld = Bld;
               if( BisR || ( BmyprocR == BcurrocR ) )
                  Bbuf = Mptr( B, BiiD, BkkR + Boff, Bld, size );
            }
            PB_Cdescset( DBUFB, N, ABrocs, Binb1D, ABrocs, BnbD, ABrocs,
                         BrocD, BcurrocR, ctxt, Bbufld );
         }
         
/*
*  Update the local indexes of sub( A ) and sub( B )
*/
         PB_CVMupdate( &VM, ABrocs, &BkkR, &AkkR );
/*
*  Replicate panels of rows or columns of sub( A ) and sub( B ) over sub( C )
*  -> WA, WB
*/
         // magma_queue_sync(myqueue);
         if(buf_tag){
            if(W_ir == 1){
               //if(!checkhipblasErrors(hipblasSetMatrixAsync(Cmp, kbb_old, sizeof(double), WA_BUF, WAd0[LLD_], d_A, Cmp, mystream))) return;
               magma_setmatrix_async(Cmp, kbb_old, sizeof(double), WA_BUF, WAd0[LLD_], d_A, Cmp, myqueue);
               W_ir++;
            }
            hipPB_CInV2( TYPE, &TrA, COLUMN, M, N, Cd0, ABrocs, Abuf, 0, 0,
                         DBUFA, &Aroc, WA, Wkbb, WAd0 );
            if(W_ir == 2){  
               //if(!checkhipblasErrors(hipblasSetMatrixAsync(kbb_old, Cnq, sizeof(double), WB_BUF, WBd0[LLD_], d_B, kb, mystream))) return;
               magma_setmatrix_async(kbb_old, Cnq, sizeof(double), WB_BUF, WBd0[LLD_], d_B, kbb_old, myqueue);
               W_ir++;
            }
            hipPB_CInV2( TYPE, &TrB, ROW,    M, N, Cd0, ABrocs, Bbuf, 0, 0,
                         DBUFB, &Broc, WB, Wkbb, WBd0 );
         }else{
            if(W_ir == 1){
               //if(!checkhipblasErrors(hipblasSetMatrixAsync(Cmp, kbb_old, sizeof(double), WA, WAd0[LLD_], d_A, Cmp, mystream))) return;
               magma_setmatrix_async(Cmp, kbb_old, sizeof(double), WA, WAd0[LLD_], d_A, Cmp, myqueue);
               W_ir++;
            }
            hipPB_CInV2( TYPE, &TrA, COLUMN, M, N, Cd0, ABrocs, Abuf, 0, 0,
                         DBUFA, &Aroc, WA_BUF, Wkbb, WAd0 );
            if(W_ir == 2){
               //if(!checkhipblasErrors(hipblasSetMatrixAsync(kbb_old, Cnq, sizeof(double), WB, WBd0[LLD_], d_B, kb, mystream))) return;
               magma_setmatrix_async(kbb_old, Cnq, sizeof(double), WB, WBd0[LLD_], d_B, kbb_old, myqueue);
               W_ir++;
            }
            hipPB_CInV2( TYPE, &TrB, ROW,    M, N, Cd0, ABrocs, Bbuf, 0, 0,
                         DBUFB, &Broc, WB_BUF, Wkbb, WBd0 );
         }
         if(W_ir == 3){ 
            // gettimeofday(&start, NULL);
            // if(!checkhipblasErrors(hipblasDgemm(handle, hipblas_trans_const(NOTRAN), hipblas_trans_const(NOTRAN),
            //                                     Cmp, Cnq, kbb_old,
            //                                     (double*)ALPHA, d_A, Cmp,
            //                                     d_B, kb,
            //                                     (double*)tbeta, d_C, Cmp)))    return;
            magma_dgemm(magma_trans_const(NOTRAN), magma_trans_const(NOTRAN),
                                             Cmp, Cnq, kbb_old,
                                             *(double*)ALPHA, d_A, Cmp,
                                             d_B, kbb_old,
                                             *(double*)tbeta, d_C, Cmp, myqueue); 
               // gettimeofday(&end, NULL);
               // timeuse = (end.tv_sec - start.tv_sec) * 1000.0 + (end.tv_usec - start.tv_usec) / 1000.0;
               // printf("calculate_time = %lf\n", timeuse);
            tbeta = one;
            W_ir++;
         }
         if( Afr & ( AisR || ( AmyprocR == AcurrocR ) ) )
            if( Abuf ) free( Abuf );
         if( Bfr & ( BisR || ( BmyprocR == BcurrocR ) ) )
            if( Bbuf ) free( Bbuf );
/*
*  ABrocs rows or columns of sub( A ) and sub( B ) have been replicated,
*  update the number of diagonals in this virtual process as well as the
*  number of rows or columns of sub( A ) and sub( B ) that are in WA, WB.
*/
         npq  -= ABrocs;
         Wkbb += ABrocs;
      }
      buf_tag = (buf_tag + 1) % 2;
      W_ir = 1;
      kbb_old = kbb;
/*
*  Perform local update
*/
         // if( Cmp > 0 && Cnq > 0 && kbb > 0)
         // {
            // for hipBLAS
              // HtoD
            // magma_setmatrix(Cmp, kbb, sizeof(double), WA, WAd0[LLD_], (double*)d_A, Cmp, myqueue);
            // magma_setmatrix(kbb, Cnq, sizeof(double), WB, WBd0[LLD_], (double*)d_B, kbb, myqueue);
            //   if(!checkhipblasErrors(hipblasSetMatrix(Cmp, kbb, sizeof(double), WA, WAd0[LLD_], (double*)d_A, Cmp))) return;
            //   if(!checkhipblasErrors(hipblasSetMatrix(kbb, Cnq, sizeof(double), WB, WBd0[LLD_], (double*)d_B, kbb))) return;
              // Run
            // magma_queue_sync(myqueue);
            // hipDeviceSynchronize();
            // magma_dgemm(magma_trans_const(NOTRAN), magma_trans_const(NOTRAN),
            //                                      Cmp, Cnq, kbb,
            //                                      *(double*)ALPHA, (double*)d_A, Cmp,
            //                                      (double*)d_B, kbb,
            //                                      *(double*)tbeta, (double*)d_C, Cmp, myqueue);
            // tbeta = one;
            // magma_queue_sync(myqueue);
            // hipDeviceSynchronize();
            //   if(!checkhipblasErrors(hipblasDgemm(handle, hipblas_trans_const(NOTRAN), hipblas_trans_const(NOTRAN),
            //                                     Cmp, Cnq, kbb,
            //                                     (double*)ALPHA, d_A, Cmp,
            //                                     d_B, kbb,
            //                                     (double*)tbeta, d_C, Cmp)))    return;
            //    tbeta = one;
                  //  gemm( C2F_CHAR( NOTRAN ), C2F_CHAR( NOTRAN ), &Cmp, &Cnq, &kbb,
                  //        ALPHA, WA, &WAd0[LLD_], WB, &WBd0[LLD_], tbeta, Cptr, &Cld );
                  //         tbeta = one;
         // }
        Wkbb = 0;
   }
   if(buf_tag){
      //if(!checkhipblasErrors(hipblasSetMatrixAsync(Cmp, kbb_old, sizeof(double), WA_BUF, WAd0[LLD_], d_A, Cmp, mystream))) return;
      //if(!checkhipblasErrors(hipblasSetMatrixAsync(kbb_old, Cnq, sizeof(double), WB_BUF, WBd0[LLD_], d_B, kb, mystream))) return;
      magma_setmatrix_async(Cmp, kbb_old, sizeof(double), WA_BUF, WAd0[LLD_], d_A, Cmp, myqueue);
      magma_setmatrix_async(kbb_old, Cnq, sizeof(double), WB_BUF, WBd0[LLD_], d_B, kbb_old, myqueue);
      //magmablas_dgeadd(Cmp, Cnq, *(double*)ALPHA, d_C + d_C_size, Cmp, d_C, Cmp, myqueue_buf);
   }else{
      //if(!checkhipblasErrors(hipblasSetMatrixAsync(Cmp, kbb_old, sizeof(double), WA, WAd0[LLD_], d_A, Cmp, mystream))) return;
      //if(!checkhipblasErrors(hipblasSetMatrixAsync(kbb_old, Cnq, sizeof(double), WB, WBd0[LLD_], d_B, kb, mystream))) return;
      magma_setmatrix_async(Cmp, kbb_old, sizeof(double), WA, WAd0[LLD_], d_A, Cmp, myqueue);
      magma_setmatrix_async(kbb_old, Cnq, sizeof(double), WB, WBd0[LLD_], d_B, kbb_old, myqueue);
      //magmablas_dgeadd(Cmp, Cnq, *(double*)ALPHA, d_C + d_C_size, Cmp, d_C, Cmp, myqueue);
   }
   magma_dgemm(magma_trans_const(NOTRAN), magma_trans_const(NOTRAN),
                                                   Cmp, Cnq, kbb_old,
                                                   *(double*)ALPHA, d_A, Cmp,
                                                   d_B, kbb_old,
                                                   *(double*)tbeta, d_C, Cmp, myqueue);
   // if(!checkhipblasErrors(hipblasDgemm(handle, hipblas_trans_const(NOTRAN), hipblas_trans_const(NOTRAN),
   //                                              Cmp, Cnq, kbb_old,
   //                                              (double*)ALPHA, d_A, Cmp,
   //                                              d_B, kb,
   //                                              (double*)tbeta, d_C, Cmp)))    return;
   // magma_dgemm(magma_trans_const(NOTRAN), magma_trans_const(NOTRAN),
   //                                                 Cmp, Cnq, kbb_old,
   //                                                 *(double*)ALPHA, d_A, Cmp,
   //                                                 d_B, kbb_old,
   //                                                 *(double*)tbeta, d_C, Cmp, myqueue);
   // // DtoH
   //if(!checkhipblasErrors(hipblasGetMatrix(Cmp, Cnq, sizeof(double), d_C, Cmp, Cptr, Cld))) return;
   magma_getmatrix(Cmp, Cnq, sizeof(double), d_C, Cmp, Cptr, Cld, myqueue);
   // Deallocate
   // hipFree(d_A);
   // hipFree(d_B);
   // hipFree(d_C);
   magma_free(d_A);
   magma_free(d_B);
   magma_free(d_C);

   hipHostUnregister(WA);
   hipHostUnregister(WB);
   hipHostUnregister(WA_BUF);
   hipHostUnregister(WB_BUF);
   free( WA_BUF );
   free( WB_BUF );
   // hipblasDestroy(handle);
   // hipStreamDestroy(mystream);
   magma_queue_destroy(myqueue);
   //magma_queue_destroy(myqueue_buf);

   if(!checkMagmaErrors(magma_finalize()))  return;
   if( WAfr ) free( WA );
   if( WBfr ) free( WB );
/*
*  End of hipPB_CpgemmAB
*/
}