/* ---------------------------------------------------------------------
*
*  -- PBLAS auxiliary routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*  ---------------------------------------------------------------------
*/
/*
*  Include files
*/
#include "../pblas.h"
#include "../PBpblas.h"
#include "../PBtools.h"
#include "../PBblacs.h"
#include "../PBblas.h"
#include "../PB_hipBLAS.h"

#ifdef __STDC__
void hipPB_CpgemmAC( PBTYP_T * TYPE, char * DIRECA, char * DIRECC,
                  char * TRANSA, char * TRANSB, int M, int N, int K,
                  char * ALPHA, char * A, int IA, int JA, int * DESCA,
                  char * B, int IB, int JB, int * DESCB, char * BETA,
                  char * C, int IC, int JC, int * DESCC )
#else
void hipPB_CpgemmAC( TYPE, DIRECA, DIRECC, TRANSA, TRANSB, M, N, K, ALPHA,
                  A, IA, JA, DESCA, B, IB, JB, DESCB, BETA, C, IC, JC,
                  DESCC )
/*
*  .. Scalar Arguments ..
*/
   char           * DIRECA, * DIRECC, * TRANSA, * TRANSB;
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
   char           Aroc, GemmTa, GemmTb, TrA, TrB, * one, * talpha, * tbeta,
                  top, * zero;
   int            Abufld, AcurrocR, Afr, Afwd, AiD, AiR, AiiD, AiiR, AinbD,
                  AinbR, Ainb1D, Ainb1R, AisR, Akk, Ald, AmyprocD, AmyprocR,
                  AnbD, AnbR, AnpD, AnpR, AnprocsD, AnprocsR, Aoff, ArocD,
                  ArocR, AsrcR, Asrc_, Bcol, Bii, Bimb1, Binb1, Bjj, Bld, Bm,
                  Bmb, Bmp, Bn, Bnb, Bnq, Brow, Cbufld, Ccol, Ccurrow, Cfr,
                  Cfwd, Cii, Cimb, Cimb1, Cinb, Cinb1, CisR, Cjj, Ckk, Cld,
                  Cmb, Cmp, Cnb, Cnq, Coff, Crow, Csrc, WAfr, WCfr, WCsum,
                  ctxt, lcmb, m, maxp, maxpm1, maxq, mb, mbb, mycol, myrow,
                  ncpq, nota, notb, npcol, npq=0, nprow, nrpq, p=0, q=0, size,
                  tmp;
   GEMM_T         gemm;
   GSUM2D_T       gsum2d;
/*
*  .. Local Arrays ..
*/
   PB_VM_T        VM;
   int            Bd0[DLEN_], DBUFA[DLEN_], DBUFC[DLEN_], WAd[DLEN_],
                  WCd[DLEN_];
   char           * Abuf = NULL, * Bptr = NULL, * Cbuf = NULL, * WA = NULL,
                  * WC   = NULL;
/* ..
*  .. Executable Statements ..
*
*/
   int WC_size;
   int A_row, B_row, A_column, B_column;
   bool first_ir = true;
   hipblasHandle_t handle;
   char* d_A = NULL;
   char* d_B = NULL;
   char* d_C = NULL;

   hipblasCreate(&handle);

   Cblacs_gridinfo( ( ctxt = DESCC[CTXT_] ), &nprow, &npcol, &myrow, &mycol );

   Afwd = ( Mupcase( DIRECA[0] ) == CFORWARD );
   Cfwd = ( Mupcase( DIRECC[0] ) == CFORWARD );
   nota = ( ( TrA = Mupcase( TRANSA[0] ) ) == CNOTRAN );
   notb = ( ( TrB = Mupcase( TRANSB[0] ) ) == CNOTRAN );

   size = TYPE->size; one  = TYPE->one; zero = TYPE->zero;
   gemm = TYPE->Fgemm; gsum2d = TYPE->Cgsum2d;
   mb   = pilaenv_( &ctxt, C2F_CHAR( &TYPE->type ) );
/*
*  Compute local information for sub( A ), sub( B ) and sub( C )
*/
   if( nota )
   {
      AiD      = JA;           AiR      = IA;
      Asrc_    = RSRC_;        Aroc     = CROW;
      AinbR    = DESCA[IMB_ ]; AinbD    = DESCA[INB_];
      AnbR     = DESCA[MB_  ]; AnbD     = DESCA[NB_ ];
      AsrcR    = DESCA[Asrc_]; Ald      = DESCA[LLD_];
      AmyprocD = mycol;        AnprocsD = npcol;
      AmyprocR = myrow;        AnprocsR = nprow;
      PB_Cinfog2l( IA, JA, DESCA, AnprocsR, AnprocsD, AmyprocR, AmyprocD,
                   &AiiR, &AiiD, &ArocR, &ArocD );
   }
   else
   {
      AiD      = IA;           AiR      = JA;
      Asrc_    = CSRC_;        Aroc     = CCOLUMN;
      AinbD    = DESCA[IMB_ ]; AinbR    = DESCA[INB_];
      AnbD     = DESCA[MB_  ]; AnbR     = DESCA[NB_ ];
      AsrcR    = DESCA[Asrc_]; Ald      = DESCA[LLD_];
      AmyprocD = myrow;        AnprocsD = nprow;
      AmyprocR = mycol;        AnprocsR = npcol;
      PB_Cinfog2l( IA, JA, DESCA, AnprocsD, AnprocsR, AmyprocD, AmyprocR,
                   &AiiD, &AiiR, &ArocD, &ArocR );
   }
   Ainb1D = PB_Cfirstnb( K, AiD, AinbD, AnbD );
   AnpD   = PB_Cnumroc( K, 0, Ainb1D, AnbD, AmyprocD, ArocD, AnprocsD );
   Ainb1R = PB_Cfirstnb( M, AiR, AinbR, AnbR );

   Cimb   = DESCC[IMB_ ]; Cinb = DESCC[INB_];
   Cmb    = DESCC[MB_  ]; Cnb  = DESCC[NB_ ];
   Csrc   = DESCC[RSRC_]; Cld  = DESCC[LLD_];
   PB_Cinfog2l( IC, JC, DESCC, nprow, npcol, myrow, mycol, &Cii, &Cjj,
                &Crow, &Ccol );
   Cimb1 = PB_Cfirstnb( M, IC, Cimb, Cmb );
   Cinb1 = PB_Cfirstnb( N, JC, Cinb, Cnb );
   Cnq   = PB_Cnumroc( N, 0, Cinb1, Cnb, mycol, Ccol, npcol );
/*
*  Retrieve the BLACS combine topology, compute conjugate of alpha for the
*  conjugate transpose case and set the transpose parameters to be passed to
*  the BLAS matrix multiply routine.
*/
   if( notb )
   {
      Bm     = K; Bn     = N;
      top    = *PB_Ctop( &ctxt, COMBINE, COLUMN, TOP_GET );
      talpha = ALPHA; GemmTa = ( nota ? CTRAN : TrA ); GemmTb = CNOTRAN;
   }
   else
   {
      Bm     = N; Bn     = K;
      top    = *PB_Ctop( &ctxt, COMBINE, ROW,    TOP_GET );
      if( TrB == CCOTRAN )
      {
         talpha = hipPB_Cmalloc( size ); PB_Cconjg( TYPE, ALPHA, talpha );
         GemmTb = ( ( TrA == CCOTRAN ) ? CTRAN : CCOTRAN );
      }
      else
      {
         talpha = ALPHA;
         GemmTb = ( ( TrA == CCOTRAN ) ? CCOTRAN : CTRAN );
      }
      GemmTa = CNOTRAN;
   }
/*
*  Compute descriptor Bd0 for sub( B )
*/
   PB_Cdescribe( Bm, Bn, IB, JB, DESCB, nprow, npcol, myrow, mycol, &Bii, &Bjj,
                 &Bld, &Bimb1, &Binb1, &Bmb, &Bnb, &Brow, &Bcol, Bd0 );

   Bmp = PB_Cnumroc( Bm, 0, Bimb1, Bmb, myrow, Brow, nprow );
   Bnq = PB_Cnumroc( Bn, 0, Binb1, Bnb, mycol, Bcol, npcol );
   if( ( Bmp > 0 ) && ( Bnq > 0 ) ) Bptr = Mptr( B, Bii, Bjj, Bld, size );
/*
*  When sub( A ) is not replicated and backward pass on sub( A ), find the
*  virtual process q owning the last row or column of sub( A ).
*/
   if( !( AisR = ( ( AsrcR < 0 ) || ( AnprocsR == 1 ) ) ) && !Afwd )
   {
      tmp = PB_Cindxg2p( M - 1, Ainb1R, AnbR, ArocR, ArocR, AnprocsR );
      q   = MModSub( tmp, ArocR, AnprocsR );
   }
/*
*  When sub( C ) is not replicated and backward pass on sub( C ), find the
*  virtual process p owning the last row or column of sub( C ).
*/
   if( !( CisR = ( ( Crow < 0 ) || ( nprow == 1 ) ) ) && !Cfwd )
   {
      tmp = PB_Cindxg2p( M - 1, Cimb1, Cmb, Crow, Crow, nprow );
      p   = MModSub( tmp, Crow, nprow );
   }
/*
*  Loop over the virtual process grid induced by the rows or columns of
*  sub( A ) and sub( C ).
*/
   lcmb   = PB_Clcm( ( maxp = ( CisR ? 1 : nprow    ) ) * Cmb,
                     ( maxq = ( AisR ? 1 : AnprocsR ) ) * AnbR );
   m      = M;
   maxpm1 = maxp - 1;
   while( m > 0 )
   {
/*
*  Initialize local virtual matrix in process (p,q)
*/
      AcurrocR = ( AisR ? -1 : MModAdd( ArocR, q, AnprocsR ) );
      Akk      = PB_Cg2lrem( AiR, AinbR, AnbR, AcurrocR, AsrcR, AnprocsR );
      AnpR     = PB_Cnumroc( M, 0, Ainb1R, AnbR, AcurrocR, ArocR, AnprocsR );

      Ccurrow  = ( CisR ? -1 : MModAdd( Crow,  p, nprow    ) );
      Ckk      = PB_Cg2lrem( IC, Cimb, Cmb, Ccurrow, Csrc, nprow );
      Cmp      = PB_Cnumroc( M, 0, Cimb1, Cmb, Ccurrow, Crow, nprow );

      PB_CVMinit( &VM, 0, Cmp, AnpR, Cimb1, Ainb1R, Cmb, AnbR, p, q,
                  maxp, maxq, lcmb );
/*
*  Find how many diagonals in this virtual process
*/
      npq = PB_CVMnpq( &VM );

      m  -= npq;
/*
*  Re-adjust the number of rows or columns to be (un)packed, in order to
*  average the message sizes.
*/
      if( npq ) mbb = npq / ( ( npq - 1 ) / mb + 1 );
      //由于mbb的大小的实际输出时都差不多一样大，所以这里只hipMalloc一次。
      while( npq )
      {
         mbb = MIN( mbb, npq );
/*
*  Find out how many rows or columns of sub( A ) and sub( C ) are contiguous
*/
         PB_CVMcontig( &VM, &nrpq, &ncpq, &Coff, &Aoff );

         if( nota )
         {
/*
*  Compute the descriptor DBUFA for the buffer that will contained the packed
*  columns of sub( A ).
*/
            if( ( Afr = ( ncpq < mbb ) ) != 0 )
            {
/*
*  If rows of sub( A ) are not contiguous, then allocate the buffer and
*  pack the mbb rows of sub( A ).
*/
               Abufld = mbb;
               if( AisR || ( AmyprocR == AcurrocR ) )
               {
                  Abuf   = hipPB_Cmalloc( AnpD * mbb * size );
                  PB_CVMpack( TYPE, &VM, COLUMN, &Aroc, PACKING, NOTRAN, mbb,
                              AnpD, one, Mptr( A, Akk, AiiD, Ald, size ), Ald,
                              zero, Abuf, Abufld );
               }
            }
            else
            {
/*
*  Otherwise, re-use sub( B ) directly.
*/
               Abufld = Ald;
               if( AisR || ( AmyprocR == AcurrocR ) )
                  Abuf = Mptr( A, Akk+Aoff, AiiD, Ald, size );
            }
            PB_Cdescset( DBUFA, mbb, K, mbb, Ainb1D, mbb, AnbD, AcurrocR,
                         ArocD, ctxt, Abufld );
         }
         else
         {
/*
*  Compute the descriptor DBUFA for the buffer that will contained the packed
*  columns of sub( A ).
*/
            if( ( Afr = ( ncpq < mbb ) ) != 0 )
            {
/*
*  If columns of sub( A ) are not contiguous, then allocate the buffer and pack
*  the mbb columns of sub( A ).
*/
               Abufld = MAX( 1, AnpD );
               if( AisR || ( AmyprocR == AcurrocR ) )
               {
                  Abuf   = hipPB_Cmalloc( AnpD * mbb * size );
                  PB_CVMpack( TYPE, &VM, COLUMN, &Aroc, PACKING, NOTRAN, mbb,
                              AnpD, one, Mptr( A, AiiD, Akk, Ald, size ), Ald,
                              zero, Abuf, Abufld );
               }
            }
            else
            {
/*
*  Otherwise, re-use sub( A ) directly.
*/
               Abufld = Ald;
               if( AisR || ( AmyprocR == AcurrocR ) )
                  Abuf = Mptr( A, AiiD, Akk+Aoff, Ald, size );
            }
            PB_Cdescset( DBUFA, K, mbb, Ainb1D, mbb, AnbD, mbb, ArocD,
                         AcurrocR, ctxt, Abufld );
         }

         if( notb )
         {
/*
*  Replicate this panel of rows or columns of sub( A ) over sub( B ) -> WA
*/
            hipPB_CInV( TYPE, NOCONJG, COLUMN, Bm, Bn, Bd0, mbb, Abuf, 0, 0,
                     DBUFA, &Aroc, &WA, WAd, &WAfr );
/*
*  Allocate space for temporary results in scope of sub( B ) -> WC
*/
            hipPB_COutV( TYPE, ROW,    INIT, Bm, Bn, Bd0, mbb, &WC, WCd, &WCfr,
                      &WCsum, &WC_size);
/*
*  Local matrix-matrix multiply iff I own some data
*/
            
            if( Bmp > 0 && Bnq > 0 && mbb > 0)
            {
               if(first_ir){
                  if(GemmTb == 'N'){
                     B_row = Bmp;
                     B_column = Bnq;
                  }else{
                     B_row = Bnq;
                     B_column = Bmp;
                  }
                  if(!checkhipErrors(hipMalloc((void**)&d_B, sizeof(double)*Bmp*Bnq)))  return;
                  if(!checkhipblasErrors(hipblasSetMatrix(B_row, B_column, sizeof(double), Bptr, Bld, (double*)d_B, B_row))) return;
                  first_ir = false;
               }
               //HtoD
               if(GemmTa == 'N'){
                  A_row = mbb;
                  A_column = Bmp;
               }else{
                  A_row = Bmp;
                  A_column = mbb;
               }
               if(!checkhipErrors(hipMalloc((void**)&d_A, sizeof(double)*Bmp*mbb)))  return;
               if(!checkhipErrors(hipMalloc((void**)&d_C, sizeof(double)*mbb*Bnq)))  return;
               if(!checkhipblasErrors(hipblasSetMatrix(A_row, A_column, sizeof(double), WA, WAd[LLD_], (double*)d_A, A_row))) return;
               //if(!checkhipblasErrors(hipblasSetMatrix(mbb  , Bnq     , sizeof(double), WC, WCd[LLD_], (double*)d_C, mbb  ))) return;
               // Run
               if(!checkhipblasErrors(hipblasDgemm(handle, hipblas_trans_const(&GemmTa), hipblas_trans_const(&GemmTb),
                                                   mbb, Bnq, Bmp,
                                                   (double*)talpha, (double*)d_A, A_row,
                                                   (double*)d_B, B_row,
                                                   (double*)zero,  (double*)d_C, mbb)))    return;
               // DtoH
               if(!checkhipblasErrors(hipblasGetMatrix(mbb, Bnq, sizeof(double), (double*)d_C, mbb, WC, WCd[LLD_]))) return;
               //gemm( C2F_CHAR( &GemmTa ), C2F_CHAR( &GemmTb ), &mbb, &Bnq, &Bmp,
                     //talpha, WA, &WAd[LLD_], Bptr, &Bld, zero, WC, &WCd[LLD_] );
               hipFree(d_A);
               hipFree(d_C);
            }
            else
            {
               gemm( C2F_CHAR( &GemmTa ), C2F_CHAR( &GemmTb ), &mbb, &Bnq, &Bmp,
                     talpha, WA, &WAd[LLD_], Bptr, &Bld, zero, WC, &WCd[LLD_] );
            }
            if( WAfr ) free( WA );
            if( Afr && ( AisR || ( AmyprocR == AcurrocR ) ) )
               if( Abuf ) free( Abuf );
/*
*  Accumulate the intermediate results in WC
*/
            if( WCsum )
            {
               WCd[RSRC_] = Ccurrow;
               if( Bnq > 0 )
               //blacs中进程间的规约操作
                  gsum2d( ctxt, COLUMN, &top, mbb, Bnq, WC, WCd[LLD_],
                          WCd[RSRC_], mycol );
            }
/*
*  Compute the descriptor DBUFC for the buffer that will contained the packed
*  rows of sub( C ). Allocate it.
*/
            if( ( Cfr = ( nrpq < mbb ) ) != 0 )
            {
/*
*  If rows of sub( C ) are not contiguous, then allocate the buffer
*/
               Cbufld = mbb; tbeta = zero;
               if( CisR || ( myrow == Ccurrow ) )
                  Cbuf = hipPB_Cmalloc( Cnq * mbb * size );
            }
            else
            {
/*
*  Otherwise re-use sub( C )
*/
               Cbufld = Cld; tbeta = BETA;
               if( CisR || ( myrow == Ccurrow ) )
                  Cbuf = Mptr( C, Ckk+Coff, Cjj, Cld, size );
            }
            PB_Cdescset( DBUFC, mbb, N, mbb, Cinb1, mbb, Cnb, Ccurrow, Ccol,
                         ctxt, Cbufld );
/*
*  Cbuf := Cbuf + WC
*/
            PB_Cpaxpby( TYPE, NOCONJG, mbb, N, one, WC, 0, 0, WCd, ROW, tbeta,
                        Cbuf, 0, 0, DBUFC, ROW );
/*
*  Unpack the mbb rows of sub( C ) and release the buffer containing them.
*/
            if( Cfr && ( CisR || ( myrow == Ccurrow ) ) )
            {
               PB_CVMpack( TYPE, &VM, ROW, ROW,    UNPACKING, NOTRAN, mbb, Cnq,
                           BETA, Mptr( C, Ckk, Cjj, Cld, size ), Cld, one, Cbuf,
                           Cbufld );
               if( Cbuf ) free( Cbuf );
            }
            if( WCfr ) free( WC );
         }
         else
         {
/*
*  Replicate this panel of rows or columns of sub( A ) over sub( B ) -> WA
*/
            hipPB_CInV( TYPE, NOCONJG, ROW,    Bm, Bn, Bd0, mbb, Abuf, 0, 0,
                     DBUFA, &Aroc, &WA, WAd, &WAfr );
/*
*  Allocate space for temporary results in scope of sub( A ) -> WC
*/
            hipPB_COutV( TYPE, COLUMN, INIT, Bm, Bn, Bd0, mbb, &WC, WCd, &WCfr,
                      &WCsum, &WC_size);
/*
*  Local matrix-matrix multiply iff I own some data
*/
            if( Bmp > 0 && Bnq > 0 )
               gemm( C2F_CHAR( &GemmTa ), C2F_CHAR( &GemmTb ), &Bmp, &mbb, &Bnq,
                     talpha, Bptr, &Bld, WA, &WAd[LLD_], zero, WC, &WCd[LLD_] );
            if( WAfr ) free( WA );
            if( Afr && ( AisR || ( AmyprocR == AcurrocR ) ) )
               if( Abuf ) free( Abuf );
/*
*  Accumulate the intermediate results in WC
*/
            if( WCsum )
            {
               WCd[CSRC_] = 0;
               if( Bmp > 0 )
                  gsum2d( ctxt, ROW,    &top, Bmp, mbb, WC, WCd[LLD_], myrow,
                          WCd[CSRC_] );
            }
/*
*  Compute the descriptor DBUFC for the buffer that will contained the packed
*  rows of sub( C ). Allocate it.
*/
            if( ( Cfr = ( nrpq < mbb ) ) != 0 )
            {
/*
*  If rows of sub( C ) are not contiguous, then allocate the buffer
*/
               Cbufld = mbb; tbeta = zero;
               if( CisR || ( myrow == Ccurrow ) )
                  Cbuf = hipPB_Cmalloc( Cnq * mbb * size );
            }
            else
            {
/*
*  Otherwise re-use sub( C )
*/
               Cbufld = Cld; tbeta = BETA;
               if( CisR || ( myrow == Ccurrow ) )
                  Cbuf = Mptr( C, Ckk+Coff, Cjj, Cld, size );
            }
            PB_Cdescset( DBUFC, mbb, N, mbb, Cinb1, mbb, Cnb, Ccurrow, Ccol,
                         ctxt, Cbufld );
/*
*  Cbuf := Cbuf + WC'
*/
            PB_Cpaxpby( TYPE, ( TrB == CCOTRAN ? CONJG : NOCONJG ), N, mbb,
                        one, WC, 0, 0, WCd, COLUMN, tbeta, Cbuf, 0, 0, DBUFC,
                        ROW    );
/*
*  Unpack the mbb rows of sub( C ) and release the buffer containing them.
*/
            if( Cfr && ( CisR || ( myrow == Ccurrow ) ) )
            {
               PB_CVMpack( TYPE, &VM, ROW, ROW,    UNPACKING, NOTRAN, mbb, Cnq,
                           BETA, Mptr( C, Ckk, Cjj, Cld, size ), Cld, one, Cbuf,
                           Cbufld );
               if( Cbuf ) free( Cbuf );
            }
            if( WCfr ) free( WC );
         }
/*
*  Update the local indexes of sub( B ) and sub( C )
*/
         PB_CVMupdate( &VM, mbb, &Ckk, &Akk );

         npq -= mbb;
      }
/*
*  Go to next or previous virtual process row or column
*/
      if( ( Cfwd      && ( p == maxpm1 ) ) ||
          ( !( Cfwd ) && ( p == 0      ) ) )
         q = ( Afwd ? MModAdd1( q, maxq ) : MModSub1( q, maxq ) );
      p = ( Cfwd ? MModAdd1( p, maxp ) : MModSub1( p, maxp ) );
   }
   if( Bmp > 0 && Bnq > 0 )
   {
      // Deallocate
      hipFree(d_B);
      
   }
   hipblasDestroy(handle);
   if( TrB == CCOTRAN ) free( talpha );

/*
*  End of PB_CpgemmAC
*/
}
