

#include "mex.h"
#include <assert.h>
#include <math.h>




//=============================================================================//
//=== Auxilliary
//=============================================================================//

//-----------------------------------------------------------------------------//
//--- memory functions
//-----------------------------------------------------------------------------//

void * operator new( size_t s ) 
{
  return mxMalloc( s );
}


void * operator new[]( size_t s ) 
{
  return mxMalloc( s );
}


void operator delete( void *p ) 
{
  mxFree( p );
}


void operator delete[]( void *p )
{
  mxFree( p );
}

/*

void* operator new( size_t s, const nothrow_t& n ) throw()
{
  return mxMalloc( s );
}


void* operator new[]( size_t s, const nothrow_t& n ) throw()
{
  return mr_malloc(s);
}


void operator delete( void *p, const nothrow_t& n ) throw()
{
  mxFree( p );
}


void operator delete[]( void *p, const nothrow_t& n ) throw()
{
  mxFree( p );
}

*/


//-----------------------------------------------------------------------------//
//--- list functions
//-----------------------------------------------------------------------------//

/* find same or just higher item in list
 * !!! to be replaced with binary search !!!
 */
inline int find( const double item, const int len, const double* const list )
{
  register int j;
  j = 0;
  while( j < len && list[j] < item ) {
    assert( j == 0 || list[j-1] < list[j] );
    ++j;
  }
  return j;
}


int insertUniques( int nofUniques, double* const list, int* const nums, const int n, const double* const items )
{
  register int i;
  for( i = 0; i < n; ++i ) {
    register const double item = items[i];
    register const int j = find( item, nofUniques, list );
    assert( 0 <= j && j <= nofUniques );
    // insert item, if not yet in list
    if( j == nofUniques ) {
      assert( j == 0 || list[j-1] < item );
      list[nofUniques] = item;
      nums[nofUniques] = 1;
      ++nofUniques;
    } else if( list[j] == item ) {
      ++nums[j];
    } else {
      assert( list[j] > item );
      assert( j == 0 || list[j-1] < item );
      register int l;
      for( l = nofUniques; l > j; --l ) {
	list[l] = list[l-1];
	nums[l] = nums[l-1];
      }
      list[j] = item;
      nums[j] = 1;
      ++nofUniques;
    }
  }
  return nofUniques;
}


void getNewSortedUniques( int& nofUniques, double*& sortedUniques, int*& uniqueCounts,
			  const int N, const int* ns, const double* const* const items )
{
  assert( N > 0 );
  int sumOfNs;
  register int i;
  sumOfNs = 0;
  for( i = 0; i < N; ++i ) {
    assert( ns[i] > 0 );
    sumOfNs += ns[ i ];
  }
  double* list = new double[ sumOfNs ];
  int* nums = new int[ sumOfNs ];
  nofUniques = 0;
  for( i = 0; i < N; ++i ) {
    nofUniques = insertUniques( nofUniques, list, nums, ns[i], items[i] );
  }
  // compact memory
  sortedUniques = new double[ nofUniques ];
  uniqueCounts = new int[ nofUniques ];
  for( i = 0; i < nofUniques; ++i ) {
    sortedUniques[i] = list[i];
    uniqueCounts[i] = nums[i];
  }
  delete[] nums;
  delete[] list;
}


/*
void getNewSortedUniques1( int& nofUniques, double*& sortedUniques, int*& uniqueCounts,
			   const int n, const double* const items )
{
  getNewSortedUniques( nofUniques, sortedUniques, uniqueCounts, 1, &n, &items );
}
*/


int* newPermutation( int const nofUniques, const int n, const int* const items )
{
  int* uniqueCounts = new int[ nofUniques ];
  int* positions = new int[ nofUniques ];
  int pos;
  int i;
  // count items of each type
  for( i = 0; i < nofUniques; ++i ) {
    uniqueCounts[ i ] = 0;
  }
  for( i = 0; i < n; ++i ) {
    ++uniqueCounts[ items[i] ];
  }
  // store positions of item types
  pos = 0;
  for( i = 0; i < nofUniques; ++i ) {
    positions[ i ] = pos;
    pos += uniqueCounts[ i ];
  }
  assert( pos == n );
  // build perutation
  int* perm = new int[ n ];
  for( i = 0; i < n; ++i ) {
    perm[i] = -1;
  }
  for( i = 0; i < n; ++i ) {
    const int item = items[ i ];
    pos = positions[ item ];
    assert( perm[pos] == -1 );
    perm[ pos ] = i;
    ++positions[ item ];
  }
  #ifdef DEBUG
  pos = 0;
  for( i = 0; i < nofUniques; ++i ) {
    pos += uniqueCounts[i];
    assert( positions[i] == pos );
  }
  #endif
  delete[] positions;
  delete[] uniqueCounts;
  return perm;
}


//-----------------------------------------------------------------------------//
//--- vector and matrix functions
//-----------------------------------------------------------------------------//

int* newIntVector( const int n, const double* const vector,
		   int const nofUniques, const double* const sortedUniques )
{
  int* result = new int[ n ];
  register int l;
  for( l = 0; l < n; ++l ) {
    const double item = vector[ l ];
    const int type = find( item, nofUniques, sortedUniques );
    assert( sortedUniques[type] == item );
    result[ l ] = type;
  }
  return result;
}


int* newPermutedVector( const int n, const int* const vector, const int* const perm )
{
  int* result = new int[ n ];
  register int i;
  for( i = 0; i < n; ++i ) {
    result[i] = vector[ perm[i] ];
  }
  return result;
}


int* newPermutedMatrix( const int n, const int* const matrix, const int* const perm )
{
  int* result = new int[ n * n ];
  register int i;
  register int j;
  int l;
  l = 0;
  for( j = 0; j < n; ++j ) {
    for( i = 0; i < n; ++i ) {
      const int k = perm[i] + n*perm[j];
      result[l] = matrix[k];
      ++l;
    }
  }
  return result;
}


//=============================================================================//
//=== Graph Kernel
//=============================================================================//

//-----------------------------------------------------------------------------//
//--- core kernel function
//-----------------------------------------------------------------------------//

double graphKernel( const int n1, const int n2,
		    const double* edges1, const double* edges2,
		    const double* nodes1, const double* nodes2,
		    const double* probStart1s, const double* probStart2s,
		    const double* probQuit1s,  const double* probQuit2s,
		    const double* probTrans1s, const double* probTrans2s,
		    const double epsilon )
{
  // obtain and sort node types
  const int nodeNums[2] = { n1, n2 };
  const double* nodeData[2] = { nodes1, nodes2 };
  double* sortedNodeTypes;
  int* nodeTypeCounts;
  int nofNodeTypes;
  getNewSortedUniques( nofNodeTypes, sortedNodeTypes, nodeTypeCounts, 2, nodeNums, nodeData );
  const int* intNodes1 = newIntVector( n1, nodes1, nofNodeTypes, sortedNodeTypes );
  const int* intNodes2 = newIntVector( n2, nodes2, nofNodeTypes, sortedNodeTypes );

  // obtain and sort edge types
  const int edgeNums[2] = { n1*n1, n2*n2 };
  const double* edgeData[2] = { edges1, edges2 };
  double* sortedEdgeTypes;
  int* edgeTypeCounts;
  int nofEdgeTypes;
  getNewSortedUniques( nofEdgeTypes, sortedEdgeTypes, edgeTypeCounts, 2, edgeNums, edgeData );
  const int* intEdges1 = newIntVector( n1*n1, edges1, nofEdgeTypes, sortedEdgeTypes );
  const int* intEdges2 = newIntVector( n2*n2, edges2, nofEdgeTypes, sortedEdgeTypes );
  delete[] edgeTypeCounts;  edgeTypeCounts = NULL;
  delete[] sortedEdgeTypes; sortedEdgeTypes = NULL;

  // permute matrices according to node types
  const int* perm1 = newPermutation( nofNodeTypes, n1, intNodes1 );
  const int* perm2 = newPermutation( nofNodeTypes, n2, intNodes2 );
  delete[] nodeTypeCounts;  nodeTypeCounts = NULL; 
  delete[] sortedNodeTypes; sortedNodeTypes = NULL;
  const int* permIntNodes1 = newPermutedVector( n1, intNodes1, perm1 );
  const int* permIntNodes2 = newPermutedVector( n2, intNodes2, perm2 );
  delete[] (void*)intNodes1; intNodes1 = NULL;
  delete[] (void*)intNodes2; intNodes2 = NULL;
  const int* permIntEdges1 = newPermutedMatrix( n1, intEdges1, perm1 );
  const int* permIntEdges2 = newPermutedMatrix( n2, intEdges2, perm2 );
  delete[] (void*)intEdges1; intEdges1 = NULL;
  delete[] (void*)intEdges2; intEdges2 = NULL;

  // setup variables
  const int n12 = n1 * n2;
  register int i1;
  register int i2;
  register int l;
  int k;

  // initialize s and q
  double* s = new double[ n12 ];
  double* q = new double[ n12 ];
  l = 0;
  for( i2 = 0; i2 < n2; ++i2 ) {
    const int i2p = perm2[ i2 ];
    for( i1 = 0; i1 < n1; ++i1 ) {
      const int i1p = perm1[ i1 ];
      s[l] = ( permIntNodes1[i1] == permIntNodes2[i2] ) ? probStart1s[i1p]*probStart2s[i2p] : 0;
      q[l] = probQuit1s[i1p]*probQuit2s[i2p];
      ++l;
    }
  }

  // determine blocks of t
  double** tBlocks = new double*[ nofNodeTypes ];
  int* tBlockSizes1 = new int[ nofNodeTypes ];
  int* tBlockSizes2 = new int[ nofNodeTypes ];
  for( k = 0; k < nofNodeTypes; ++k ) {
    tBlockSizes1[k] = 0;
    tBlockSizes2[k] = 0;
  }
  for( i1 = 0; i1 < n1; ++i1 ) {
    assert( 0 <= permIntNodes1[ i1 ] && permIntNodes1[ i1 ] < nofNodeTypes );
    ++tBlockSizes1[ permIntNodes1[ i1 ] ];
  }
  for( i2 = 0; i2 < n2; ++i2 ) {
    assert( 0 <= permIntNodes2[ i2 ] && permIntNodes2[ i2 ] < nofNodeTypes );
    ++tBlockSizes2[ permIntNodes2[ i2 ] ];
  }
  for( k = 0; k < nofNodeTypes; ++k ) {
    //tBlocks[k] = new double[ n1*tBlockSizes1[k] * n2*tBlockSizes2[k] ];
    //tBlocks[k] = (double*) mxCalloc( n12 * tBlockSizes1[k] * tBlockSizes2[k], sizeof(double) );
    if( sizeof(double) * n12 * tBlockSizes1[k] * tBlockSizes2[k] < 256 * 1024 * 1024 ) {
      //printf( "allocating %d bytes.\n", sizeof(double) * n12 * tBlockSizes1[k] * tBlockSizes2[k] );
      tBlocks[k] = new double[ n12 * tBlockSizes1[k] * tBlockSizes2[k] ];
    } else {
      //printf( "NOT allocating %d bytes.\n", sizeof(double) * n12 * tBlockSizes1[k] * tBlockSizes2[k] );
      tBlocks[k] = NULL;
    }
  }
  //delete[] permIntNodes1; permIntNodes1 = NULL;
  //delete[] permIntNodes2; permIntNodes2 = NULL;

  // initialize t
  int pos1;
  int pos2;
  int j1;
  int j2;
  pos1 = 0;
  pos2 = 0;
  for( k = 0; k < nofNodeTypes; ++k ) {
    const int l1 = tBlockSizes1[ k ];
    const int l2 = tBlockSizes2[ k ];
    double* const tBlock = tBlocks[ k ];
    if( tBlock == NULL ) {
      pos1 += l1;
      pos2 += l2;
      continue;
    }
    if( l1 > 0 && l2 > 0 ) {
      for( j1 = 0; j1 < l1; ++j1 ) {
	const int j1p = perm1[ pos1+j1 ];
	assert( permIntNodes1[ pos1+j1 ] == k );
	for( j2 = 0; j2 < l2; ++j2 ) {
	  const int j2p = perm2[ pos2+j2 ];
	  assert( permIntNodes2[ pos2+j2 ] == k );
	  const int blockIndexBase = n12 * ( j1 + l1*j2 );
	  for( i1 = 0; i1 < n1; ++i1 ) {
	    const int edgeType1 = permIntEdges1[ i1 + n1*(pos1+j1) ];
	    if( edgeType1 == 0 ) {
	      // no continue, since matrix entry must be set to 0
	      //continue;
	    }
	    const int i1p = perm1[ i1 ];
	    for( i2 = 0; i2 < n2; ++i2 ) {
	      const int edgeType2 = permIntEdges2[ i2 + n2*(pos2+j2) ];
	      const int i2p = perm2[ i2 ];
	      const int ii = i1 + n1*i2;
	      assert( ii + blockIndexBase < n12 * l1 * l2 );
	      tBlock[ ii + blockIndexBase ] =
		( edgeType1 != 0 && edgeType1 == edgeType2 )
		? probTrans1s[ i1p+n1*j1p ] * probTrans2s[ i2p+n2*j2p ]
		: 0;
	    }
	  }
	}
      }
    }
    pos1 += l1;
    pos2 += l2;
  }
  assert( pos1 == n1 );
  assert( pos2 == n2 );
  //delete[] perm1; perm1 = NULL; // !!!
  //delete[] perm2; perm2 = NULL; // !!!
  delete[] (void*)permIntNodes1; permIntNodes1 = NULL; // !!!
  delete[] (void*)permIntNodes2; permIntNodes2 = NULL; // !!!

  // do iterations until convergence
  const double epsSquare = epsilon * epsilon;
  register double* R_inf  = new double[ n12 ];
  register double* R_last = new double[ n12 ];
  double relDeltaSquare = epsSquare + 1;
  int nofIters;
  for( l = 0; l < n12; ++l ) {
    R_inf[l] = q[l];
  }
  nofIters = 0;
  while( relDeltaSquare > epsSquare ) {
    // - swap R_inf and R_last
    double* const temp = R_last;
    R_last = R_inf;
    R_inf = temp;
    // - compute: R_inf = q
    for( l = 0; l < n12; ++l ) {
      R_inf[l] = q[l];
    }
    // - compute: R_inf += t' * R_last;
    pos1 = 0;
    pos2 = 0;
    for( k = 0; k < nofNodeTypes; ++k ) {
      const int l1 = tBlockSizes1[ k ];
      const int l2 = tBlockSizes2[ k ];
      const double* const tBlock = tBlocks[ k ];
      if( l1 > 0 && l2 > 0 ) {
        for( j1 = 0; j1 < l1; ++j1 ) {
          for( j2 = 0; j2 < l2; ++j2 ) {
	    const double R_last_jj = R_last[ (pos1+j1) + n1*(pos2+j2) ];
            const int blockIndexBase = n12 * ( j1 + l1*j2 );
            const int njpos1 = n1 * ( pos1 + j1 );
            const int njpos2 = n2 * ( pos2 + j2 );
            for( i1 = 0; i1 < n1; ++i1 ) {
              const int edgeType1 = permIntEdges1[ i1 + njpos1 ];
              if( edgeType1 == 0 ) {
                // continue, since all corresponding matrix entries must be 0
                continue;
              }
              for( i2 = 0; i2 < n2; ++i2 ) {
                const int edgeType2 = permIntEdges2[ i2 + njpos2 ];
		if( edgeType2 != edgeType1 ) {
		  // continue, since all corresponding matrix entries must be 0
		  continue;
		}
		const int ii = i1 + n1*i2;
		assert( ii + blockIndexBase < n12 * l1 * l2 );
		if( tBlock == NULL ) {
		  const int j1p = perm1[ pos1+j1 ];
		  const int j2p = perm2[ pos2+j2 ];
		  const int i1p = perm1[ i1 ];
		  const int i2p = perm2[ i2 ];
		  R_inf[ ii ] += R_last_jj * probTrans1s[ i1p+n1*j1p ] * probTrans2s[ i2p+n2*j2p ];
		  continue;
		}
		R_inf[ ii ] += R_last_jj * tBlock[ ii + blockIndexBase ];
              }
	    }
          }
        }
      }
      pos1 += l1;
      pos2 += l2;
    }
    assert( pos1 == n1 );
    assert( pos2 == n2 );
    // - compute: relDeltaSquare = || R_inf - R_last ||^2 / || R_inf ||^2
    double deltaSquare;
    double normSquare;
    deltaSquare = 0;
    normSquare = 0;
    for( l = 0; l < n12; ++l ) {
      const double value = R_inf[l];
      const double diff = value - R_last[l];
      normSquare += value * value;
      deltaSquare += diff * diff;
    }
    relDeltaSquare = deltaSquare > 0 ? deltaSquare / normSquare : 0;
    //printf( "Delta: rel diff = %f = %f / %f\n", sqrt( relDeltaSquare ), sqrt( deltaSquare ), sqrt( normSquare ) );
    ++nofIters;
  }
  //printf( "Delta: number of iterations = %d\n", nofIters );
  delete[] (void*)perm1; perm1 = NULL; // !!!
  delete[] (void*)perm2; perm2 = NULL; // !!!

  // compute: result = s'*Rinf;
  double result;
  result = 0;
  for( l = 0; l < n12; ++l ) {
    result += s[l] * R_inf[l];
  }
  delete[] R_inf;
  delete[] R_last;
  delete[] tBlockSizes1;
  delete[] tBlockSizes2;
  for( k = 0; k < nofNodeTypes; ++k ) {
    if( tBlocks[k] != NULL ) {
      delete[] tBlocks[k];
      //mxFree( tBlocks[k] );
    }
  }
  delete[] tBlocks;
  delete[] q;
  delete[] s;
  delete[] (void*)permIntEdges1;
  delete[] (void*)permIntEdges2;
  return result;
}


//-----------------------------------------------------------------------------//
//--- matlab interface
//-----------------------------------------------------------------------------//

static const char* const copyright = "written by Alexander Zien, MPI-Kyb, 2003";
static const char* const syntax = "USAGE: res = graphKernelDelta( edges1, edges2, nodes1, nodes2, lambda, epsilon )";

/*
 * function [ res ] = graphKernelDelta( edges1, edges2, nodes1, nodes2, lambda, epsilon );
 */
void mexFunction( int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[] )
{
  // check number of output params
  if( nlhs == 0 ) {
    printf( "%s\n", copyright );
    mexErrMsgTxt( syntax );
    return;
  }
  if( !( nlhs == 1 ) ) {
    int i;
    for( i = 1; i < nlhs; ++i ) {
      plhs[i] = mxCreateDoubleMatrix( 1, 1, mxREAL );
      double* const matPtr = mxGetPr( plhs[i] );
      *matPtr = mxGetNaN();
    }
  }

  // prepare for unsuccessful return
  plhs[0] = mxCreateDoubleMatrix( 1, 1, mxREAL );
  double* const resPtr = mxGetPr( plhs[0] );
  *resPtr = mxGetNaN();

  // check number of input params
  if( !( nrhs == 6 ) ) {
    printf( "found %d parameters\n", nrhs );
    mexErrMsgTxt( syntax );
    return;
  }

  // setup constants
  #define edges1_ ( prhs[ 0 ] )
  #define edges2_ ( prhs[ 1 ] )
  #define nodes1_ ( prhs[ 2 ] )
  #define nodes2_ ( prhs[ 3 ] )
  #define lambda_ ( prhs[ 4 ] )
  #define epsilon_ ( prhs[ 5 ] )
  const double* const edges1 = (double *) mxGetPr( edges1_ );
  const double* const edges2 = (double *) mxGetPr( edges2_ );
  const double* const nodes1 = (double *) mxGetPr( nodes1_ );
  const double* const nodes2 = (double *) mxGetPr( nodes2_ );
  const double lambda = *( (double *) mxGetPr( lambda_ ) );
  const double epsilon = *( (double *) mxGetPr( epsilon_ ) );
  const int n1 = (int)( mxGetM( nodes1_ ) );
  const int n2 = (int)( mxGetM( nodes2_ ) );

  // check input
  if( !( mxGetN( nodes1_ ) == 1 && mxGetN( nodes2_ ) == 1 ) ) {
    printf( "second dimension of nodes must be 1\n" );
    mexErrMsgTxt( syntax );
    return;
  }
  if( !( mxGetM( edges1_ ) == n1 && mxGetN( edges1_ ) == n1 ) ) {
    printf( "edges1 must be square\n" );
    mexErrMsgTxt( syntax );
    return;
  }
  if( !( mxGetM( edges2_ ) == n2 && mxGetN( edges2_ ) == n2 ) ) {
    printf( "edges2 must be square\n" );
    mexErrMsgTxt( syntax );
    return;
  }
  if( !( 0 <= lambda && lambda <= 1 ) ) {
    printf( "lambda must be between 0 and 1\n" );
    mexErrMsgTxt( syntax );
    return;
  }

  // setup variables
  double* probStart1s = new double[ n1 ];
  double* probStart2s = new double[ n2 ];
  double* probQuit1s = new double[ n1 ];
  double* probQuit2s = new double[ n2 ];
  double* probTrans1s = new double[ n1*n1 ];
  double* probTrans2s = new double[ n2*n2 ];

  // initialize probability distributions
  const double probStart1 = 1.0 / n1;
  const double probStart2 = 1.0 / n2;
  const double lambdaNeg = 1 - lambda;
  register int i;
  register int j;
  register int l;
  // - start and quit probabilities, graph 1
  for( i = 0; i < n1; ++i ) {
    probStart1s[ i ] = probStart1;
    probQuit1s[ i ] = lambda;
  }
  // - start and quit probabilities, graph 2
  for( i = 0; i < n2; ++i ) {
    probStart2s[ i ] = probStart2;
    probQuit2s[ i ] = lambda;
  }
  // - transition probabilities, graph 1
  l = 0;
  for( j = 0; j < n1; ++j ) {
    register int outDegree = 0;
    for( i = 0; i < n1; ++i ) {
      outDegree += ( edges1[l] != 0 );
      ++l;
    }
    l -= n1;
    register const double probTrans = lambdaNeg / outDegree;
    for( i = 0; i < n1; ++i ) {
      probTrans1s[l] = ( edges1[l] != 0 ) ? probTrans : 0;
      ++l;
    }
  }
  // - transition probabilities, graph 2
  l = 0;
  for( j = 0; j < n2; ++j ) {
    register int outDegree = 0;
    for( i = 0; i < n2; ++i ) {
      outDegree += ( edges2[l] != 0 );
      ++l;
    }
    l -= n2;
    register const double probTrans = lambdaNeg / outDegree;
    for( i = 0; i < n2; ++i ) {
      probTrans2s[l] = ( edges2[l] != 0 ) ? probTrans : 0;
      ++l;
    }
  }

  // compute
  const double result = graphKernel( n1, n2, edges1, edges2, nodes1, nodes2, probStart1s, probStart2s, probQuit1s, probQuit2s, probTrans1s, probTrans2s, epsilon );

  // return the final result
  *resPtr = result;
  delete[] probStart1s;
  delete[] probStart2s;
  delete[] probQuit1s;
  delete[] probQuit2s;
  delete[] probTrans1s;
  delete[] probTrans2s;
}



