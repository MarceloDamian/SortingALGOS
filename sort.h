
// Code from Mark Allen Weiss textbook
// Modified by: Erick Cabrera

#ifndef SORT_H
#define SORT_H

/**
 * Several sorting routines.
 * Arrays are rearranged with smallest item first.
 */

#include <vector>
#include <functional>
using namespace std;

/**
 * Simple insertion sort.
 */
template <typename Comparable>
void insertionSort( vector<Comparable> & a )
{
    for( int p = 1; p < a.size( ); ++p )
    {
        Comparable tmp = std::move( a[ p ] );

        int j;
        for( j = p; j > 0 && tmp < a[ j - 1 ]; --j )
            a[ j ] = std::move( a[ j - 1 ] );
        a[ j ] = std::move( tmp );
    }
}


/**
 * Internal insertion sort routine for subarrays
 * that is used by quicksort.
 * a is an array of Comparable items.
 * left is the left-most index of the subarray.
 * right is the right-most index of the subarray.
 */
template <typename Comparable>
void insertionSort( vector<Comparable> & a, int left, int right )
{
    for( int p = left + 1; p <= right; ++p )
    {
        Comparable tmp = std::move( a[ p ] );
        int j;

        for( j = p; j > left && tmp < a[ j - 1 ]; --j )
            a[ j ] = std::move( a[ j - 1 ] );
        a[ j ] = std::move( tmp );
    }
}


template <typename Comparable, typename Comparator>
void insertionSort(vector<Comparable>& a, int left, int right, Comparator lessthan){
	for (int p = left + 1; p <= right; ++p) // Alternate helper function derived from the last helper function
	{
		Comparable tmp = std::move(a[p]);
		int j;
		for (j = p; j > left && lessthan(tmp, a[j - 1]); --j)  // changed condition to account for lessthan and compare j and j-1    
			a[j] = std::move(a[j - 1]);
		a[j] = std::move(tmp);
	}
}


/**
 * Shellsort, using Shell's (poor) increments.
 */
template <typename Comparable>
void shellsort( vector<Comparable> & a )
{
    for( int gap = a.size( ) / 2; gap > 0; gap /= 2 )
        for( int i = gap; i < a.size( ); ++i )
        {
            Comparable tmp = std::move( a[ i ] );
            
            int j = i;

            for( ; j >= gap && tmp < a[ j - gap ]; j -= gap )
                a[ j ] = std::move( a[ j - gap ] );
            a[ j ] = std::move( tmp );
        }
}


/*  Standard heapsort. */

template <typename Comparable>
void heapsort( vector<Comparable> & a )
{
    for( int i = a.size( ) / 2 - 1; i >= 0; --i )  // Builds heap
    {  
        percDown( a, i, a.size( ) );
    }
    for( int k = a.size( ) - 1; k > 0; --k )
    {
        std::swap( a[ 0 ], a[ k ] );              // Deletes max
        percDown( a, 0, k );
    }
}


/**
 * Internal method for heapsort.
 * i is the index of an item in the heap.
 * Returns the index of the left child.
 */
inline int leftChild( int i )
{
    return 2 * i + 1;
}





/**
 * Internal method for heapsort that is used in
 * deleteMax and buildHeap.
 * i is the position from which to percolate down.
 * n is the logical size of the binary heap.
 */
template <typename Comparable>
void percDown( vector<Comparable> & a, int i, int n )
{
    int child;
    Comparable tmp;

    for( tmp = std::move( a[ i ] ); leftChild( i ) < n; i = child )
    {
        child = leftChild( i );
        if( child != n - 1 && a[ child ] < a[ child + 1 ] )
            ++child;
        if( tmp < a[ child ] )
            a[ i ] = std::move( a[ child ] );
        else
            break;
    }

    a[ i ] = std::move( tmp );
}

template <typename Comparable, typename Comparator>
void percDown(vector<Comparable>& a, int i, int n, Comparator lessthan)    // Added a compartor datatype
{
    int child;  
    Comparable tmp;
    for (tmp = std::move(a[i]); leftChild(i) < n; i = child)
    {
        child = leftChild(i);
        if (child != n - 1 && lessthan(a[child], a[child + 1]))    // Modified to account for lessthan and not account for the < operator
        {
            ++child;
        }
        if (lessthan(tmp, a[child]))   // replaced < operator with compartor
          {
              a[i] = std::move(a[child]);
          }
        else
        {
            break; 
        }
        
    }
    a[i] = std::move(tmp);
}

template <typename Comparable, typename Comparator>
void HeapSort(vector<Comparable> &a, Comparator lessthan) // copied but modified helper function from above (heapsort)
{
    for (int i = a.size() / 2 - 1; i >= 0; --i) 
    {
        percDown(a, i, a.size(), lessthan); // Extra comparator parameter in this function call to call the new percDown
    }
    for (int k = a.size() - 1; k > 0; --k)
    {
        std::swap(a[0], a[k]);       // swaps first positon for kth position       
        percDown(a, 0, k, lessthan); 
    }	
}	

template <typename Comparable>
void mergeSort( vector<Comparable> & a,
                vector<Comparable> & tmpArray, int left, int right )
{
    if( left < right ) // if right item is greater
    {
        int center = ( left + right ) / 2; // average

        mergeSort( a, tmpArray, left, center ); // rearrangment

        mergeSort( a, tmpArray, center + 1, right ); 

        merge( a, tmpArray, left, center + 1, right ); // placed in right order
    }

}

/**
 * Mergesort algorithm    (driver).
 */
template <typename Comparable>
void mergeSort( vector<Comparable> & a )
{
    vector<Comparable> tmpArray( a.size( ) ); // dynamically created vector with given size

    mergeSort( a, tmpArray, 0, a.size( ) - 1 ); // call to actual mergesort

}


/**
 * Internal method that merges two sorted halves of a subarray.
 * a is an array of Comparable items.
 * tmpArray is an array to place the merged result.
 * leftPos is the left-most index of the subarray.
 * rightPos is the index of the start of the second half.
 * rightEnd is the right-most index of the subarray.
 */
template <typename Comparable>
void merge( vector<Comparable> & a, vector<Comparable> & tmpArray, int leftPos, int rightPos, int rightEnd )
{
    int leftEnd = rightPos - 1;
    int tmpPos = leftPos;
    int numElements = rightEnd - leftPos + 1;

    // Main loop
    while( leftPos <= leftEnd && rightPos <= rightEnd )
        if( a[ leftPos ] <= a[ rightPos ] )
            tmpArray[ tmpPos++ ] = std::move( a[ leftPos++ ] );
        else
            tmpArray[ tmpPos++ ] = std::move( a[ rightPos++ ] );

    while( leftPos <= leftEnd )    // Copy rest of first half
        tmpArray[ tmpPos++ ] = std::move( a[ leftPos++ ] );

    while( rightPos <= rightEnd )  // Copy rest of right half
        tmpArray[ tmpPos++ ] = std::move( a[ rightPos++ ] );

    // Copy tmpArray back
    for( int i = 0; i < numElements; ++i, --rightEnd )
        a[ rightEnd ] = std::move( tmpArray[ rightEnd ] );
}


template <typename Comparable, typename Comparator>
void Merge(vector<Comparable>& a, vector<Comparable>& tmpArray, 
            int leftPos, int rightPos, int rightEnd, Comparator lessthan)  // Added comparator
{
    int leftEnd = rightPos - 1;  // logically its one index less
    int tmpPos = leftPos;   // stores left pos
    int numElements = rightEnd - leftPos + 1; // amount of elements
    // Main loop
    while (leftPos <= leftEnd && rightPos <= rightEnd)  //used to iterate though array 
    {
        if (!lessthan(a[rightPos], a[leftPos]))
        {
            tmpArray[tmpPos++] = std::move(a[leftPos++]);
        }
        else
        {    
            tmpArray[tmpPos++] = std::move(a[rightPos++]);
        }
    }
    while (leftPos <= leftEnd)    
        tmpArray[tmpPos++] = std::move(a[leftPos++]);
    while (rightPos <= rightEnd)  
        tmpArray[tmpPos++] = std::move(a[rightPos++]);
    for (int i = 0; i < numElements; ++i, --rightEnd) // iterating backwards
        a[rightEnd] = std::move(tmpArray[rightEnd]); // rearranging elements
}

template <typename Comparable, typename Comparator>
void MergeSort(vector<Comparable>& a, vector<Comparable>& tmpArray, int left, int right, Comparator lessthan)  //added comparator
{
    if (left < right) //used from previous mergesort
    {
        int center = (left + right) / 2;
        MergeSort(a, tmpArray, left, center, lessthan);    //function call now includes comparator and rearranges
        MergeSort(a, tmpArray, center + 1, right,lessthan); //function call now includes comparator and rearranges
        Merge(a, tmpArray, left, center + 1, right, lessthan); // merges final product
    }
}


template <typename Comparable, typename Comparator>
void MergeSort(vector<Comparable> &a, Comparator lessthan) // added comparator
{
    vector<Comparable> tmpArray(a.size());  // created vector with given size

    MergeSort(a, tmpArray, 0, a.size() - 1, lessthan); //added comparator as a parameter
}

/**
 * Return median of left, center, and right.
 * Order these and hide the pivot.
 */
template <typename Comparable>
const Comparable & median3( vector<Comparable> & a, int left, int right )
{
    int center = ( left + right ) / 2;
    
    if( a[ center ] < a[ left ] )
        std::swap( a[ left ], a[ center ] );
    if( a[ right ] < a[ left ] )
        std::swap( a[ left ], a[ right ] );
    if( a[ right ] < a[ center ] )
        std::swap( a[ center ], a[ right ] );

        // Place pivot at position right - 1
    std::swap( a[ center ], a[ right - 1 ] );
    return a[ right - 1 ];
}

/**
 * Internal quicksort method that makes recursive calls.
 * Uses median-of-three partitioning and a cutoff of 10.
 * a is an array of Comparable items.
 * left is the left-most index of the subarray.
 * right is the right-most index of the subarray.
 */
template <typename Comparable>
void quicksort( vector<Comparable> & a, int left, int right )
{
    if( left + 10 <= right )
    {
        const Comparable & pivot = median3( a, left, right );

            // Begin partitioning
        int i = left, j = right - 1;
        for( ; ; )
        {
            while( a[ ++i ] < pivot ) { }
            while( pivot < a[ --j ] ) { }
            if( i < j )
                std::swap( a[ i ], a[ j ] );
            else
                break;
        }

        std::swap( a[ i ], a[ right - 1 ] );  // Restore pivot

        quicksort( a, left, i - 1 );     // Sort small elements
        quicksort( a, i + 1, right );    // Sort large elements
    }
    else  // Do an insertion sort on the subarray
        insertionSort( a, left, right );
}

/**
 * Quicksort algorithm (driver).
 */
template <typename Comparable>
void quicksort( vector<Comparable> & a )
{
    quicksort( a, 0, a.size( ) - 1 );
}

template <typename Comparable, typename Comparator>
const Comparable & median3(vector<Comparable>& a, int left, int right, Comparator lessthan) // aded comparator 
{
    int center = (left + right) / 2;    //average

    if (lessthan(a[center], a[left]))  // comapres center and left and swaps
        std::swap(a[left], a[center]);
    if (lessthan(a[right], a[left]))   // compares right and left and swaps
        std::swap(a[left], a[right]);
    if (lessthan(a[right], a[center])) // compares right and center and swaps
        std::swap(a[center], a[right]);
    // Place pivot at position right - 1
    std::swap(a[center], a[right - 1]);
    return a[right - 1]; // returns left
}

template <typename Comparable, typename Comparator>
void quicksort( vector<Comparable> & a, int left, int right, Comparator lessthan) // added comparator to quicksort
{
    if( left + 10 <= right )    //in order to iterate
    {
        const Comparable & pivot = median3( a, left, right, lessthan); //called median3 with lessthan parameter            // Begin partitioning
        int i = left;
        int j = right - 1;
        for( ; ; )  // a bit confused from this
        {
            while(lessthan(a[++i], pivot)) { } //used comparator 
            while(lessthan(pivot, a[--j])) { } 
            if( i < j ) // condition to hit every element
                std::swap( a[ i ], a[ j ] );    // swaps elements
            else
                break;
        }
        std::swap( a[ i ], a[ right - 1 ] );  // this restores the pivot

        quicksort( a, left, i - 1, lessthan);     
        quicksort( a, i + 1, right, lessthan);   
    }
    else  
        insertionSort( a, left, right, lessthan); // An insertion sort on the subarray.
}


template <typename Comparable, typename Comparator>
void QuickSort(vector<Comparable> &a, Comparator lessthan) // different quicksort with capitalized Q, but also with compartor
{
    quicksort(a, 0, a.size() - 1, lessthan);
}


template <typename Comparable, typename Comparator>
void quicksort2( vector<Comparable> & a, int left, int right, Comparator lessthan) // compartor added to similar quicksort helper function
{
    if( left + 10 <= right )   
    {
        const Comparable & pivot = a[(left+right)/2]; //finds the average and middle element ultimately
        int i = left;
        int j = right - 1;
        for( ; ; )  
        {
            while(lessthan(a[++i], pivot)) { } //used comparator
            while(lessthan(pivot, a[--j])) { }  //decremented j
            if( i < j )
                std::swap( a[ i ], a[ j ] );    //swaps element
            else
                break;
        }
        std::swap( a[ i ], a[(right+left)/2] );  // This swaps the middle element restoring it. This is one of the only things that are different from the previous quicksort.
        quicksort2( a, left, i - 1, lessthan);  
        quicksort2( a, i + 1, right, lessthan);    
    }
    else  // Do an insertion sort on the subarray
        insertionSort( a, left, right, lessthan);
}


template <typename Comparable, typename Comparator>
void QuickSort2(vector<Comparable> &a, Comparator lessthan)
{
    quicksort2(a, 0, a.size() - 1, lessthan); // using compartor this quicksorts the second vector using the pivot.
}

template <typename Comparable, typename Comparator>
void quicksort3( vector<Comparable> & a, int left, int right, Comparator lessthan) // comparator added and different pivot is included
{
    if( left + 10 <= right )   
    {
        const Comparable & pivot = a[left]; // This is one of the only differences from the previous quicksorts it looks at the left element.
        int i = left;
        int j = right - 1;
        for( ; ; ) 
        {
            while(lessthan(a[++i], pivot)) { } // used lessthan comparator
            while(lessthan(pivot, a[--j])) { } 
            if( i < j )
                std::swap( a[ i ], a[ j ] );    // swaps elements
            else
                break;
        }
        
        quicksort3( a, left, i - 1, lessthan);     
        quicksort3( a, i + 1, right, lessthan);   
    }
    else  // Do an insertion sort on the subarray
        insertionSort( a, left, right, lessthan);
}

//second additional quicksort: First. same as the second, all it does is call the helper function
template <typename Comparable, typename Comparator>
void QuickSort3(vector<Comparable> &a, Comparator lessthan)
{
    quicksort3(a, 0, a.size() - 1, lessthan); // uses comparator to get quicksort3.
}


/**
 * Internal selection method that makes recursive calls.
 * Uses median-of-three partitioning and a cutoff of 10.
 * Places the kth smallest item in a[k-1].
 * a is an array of Comparable items.
 * left is the left-most index of the subarray.
 * right is the right-most index of the subarray.
 * k is the desired rank (1 is minimum) in the entire array.
 */
template <typename Comparable>
void quickSelect( vector<Comparable> & a, int left, int right, int k )
{
    if( left + 10 <= right )
    {
        const Comparable & pivot = median3( a, left, right );

            // Begin partitioning
        int i = left, j = right - 1;
        for( ; ; )
        {
            while( a[ ++i ] < pivot ) { }
            while( pivot < a[ --j ] ) { }
            if( i < j )
                std::swap( a[ i ], a[ j ] );
            else
                break;
        }

        std::swap( a[ i ], a[ right - 1 ] );  // Restore pivot

            // Recurse; only this part changes
        if( k <= i )
            quickSelect( a, left, i - 1, k );
        else if( k > i + 1 )
            quickSelect( a, i + 1, right, k );
    }
    else  // Do an insertion sort on the subarray
        insertionSort( a, left, right );
}

/**
 * Quick selection algorithm.
 * Places the kth smallest item in a[k-1].
 * a is an array of Comparable items.
 * k is the desired rank (1 is minimum) in the entire array.
 */
template <typename Comparable>
void quickSelect( vector<Comparable> & a, int k )
{
    quickSelect( a, 0, a.size( ) - 1, k );
}


template <typename Comparable>
void SORT( vector<Comparable> & items )
{
    if( items.size( ) > 1 )
    {
        vector<Comparable> smaller;
        vector<Comparable> same;
        vector<Comparable> larger;
        
        auto chosenItem = items[ items.size( ) / 2 ];
        
        for( auto & i : items )
        {
            if( i < chosenItem )
                smaller.push_back( std::move( i ) );
            else if( chosenItem < i )
                larger.push_back( std::move( i ) );
            else
                same.push_back( std::move( i ) );
        }
        
        SORT( smaller );     // Recursive call!
        SORT( larger );      // Recursive call!
        
        std::move( begin( smaller ), end( smaller ), begin( items ) );
        std::move( begin( same ), end( same ), begin( items ) + smaller.size( ) );
        std::move( begin( larger ), end( larger ), end( items ) - larger.size( ) );

/*
        items.clear( );
        items.insert( end( items ), begin( smaller ), end( smaller ) );
        items.insert( end( items ), begin( same ), end( same ) );
        items.insert( end( items ), begin( larger ), end( larger ) );
*/
    }
}

/*
 * This is the more public version of insertion sort.
 * It requires a pair of iterators and a comparison
 * function object.
 */
template <typename RandomIterator, typename Comparator>
void insertionSort( const RandomIterator & begin, const RandomIterator & end, Comparator lessThan )
{
    if( begin == end )
        return;
        
    RandomIterator j;

    for( RandomIterator p = begin+1; p != end; ++p )
    {
        auto tmp = std::move( *p );
        for( j = p; j != begin && lessThan( tmp, *( j-1 ) ); --j )
            *j = std::move( *(j-1) );
        *j = std::move( tmp );
    }
}

/*
 * The two-parameter version calls the three parameter version, using C++11 decltype
 */
template <typename RandomIterator>
void insertionSort( const RandomIterator & begin,
                    const RandomIterator & end )
{
    insertionSort( begin, end, less<decltype(*begin )>{ } );
}



#endif
