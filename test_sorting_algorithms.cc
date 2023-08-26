
#include "Sort.h"
#include <chrono>
#include <iostream>
#include <fstream>
#include <functional>
#include <string>
#include <vector>
#include <time.h>  // This is neccsary to mark time stamps


using namespace std;

// Test function that shows how you can time a piece of code.
// Just times a simple loop.

void TestingTiming() 
{
  cout << "Testing Timing" << endl;
  const auto begin = chrono::high_resolution_clock::now();
  // Time this piece of code.
  int sum = 0;
  for (int i = 1; i < 10000; ++i) sum ++;
  // End of piece of code to time.
  const auto end = chrono::high_resolution_clock::now();    
  cout << chrono::duration_cast<chrono::nanoseconds>(end - begin).count() << "ns" << endl;
  cout << chrono::duration_cast<chrono::milliseconds>(end - begin).count() << "ms" << endl;

}

// Generates and returns random vector of size @size_of_vector.
vector<int> GenerateRandomVector(size_t size_of_vector) 
{
  // Add code
    srand(static_cast<unsigned>(time(0)));    //creates a time stamp
    
    vector<int> ranoutput; //creates vector to output. Datatype integers.
    
    for (int i=0; i < size_of_vector; i++) // reads vector
    {          
        ranoutput.push_back(rand() % size_of_vector);  //creates a new random value making sure to not overlap.
    }
    return ranoutput;  //returns vector
}

// Generate and returns sorted vector of size @size_of_vector
// If smaller_to_larger is true, returns vector sorted from small to large
// Otherwise returns vector sorted from large to small
vector<int> GenerateSortedVector(size_t size_of_vector, bool smaller_to_larger) 
{
  // Add code
  vector<int> sortoutput; //creates vector to output

  if(smaller_to_larger)
  {    
    for (int j=0; j < size_of_vector; j++) // reads vector
    {    //goes through the vector  
        sortoutput.push_back(j);   // pushes iteration of for loop in vector. Starting at index 0.
    }
  }
  else
  {    
    for (int i=0; i < size_of_vector; i++)
    {    //goes through the vector 
        sortoutput.push_back(size_of_vector-i-1);    // instead of incrementing. This actually decrements and pushes them into the vector.
    }
  }
  return sortoutput;    //returns vector
}

// Verifies that a vector is sorted given a comparator
template <typename Comparable, typename Comparator>
bool VerifyOrder(const vector<Comparable> &input, Comparator less_than) 
{
  // Add code
      for (int i=1; i<input.size(); i++) // reads vector input.
      { 
        if(less_than(input[i], input[i-1])) // compares vector to element to the left.
        {    
            return false;       // this it is not sorted if it is not less than. 
        }
    }
    return true;
}

// Computes duration given a start time and a stop time in nano seconds
long long ComputeDuration(chrono::high_resolution_clock::time_point start_time, chrono::high_resolution_clock::time_point end_time) 
{
  // Add code

    long long ns = chrono::duration_cast<chrono::nanoseconds>(end_time - start_time).count(); // Required to store it as nanoseconds.
    return ns;    // returns it as a long long so at least 64 bits.
}

void testSortingWrapper(int argc, char **argv) 
{
  const string inputtype = string(argv[1]); // stores type as string
  const int inputsize = stoi(string(argv[2])); // converts string into integer.
  const string comparisontype = string(argv[3]); // needed as the command line requires it to be a string.
  
  if (inputtype != "random" && inputtype != "sorted_small_to_large" && inputtype != "sorted_large_to_small")  // if the command line is not one of these 3.
  {
    cout << "Invalid input type" << endl;
    
    return; // exit (0) also works. Return 0 doesnt since it is a void function.
  }
  if (inputsize <= 0) // there has to be something greater that 0 for it to be a valid size.
  {
    cout << "Invalid size" << endl;
    return;
  }
  if (comparisontype != "less" && comparisontype != "greater")  // The only comparison types are either greater or lesser. If the command line doesnt read either of those strings then exit.
  {
    cout << "Invalid comparison type" << endl;
    return;
  }

  cout << "Running sorting algorithms: " << inputtype << " " << inputsize << " numbers " << comparisontype << endl;
  
  vector<int> inputvector;
  
  
  if (inputtype == "random") // in the case that the string from the command line is random then continue. 
  {
    // Generate random vector
    inputvector = GenerateRandomVector(inputsize);    // call generate random vector.
  } 
  else  // it is sorted
  {  
  
    if(inputtype== "sorted_small_to_large")
    {   
        inputvector = GenerateSortedVector(inputsize, true); // if it is sorted from left to right 
    } 
    else // in the case it isnt
    { 
        inputvector = GenerateSortedVector(inputsize, false); // switch it to false thus suggesting it isnt sorted left to right. 
    }
  }

    using Clock = std::chrono::high_resolution_clock; // creates clock.
    
    if (comparisontype == "less")  // if string from the command line is less.
    {    
      auto begin = Clock::now();  //starts the clock timer.
      HeapSort(inputvector, less<int>{}); // calls heapsort
      auto end = Clock::now();    //stops clock and stores the ns time stopped.
      
      cout << "Heapsort" << endl << endl;     //declares that heapsort ran, the ns it took to run and that it is ordered
      cout << "Runtime: " << ComputeDuration(begin, end) << "ns" << endl; // expresses complete time. 
      cout << "Verified: " << VerifyOrder(inputvector, less<int>{}) << endl << endl;

      begin = Clock::now(); // new time 
      MergeSort(inputvector, less<int>{});   //same as above except for merge
      end = Clock::now(); // new stop time 
      
      cout << "MergeSort" << endl << endl;
      cout << "Runtime: " << ComputeDuration(begin, end) << "ns" << endl; // expresses complete time. 
      cout << "Verified: " << VerifyOrder(inputvector, less<int>{}) << endl << endl;

      begin = Clock::now();   // repeat process from previous but now for quicksort.
      QuickSort(inputvector, less<int>{});
      end = Clock::now();
      
      cout << "QuickSort" << endl << endl;
      cout << "Runtime: " << ComputeDuration(begin, end) << "ns" << endl;
      cout << "Verified: " << VerifyOrder(inputvector, less<int>{}) << endl << endl;

      //Q2
      cout << "Testing Quicksort Pivet Implementations" << endl << endl;  
      begin = Clock::now();   
      QuickSort(inputvector, less<int>{});   // this is a different quicksort. 
      end = Clock::now();
      
      cout << "Median of Three" << endl << endl;  //PDF designated structure
      cout << "Runtime: " << ComputeDuration(begin, end) << "ns" << endl;
      cout << "Verified: " << VerifyOrder(inputvector, less<int>{}) << endl << endl;

      begin = Clock::now();   
      QuickSort2(inputvector, less<int>{});   
      end = Clock::now();
      
      cout << "Middle" << endl << endl;  // semantics to print it out in order. 
      cout << "Runtime: " << ComputeDuration(begin, end) << "ns" << endl;
      cout << "Verified: " << VerifyOrder(inputvector, less<int>{}) << endl << endl;

      begin = Clock::now();   
      QuickSort3(inputvector, less<int>{});   
      end = Clock::now();
      
      cout << "First" << endl << endl;  
      cout << "Runtime: " << ComputeDuration(begin, end) << "ns" << endl;
      cout << "Verified: " << VerifyOrder(inputvector, less<int>{}) << endl << endl;
    }
  else if (comparisontype == "greater") 
    {  //case if comparison type is greater
      auto begin = Clock::now();
      HeapSort(inputvector, greater<int>{});
      auto end = Clock::now();
      
      cout << "Heapsort" << endl << endl;
      cout << "Runtime: " << ComputeDuration(begin, end) << "ns" << endl;
      cout << "Verified: " << VerifyOrder(inputvector, greater<int>{}) << endl << endl;

      begin = Clock::now();
      MergeSort(inputvector, greater<int>{});
      end = Clock::now();
     
      cout << "MergeSort" << endl << endl;
      cout << "Runtime: " << ComputeDuration(begin, end) << "ns" << endl;
      cout << "Verified: " << VerifyOrder(inputvector, greater<int>{}) << endl << endl;

      begin = Clock::now();
      QuickSort(inputvector, greater<int>{});
      end = Clock::now();
      
      cout << "QuickSort" << endl << endl;
      cout << "Runtime: " << ComputeDuration(begin, end) << "ns" << endl;
      cout << "Verified: " << VerifyOrder(inputvector, greater<int>{}) << endl << endl;

      //Q2
      cout << "Testing Quicksort Pivet Implementations" << endl << endl;
      begin = Clock::now();
      QuickSort(inputvector, greater<int>{});
      end = Clock::now();
      
      cout << "Median of Three" << endl << endl;  //PDF designated structure
      cout << "Runtime: " << ComputeDuration(begin, end) << "ns" << endl;
      cout << "Verified: " << VerifyOrder(inputvector, greater<int>{}) << endl << endl;

      begin = Clock::now();
      QuickSort2(inputvector, greater<int>{});
      end = Clock::now();
      
      cout << "Middle" << endl << endl;  //PDF designated structure
      cout << "Runtime: " << ComputeDuration(begin, end) << "ns" << endl;
      cout << "Verified: " << VerifyOrder(inputvector, greater<int>{}) << endl << endl;

      begin = Clock::now();
      QuickSort3(inputvector, greater<int>{});
      end = Clock::now();
      
      cout << "First" << endl << endl;  //PDF designated structure
      cout << "Runtime: " << ComputeDuration(begin, end) << "ns" << endl;
      cout << "Verified: " << VerifyOrder(inputvector, greater<int>{}) << endl << endl;
    }
}


// Do not change anything below

int main(int argc, char **argv) 
{
  if (argc != 4) 
  {
    cout << "Usage: " << argv[0] << "<input_type> <input_size> <comparison_type>" << endl;
    return 0;
  }
  
  testSortingWrapper(argc, argv);

  return 0;
}
