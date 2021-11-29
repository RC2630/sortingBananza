#include "vectorUtility.h"

#include <iostream>
#include <numeric>
#include <utility>
#include <algorithm>
#include <ctime>
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <climits>

using namespace std;
using vecUtil::operator<<;

struct Number {

  double value = INT_MAX;
  int numOps = 0;

  bool operator < (const Number& other) {
    numOps++;
    return value < other.value;
  }

};

struct MinHeapQueue {

private:

  vector<Number> heap = {Number()}; // initialized with useless element because index 0 is not used

  void heapifyUp(int cIndex) {
    if (cIndex > 1 && heap.at(cIndex) < heap.at(parent(cIndex))) {
      swap(heap.at(cIndex), heap.at(parent(cIndex)));
      heapifyUp(parent(cIndex));
    }
  }

  void heapifyDown(int cIndex) {
    if (hasChild(cIndex)) {
      int minChildIndex = minimumChildIndex(cIndex);
      if (heap.at(minChildIndex) < heap.at(cIndex)) {
        swap(heap.at(minChildIndex), heap.at(cIndex));
      }
      heapifyDown(minChildIndex);
    }
  }

  static int parent(int cIndex) {
    return cIndex / 2; // integer division
  }

  bool hasChild(int cIndex) const {
    return 2 * cIndex <= size();
  }

  int minimumChildIndex(int cIndex) {
    int leftIndex = 2 * cIndex; // we know for sure that this left index is within bounds
    int rightIndex = leftIndex + 1; // but not for the right one
    if (rightIndex > size()) {
      return leftIndex;
    }
    return (heap.at(leftIndex) < heap.at(rightIndex)) ? leftIndex : rightIndex;
  }

  static MinHeapQueue makeHeap(const vector<Number>& v) {
    MinHeapQueue minHeap;
    minHeap.heap.insert(minHeap.heap.end(), v.begin(), v.end());
    for (int i = parent(minHeap.size()); i > 0; i--) {
      minHeap.heapifyDown(i);
    }
    return minHeap;
  }

public:

  void insert(const Number& data) {
    heap.push_back(data);
    heapifyUp(size());
  }

  Number removeMin() {
    Number minVal = heap.at(1);
    heap.at(1) = heap.back();
    heap.pop_back();
    heapifyDown(1);
    return minVal;
  }

  int size() const {
    return heap.size() - 1;
  }

  static void heapsort(vector<Number>& v) {
    int initialSize = v.size();
    MinHeapQueue minHeap = makeHeap(v);
    v.clear();
    for (int i = 0; i < initialSize; i++) {
      v.push_back(minHeap.removeMin());
    }
  }

};

double logarithm(int n, int b) {
  return log(n) / log(b);
}

int operator + (int n1, const Number& n2) {
  return n1 + n2.numOps;
}

ostream& operator << (ostream& out, const Number& n) {
  out << n.value;
  return out;
}

// ------------------------------ sorting algorithms begin here ------------------------------

void slide(vector<Number>& v, int pos) {
  Number n = v.at(pos);
  int i = pos;
  while (i > 0 && n < v.at(i - 1)) {
    v.at(i) = v.at(i - 1);
    i--;
  }
  v.at(i) = n;
}

void insertionSort(vector<Number>& v) {
  for (int i = 1; i < v.size(); i++) {
    slide(v, i);
  }
}

void bubbleSort(vector<Number>& v) {
  for (int i = v.size() - 1; i > 0; i--) {
    for (int j = 1; j <= i; j++) {
      if (v.at(j) < v.at(j - 1)) {
        swap(v.at(j), v.at(j - 1));
      }
    }
  }
}

void standardSort(vector<Number>& v) {
  sort(v.begin(), v.end());
}

int indexMin(vector<Number>& v, int initialPos) {
  int retloc = initialPos;
  for (int i = initialPos + 1; i < v.size(); i++) {
    if (v.at(i) < v.at(retloc)) {
      retloc = i;
    }
  }
  return retloc;
}

void selectionSort(vector<Number>& v) {
  for (int i = 0; i < v.size(); i++) {
    swap(v.at(i), v.at(indexMin(v, i)));
  }
}

template <typename T>
void merge(vector<T>& v, int low, int mid, int high) {
	int i = low;
	int j = mid + 1;
	vector<T> nv;
	while (true) {
		if (i > mid && j <= high) {
			nv.push_back(v.at(j));
			j++;
		} else if (i <= mid && j > high) {
			nv.push_back(v.at(i));
			i++;
		} else if (i > mid && j > high) {
			break;
		} else if (v.at(i) < v.at(j)) { // by this point, i <= mid && j <= high
			nv.push_back(v.at(i));
			i++;
		} else {
			nv.push_back(v.at(j));
			j++;
		}
	}
  for (int i = low; i <= high; i++) {
    v.at(i) = nv.at(i - low);
  }
	//vector<vector<T>> frontMergeBack = {vecUtil::subvector(v, 0, low - 1), nv, vecUtil::subvector(v, high + 1, v.size() - 1)};
	//v = vecUtil::concatenate(frontMergeBack);
}

template <typename T>
void mergeSort(vector<T>& v, int low, int high) {
	if (high > low) {
		int mid = (low + high) / 2;
		mergeSort(v, low, mid);
		mergeSort(v, mid + 1, high);
		merge(v, low, mid, high);
	}
}

template <typename T>
void mergeSortAll(vector<T>& v) {
	if (v.size() >= 2) { // vectors of size 0 or 1 are vacuously/trivially sorted, so this function just does nothing
		mergeSort(v, 0, v.size() - 1);
	}
}

int partition(vector<Number>& v, int low, int high) {
  int pivot = low;
  for (int i = low + 1; i <= high; i++) {
    if (v.at(i) < v.at(low)) {
      pivot++;
      swap(v.at(pivot), v.at(i));
    }
  }
  swap(v.at(low), v.at(pivot));
  return pivot;
}

void quicksort(vector<Number>& v, int low, int high) {
  if (low < high) {
    int pivot = partition(v, low, high);
    quicksort(v, low, pivot - 1);
    quicksort(v, pivot + 1, high);
  }
}

void quicksortAll(vector<Number>& v) {
  if (v.size() >= 2) {
    quicksort(v, 0, v.size() - 1);
  }
}

int comparePointersToNumbers(const void* a, const void* b) {
  Number* n1 = (Number*) a;
  Number* n2 = (Number*) b;
  return (*n1 < *n2) ? -1 : 1;
}

void standardQuicksort(vector<Number>& v) {
  qsort(&v.at(0), v.size(), sizeof(Number), comparePointersToNumbers);
}

// ------------------------------ sorting algorithms end here ------------------------------

vector<Number> numberify(const vector<double>& v) {
  vector<Number> nums;
  for (double n : v) {
    nums.push_back({n, 0});
  }
  return nums;
}

void clear(vector<Number>& v) {
  for (Number& n : v) {
    n.numOps = 0;
  }
}

template <typename Function>
void sortAndDisplay(const vector<Number>& og, const Function& sortFunc, const string& sortAlgorithm) {
  vector<Number> v = og; // making a copy
  sortFunc(v);
  //cout << "\nv (" << sortAlgorithm << ") = " << v << "\n";
  int totalNumOps = accumulate(v.begin(), v.end(), 0);
  cout << "total # of comparison operations (<) for " << sortAlgorithm << " is on the order of 10^" << logarithm(totalNumOps, 10) << "\n";
}

vector<double> randomVector(int size) {
  vector<double> v;
  for (int i = 1; i <= size; i++) {
    v.push_back(rand() % 201 - 100); // range: from -100 to 100
  }
  return v;
}

int main() {

  cout << "\n" << fixed << setprecision(1);
  srand(time(nullptr));

  int size = 20000;
  vector<double> raw = randomVector(size);
  vector<Number> v = numberify(raw);
  //cout << "v (before sorting) = " << v << "\n";
  cout << "size of vector to sort (n) = " << size << "\n\n";

  sortAndDisplay(v, insertionSort, "insertion sort");
  sortAndDisplay(v, bubbleSort, "bubble sort");
  sortAndDisplay(v, standardSort, "std::sort");
  sortAndDisplay(v, selectionSort, "selection sort");
  sortAndDisplay(v, mergeSortAll<Number>, "merge sort");
  sortAndDisplay(v, MinHeapQueue::heapsort, "heap sort");
  sortAndDisplay(v, quicksortAll, "quick sort");
  sortAndDisplay(v, standardQuicksort, "std::qsort");
  cout << "\n";
  
}