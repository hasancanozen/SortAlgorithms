import java.time.Instant;
import java.util.Arrays;

import static java.time.temporal.ChronoUnit.*;

public class Main {

    public static void main(String[] args) {

        int[] array = new int[10000];

        MedianSort medianSort= new MedianSort(array.length);

        for (int i=0; i<array.length; i++){

            int n = (int)(java.lang.Math.random()*10000);
            array[i] = (int) n;
            medianSort.insert(n);
        }
        //temporary array that will use in methods
        int[] temp = new int[10000];
        //copy our array into temp array
        copyArray(temp, array);

        //Counting Sort Algorithm
        Instant countingStart = Instant.now();
        countingSort(temp, temp.length);
        Instant countingEnd = Instant.now();
        long countingDiff = MICROS.between(countingStart,countingEnd);
        System.out.println("Counting Sort \nDuration Time : " + countingDiff +  " µs"+"\n");

        //reset sorted temp array elements
        copyArray(temp, array);

        //Binary Insertion Sort Algorithm
        Instant binaryStart = Instant.now();
        binaryInsSort(temp);
        Instant binaryEnd = Instant.now();
        long binaryDiff = MICROS.between(binaryStart,binaryEnd);
        System.out.println("Binary Insertion Sort \nDuration Time : " + binaryDiff +  " µs"+"\n");

        copyArray(temp, array);

        //Insertion Sort Algorithm
        Instant insertionStart = Instant.now();
        insertionSort(temp);
        Instant insertionEnd = Instant.now();
        long insertionDiff = MICROS.between(insertionStart, insertionEnd);
        System.out.println("Insertion Sort \nDuration Time : " + insertionDiff +  " µs"+"\n");

        copyArray(temp, array);

        //Heap Sort Algorithm
        Instant heapStart = Instant.now();
        heapSort(temp);
        Instant heapEnd = Instant.now();
        long heapDiff = MICROS.between(heapStart, heapEnd);
        System.out.println("Heap Sort \nDuration Time : " + heapDiff +  " µs"+"\n");

        copyArray(temp, array);

        //Merge Sort Algorithm
        Instant mergeStart = Instant.now();
        mergeSort(temp, 0, temp.length-1);
        Instant mergeEnd = Instant.now();
        long mergeDiff = MICROS.between(mergeStart, mergeEnd);
        System.out.println("Merge Sort \nDuration Time : " + mergeDiff + " µs"+"\n");

        copyArray(temp, array);

        //Quick Sort Algorithm (Pivot First)
        Instant quickFirstStart = Instant.now();
        quickSortFirst(temp, 0, temp.length-1);
        Instant quickFirstEnd = Instant.now();
        long quickFirstDiff = MICROS.between(quickFirstStart, quickFirstEnd);
        System.out.println("Quick Sort First \nDuration Time : " + quickFirstDiff + " µs"+"\n");

        copyArray(temp, array);

        //Quick Sort Algorithm (Pivot Median)
        Instant quickMedianStart = Instant.now();
        medianSort.quickSort();
        Instant quickMedianEnd = Instant.now();
        long quickMedianDiff = MICROS.between(quickMedianStart, quickMedianEnd);
        System.out.println("Quick Sort Median \nDuration Time : " + quickMedianDiff + " µs"+"\n");

    }

    private static void copyArray(int[] temp, int[] array){
        System.arraycopy(array, 0, temp, 0, array.length);
    }

    //Method for displaying desired array elements
    private static void printArray(int[] array) {
        int n = array.length;
        System.out.println();
        for (int i = 0; i < n; ++i)
            System.out.print(array[i] + " ");
        System.out.println();
    }

    //Method for binary insertion algorithm
    public static void binaryInsSort(int[] array){

        int index = array.length;

        for (int i=1; i<index; i++){
            int x = array[i];
            //eklediğimiz yeri belirleme ve binarysearch kullanımı
            int number = Math.abs(Arrays.binarySearch(array, 0, i ,x) + 1);
            //shifting array
            System.arraycopy(array, number, array, number+1, i-number);
            array[number] =  x;
        }
    }

    //Method for counting sort algorithm
    public static void countingSort(int[] array, int size){

        int[] output = new int[size + 1];

        int max = array[0];
        for (int i = 1; i < size; i++) {
            if (array[i] > max)
                max = array[i];
        }
        int[] count = new int[max + 1];

        //first for loop
        for (int i = 0; i < max; ++i) {
            count[i] = 0;
        }
        //second for loop
        for (int i = 0; i < size; i++) {
            count[array[i]]++;
        }
        //third for loop
        for (int i = 1; i <= max; i++) {
            count[i] += count[i - 1];
        }
        //fourth for loop
        for (int i = size - 1; i >= 0; i--) {
            output[count[array[i]] - 1] = array[i];
            count[array[i]]--;
        }

        if (size >= 0) System.arraycopy(output, 0, array, 0, size);

    }

    //Method for the main heap sort algorithm
    public static void heapSort(int[] array){


        //arrayi heap'e göre şekillendirme
        for (int i=array.length/2 - 1; i>=0; i--){
            heapify(array, array.length, i);
        }

        //heapteki elemanları ayıklıyoruz
        for (int i = array.length - 1; i > 0; i--) {
            int temp = array[0];
            array[0] = array[i];
            array[i] = temp;

            heapify(array, i, 0);
        }
    }

    //Method needed to complete heap sort algorithm
    private static void heapify(int[] array, int index, int i) {

        int largest = i;
        int left = 2 * i + 1;
        int right = 2 * i + 2;

        //sol çocuk roottan büyükse
        if (left < index && array[left] > array[largest])
            largest = left;

        //sağ çocuk roottan büyükse
        if (right < index && array[right] > array[largest]){
            largest = right;
        }

        //largest eğer root değilse
        if (largest != i){

            int swap = array[i];
            array[i] = array[largest];
            array[largest] = swap;

            heapify(array, index, largest);
        }
    }

    //Method for insertion sort algorithm
    public static void insertionSort(int[] array){

        int index = array.length;

        for (int i=1; i<index; i++){
            int key = array[i];
            int j = i-1;

            while (j >= 0 && array[j] > key){
                array[j+1] = array[j];
                j = j-1;
            }

            array[j+1] = key;
        }
    }

    //Method needed to complete merge sort algorithm
    private static void merge(int[] array, int p, int q, int r){
        int n1 = q - p + 1;
        int n2 = r - q;

        int[] L = new int[n1];
        int[] M = new int[n2];

        System.arraycopy(array, p, L, 0, n1);
        for (int j = 0; j < n2; j++)
            M[j] = array[q + 1 + j];

        int i, j, k;
        i = 0;
        j = 0;
        k = p;

        while (i < n1 && j < n2) {
            if (L[i] <= M[j]) {
                array[k] = L[i];
                i++;
            } else {
                array[k] = M[j];
                j++;
            }
            k++;
        }

        while (i < n1) {
            array[k] = L[i];
            i++;
            k++;
        }

        while (j < n2) {
            array[k] = M[j];
            j++;
            k++;
        }
    }

    //Method for main merge sort algorithm
    public static void mergeSort(int[] array, int left, int right) {
        if (left < right) {
            int mid = (left + right) / 2;
            mergeSort(array, left, mid);
            mergeSort(array, mid + 1, right);
            merge(array, left, mid, right);
        }
    }

    //method for choosing pivot element of quick sort algorithm
    private static int partitionFirst(int[] array, int low, int high) {

        //pivot first element
        int pivot = array[low];
        int i = (low + 1);

        for (int j = low+1; j <= high; j++) {
            if (array[j] < pivot) {
                if (j != i){
                    int temp = array[i];
                    array[i] = array[j];
                    array[j] = temp;
                }
                i++;
            }
        }
        array[low] = array[i-1];
        array[i-1] = pivot;

        return (i - 1);
    }

    //Method for main quick sort algorithm (pivot first)
    private static void quickSortFirst(int[] array, int low, int high) {
        if (low < high) {
            int pi = partitionFirst(array, low, high);
            quickSortFirst(array, low, pi - 1);
            quickSortFirst(array, pi + 1, high);
        }
    }

    //Implementation of quick sort algorithm (pivot median)
    private static class MedianSort {

        private int[] theArray;
        private int nElems;

        public MedianSort(int max) {
            theArray = new int[max];
            nElems = 0;
        }
        public void insert(int value) {
            theArray[nElems] = value;
            nElems++;
        }

        public void quickSort() {
            recQuickSort(0, nElems-1);
        }

        public void recQuickSort(int left, int right) {
            int size = right-left+1;
            if(size <= 3)
                manualSort(left, right);
            else {
                long median = medianOf3(left, right);
                int partition = partitionIt(left, right, median);
                recQuickSort(left, partition-1);
                recQuickSort(partition+1, right);
            }
        }

        public long medianOf3(int left, int right) {
            int center = (left+right)/2;

            if( theArray[left] > theArray[center] )
                swap(left, center);

            if( theArray[left] > theArray[right] )
                swap(left, right);

            if( theArray[center] > theArray[right] )
                swap(center, right);

            swap(center, right-1);
            return theArray[right-1];
        }
        public void display(){
            for(int j=0; j<nElems; j++)
                System.out.print(theArray[j] + " ");
            System.out.println();
        }

        public void swap(int dex1, int dex2) {
            int temp = theArray[dex1];
            theArray[dex1] = theArray[dex2];
            theArray[dex2] = temp;
        }

        public int partitionIt(int left, int right, long pivot) {
            int leftPtr = left;
            int rightPtr = right - 1;

            while(true)
            {
                while( theArray[++leftPtr] < pivot ) ;
                while( theArray[--rightPtr] > pivot ) ;
                if(leftPtr >= rightPtr)
                    break;
                else
                    swap(leftPtr, rightPtr);
            }
            swap(leftPtr, right-1);
            return leftPtr;
        }

        public void manualSort(int left, int right) {
            int size = right-left+1;
            if(size <= 1)
                return;
            if(size == 2)
            {
                if( theArray[left] > theArray[right] )
                    swap(left, right);
                return;
            }
            else
            {
                if( theArray[left] > theArray[right-1] )
                    swap(left, right-1);
                if( theArray[left] > theArray[right] )
                    swap(left, right);
                if( theArray[right-1] > theArray[right] )
                    swap(right-1, right);
            }
        }
    }
}
