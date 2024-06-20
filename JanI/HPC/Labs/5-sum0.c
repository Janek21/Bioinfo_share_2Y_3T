#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
 
#define N 2000

//int sum[N]={1, 4, -2, -2, 5, -4, 3, 8};
int sum[N];

int sumCounter(int arr[], int n){
    int c=0;
    
    for (int i=0;i<n;i++){

        int sum = arr[i]; //if element is 0 -> +1 to counter
        if (sum == 0){c++;}
            
        for (int j = i + 1; j < n; j++) { //setarting from current position and forward   
            sum += arr[j]; //every time you advance a position add the corresponding number
            if (sum == 0) // if the sum of the numers from i so far is 0, add 1 to the couter
                c++;
        }
    }
    return c;
}

int PsumCounter(int arr[], int n){ //same reasonings as sumCounter
    int i, j, c=0; //i, j defined here for parallels later

    #pragma omp parallel for reduction(+:c) //distributes an i for each process, makes c private and adds them up after the for loop is finished to ensure its correct
    for (i=0;i<n;i++){
        int sum = arr[i];
        if (sum == 0){c++;}

        #pragma omp parallel for reduction(+:c) // distributes a j for each process, makes c private and adds them up when the loop is finished to ensure its correct
        for (j = i + 1; j < n; j++){
            sum+=arr[j];
            if (sum==0){c++;}
        }
    }
    return c;
}

bool subArrayExists(int arr[], int n)
{
    for (int i = 0; i < n; i++) {
       
        // starting point of the sub arrray and sum is initialized 
        // with first element of subarray a[i]
        int sum = arr[i];
        if (sum == 0)
            return true;
            
        for (int j = i + 1; j < n; j++) {        
        // we are finding the sum till jth index starting from ith index
            sum += arr[j];
            if (sum == 0)
                return true; 
        }
    }
    return false;
}  
    
int main ()  {

    int seed=2; //2
    int x, i, upper=10, lower=-10;
    
    srand(seed);
    printf ("\n"); 
    for (i=0; i < N; i++){
        x = (rand() % (upper - lower + 1)) + lower;
        sum[i] = x;
        //printf ("%d ", x);
    }
    //int sum[8]={-9, -6, -1, 9, -2, 0, 0, -1};

    if (subArrayExists(sum, N))
            printf ("\n \n ***********  Found a subarray with 0 sum ************ \n \n"); 
    else
            printf ("\n \n ***********  No Such Sub Array Exists! ************ \n \n"); 

    int length_sum = sizeof(sum) / sizeof(sum[0]);
    int Seq_c=sumCounter(sum, length_sum);
    printf ("\n \n ****  Found %d subarrays with 0 sum using a sequential program  **** \n \n", Seq_c);

    int P_c=PsumCounter(sum, length_sum);
    printf ("\n \n *****  Found %d subarrays with 0 sum using a parallel program  ***** \n \n", P_c);
    return 0;
}

