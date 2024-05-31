#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <omp.h>

#define N 15

//int sum[N]={1, 4, -2, -2, 5, -4, 3, 8}; 
int sum[N];


int contando(int arr[], int n)
{
    int cont=0;
    for (int i = 0; i < n; i++) {

        // starting point of the sub arrray and sum is initialized 
        // with first element of subarray a[i]
        int sum = arr[i];
        if (sum == 0)
            cont+=1;

        for (int j = i + 1; j < n; j++) {
        // we are finding the sum till jth index starting from ith index
            sum += arr[j];
            if (sum == 0)
                cont+=1;
        }
    }
    return cont;
}
int contandoenparalelo(int arr[], int n)
{
    int cont = 0;

    #pragma omp parallel
    {
        int local_cont = 0;
        #pragma omp for
        for (int i = 0; i < n; i++) {
            // starting point of the subarray and sum is initialized
            // with the first element of subarray a[i]
            int sum = arr[i];
            if (sum == 0)
                local_cont += 1;

            for (int j = i + 1; j < n; j++) {
                // we are finding the sum till jth index starting from ith index
                sum += arr[j];
                if (sum == 0)
                    local_cont += 1;
            }
        }

                                                                                                                                                                                          1,18          Top

        #pragma omp atomic
        cont += local_cont;
    }

    return cont;
}
int contandoenparalel(int arr[], int n)
{
    int cont=0;
    #pragma omp for
    for (int i = 0; i < n; i++) {

        // starting point of the sub arrray and sum is initialized 
        // with first element of subarray a[i]
        int sum = arr[i];
        if (sum == 0)
            cont+=1;

        for (int j = i + 1; j < n; j++) {
        // we are finding the sum till jth index starting from ith index
            sum += arr[j];
            if (sum == 0)
                cont+=1;
        }
    }
    return cont;
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

int seed=2;
int x, i, upper=10, lower=-10;

  srand(seed);
  printf ("\n");
  for (i=0; i < N; i++){
     x = (rand() % (upper - lower + 1)) + lower;
     sum[i] = x;
     printf ("%d ", x);
   }
  int sume[7]={1, -1, 0, 0, -3, 5, -2};
  printf("\n %d %d %d %d %d %d %d ", sume[0], sume[1], sume[2], sume[3], sume[4], sume[5], sume[6]);
  int t=contando(sume,7);
  int d=contandoenparalelo(sume,7);
  printf("\n \n %d \n \n",t);
  printf("\n \n %d \n \n",d);

  if (subArrayExists(sum, N))
        printf ("\n \n ***********  Found a subarray with 0 sum ************ \n \n");
  else
        printf ("\n \n ***********  No Such Sub Array Exists! ************ \n \n");

  return (0);
}
