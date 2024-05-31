#include <stdio.h>

int main() {

  // Declare an integer variable
  int num = 23;

  // Declare a pointer to an integer (int*)
  // This variable will store memory addresses
  int *ptr;

  // Print the value of the variable 'num'
  printf("Value of the variable 'num': %d\n", num);

  // Get the memory address of 'num' using the '&' operator
  // and store it in the pointer 'ptr'
  ptr = &num;

  // Print the memory address stored in the pointer 'ptr'
  printf("Memory address stored in 'ptr': %p\n", ptr);

  // Access the value stored at the memory address pointed to by 'ptr'
  // We use the '*' operator (dereference) to get the value
  printf("Value accessed through the pointer 'ptr': %d\n", *ptr);

  // Modify the value at the memory address pointed to by 'ptr'
  *ptr = 87;

  // Print the value of the variable 'num' again (it has changed)
  printf("Value of the variable 'num' after modification through pointer: %d\n", num);

  return 0;
}
