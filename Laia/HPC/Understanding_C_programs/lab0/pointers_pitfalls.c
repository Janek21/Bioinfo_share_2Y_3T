#include <stdio.h>

int main() {

  // Declare an array of characters  
  char num[10] = "123456789";

  // Declare a pointer to an integer (int*)
  // This variable will store memory addresses
  int *ptr;

  // Print the value of the variable 'num'
  printf("Value of the variable 'num': %s\n", num);

  // Get the memory address of 'num' using the '&' operator
  // and store it in the pointer 'ptr'
  ptr = &num;

  // Print the memory address stored in the pointer 'ptr'
  printf("Memory address stored in 'ptr': %p\n", ptr);

  // Access the value stored at the memory address pointed to by 'ptr'
  // We use the '*' operator (dereference) to get the value
  printf("Value accessed through the pointer 'ptr': %d\n", *ptr);

  // Modify the value at the memory address pointed to by 'ptr'
  *ptr = -43210;

  // Print the value of the variable 'num' again (it has changed)
  printf("Changing the value through the pointer 'ptr' to: %d\n", *ptr);

  // Print the value of the variable 'num' again (it has changed)
  printf("Value of the variable 'num' after modification through pointer: %s\n", num);

  return 0;
}
