#include <stdio.h>

int main() {

  // ---------- Integer Array and Pointer Arithmetic ----------

  // Declare an integer array `int_arr` of size 5
  int int_arr[5] = {10, 20, 30, 40, 50};

  // Declare a pointer `int_ptr` of type `int*` to point to integers
  // and initialize it to point to the beginning of the integer array.
  int *int_ptr = int_arr;

  // Print the initial address of the array (address of int_arr[0])
  printf("Address of int_arr[0] using int_arr: %p\n",      int_arr   );
  printf("Address of int_arr[0] using &int_arr: %p\n",    &int_arr   );
  printf("Address of int_arr[0] using &int_arr[0]: %p\n", &int_arr[0]);

  // Print the address stored in the pointer (also points to int_arr[0])
  printf("Address stored in int_ptr: %p\n", int_ptr);

  // Access and print the value at the memory location pointed to by int_ptr
  printf("Value at int_arr[0] using dereference (*int_ptr): %d\n", *int_ptr);

  // Move the pointer 2 steps forward (points to int_arr[2]) using pointer arithmetic
  int_ptr += 2;

  // Print the new address stored in the pointer (points to int_arr[2])
  printf("Address stored in int_ptr after arithmetic (points to int_arr[2]): %p\n", int_ptr);

  // Access and print the value at the new memory location pointed to by int_ptr
  printf("Value at int_arr[2] using dereference (*int_ptr): %d\n", *int_ptr);

  // Calculate and print the difference in addresses for integers (typically 4 bytes)
  long int_address_diff = (long)int_ptr - (long)&int_arr;
  printf("Difference in addresses for integers (should be 2 * sizeof(int)): %ld\n\n", int_address_diff);


  // ---------- Character Array and Pointer Arithmetic ----------

  // Declare a character array `char_arr` of size 6 (including null terminator)
  char char_arr[] = "Hello";

  // Declare a pointer `char_ptr` of type `char*` to point to characters
  // and initialize it to point to the beginning of the character array.
  char *char_ptr = char_arr;

  // Print the address of the first character in the array (char_arr[0])
  printf("Address of char_arr[0] using &char_arr: %p\n", &char_arr);

  // Print the address stored in the pointer (also points to char_arr[0])
  printf("Address stored in char_ptr: %p\n", char_ptr);

  // Access and print the character value at the memory location pointed to by char_ptr
  printf("Value at char_arr[0] using dereference (*char_ptr): %c\n", *char_ptr);

  // Move the pointer 2 steps forward (points to char_arr[2]) using pointer arithmetic
  // Since characters are typically 1 byte, this moves 2 bytes forward
  char_ptr += 2;

  // Print the new address stored in the pointer (points to char_arr[2])
  printf("Address stored in char_ptr after arithmetic (points to char_arr[2]): %p\n", char_ptr);

  // Access and print the character value at the new memory location pointed to by char_ptr
  printf("Value at char_arr[2] using dereference (*char_ptr): %c\n", *char_ptr);

  // Calculate and print the difference in addresses for characters (typically 1 byte)
  long char_address_diff = (long)char_ptr - (long)&char_arr;
  printf("Difference in addresses for characters (should be 2 * sizeof(char)): %ld\n", char_address_diff);

  return 0;
}

