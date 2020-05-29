#include <iostream>
#include <stdlib.h>
using namespace std;

int main(int argc, char* argv[]) {
  int cores = atoi(argv[1]);
  int times = atoi(argv[2]);
  for (int j = 4; j <= times; j++) {
    int n = cores * j;
    cout << j << " " << n << '\n';
    cout << "1 ";
    for (int i = 2; i <= cores; ++i)
      if (n % i == 0 && (i & 1) == 0) {
        cout << i << " ";
      }
    cout << '\n';
  }

    return 0;
}