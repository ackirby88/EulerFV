#ifndef C_PRINTER_CXX
#define C_PRINTER_CXX

#include <iostream>
#include <iomanip>
#include <fstream>

extern "C" {
    void c_printer(int *mode, double *arr, int *len, int *step);
}

void c_printer(int *mode, double *arr, int *len, int *step) {
  std::string filename;

  if (*mode == 1) {
    filename = "res_" + std::to_string(*step) + ".txt";
  } else if (*mode == 2) {
    filename = "soln_" + std::to_string(*step) + ".txt";
  } else {
    return;  // Handle unexpected mode_id (optional)
  }
  std::ofstream outfile(filename);

  for (int i = 0; i < *len; ++i) {
    outfile << std::uppercase
            << std::scientific     
            << std::setprecision(16)
            << i << " " << arr[i] << std::endl;
  }
}
#endif //C_PRINTER_CXX
