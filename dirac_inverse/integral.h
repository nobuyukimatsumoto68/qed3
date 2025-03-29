#include <gsl/gsl_integration.h>


// https://martin-ueding.de/posts/passing-a-c-lambda-to-a-gsl-function/
// double unwrap(double x, void *p) {
//   auto fp = static_cast<std::function<double(double)> *>(p);
//   return (*fp)(x);
// }

double unwrap(double x, void *p) {
  auto fp = static_cast<std::function<double(double)> *>(p);
  return (*fp)(x);
}

