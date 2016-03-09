#include <stdint.h>
#include <Float64.h>

extern const f64 pio2;
#define pi (pio2*f64(2))

/* basic functions */
f64 sin64(f64 x);
f64 atan64(f64 z);
f64 log64(f64 z);
f64 exp64(f64 z);

/* derived functions */
f64 cos64(f64 z);
f64 tan64(f64 z);
f64 asin64(f64 z);
f64 acos64(f64 z);
f64 atan264(f64 y, f64 x);
f64 fact64(int16_t n);
f64 abs64(f64 z);
f64 cosh64(f64 z);
f64 sinh64(f64 z);
f64 tanh64(f64 z);
f64 acosh64(f64 z);
f64 asinh64(f64 z);
f64 atanh64(f64 z);

f64 sqrt64(f64 z);
