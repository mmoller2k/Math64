#include <stdio.h>
#include <stdint.h>
#include <Math64.h>

const f64 pio2 = f64( 0x3FF921FB, 0x54442D18); //3.1415926535897932/2.0;
//const f64 rootThree = 1.7320508075688772;

f64 sin64(f64 z)
{
  int32_t n;
  f64 v,s,e;
  int16_t i,j;
  bool m = false; // m= false => z positive, m= true => z negative)
  
  if(z.isNaN())return z;
  if (z.isNegative()) {
    m = true;
    z = -z;
  }
  
  // range reduction to (-Pi/2,Pi/2)
  n = ((z/pio2).ipart()+1)/2;
  z -= pio2*f64(n*2);
  if (n%2) z = -z;

  // taylor series
  v = z;
  e = z;
  s = -(z*z);
  for (i=3;;i+=2) {
    j=i-1;
    e /= f64(i*j);
    e *= s;
    if (e.isZero()) {
      if (m) return -v;
      return v;
    }
    v += e;
  }
}

f64 atan64(f64 z)
{
  f64 f,n,s,e,v;
  bool m = false;   // m= false => z positive, m= true => z negative)
  f64 rootThree = sqrt64(3);
  int16_t i;
  
  if(z.isNaN()) return z;

  // negative z
  if (z.isNegative()) {
    m = true;
    z = -z;
  }

  // special case and for fast answers
  if (z==1) { // ArcTan[1] = Pi/4
    if(m) return -(pio2/f64(2));
    return pio2/f64(2);
  }

  /*
    if (z==recipRootThree) { // ArcTan[1/Sqrt[3]] = pi/6.0
    if(m) return -pi/6.0;
    return pi/6.0;
    }
  */  
  
  // range reduction to (-2+Sqrt[3],2-Sqrt[3])
  // optimised repeated application of (z-1/Sqrt[3])/(1+z*1/Sqrt[3])
  if (z > rootThree+f64(2)) {
    f = f64(3);
    z = f64(-1)/z;
  }
  else if (z > 1) {
    f = f64(2);
    z = (z-rootThree)/(rootThree*z+f64(1));
  } 
  else if (z > f64(2)-rootThree) {
    f = f64(1);
    z = (z*f64(3)-rootThree)/(rootThree*z+f64(3));
  }
  else { // else no range reduction necessary 
    f = f64(0);
  }
	
  // taylor series
  v = z;
  n = z;
  s = -(z*z);
  for (i=3;;i+=2) {
    n *= s;
    e = n/f64(i);
    if (e.isZero()) {
      if(m) return (-(f*pio2/f64(3))+v);
      return f*pio2/f64(3)+v;
    }
    v += e;
  }
}

f64 exp64(f64 z)
{
  int i,n=0;
  f64 v,e;
  f64 f=f64(0x3e7e3139, 0xb54f3b84); //e^-16

  if(z.isNaN())return z;
  if(z.isZero())return 1;

  if(z<-800)return 0; //too small
  // range reduction to (-16,inf)
  if(z<-16){
    n = -z.ipart()/16;
    z += f64(n*16);
  }
  
  // taylor series
  // this fails if z<-18 or so
  v = 1;
  e = 1;
  for (i=1;;i++) {
    e *= z/f64(i);
    if(e.isZero()){
      for(i=0;i<n;i++){
	v *= f;
      }
      break;
    }
    v += e;
  }
  return v;
}

f64 cos64(f64 z)
{
    z = sin64(pio2-z);
    return z;
}

f64 tan64(f64 z)
{
    z = sin64(z)/cos64(z);
    return z;
}

f64 asin64(f64 z)
{
  // Complex numbers not implemented
  if((z>1) || (z<-1)) return z.setNaN();

  // Special cases
  if(z==1) return pio2;
  if(z==-1) return pio2*f64(-1);
	
  z=z/sqrt64(f64(1)-z*z);
	
  return atan64(z);
}

f64 acos64(f64 z)
{
  return pio2 - asin64(z);
}

f64 atan264(f64 y, f64 x)
{
  // Special cases
  if(y.isZero() && x.isZero()) return 0; // undefined
  if(x.isZero() && y.isNegative()) return pio2*f64(-1);
  if(x.isZero() && !y.isNegative()) return pio2;	
  if(y.isNegative() && x.isNegative()) return atan64(y/x)-pio2*f64(2);
  if(x.isNegative()) return atan64(y/x)+pio2*f64(2); // y>=zero
  if(y.isNegative() && !x.isNegative()) return -atan64(y/x);
  return atan64(y/x); // x>zero
}

f64 fact64(int16_t n)
{
  int16_t i;
  f64 f=1;
  for(i=2;i<=n;i++){
    f *= i;
  }
  return f;
}

f64 abs64(f64 z)
{
  if(z.isNegative())return -z;
  return z;
}

#if 0
f64 log64(f64 z) //taylor series - slow convergence for z<<1
{
  f64 v,e;
  
  v = z-f64(1);
  e = v;
  for (int i=2;;i++) {
    e *= -(z-f64(1))*(i-f64(1))/f64(i);
    if (abs64(e) < 1e-199) { /* doesn't get to zero! */
      printf("%d]",i);
      return v;
    }
    v += e;
  }
}
#endif

//#define LOG_TAYLOR
#ifdef LOG_TAYLOR
//twice as fast
f64 log64(f64 z) //modified taylor series - better convergence
{
  f64 v,s,e,y;
  bool m = false;
  int16_t i;

  if(z==1)return 0;
  if(z.isNegative() || z.isZero())return z.setNaN();

  // range reduction (0<z<1)
  if(z>1){
    m=true;
    z=f64(1)/z;
  }
  
  y = (z-f64(1))/(z+f64(1));
  v = 1;
  e = v;
  for (i=2;;i+=2) {
    e *= (y*y)*(f64(i)-f64(1))/(f64(i)+f64(1));
    //if(e == 0.0){
    if (abs64(e) < f64(0xe940b8e0,0xacac4ec6)) { /* e < 1e-199: doesn't get to zero! */
      printf("%d]",i);
      if(m)return -(v*y*f64(2));
      return v*y*f64(2);
    }
    v += e;
  }
}

#else
f64 log64(f64 z) //Newton's method - very fast
{
  f64 e,v;
  bool m = false;

  if(z==1)return 0;
  if(z.isZero() || z.isNegative() || z.isNaN())return z.setNaN();

  // range reduction (0<z<1)
  if(z>1){
    m=true;
    z=f64(1)/z;
  }
  
  v = z-f64(1);
  e = f64(0);
  for (;;) {
    e = (z-exp64(v)) / (z+exp64(v));
    if(!e.isNegative()){ /* stop error oscillating around answer */
      if(m)return -v;
      return v;
    }
    v+=e*f64(2);
  }
}

#endif

f64 sqrt64(f64 z)
{
  f64 result;
  if(z.isNegative() || z.isNaN())return result.setNaN();
  return exp64(log64(z)/f64(2));
}

f64 sinh64(f64 z)
{
  return (exp64(z) - exp64(-z))/f64(2);
}

f64 cosh64(f64 z)
{
  return (exp64(z) + exp64(-z))/f64(2);
}

f64 tanh64(f64 z)
{
  return sinh64(z)/cosh64(z);
}

f64 asinh64(f64 z)
{
  return log64(z+sqrt64(z*z+f64(1)));
}
f64 acosh64(f64 z)
{
  return log64(z+sqrt64(z*z-f64(1)));
}
f64 atanh64(f64 z)
{
  return log64((f64(1)+z)/(f64(1)-z))/f64(2);
}

