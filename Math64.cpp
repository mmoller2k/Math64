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
  
  if (z < 0.0) {
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
    if (e == 0.0) {
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
  f64 rootThree = sqrt64(3.0);
  int16_t i;
  
  // negative z
  if (z<0.0) {
    m = true;
    z = -z;
  }

  // special case and for fast answers
  if (z==1.0) { // ArcTan[1] = Pi/4
    if(m) return -(pio2/2.0);
    return pio2/2.0;
  }

  /*
    if (z==recipRootThree) { // ArcTan[1/Sqrt[3]] = pi/6.0
    if(m) return -pi/6.0;
    return pi/6.0;
    }
  */  
  
  // range reduction to (-2+Sqrt[3],2-Sqrt[3])
  // optimised repeated application of (z-1/Sqrt[3])/(1+z*1/Sqrt[3])
  if (z > rootThree+2.0) {
    f = 3.0;
    z = f64(-1)/z;
  }
  else if (z > 1.0) {
    f = 2.0;
    z = (z-rootThree)/(rootThree*z+1.0);
  } 
  else if (z > f64(2.0)-rootThree) {
    f = 1.0;
    z = (z*3.0-rootThree)/(rootThree*z+3.0);
  }
  else { // else no range reduction necessary 
    f = 0.0;
  }
	
  // taylor series
  v = z;
  n = z;
  s = -(z*z);
  for (i=3;;i+=2) {
    n *= s;
    e = n/f64(i);
    if (e == 0.0) {
      if(m) return (-(f*pio2/3.0)+v);
      return f*pio2/3.0+v;
    }
    v += e;
  }
}

f64 exp64(f64 z)
{
  int16_t n=0;
  f64 v,e;
  f64 f=1.0;
  int16_t i;

  // range reduction to (-16,inf)
  if(z<-800.0)return 0.0; //too small
  if(z<-16.0){
    n = ((-z)/16.0).ipart();
    z += n*16.0;
    f = exp64(-16.0);
  }
  
  // taylor series
  // this fails if z<-18 or so
  v = 1;
  e = 1;
  for (i=1;;i++) {
    e *= z/f64(i);
    if (e == 0.0) {
      for(i=0;i<n;i++){
	v *= f;
      }
      return v;
    }
    v += e;
  }
}

f64 cos64(f64 z)
{
    z = sin64(pio2-z);
    return z;
}

f64 asin64(f64 z)
{
  // Complex BigNumbers not implemented
  if(z>1.0) return 0.0;
  if(z<-1.0) return 0.0;
  // Special cases
  if(z==1.0) return pio2;
  if(z==-1.0) return pio2*f64(-1);
	
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
  if(y==0.0 && x==0.0) return 0.0; // undefined
  if(x==0.0 && y==0.0) return pio2*f64(-1);
  if(x==0.0 && y>0.0) return pio2;	
  if(y<0.0 && x<0.0) return atan64(y/x)-pio2*2.0;
  if(x<0.0) return atan64(y/x)+pio2*2.0; // y>=zero

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
  if(z<0.0)return -z;
  return z;
}

#if 0
f64 log(f64 z) //taylor series - slow convergence for z<<1
{
  f64 v,e;
  
  v = z-1.0;
  e = v;
  for (int i=2;;i++) {
    e *= -(z-1.0)*(i-1.0)/(i*1.0);
    if (abs64(e) < 1e-199) { /* doesn't get to zero! */
      printf("%d]",i);
      return v;
    }
    v += e;
  }
}
#endif

#ifdef LOG_TAYLOR
//twice as fast
f64 log(f64 z) //modified taylor series - better convergence
{
  f64 v,s,e,y;
  bool m = false;
  int16_t i;

  if(z==1.0)return 0;
  if(z<=0.0)return z.NaN();

  // range reduction (0<z<1)
  if(z>1.0){
    m=true;
    z=f64(1)/z;
  }
  
  y = (z-1)/(z+1);
  v = 1.0;
  e = v;
  for (i=2;;i+=2) {
    e *= (y*y)*(i-1.0)/(i+1.0);
    //if(e == 0.0){
    if (abs64(e) < 1e-199) { /* doesn't get to zero! */
      printf("%d]",i);
      if(m)return -(v*y*2.0);
      return v*y*2.0;
    }
    v += e;
  }
}

#else
f64 log64(f64 z) //Newton's method - very fast
{
  f64 e,v;
  bool m = false;

  if(z==1.0)return 0;
  if(z<=0.0)return z.NaN();

  // range reduction (0<z<1)
  if(z>1.0){
    m=true;
    z=f64(1)/z;
  }
  
  v = z-1.0;
  e = 0.0;
  for (;;) {
    e = (z-exp64(v)) / (z+exp64(v));
    if(e>=0.0){ /* stop error oscillating around answer */
      if(m)return -v;
      return v;
    }
    v+=e*2.0;
  }
}

#endif

f64 sqrt64(f64 z)
{
  return exp64(log64(z)*0.5);
}
