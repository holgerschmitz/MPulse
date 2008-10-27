#include "../src/hdfstream.h"

#include <schnek/matrix.h>
#include <schnek/functions.h>

#include <cmath>

int main()
{
  typedef schnek::Matrix<double, 2> Matrix2;
  typedef schnek::Matrix<double, 3> Matrix3;
  typedef Matrix2::IndexType Index2;
  typedef Matrix3::IndexType Index3;
  
  Index2 lowA(2,3), highA(6,8);
  Index3 lowB(4,5,6), highB(8,9,10);
  
  Matrix2 A(lowA, highA);
  Matrix3 B(lowB, highB);
  
  for (int i=lowA[0]; i<=highA[0]; ++i)
    for (int j=lowA[1]; j<=highA[1]; ++j)
      A(i,j) = exp(j-i)*sqrt(j);
    
  for (int i=lowB[0]; i<=highB[0]; ++i)
    for (int j=lowB[1]; j<=highB[1]; ++j)
      for (int k=lowB[2]; k<=highB[2]; ++k)
        B(i,j,k) = A(i,j)*A(j,k)*sin(k);

  
  HDFostream hdfout("test.h5");
  hdfout << A << B;
  hdfout.close();
  
  Matrix2 C;
  Matrix3 D;
  
  HDFistream hdfin("test.h5");
  hdfin >> C >> D;
  hdfin.close();
  
  Index2 lowC = C.getLow();
  Index2 highC = C.getHigh();
  
  Index3 lowD = D.getLow();
  Index3 highD = D.getHigh();

  if (lowC != lowA) {
    std::cerr << "lowC not restored!  (" << lowC[0] << ", " << lowC[1] << ")\n";
    exit(-1);
  }
  
  if (lowD != lowB) {
    std::cerr << "lowD not restored!  (" << lowD[0] << ", " << lowD[1] << ", " << lowD[2] << ")\n";
    exit(-1);
  }
  
  for (int i=lowC[0]; i<=highC[0]; ++i)
    for (int j=lowC[1]; j<=highC[1]; ++j)
      if (C(i,j) != A(i,j)) {
        std::cerr << "Matrix C not restored!  C("<<i<<", "<<j<<") = "<< C(i,j) << " != " << A(i,j) << "\n";
        exit(-1);
      }

  for (int i=lowD[0]; i<=highD[0]; ++i)
    for (int j=lowD[1]; j<=highD[1]; ++j)
      for (int k=lowD[2]; k<=highD[2]; ++k)
        if (D(i,j,k) != B(i,j,k)) {
          std::cerr << "Matrix D not restored!  D("<<i<<", "<<j<<", "<<k<<") = "<< D(i,j,k) << " != " << B(i,j,k) << "\n";
          exit(-1);
        }
  

}
