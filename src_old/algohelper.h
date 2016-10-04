#ifndef ALGOHELPER_H
#define ALGOHELPER_H

template<typename T>
class calcsum
{
  private:
    T sum;
  public:
    calcsum() : sum(0) {}
    calcsum(const calcsum &csum) : sum(csum.sum) {}
    
    void operator()(const T& x) { sum += x; }
    T result() { return sum; }
};

#endif
