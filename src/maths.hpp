#include "../huerto/maths/random.hpp"

namespace parserFunc {

  double random(double) {
    return random_unit_interval(rng);
  }

}