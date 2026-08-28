#include <mrock/symbolic_operators/Exceptions.hpp>

#include <stdexcept>
#include <string>

namespace mrock::symbolic_operators {

momentum_replacement_error::momentum_replacement_error(const char replaceWhat, const std::string& replaceWith)
    : std::invalid_argument(std::string("You are trying to replace a momentum with itself. The task was: replace <") +
                            replaceWhat + std::string("> with <") + replaceWith + ">.") {}

redistribution_error::redistribution_error(const std::string which)
    : std::invalid_argument("There is no summation that would allow the desired " + which + " transformation!") {}

renaming_error::renaming_error(const std::string which)
    : std::invalid_argument("You are replacing a " + which + " sum with an index that already exists!") {}

sum_transformation_error::sum_transformation_error()
    : std::invalid_argument(
          "You are trying to perform a sum transformation on a momentum that is not being summed over!") {}

}  // namespace mrock::symbolic_operators