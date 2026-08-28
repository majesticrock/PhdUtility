#ifndef MROCK_SYMBOLIC_OPERATORS_INCLUDE_MROCK_SYMBOLIC_OPERATORS_EXCEPTIONS_HPP
#define MROCK_SYMBOLIC_OPERATORS_INCLUDE_MROCK_SYMBOLIC_OPERATORS_EXCEPTIONS_HPP

#include <stdexcept>
#include <string>

namespace mrock::symbolic_operators {
class momentum_replacement_error : public std::invalid_argument {
public:
    momentum_replacement_error(const char replaceWhat, const std::string& replaceWith);
};

class redistribution_error : public std::invalid_argument {
public:
    redistribution_error(const std::string which);
};

class renaming_error : public std::invalid_argument {
public:
    renaming_error(const std::string which);
};

class sum_transformation_error : public std::invalid_argument {
public:
    sum_transformation_error();
};

}  // namespace mrock::symbolic_operators
#endif  // MROCK_SYMBOLIC_OPERATORS_INCLUDE_MROCK_SYMBOLIC_OPERATORS_EXCEPTIONS_HPP
