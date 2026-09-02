#ifndef MROCK_SYMBOLIC_OPERATORS_INCLUDE_MROCK_SYMBOLIC_OPERATORS_WICKTERMCOLLECTOR_HPP
#define MROCK_SYMBOLIC_OPERATORS_INCLUDE_MROCK_SYMBOLIC_OPERATORS_WICKTERMCOLLECTOR_HPP

#include "AbstractCollector.hpp"
#include "WickTerm.hpp"

#include <memory>
#include <vector>

namespace mrock::symbolic_operators {

struct WickSymmetry;

/**
 * @class WickTermCollector
 * @brief A wrapper for a vector of WickTerm objects.
 */
struct WickTermCollector : public AbstractCollector<WickTerm> {
    /**
     * @brief Serializes the WickTermCollector object.
     *
     * @tparam Archive The type of the archive.
     * @param ar The archive object.
     * @param version The version of the serialization.
     */
    template <class Archive>
    void serialize(Archive& ar, [[maybe_unused]] const unsigned int version) {
        ar & terms;
    };

    MROCK_FORWARD_CONSTRUCTORS(WickTermCollector, AbstractCollector<WickTerm>);

    MROCK_FORWARD_ASSIGNMENT(WickTermCollector, terms);

    /**
     * @brief Clears eta terms from the WickTermCollector. Intended for use if <eta>=0
     */
    void clear_etas();

    /**
     * @brief Cleans Wick terms using the provided symmetries.
     *
     * @param symmetries The vector of unique pointers to WickSymmetry objects.
     */
    void clean_up();
    void clean_up(const std::vector<std::unique_ptr<WickSymmetry>>& symmetries);
};

/**
 * @brief Addition assignment operator for WickTermCollector and WickTerm.
 *
 * @param lhs The left-hand side WickTermCollector.
 * @param rhs The right-hand side WickTerm.
 * @return WickTermCollector& The updated WickTermCollector.
 */
WickTermCollector& operator+=(WickTermCollector& lhs, const WickTerm& rhs);

/**
 * @brief Subtraction assignment operator for WickTermCollector and WickTerm.
 *
 * @param lhs The left-hand side WickTermCollector.
 * @param rhs The right-hand side WickTerm.
 * @return WickTermCollector& The updated WickTermCollector.
 */
WickTermCollector& operator-=(WickTermCollector& lhs, const WickTerm& rhs);

/**
 * @brief Addition assignment operator for two WickTermCollector objects.
 *
 * @param lhs The left-hand side WickTermCollector.
 * @param rhs The right-hand side WickTermCollector.
 * @return WickTermCollector& The updated WickTermCollector.
 */
WickTermCollector& operator+=(WickTermCollector& lhs, const WickTermCollector& rhs);

/**
 * @brief Subtraction assignment operator for two WickTermCollector objects.
 *
 * @param lhs The left-hand side WickTermCollector.
 * @param rhs The right-hand side WickTermCollector.
 * @return WickTermCollector& The updated WickTermCollector.
 */
WickTermCollector& operator-=(WickTermCollector& lhs, const WickTermCollector& rhs);

/**
 * @brief Addition operator for WickTermCollector and WickTerm.
 *
 * @param lhs The left-hand side WickTermCollector.
 * @param rhs The right-hand side WickTerm.
 * @return WickTermCollector The resulting WickTermCollector.
 */
inline WickTermCollector operator+(WickTermCollector lhs, const WickTerm& rhs) {
    lhs += rhs;
    return lhs;
};

/**
 * @brief Subtraction operator for WickTermCollector and WickTerm.
 *
 * @param lhs The left-hand side WickTermCollector.
 * @param rhs The right-hand side WickTerm.
 * @return WickTermCollector The resulting WickTermCollector.
 */
inline WickTermCollector operator-(WickTermCollector lhs, const WickTerm& rhs) {
    lhs -= rhs;
    return lhs;
};

/**
 * @brief Addition operator for WickTerm and WickTermCollector.
 *
 * @param lhs The left-hand side WickTerm.
 * @param rhs The right-hand side WickTermCollector.
 * @return WickTermCollector The resulting WickTermCollector.
 */
inline WickTermCollector operator+(const WickTerm& lhs, WickTermCollector rhs) {
    rhs += lhs;
    return rhs;
};

/**
 * @brief Subtraction operator for WickTerm and WickTermCollector.
 *
 * @param lhs The left-hand side WickTerm.
 * @param rhs The right-hand side WickTermCollector.
 * @return WickTermCollector The resulting WickTermCollector.
 */
inline WickTermCollector operator-(const WickTerm& lhs, WickTermCollector rhs) {
    rhs -= lhs;
    return rhs;
};

/**
 * @brief Addition operator for two WickTermCollector objects.
 *
 * @param lhs The left-hand side WickTermCollector.
 * @param rhs The right-hand side WickTermCollector.
 * @return WickTermCollector The resulting WickTermCollector.
 */
inline WickTermCollector operator+(WickTermCollector lhs, const WickTermCollector& rhs) {
    lhs += rhs;
    return lhs;
};

/**
 * @brief Subtraction operator for two WickTermCollector objects.
 *
 * @param lhs The left-hand side WickTermCollector.
 * @param rhs The right-hand side WickTermCollector.
 * @return WickTermCollector The resulting WickTermCollector.
 */
inline WickTermCollector operator-(WickTermCollector lhs, const WickTermCollector& rhs) {
    lhs -= rhs;
    return lhs;
};

/**
 * @brief Stream insertion operator for WickTermCollector.
 *
 * @param os The output stream.
 * @param terms The WickTermCollector object.
 * @return std::ostream& The updated output stream.
 */
std::ostream& operator<<(std::ostream& os, const WickTermCollector& terms);

}  // namespace mrock::symbolic_operators

#endif  // MROCK_SYMBOLIC_OPERATORS_INCLUDE_MROCK_SYMBOLIC_OPERATORS_WICKTERMCOLLECTOR_HPP