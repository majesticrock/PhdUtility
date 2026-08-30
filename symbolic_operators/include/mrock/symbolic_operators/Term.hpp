#ifndef MROCK_SYMBOLIC_OPERATORS_INCLUDE_MROCK_SYMBOLIC_OPERATORS_TERM_HPP
#define MROCK_SYMBOLIC_OPERATORS_INCLUDE_MROCK_SYMBOLIC_OPERATORS_TERM_HPP
/**
 * @file Term.hpp
 * @brief Defines the Term class and related functions for symbolic operators.
 */

#include "AbstractTerm.hpp"
#include "Coefficient.hpp"
#include "Fractional.hpp"
#include "IndexWrapper.hpp"
#include "Momentum.hpp"
#include "MomentumSymbol.hpp"
#include "Operator.hpp"
#include "SumContainer.hpp"

#include <algorithm>
#include <iosfwd>
#include <string>
#include <utility>
#include <vector>

namespace mrock::symbolic_operators {
/**
 * @class Term
 * @brief Represents a term in symbolic operator expressions.
 *
 * This class represents a Term. It has various kind of constructors that allow setting coefficient(s), sums, operators
 * and deltas. Using IntFractional, the term can have rational prefactors, e.g., 1/2.
 *
 * A Hamiltonian (or any other summation of operators) is characterized as \c TermCollector.
 * It can consist of any number of individual terms.
 * For a few practical examples, see the files in the tests folder.
 * See bosons.cpp, continuum.cpp, and compare_test.hpp.
 * My own projects using this library are, e.g.,
 * https://github.com/majesticrock/FermionCommute and https://github.com/majesticrock/FlowCommutators.
 *
 * After creating atleast two \c Terms (or \c TermCollector), you may commute them by calling
 * \code
 * TermCollector result = commutator(A, B);
 * result.clean_up();
 * \endcode
 * After calling the commutator, you should pretty much always call
 * TermCollector::clean_up() because commutator performs the normal ordering procedure,
 * however, does not attempt to beautify the result. clean_up then sorts the terms, adds identical ones togeFther and
 * removes those that are equal to 0.
 *
 * Similarly, a double commutator \f$ [C, [A, B]] \f$ can be evaluated by
 * \code
 * TermCollector inner_result = commutator(A, B);
 * inner_result.clean_up();
 * TermCollector result = commutator(C, inner_result);
 * result.clean_up();
 * \endcode
 *
 * To output the results, an overload of \c operator<< is provided for both \c Term and \c TermCollector.
 * The out put is formatted so that it can be used within an align-environment within LaTeX.
 *
 * @sa Coefficient, SumContainer, Operator, KroneckerDelta
 */
class Term : public AbstractTerm<Operator> {
public:
    /**
     * @brief Serializes the term.
     * @tparam Archive The archive type.
     * @param ar The archive.
     * @param version The version.
     */
    template <class Archive>
    void serialize(Archive& ar, [[maybe_unused]] const unsigned int version) {
        ar & coefficients;
        ar & sums;
        ar & operators;
        ar & delta_momenta;
        ar & delta_indices;
        ar & multiplicity;
    }

    /**
     * @brief Constructs a Term with a summation over momenta and spins and multiple coefficients
     *
     * @param _multiplicity The _multiplicity of the term
     * @param _coefficients The coefficients
     * @param _sums The sums
     * @param _operators The operators of the term
     */
    Term(IntFractional _multiplicity,
         std::vector<Coefficient> _coefficients,
         const SumContainer& _sums,
         const std::vector<Operator>& _operators = std::vector<Operator>());

    /**
     * @brief Constructs a Term with a summation over momenta and spins and a coefficient
     *
     * @param _multiplicity The _multiplicity of the term
     * @param _coefficient The coefficient
     * @param _sums The sums
     * @param _operators The operators of the term
     */
    Term(IntFractional _multiplicity,
         Coefficient _coefficient,
         const SumContainer& _sums,
         const std::vector<Operator>& _operators = std::vector<Operator>());

    /**
     * @brief Constructs a Term with a summation over momenta and a coefficient
     *
     * @param _multiplicity The _multiplicity of the term
     * @param _coefficient The coefficient
     * @param _sum_momenta Sum over momenta
     * @param _operators The operators of the term
     */
    Term(IntFractional _multiplicity,
         Coefficient _coefficient,
         const MomentumSum& _sum_momenta,
         const std::vector<Operator>& _operators = std::vector<Operator>());

    /**
     * @brief Constructs a Term with a summation over spins (or other indices) and a coefficient
     *
     * @param _multiplicity The _multiplicity of the term
     * @param _coefficient The coefficient
     * @param _sum_spins Sum over spins (or other indices)
     * @param _operators The operators of the term
     */
    Term(IntFractional _multiplicity,
         Coefficient _coefficient,
         const IndexSum& _sum_spins,
         const std::vector<Operator>& _operators = std::vector<Operator>());

    /**
     * @brief Constructs a Term with a coefficient
     *
     * @param _multiplicity The _multiplicity of the term
     * @param _coefficient The coefficient
     * @param _operators The operators of the term
     */
    Term(IntFractional _multiplicity,
         Coefficient _coefficient,
         const std::vector<Operator>& _operators = std::vector<Operator>());

    /**
     * @brief Constructs a Term with a summation over momenta and indices
     *
     * @param _multiplicity The _multiplicity of the term
     * @param _sums Sums
     * @param _operators The operators of the term
     */
    Term(IntFractional _multiplicity,
         const SumContainer& _sums,
         const std::vector<Operator>& _operators = std::vector<Operator>());

    /**
     * @brief Constructs a Term with a summation over momenta
     *
     * @param _multiplicity The _multiplicity of the term
     * @param _sum_momenta Sum over momenta
     * @param _operators The operators of the term
     */
    Term(IntFractional _multiplicity,
         const MomentumSum& _sum_momenta,
         const std::vector<Operator>& _operators = std::vector<Operator>());

    /**
     * @brief Constructs a Term with a summation over spins (or other indices)
     *
     * @param _multiplicity The _multiplicity of the term
     * @param _sum_spins Sum over spins (or other indices)
     * @param _operators The operators of the term
     */
    Term(IntFractional _multiplicity,
         const IndexSum& _sum_spins,
         const std::vector<Operator>& _operators = std::vector<Operator>());

    /**
     * @brief Constructs a Term with only a multiplicity
     *
     * @param _multiplicity The _multiplicity of the term
     * @param _operators The operators of the term
     */
    explicit Term(IntFractional _multiplicity, const std::vector<Operator>& _operators = std::vector<Operator>());

    /**
     * @brief Default constructor.
     */
    Term() = default;

    /**
     * @brief Checks if the term contains a boson.
     * @return True if the term contains a boson, false otherwise.
     */
    inline bool contains_boson() const;

    /**
     * @brief Checks if the term contains a fermion.
     * @return True if the term contains a fermion, false otherwise.
     */
    inline bool contains_fermion() const;

    /**
     * @brief Counts the number of bosons in the term.
     * @return The number of bosons.
     */
    inline int count_bosons() const;

    /**
     * @brief Counts the number of fermions in the term.
     * @return The number of fermions.
     */
    inline int count_fermions() const;

    /**
     * @brief Prints the term.
     */
    void print() const;

    /**
     * @brief Swaps two operators in the term. Does NOT consider possible additional terms spawned by this operation due
     * to non-commutivity!
     * @param i The index of the first operator.
     * @param j The index of the second operator.
     */
    inline void perform_operator_swap(const std::size_t i, const std::size_t j);

    /**
     * @brief Resolves the Kronecker deltas in the term ( calls \c resolve_momentum_deltas() and \c
     * resolve_index_deltas() )
     * @return True if successful, false otherwise.
     */
    bool resolve_deltas();

    /**
     * @brief Puts the operators in the term in a specific order based on their indices.
     */
    void sort_operators_by_indices();

    /**
     * @brief Tries to bring the momentum dependencies of the operators in \c *this to a fixed notation
     */
    void structure_momentum_dependencies();

    /**
     * @brief Tries to bring \c *this to a fixed notation
     */
    void structure();

    /**
     * @brief Checks if the term is equal to another term (excluding multiplicity).
     * @param other The other term.
     * @return True if equal, false otherwise.
     */
    bool is_equal(const Term& other) const;

    /**
     * @brief Checks if the term is in normal order.
     * @return True if in normal order, false otherwise.
     */
    bool is_normal_ordered() const;

    /**
     * @brief Converts the term to a string without the prefactor.
     * @return The string representation.
     */
    std::string to_string_without_prefactor() const;

    /**
     * @brief Applies the Hermitian conjugate to the term.
     * @return A reference to *this
     */
    Term& hermitian_conjugate_inplace();

    /**
     * @brief Creates hermitian conjugate of this as a new object.
     * @return Returns the new object.
     */
    Term hermitian_conjugate() const;

    /**
     * @brief Renames indices in the term.
     * @param what The index to rename.
     * @param to The new index.
     */
    void rename_indices(const Index what, const Index to);

    /**
     * @brief Compute other * this
     * IMPOARTANT: The result will not be normal ordered! If you require
     * a normal ordered expression, please call normal_order!
     * Note, that doing so may create additional terms, and the result must therefore
     * be \c TermCollector and can therefore not return a reference to \c *this
     *
     * @param other the other Term object.
     * @return Reference to \c *this containing the result.
     */
    Term& multiply_from_the_left(const Term& other);

    /**
     * @brief Compute this * other
     * IMPOARTANT: The result will not be normal ordered! If you require
     * a normal ordered expression, please call normal_order!
     * Note, that doing so may create additional terms, and the result must therefore
     * be \c TermCollector and can therefore not return a reference to \c *this
     *
     * @param other the other Term object.
     * @return Reference to \c *this containing the result.
     */
    Term& multiply_from_the_right(const Term& other);

    /**
     * @brief Multiplies this by rhs
     * IMPOARTANT: The result will not be normal ordered! If you require
     * a normal ordered expression, please call normal_order!
     * Note, that doing so may create additional terms, and the result must therefore
     * be \c TermCollector and can therefore not return a reference to \c *this
     *
     * @param other the other Term object.
     * @return Reference to \c *this containing the result.
     */
    Term& operator*=(const Term& rhs);
};

/**
 * @brief Checks if two terms are equal.
 * @param lhs The left-hand side term.
 * @param rhs The right-hand side term.
 * @return True if equal, false otherwise.
 */
inline bool operator==(const Term& lhs, const Term& rhs) {
    return lhs.is_equal(rhs);
}

/**
 * @brief Checks if two terms are not equal.
 * @param lhs The left-hand side term.
 * @param rhs The right-hand side term.
 * @return True if not equal, false otherwise.
 */
inline bool operator!=(const Term& lhs, const Term& rhs) {
    return !(lhs == rhs);
}

/**
 * @brief Multiplies a Term by another Term.
 * IMPOARTANT: The result will not be normal ordered! If you require
 * a normal ordered expression, please call normal_order!
 * See Term::operator*=(const Term& rhs)
 *
 * @param lhs The left-hand side Term.
 * @param rhs The right-hand side Term.
 * @return The result of the multiplication.
 */
inline Term operator*(Term lhs, const Term& rhs) {
    lhs *= rhs;
    return lhs;
}

// Inline definitions
bool Term::contains_boson() const {
    return std::any_of(operators.begin(), operators.end(), [](Operator const& op) { return (!op.is_fermion); });
}
bool Term::contains_fermion() const {
    return std::any_of(operators.begin(), operators.end(), [](Operator const& op) { return op.is_fermion; });
}
int Term::count_bosons() const {
    return std::count_if(operators.begin(), operators.end(), [](Operator const& op) { return (!op.is_fermion); });
}
int Term::count_fermions() const {
    return std::count_if(operators.begin(), operators.end(), [](Operator const& op) { return op.is_fermion; });
}
void Term::perform_operator_swap(const std::size_t i, const std::size_t j) {
    std::size_t min = std::min(i, j);
    const std::size_t max = std::max(i, j);
    bool parity = operators[min].is_fermion && operators[max].is_fermion;
    while (++min < max) {
        if (operators[min].is_fermion && operators[max].is_fermion) {
            parity = !parity;
        }
    }
    if (parity) {
        this->multiplicity *= -1;
    }
    std::swap(operators[i], operators[j]);
}
}  // namespace mrock::symbolic_operators
#endif  // MROCK_SYMBOLIC_OPERATORS_INCLUDE_MROCK_SYMBOLIC_OPERATORS_TERM_HPP
