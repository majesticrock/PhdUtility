/**
 * @file ElementwiseVector.hpp
 * @brief Provides the ElementwiseVector template class and reduction function objects.
 */

#pragma once
#include "detail/vector_macro.h"
#include <numeric>
#include <algorithm>
#include <iostream>
#include <functional>
#include <cassert>
#include <concepts>
#include <type_traits>

namespace mrock::utility {
    /**
     * @brief Functor for computing the L2 squared norm of a vector.
     */
    struct L2SquaredNorm {
        /**
         * @brief Computes the L2 squared norm of the elements in the vector.
         * 
         * @tparam vector_type The type of the vector.
         * @param elements The vector whose L2 squared norm is to be computed.
         * @return The L2 squared norm of the elements.
         */
        template<class vector_type>
        auto operator()(const vector_type& elements) const {
            using value_type = typename vector_type::value_type;
            return std::transform_reduce(elements.begin(), elements.end(), elements.begin(), value_type{});
        }
    };

    /**
     * @brief Functor for retrieving the first element of a vector.
     */
    struct FirstElement {
        /**
         * @brief Retrieves the first element of the vector.
         * 
         * @tparam vector_type The type of the vector.
         * @param elements The vector from which the first element is to be retrieved.
         * @return The first element of the vector.
         */
        template<class vector_type>
        auto operator()(const vector_type& elements) const {
            return elements.front();
        }
    };

    /**
     * @brief Template class for performing element-wise operations on a vector with a reduction function.
     * 
     * @tparam vector_type The type of the vector.
     * @tparam reduction_type The type of the reduction function. Defaults to L2SquaredNorm.
     */
    template<class vector_type>
    struct ElementwiseVector {
    private:
        /**
         * @brief checks if \c reduction holds a valid function
         * 
         * @return returns \c true if \c reduction is valid.
         */
        bool has_valid_reduction() const {
            return static_cast<bool>(reduction);
        }
    public:
        MROCK_VECTOR_WRAPPER_TYPEDEFS(vector_type);

        typedef std::add_lvalue_reference_t<std::add_const_t<vector_type>> const_reference_vector;
        typedef std::function<value_type(const_reference_vector)> reduction_type;

        reduction_type reduction; ///< The reduction function.
        vector_type elements{}; ///< The vector of elements.
        value_type compare_value{}; ///< The value obtained by applying the reduction function to the elements.
        
        MROCK_VECTOR_WRAPPER_MEMBERS(elements);

        /**
         * @brief Constructs an ElementwiseVector with the given reduction function.
         * 
         * @param reduction The reduction function. Defaults to L2SquaredNorm().
         */
        ElementwiseVector(reduction_type reduction = L2SquaredNorm()) 
            : reduction(reduction), elements(), compare_value(reduction(elements)) {}

        /**
         * @brief Constructs an ElementwiseVector with the given allocator and reduction function.
         * 
         * @param alloc The allocator.
         * @param reduction The reduction function. Defaults to L2SquaredNorm().
         */
        explicit ElementwiseVector(const allocator_type& alloc, reduction_type reduction = L2SquaredNorm()) 
            : reduction(reduction), elements(alloc), compare_value(reduction(elements)) {}

        /**
         * @brief Constructs an ElementwiseVector with the given count, value, allocator, and reduction function.
         * 
         * @param count The number of elements.
         * @param value The value to initialize the elements with.
         * @param alloc The allocator. Defaults to allocator_type().
         * @param reduction The reduction function. Defaults to L2SquaredNorm().
         */
        ElementwiseVector(size_type count, const value_type& value, const allocator_type& alloc = allocator_type(), reduction_type reduction = L2SquaredNorm()) 
            : reduction(reduction), elements(count, value, alloc), compare_value(reduction(elements)) {}

        /**
         * @brief Constructs an ElementwiseVector with the given count and reduction function.
         * 
         * @param count The number of elements.
         * @param reduction The reduction function - has no default because there would be ambiguity otherwise
         */
        explicit ElementwiseVector(size_type count, reduction_type reduction) 
            : reduction(reduction), elements(count), compare_value(reduction(elements)) {}

        /**
         * @brief Constructs an ElementwiseVector with the elements in the range [first, last), allocator, and reduction function.
         * 
         * @tparam InputIt The type of the input iterator.
         * @param first The beginning of the range.
         * @param last The end of the range.
         * @param alloc The allocator. Defaults to allocator_type().
         * @param reduction The reduction function. Defaults to L2SquaredNorm().
         */
        template<class InputIt>
        ElementwiseVector(InputIt first, InputIt last, const allocator_type& alloc = allocator_type(), reduction_type reduction = L2SquaredNorm()) 
            : reduction(reduction), elements(first, last, alloc), compare_value(reduction(elements)) {}

        /**
         * @brief Copy constructor with the given reduction function.
         * 
         * @param other The other ElementwiseVector to copy from.
         * @param reduction The reduction function. Defaults to L2SquaredNorm().
         */
        ElementwiseVector(const ElementwiseVector& other, reduction_type reduction = L2SquaredNorm()) 
            : reduction(reduction), elements(other.elements), compare_value(reduction(elements)) {}

        /**
         * @brief Move constructor with the given reduction function.
         * 
         * @param other The other ElementwiseVector to move from.
         * @param reduction The reduction function. Defaults to L2SquaredNorm().
         */
        ElementwiseVector(ElementwiseVector&& other, reduction_type reduction = L2SquaredNorm()) noexcept 
            : reduction(reduction), elements(std::move(other.elements)), compare_value(reduction(elements)) {}

        /**
         * @brief Constructs an ElementwiseVector with the elements in the initializer list, allocator, and reduction function.
         * 
         * @param init The initializer list.
         * @param alloc The allocator. Defaults to allocator_type().
         * @param reduction The reduction function. Defaults to L2SquaredNorm().
         */
        ElementwiseVector(std::initializer_list<value_type> init, const allocator_type& alloc = allocator_type(), reduction_type reduction = L2SquaredNorm()) 
            : reduction(reduction), elements(init, alloc), compare_value(reduction(elements)) {}

        /**
         * @brief Copy constructor with the given allocator and reduction function.
         * 
         * @param other The other ElementwiseVector to copy from.
         * @param alloc The allocator.
         * @param reduction The reduction function. Defaults to L2SquaredNorm().
         */
        ElementwiseVector(const ElementwiseVector& other, const allocator_type& alloc, reduction_type reduction = L2SquaredNorm()) 
            : reduction(reduction), elements(other.elements, alloc), compare_value(reduction(elements)) {}

        /**
         * @brief Move constructor with the given allocator and reduction function.
         * 
         * @param other The other ElementwiseVector to move from.
         * @param alloc The allocator.
         * @param reduction The reduction function. Defaults to L2SquaredNorm().
         */
        ElementwiseVector(ElementwiseVector&& other, const allocator_type& alloc, reduction_type reduction = L2SquaredNorm()) 
            : reduction(reduction), elements(std::move(other.elements), alloc), compare_value(reduction(elements)) {}

        /**
         * @brief A 0-constructor so that something like vec = 0 will work.
         * This is required for certain boost functionality
         */
        template<class T> requires std::convertible_to<T, value_type>
        ElementwiseVector(T val)
            : reduction(nullptr), elements(), compare_value(val) 
        {
            assert(compare_value == value_type{});
        }

        template<class T> requires std::convertible_to<T, value_type>
        ElementwiseVector& operator=(T val) {
            return ((*this) = ElementwiseVector(val));
        }

        /**
         * @brief Copy assignment operator.
         * 
         * @param other The other ElementwiseVector to copy from.
         * @return Reference to this ElementwiseVector.
         */
        ElementwiseVector& operator=(const ElementwiseVector& other) = default;

        /**
         * @brief Move assignment operator.
         * 
         * @param other The other ElementwiseVector to move from.
         * @return Reference to this ElementwiseVector.
         */
        ElementwiseVector& operator=(ElementwiseVector&& other) noexcept = default;

        /**
         * @brief Adds the elements of another ElementwiseVector to this one.
         * 
         * @param other The other ElementwiseVector.
         * @return Reference to this ElementwiseVector.
         */
        ElementwiseVector& operator+=(const ElementwiseVector& other) {
            assert(other.has_valid_reduction() || this->has_valid_reduction());
            if (!this->has_valid_reduction()) { 
                this->reduction = other.reduction;
                this->resize(other.size());
            }
            if (!other.has_valid_reduction() || other.empty()) {
                return *this;
            }
            std::transform(elements.begin(), elements.end(), other.elements.begin(), elements.begin(), std::plus<value_type>());
            compare_value = reduction(elements);
            return *this;
        }

        /**
         * @brief Subtracts the elements of another ElementwiseVector from this one.
         * 
         * @param other The other ElementwiseVector.
         * @return Reference to this ElementwiseVector.
         */
        ElementwiseVector& operator-=(const ElementwiseVector& other) {
            assert(other.has_valid_reduction() || this->has_valid_reduction());
            if (!this->has_valid_reduction()) { 
                this->reduction = other.reduction;
                this->resize(other.size());
            }
            if (!other.has_valid_reduction() || other.empty()) {
                return *this;
            }
            std::transform(elements.begin(), elements.end(), other.elements.begin(), elements.begin(), std::minus<value_type>());
            compare_value = reduction(elements);
            return *this;
        }

        /**
         * @brief Multiplies the elements of this ElementwiseVector by the elements of another.
         * 
         * @param other The other ElementwiseVector.
         * @return Reference to this ElementwiseVector.
         */
        ElementwiseVector& operator*=(const ElementwiseVector& other) {
            assert(other.has_valid_reduction() || this->has_valid_reduction());
            if (!this->has_valid_reduction()) { 
                this->reduction = other.reduction;
                this->resize(other.size());
            }
            if (!other.has_valid_reduction() || other.empty()) {
                std::ranges::fill(this->begin(), this->end(), value_type{});
                return *this;
            }
            std::transform(elements.begin(), elements.end(), other.elements.begin(), elements.begin(), std::multiplies<value_type>());
            compare_value = reduction(elements);
            return *this;
        }

        /**
         * @brief Divides the elements of this ElementwiseVector by the elements of another.
         * 
         * @param other The other ElementwiseVector.
         * @return Reference to this ElementwiseVector.
         */
        ElementwiseVector& operator/=(const ElementwiseVector& other) {
            assert((other.has_valid_reduction() && !other.empty()) && "Other is invalid! Are you dividing by 0?");
            if (!this->has_valid_reduction()) { 
                this->reduction = other.reduction;
                this->resize(other.size());
            }
            std::transform(elements.begin(), elements.end(), other.elements.begin(), elements.begin(), std::divides<value_type>());
            compare_value = reduction(elements);
            return *this;
        }

        /**
         * @brief Adds a number to the elements of this ElementwiseVector.
         * 
         * @param other The number to add.
         * @return Reference to this ElementwiseVector.
         */
        template<typename T>
        ElementwiseVector& operator+=(const T& other) {
            std::transform(elements.begin(), elements.end(), elements.begin(), [&other](const value_type& elem) { return elem + other; });
            compare_value = reduction(elements);
            return *this;
        }

        /**
         * @brief Subtracts a number from the elements of this ElementwiseVector.
         * 
         * @param other The number to subtract.
         * @return Reference to this ElementwiseVector.
         */
        template<typename T>
        ElementwiseVector& operator-=(const T& other) {
            std::transform(elements.begin(), elements.end(), elements.begin(), [&other](const value_type& elem) { return elem - other; });
            compare_value = reduction(elements);
            return *this;
        }

        /**
         * @brief Multiplies the elements of this ElementwiseVector by a number.
         * 
         * @param other The number to multiply by.
         * @return Reference to this ElementwiseVector.
         */
        template<typename T>
        ElementwiseVector& operator*=(const T& other) {
            std::transform(elements.begin(), elements.end(), elements.begin(), [&other](const value_type& elem) { return elem * other; });
            compare_value = reduction(elements);
            return *this;
        }

        /**
         * @brief Divides the elements of this ElementwiseVector by a number.
         * 
         * @param other The number to divide by.
         * @return Reference to this ElementwiseVector.
         */
        template<typename T>
        ElementwiseVector& operator/=(const T& other) {
            std::transform(elements.begin(), elements.end(), elements.begin(), [&other](const value_type& elem) { return elem / other; });
            compare_value = reduction(elements);
            return *this;
        }

        /**
         * @brief Unary minus operator.
         * 
         * @return A new ElementwiseVector with negated elements.
         */
        ElementwiseVector operator-() const {
            ElementwiseVector result(*this);
            std::transform(result.elements.begin(), result.elements.end(), result.elements.begin(), std::negate<value_type>());
            result.compare_value = result.reduction(result.elements);
            return result;
        }

        /**
         * @brief Adds a number to the elements of an ElementwiseVector.
         * 
         * @param lhs The ElementwiseVector.
         * @param rhs The number.
         * @return A new ElementwiseVector containing the result.
         */
        template<typename T>
        friend ElementwiseVector operator+(ElementwiseVector lhs, const T& rhs) {
            lhs += rhs;
            return lhs;
        }

        /**
         * @brief Subtracts a number from the elements of an ElementwiseVector.
         * 
         * @param lhs The ElementwiseVector.
         * @param rhs The number.
         * @return A new ElementwiseVector containing the result.
         */
        template<typename T>
        friend ElementwiseVector operator-(ElementwiseVector lhs, const T& rhs) {
            lhs -= rhs;
            return lhs;
        }

        /**
         * @brief Multiplies the elements of an ElementwiseVector by a number.
         * 
         * @param lhs The ElementwiseVector.
         * @param rhs The number.
         * @return A new ElementwiseVector containing the result.
         */
        template<typename T>
        friend ElementwiseVector operator*(ElementwiseVector lhs, const T& rhs) {
            lhs *= rhs;
            return lhs;
        }

        /**
         * @brief Divides the elements of an ElementwiseVector by a number.
         * 
         * @param lhs The ElementwiseVector.
         * @param rhs The number.
         * @return A new ElementwiseVector containing the result.
         */
        template<typename T>
        friend ElementwiseVector operator/(ElementwiseVector lhs, const T& rhs) {
            lhs /= rhs;
            return lhs;
        }

        /**
         * @brief Adds a number to the elements of an ElementwiseVector.
         * 
         * @param rhs The ElementwiseVector.
         * @param lhs The number.
         * @return A new ElementwiseVector containing the result.
         */
        template<typename T>
        friend ElementwiseVector operator+(const T& lhs, ElementwiseVector rhs) {
            rhs += lhs;
            return lhs;
        }

        /**
         * @brief Subtracts a number from the elements of an ElementwiseVector.
         * 
         * @param rhs The ElementwiseVector.
         * @param lhs The number.
         * @return A new ElementwiseVector containing the result.
         */
        template<typename T>
        friend ElementwiseVector operator-(const T& lhs, ElementwiseVector rhs) {
            rhs -= lhs;
            return -rhs;
        }

        /**
         * @brief Multiplies the elements of an ElementwiseVector by a number.
         * 
         * @param rhs The ElementwiseVector.
         * @param lhs The number.
         * @return A new ElementwiseVector containing the result.
         */
        template<typename T>
        friend ElementwiseVector operator*(const T& lhs, ElementwiseVector rhs) {
            rhs *= lhs;
            return rhs;
        }

        /**
         * @brief Divides the elements of an ElementwiseVector by a number.
         * 
         * @param rhs The ElementwiseVector.
         * @param lhs The number.
         * @return A new ElementwiseVector containing the result.
         */
        template<typename T>
        friend ElementwiseVector operator/(const T& lhs, ElementwiseVector rhs) {
            for (auto& element : rhs.elements) {
                element = lhs / element;
            }
            rhs.compare_value = rhs.reduction(rhs.elements);
            return rhs;
        }

        /**
         * @brief Adds the elements of two ElementwiseVectors.
         * 
         * @param lhs The left-hand side ElementwiseVector.
         * @param rhs The right-hand side ElementwiseVector.
         * @return A new ElementwiseVector containing the result.
         */
        friend ElementwiseVector operator+(ElementwiseVector lhs, const ElementwiseVector& rhs) {
            lhs += rhs;
            return lhs;
        }

        /**
         * @brief Subtracts the elements of two ElementwiseVectors.
         * 
         * @param lhs The left-hand side ElementwiseVector.
         * @param rhs The right-hand side ElementwiseVector.
         * @return A new ElementwiseVector containing the result.
         */
        friend ElementwiseVector operator-(ElementwiseVector lhs, const ElementwiseVector& rhs) {
            lhs -= rhs;
            return lhs;
        }

        /**
         * @brief Multiplies the elements of two ElementwiseVectors.
         * 
         * @param lhs The left-hand side ElementwiseVector.
         * @param rhs The right-hand side ElementwiseVector.
         * @return A new ElementwiseVector containing the result.
         */
        friend ElementwiseVector operator*(ElementwiseVector lhs, const ElementwiseVector& rhs) {
            lhs *= rhs;
            return lhs;
        }

        /**
         * @brief Divides the elements of two ElementwiseVectors.
         * 
         * @param lhs The left-hand side ElementwiseVector.
         * @param rhs The right-hand side ElementwiseVector.
         * @return A new ElementwiseVector containing the result.
         */
        friend ElementwiseVector operator/(ElementwiseVector lhs, const ElementwiseVector& rhs) {
            lhs /= rhs;
            return lhs;
        }

        /**
         * @brief Three-way comparison operator.
         * 
         * @param other The other ElementwiseVector to compare with.
         * @return The result of the comparison.
         */
        auto operator<=>(const ElementwiseVector& other) const {
            return compare_value <=> other.compare_value;
        }

        /**
         * @brief Output stream operator.
         * If \c MROCK_ELEMENTWISE_VECTOR_PRINT_ELEMENTS is defined, the entire vector content will be printed.
         * Otherwise only \c vec.compare_value will be printed.
         * @param os The output stream.
         * @param vec The ElementwiseVector to output.
         * @return The output stream.
         */
        friend std::ostream& operator<<(std::ostream& os, const ElementwiseVector& vec) {
#ifdef MROCK_ELEMENTWISE_VECTOR_PRINT_ELEMENTS
            os << "[";
            for (size_t i = 0; i < vec.elements.size(); ++i) {
                os << vec.elements[i];
                if (i < vec.elements.size() - 1) {
                    os << ", ";
                }
            }
            os << "] ";
#endif
            os << vec.compare_value;
            return os;
        }

        friend value_type abs(const ElementwiseVector& vec) {
            return std::abs(vec.compare_value);
        }
    };
}