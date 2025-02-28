/**
 * @file ThrowException.hpp
 * @brief Provides a utility function for throwing exceptions conditionally based on the build configuration.
 */

#pragma once
#include <utility>

#ifndef NDEBUG
#define MROCK_NDEBUG_CONDITION(condition) condition
#else
#define MROCK_NDEBUG_CONDITION(condition) false
#endif

namespace mrock::utility {
    /**
    * @brief Conditionally throws an exception based on the build configuration.
    *
    * In debug builds (when \c NDEBUG is not defined), this function throws an exception of the specified type,
    * constructed with the provided arguments. In release builds (when \c NDEBUG is defined), this function does nothing.
    *
    * @tparam exception The type of the exception to be thrown.
    * @tparam Args The types of the arguments to be passed to the exception's constructor.
    * @param condition If true the exception is thrown. May be used with the macro \c MROCK_NDEBUG_CONDITION(condition) which is only evaluated if \c NDEBUG is not defined.
    * @param args The arguments to be passed to the exception's constructor.
    */
    template<class exception, class... Args>
    void throw_exception(bool condition, Args&&... args) {
#ifndef NDEBUG
        if (condition) {
            throw exception(std::forward<Args>(args)...);
        }
#endif
    }
} // namespace mrock::utility