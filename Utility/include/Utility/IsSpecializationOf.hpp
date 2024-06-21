#pragma once
#include <type_traits>

namespace Utility {
	template<typename Test, template<typename...> class Ref>
    struct IsSpecializationOf_t : std::false_type {};

    template<template<typename...> class Ref, typename... Args>
    struct IsSpecializationOf_t<Ref<Args...>, Ref>: std::true_type {};

    template <class Test, template<typename...> class Ref>
	constexpr bool is_specialization_of() {
		return IsSpecializationOf_t<Test, Ref>::value;
	};
};