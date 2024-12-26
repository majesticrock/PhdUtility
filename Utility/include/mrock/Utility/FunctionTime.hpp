#pragma once
#include <chrono>
#include <iostream>
#include <utility>
#include <functional>
#include <type_traits>

namespace mrock::Utility {
	using clock = std::chrono::high_resolution_clock;

	// Prints the runtime of a function to the console and returns the function's return value
	template <class FunctionTypes, class... Args>
	auto function_time_ms(const FunctionTypes& func, Args&&... args)
		-> decltype(func(std::forward<Args>(args)...))
	{
		clock::time_point begin = clock::now();
		if constexpr (!(std::is_same_v<decltype(func(std::forward<Args>(args)...)), void>)) {
			auto return_value = func(std::forward<Args>(args)...);
			clock::time_point end = clock::now();
			std::cout << "Runtime = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
			return return_value;
		}
		else {
			auto return_value = func(std::forward<Args>(args)...);
			clock::time_point end = clock::now();
			std::cout << "Runtime = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
			return return_value;
		}
	}

	// Prints the runtime of a function to the console and returns the function's return value
	template <class FunctionTypes, class... Args>
	auto function_time_micro(const FunctionTypes& func, Args&&... args)
		-> decltype(func(std::forward<Args>(args)...))
	{
		clock::time_point begin = clock::now();
		if constexpr (!(std::is_same_v<decltype(func(std::forward<Args>(args)...)), void>)) {
			auto return_value = func(std::forward<Args>(args)...);
			clock::time_point end = clock::now();
			std::cout << "Runtime = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[microseconds]" << std::endl;
			return return_value;
		}
		else {
			auto return_value = func(std::forward<Args>(args)...);
			clock::time_point end = clock::now();
			std::cout << "Runtime = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[microseconds]" << std::endl;
			return return_value;
		}
	}

	// Prints the runtime of a member function of ObjectType to the console and returns the function's return value
	template <class ObjectType, class ReturnType, class... Args, class... Args2>
	ReturnType member_function_time_ms(ObjectType& obj, ReturnType(ObjectType::* function)(Args...), Args2&&... args)
	{
		clock::time_point begin = clock::now();
		if constexpr (!std::is_same_v<decltype((obj.*function)(std::forward<Args2>(args)...)), void>) {
			auto return_value = (obj.*function)(std::forward<Args2>(args)...);
			clock::time_point end = clock::now();
			std::cout << "Runtime = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
			return return_value;
		}
		else {
			auto return_value = (obj.*function)(std::forward<Args2>(args)...);
			clock::time_point end = clock::now();
			std::cout << "Runtime = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
			return return_value;
		}
	}

	// Prints the runtime of a member function of ObjectType to the console and returns the function's return value
	template <class ObjectType, class ReturnType, class... Args, class... Args2>
	auto member_function_time_micro(ObjectType& obj, ReturnType(ObjectType::* function)(Args...), Args2&&... args)
	{
		clock::time_point begin = clock::now();
		if constexpr (!std::is_same_v<decltype((obj.*function)(std::forward<Args2>(args)...)), void>) {
			auto return_value = (obj.*function)(std::forward<Args2>(args)...);
			clock::time_point end = clock::now();
			std::cout << "Runtime = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[microseconds]" << std::endl;
			return return_value;
		}
		else {
			auto return_value = (obj.*function)(std::forward<Args2>(args)...);
			clock::time_point end = clock::now();
			std::cout << "Runtime = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[microseconds]" << std::endl;
			return return_value;
		}
	}
}