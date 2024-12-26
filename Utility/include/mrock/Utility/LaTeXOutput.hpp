#pragma once
#include <string>

namespace mrock::Utility {
	template <class T, class OutputStream>
	void surround_with_env(const T& toPrint, OutputStream& os, const std::string& env) {
		os << "\\begin{" << env << "}\n" << toPrint << " " << "\\end{" << env << "}";
	}

	template <class T>
	struct LaTeXWrapper {
		const T& toPrint;
		const std::string env;
	};
	template <class T>
	LaTeXWrapper<T> as_LaTeX(const T& toPrint, const std::string& env) {
		return LaTeXWrapper<T>{toPrint, env};
	}

	template <class T>
	inline std::ostream& operator<<(std::ostream& os, const LaTeXWrapper<T>& toPrint) {
		surround_with_env(toPrint.toPrint, os, toPrint.env);
		return os;
	}
}