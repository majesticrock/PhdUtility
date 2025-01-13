#pragma once
#include <iomanip>
#include <fstream>
#include <vector>
#include <ctime>
#include <cmath>
#include <iostream>

namespace mrock::utility {
	// Casts a floating point number to a std::string with desired precision n
	template <typename T>
	std::string to_string_with_precision(const T a_value, const int n)
	{
		std::ostringstream out;
		out.precision(n);
		out << std::fixed << a_value;
		return out.str();
	}

	// Returns the local time - the preprocessor statements deal with different OSs
	inline std::tm localtime_xp(std::time_t timer)
	{
		std::tm bt{};
#if defined(__unix__)
		localtime_r(&timer, &bt);
#elif defined(_MSC_VER)
		localtime_s(&bt, &timer);
#else
		static std::mutex mtx;
		std::lock_guard<std::mutex> lock(mtx);
		bt = *std::localtime(&timer);
#endif
		return bt;
	}

	// Returns the current time in a readable format
	// default = "DD-MM-YYYY HH:MM:SS"
	inline std::string time_stamp(const std::string& fmt = "%d-%m-%Y %T")
	{
		auto bt = localtime_xp(std::time(0));
		char buf[64];
		return { buf, std::strftime(buf, sizeof(buf), fmt.c_str(), &bt) };
	}

	// Checks whether the vector data contains any NaNs. Only defines if data_type is float, double or long double
	template <typename vector_type, typename data_type = typename vector_type::value_type>
	inline bool check_data_for_NaN(const vector_type& data) {
		static_assert(std::is_floating_point_v<data_type>);
		for (const auto& value : data) {
			if (std::isnan(value)) {
				std::cerr << "Atleast one of your data points is NaN!" << std::endl;
				return true;
			}
		}
		return false;
	};

	// Writes the current time stamp and comments to the file
	// If the latter is not provided only the time stamp is written
	template <typename outstream_type>
	void write_comments(outstream_type& out, const std::vector<std::string>& comments = std::vector<std::string>())
	{
		out << "# " << time_stamp() << "\n#\n";
		for (const auto& comment : comments) {
			out << "# " << comment << "\n";
		}
	};

	// data_type is a type that overloads the << operator - e.g. float or double, but custom types are allowed too
	// outstream_type is a type of output stream, e.g. std::ofstream or std::ostringstream
	template <typename outstream_type, typename vector_type, typename data_type = typename vector_type::value_type>
	class OutputWriter {
	public:
		// Appends a line consisting of <data> to <out>
		void append_line(const vector_type& data, outstream_type& out) const
		{
			if constexpr (std::is_floating_point_v<data_type>) {
				check_data_for_NaN(data);
			}
			for (size_t i = 0U; i < data.size(); ++i)
			{
				out << data[i];
				if (i < data.size() - 1) {
					out << " ";
				}
			}
			out << "\n";
		}

		void save_data(const vector_type& data, outstream_type& out,
			const std::vector<std::string>& comments = std::vector<std::string>()) const
		{
			if constexpr (std::is_floating_point_v<data_type>) {
				check_data_for_NaN(data);
			}
			write_comments(out, comments);
			out << std::scientific << std::setprecision(10);
			for (const auto& dat : data) {
				out << dat << "\n";
			}
		};

		template<typename Allocator>
		void save_data(const std::vector<vector_type, Allocator>& data, outstream_type& out,
			const std::vector<std::string>& comments = std::vector<std::string>()) const
		{
			write_comments(out, comments);
			out << std::scientific << std::setprecision(10);
			for (const auto& data_line : data) {
				append_line(data_line, out);
			}
		};
	};
}