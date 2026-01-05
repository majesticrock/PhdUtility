#pragma once
#include <vector>
#include <array>
#include <string>
#include <iostream>
#include <nlohmann/json.hpp>

#include "../RangeUtility.hpp"
#include "../OutputConvenience.hpp"

namespace mrock::utility::Numerics::resolvent_details {
    template <class RealType>
	struct ResolventData {
		std::vector<RealType> a_i;
		std::vector<RealType> b_i;
	};

	template <class RealType>
	struct ResolventDataWrapper {
		std::vector<ResolventData<RealType>> lanczos;
		std::string name;

		ResolventDataWrapper() = default;
		ResolventDataWrapper(const std::string& _name) 
			: name(_name) {};

		void push_back(ResolventData<RealType>&& data_point) {
			lanczos.push_back(std::move(data_point));
		};
		void push_back(const ResolventData<RealType>& data_point) {
			lanczos.push_back(data_point);
		};
		// Prints the computed data to <filename>
		// Asummes that the data has been computed before...
		void write_data_to_file(const std::string& filename, const std::vector<std::string>& comments = {}) const
		{
			for (const auto& res_data : lanczos) {
				if (check_data_for_NaN(res_data.a_i)) std::cerr << "Resolvent a_i" << std::endl;
				if (check_data_for_NaN(res_data.b_i)) std::cerr << "Resolvent b_i" << std::endl;
			}
			save_data(lanczos, filename + ".dat.gz");
		};
	};

	template <class RealType>
	void join_data_wrapper(std::vector<ResolventDataWrapper<RealType>>& target, ResolventDataWrapper<RealType> const& new_data)
	{
		for (auto& sub_target : target) {
			if(sub_target.name == new_data.name) {
				append_vector(sub_target.lanczos, new_data.lanczos);
				return;
			}
		}
		target.push_back(new_data);
	}

	template <class RealType>
	void join_data_wrapper(std::vector<ResolventDataWrapper<RealType>>& target, std::vector<ResolventDataWrapper<RealType>> const& new_data)
	{
		if(target.empty()) {
			target = new_data;
			return;
		}
		for (const auto& _new : new_data) {
			join_data_wrapper(target, _new);
		}
	}

    template <typename T>
	inline std::ostream& operator<<(std::ostream& os, const ResolventData<T>& data)
	{
		for (const auto& elem : data.a_i) {
			os << elem << " ";
		}
		os << "\n";
		for (const auto& elem : data.b_i) {
			os << elem << " ";
		}
		os << "\n";
		return os;
	}

	template<class RealType>
	void to_json(nlohmann::json& j, const ResolventData<RealType>& res_data) {
		j = nlohmann::json{
			{"a_i", res_data.a_i}, {"b_i", res_data.b_i}
		};
	}

	template<class RealType>
	void to_json(nlohmann::json& j, const std::vector<ResolventDataWrapper<RealType>>& vec_resolvent_data) {
		for (const auto& res : vec_resolvent_data) {
			j[res.name] = res.lanczos;
		}
	}



	template<class RealType, int n_residuals>
	struct ResidualInformation {
		std::array<RealType, n_residuals> eigenvalues{};
		std::array<std::vector<RealType>, n_residuals> eigenvectors{};
		std::array<RealType, n_residuals> weights{};
		std::array<RealType, n_residuals> residuals{};
		std::array<bool, n_residuals> converged{};
	};
}