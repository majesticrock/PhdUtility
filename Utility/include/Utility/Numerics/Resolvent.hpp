#pragma once
// Use (void) to silence unused warnings.
#define assertm(exp, msg) assert(((void)msg, exp))

#include "../OutputConvenience.hpp"
#include "GramSchmidt.hpp"
#include <type_traits>
#include <Eigen/Dense>
#include <cmath>
#include <nlohmann/json.hpp>

namespace Utility::Numerics {
	using std::abs;

	template <class RealType>
	struct ResolventData {
		std::vector<RealType> a_i;
		std::vector<RealType> b_i;
	};

	template<class RealType>
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
		void writeDataToFile(const std::string& filename, const std::vector<std::string>& comments = {}) const
		{
			for (const auto& res_data : lanczos) {
				if (checkDataForNaN(res_data.a_i)) std::cerr << "Resolvent a_i" << std::endl;
				if (checkDataForNaN(res_data.b_i)) std::cerr << "Resolvent b_i" << std::endl;
			}
			saveData(lanczos, filename + ".dat.gz");
		};
	};

	// choose the floating point precision, i.e. float, double or long double
	template <class RealType, bool isComplex>
	class Resolvent
	{
	private:
		using ComputationType = std::conditional_t<isComplex, std::complex<RealType>, RealType>;
		using matrix_t = Eigen::Matrix<ComputationType, Eigen::Dynamic, Eigen::Dynamic>;
		using vector_t = Eigen::Vector<ComputationType, Eigen::Dynamic>;
		using resolvent_data = ResolventData<RealType>;

	public:
		vector_t startingState;
		ResolventDataWrapper<RealType> data;
		// Sets the starting state
		inline void setStartingState(const vector_t& state) {
			this->startingState = state;
		};
		const vector_t& getStartingState() const {
			return this->startingState;
		}
		Resolvent(const vector_t& _StargingState, const std::string& name="") 
			: startingState(_StargingState), data(name) { };
		Resolvent(const std::string name)
			: data(name) { };
		Resolvent() = default;

		// Computes the resolvent's parameters a_i and b_i
		// Symplectic needs to be atleast positive semidefinite!
		void compute(const matrix_t& toSolve, const matrix_t& symplectic, int maxIter)
		{
			size_t matrixSize = toSolve.rows();

			if (toSolve.rows() != toSolve.cols()) {
				std::cerr << "Matrix is not square!" << std::endl;
				throw;
			}

			vector_t currentSolution(matrixSize); // corresponds to |q_(i+1)>
			// First filling
			std::vector<vector_t> basisVectors;
			vector_t first = vector_t::Zero(matrixSize); // corresponds to |q_0>
			vector_t second = this->startingState; // corresponds to |q_1>
			resolvent_data res;
			ComputationType norm_buffer = second.dot(symplectic * second);
			if constexpr (isComplex) {
				assertm(abs(norm_buffer.imag()) < 1e-6, "First norm is complex! ");
			}
			res.b_i.push_back(abs(norm_buffer));

			second /= sqrt(res.b_i.back());
			basisVectors.push_back(first);
			basisVectors.push_back(second);

			std::vector<RealType> deltas, gammas;
			gammas.push_back(1);
			size_t iterNum{};
			bool goOn = true;
			vector_t buffer;
			while (goOn) {
				// algorithm
				buffer = toSolve * basisVectors.back();
				norm_buffer = basisVectors.back().dot(symplectic * buffer);
				if constexpr (isComplex) {
					assertm(abs(norm_buffer.imag()) < 1e-6, "First norm in loop is complex!");
					deltas.push_back(norm_buffer.real());
				}
				else {
					deltas.push_back(norm_buffer);
				}

				currentSolution = (buffer - (deltas.back() * basisVectors.back())) - (gammas.back() * basisVectors[iterNum]);
				norm_buffer = sqrt(currentSolution.dot(symplectic * currentSolution));
				if constexpr (isComplex) {
					assertm(abs(norm_buffer.imag()) < 1e-6, "Second norm in loop is complex!");
				}
				gammas.push_back(abs(norm_buffer));
				basisVectors.push_back(currentSolution / gammas.back());
				++iterNum;

				// breaking conditions
				if (iterNum >= maxIter) {
					goOn = false;
				}
				if (iterNum >= toSolve.rows()) {
					goOn = false;
				}
				if (abs(gammas.back()) < 1e-10) {
					goOn = false;
				}
			}
			for (long i = 0; i < deltas.size(); i++)
			{
				res.a_i.push_back(deltas[i]);
				res.b_i.push_back(gammas[i + 1]);
			}
			data.push_back(std::move(res));
		};

		// Computes the resolvent for a Hermitian problem (i.e. the symplectic matrix is the identity)
		void compute(const matrix_t& toSolve, int maxIter)
		{
			const size_t matrixSize = toSolve.rows();

			if (toSolve.rows() != toSolve.cols()) {
				std::cerr << "Matrix is not square!" << std::endl;
				throw;
			}

			vector_t currentSolution(matrixSize); // corresponds to |q_(i+1)>
			// First filling
			std::vector<vector_t> basisVectors;
			vector_t first = vector_t::Zero(matrixSize); // corresponds to |q_0>
			vector_t second = this->startingState; // corresponds to |q_1>
			resolvent_data res;
			res.b_i.push_back(second.squaredNorm());

			second /= sqrt(res.b_i.back());
			basisVectors.push_back(first);
			basisVectors.push_back(second);

			std::vector<RealType> deltas, gammas;
			gammas.push_back(1);
			size_t iterNum{};
			bool goOn = true;
			vector_t buffer;

			while (goOn) {
				// algorithm
				buffer = toSolve * basisVectors.back();
				if constexpr (isComplex) {
					// This has to be real, as <x|H|x> is always real if H=H^+
					deltas.push_back(basisVectors.back().dot(buffer).real());
				}
				else {
					deltas.push_back(basisVectors.back().dot(buffer));
				}
				currentSolution = (buffer - (deltas.back() * basisVectors.back())) - (gammas.back() * basisVectors[iterNum]);
				gammas.push_back(currentSolution.norm());
				basisVectors.push_back(currentSolution / gammas.back());
				++iterNum;

				// breaking conditions
				if (iterNum >= maxIter) {
					goOn = false;
				}
				if (iterNum >= toSolve.rows()) {
					goOn = false;
				}
				if (abs(gammas.back()) < 1e-10) {
					goOn = false;
				}
			}
			for (size_t i = 0U; i < deltas.size(); ++i)
			{
				res.a_i.push_back(deltas[i]);
				res.b_i.push_back(gammas[i + 1] * gammas[i + 1]);
			}
			// The last b is irrelevant, it does not really exist; it's an artifact of the algorithm
			res.b_i.pop_back();
			data.push_back(std::move(res));
		};

		// Computes the resolvent directly from M and N. This might be more stable for complex matrices
		void computeFromNM(const matrix_t& toSolve, const matrix_t& symplectic, const matrix_t& N, int maxIter)
		{
			auto matrixSize = toSolve.rows();

			if (toSolve.rows() != toSolve.cols()) {
				std::cerr << "Matrix is not square!" << std::endl;
				throw;
			}

			vector_t currentSolution(matrixSize); // corresponds to |q_(i+1)>
			// First filling
			std::vector<vector_t> basisVectors;
			vector_t first = vector_t::Zero(matrixSize); // corresponds to |q_0>
			vector_t second = this->startingState; // corresponds to |q_1>
			resolvent_data res;
			ComputationType norm_buffer = second.dot(symplectic * second);
			if constexpr (isComplex) {
				assertm(abs(norm_buffer.imag()) < 1e-6, "First norm is complex! ");
			}
			res.b_i.push_back(abs(norm_buffer));

			second /= sqrt(res.b_i.back());
			basisVectors.push_back(first);
			basisVectors.push_back(second);

			std::vector<RealType> deltas, gammas;
			gammas.push_back(1);
			size_t  iterNum{};
			bool goOn = true;
			vector_t buffer;
			while (goOn) {
				// algorithm
				buffer = toSolve * basisVectors.back();
				norm_buffer = basisVectors.back().dot(N * basisVectors.back());
				if constexpr (isComplex) {
					assertm(abs(norm_buffer.imag()) < 1e-6, "First norm in loop is complex!");
					deltas.push_back(norm_buffer.real());
				}
				else {
					deltas.push_back(norm_buffer);
				}

				currentSolution = (buffer - (deltas.back() * basisVectors.back())) - (gammas.back() * basisVectors[iterNum]);
				norm_buffer = sqrt(currentSolution.dot(symplectic * currentSolution));
				if constexpr (isComplex) {
					assertm(abs(norm_buffer.imag()) < 1e-6, "Second norm in loop is complex!");
				}
				gammas.push_back(abs(norm_buffer));
				basisVectors.push_back(currentSolution / gammas.back());
				++iterNum;

				// breaking conditions
				if (iterNum >= maxIter) {
					goOn = false;
				}
				if (iterNum >= toSolve.rows()) {
					goOn = false;
				}
				if (abs(gammas.back()) < 1e-10) {
					goOn = false;
				}
			}
			for (size_t i = 0U; i < deltas.size(); ++i)
			{
				res.a_i.push_back(deltas[i]);
				res.b_i.push_back(gammas[i + 1]);
			}
			data.push_back(std::move(res));
		};

		// Computes the resolvent for a Hermitian problem (i.e. the symplectic matrix is the identity)
		// Additionally, this function orthogonalizes the Krylov basis each step
		void computeWithReorthogonalization(const matrix_t& toSolve, int maxIter)
		{
			const size_t matrixSize = toSolve.rows();

			if (toSolve.rows() != toSolve.cols()) {
				std::cerr << "Matrix is not square!" << std::endl;
				throw;
			}

			vector_t currentSolution(matrixSize); // corresponds to |q_(i+1)>
			// First filling
			std::vector<vector_t> basisVectors;
			vector_t second = this->startingState; // corresponds to |q_1>
			resolvent_data res;
			res.b_i.push_back(second.squaredNorm());

			second /= sqrt(res.b_i.back());
			basisVectors.push_back(second);

			std::vector<RealType> deltas, gammas;
			gammas.push_back(1);
			size_t iterNum{};
			bool goOn = true;
			vector_t buffer;
			while (goOn) {
				// algorithm
				buffer = toSolve * basisVectors.back();
				if constexpr (isComplex) {
					// This has to be real, as <x|H|x> is always real if H=H^+
					deltas.push_back(basisVectors.back().dot(buffer).real());
				}
				else {
					deltas.push_back(basisVectors.back().dot(buffer));
				}
				if (iterNum > 0U) {
					currentSolution = (buffer - (deltas.back() * basisVectors.back())) - (gammas.back() * basisVectors[iterNum]);
				}
				else {
					currentSolution = (buffer - (deltas.back() * basisVectors.back()));
				}
				GramSchmidt<ComputationType>::orthogonalizeSingleVector(currentSolution, basisVectors);
				gammas.push_back(currentSolution.norm());
				basisVectors.push_back(currentSolution / gammas.back());
				++iterNum;

				// breaking conditions
				if (iterNum >= maxIter) {
					goOn = false;
				}
				if (iterNum >= toSolve.rows()) {
					goOn = false;
				}
				if (abs(gammas.back()) < 1e-10) {
					goOn = false;
				}
			}
			for (size_t i = 0U; i < deltas.size(); ++i)
			{
				res.a_i.push_back(deltas[i]);
				res.b_i.push_back(gammas[i + 1] * gammas[i + 1]);
			}
			// The last b is irrelevant, it does not really exist; it's an artifact of the algorithm
			res.b_i.pop_back();
			data.push_back(std::move(res));
		};

		const ResolventDataWrapper<RealType>& getData() const {
			return data;
		};

		// Prints the computed data to <filename>
		// Asummes that the data has been computed before...
		void writeDataToFile(const std::string& filename) const
		{
			data.writeDataToFile(filename);
		};
	};

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

	//template<class RealType>
	//void to_json(nlohmann::json& j, const ResolventDataWrapper<RealType>& res_data) {
	//	j = nlohmann::json{
	//		{res_data.name, res_data.lanczos}
	//	};
	//}

	template<class RealType>
	void to_json(nlohmann::json& j, const std::vector<ResolventDataWrapper<RealType>>& vec_resolvent_data) {
		for(const auto& res : vec_resolvent_data) {
			j[res.name] = res.lanczos;
		}
	}
}