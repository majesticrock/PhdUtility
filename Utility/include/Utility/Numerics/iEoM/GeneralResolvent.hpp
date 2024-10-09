#pragma once
#include "../Resolvent.hpp"
#include "_internal_functions.hpp"
#include "../../IsComplex.hpp"
#include "../../ConstexprPower.hpp"
#include "../PivotToBlockStructure.hpp"
#include <chrono>
#include <string>

namespace Utility::Numerics::iEoM {
	template<class Derived, class NumberType>
	struct GeneralResolvent {
	public:
		using RealType = UnderlyingFloatingPoint_t<NumberType>;
		using Matrix = Eigen::Matrix<NumberType, Eigen::Dynamic, Eigen::Dynamic>;
		using Vector = Eigen::Vector<NumberType, Eigen::Dynamic>;
		using RealVector = Eigen::Vector<RealType, Eigen::Dynamic>;
		using ComplexVector = Eigen::Vector<std::complex<RealType>, Eigen::Dynamic>;

	protected:
		Matrix M, N;
		std::vector<std::string> resolvent_names;
		std::vector<ComplexVector> starting_states;

	public:
		GeneralResolvent(Derived* derived_ptr, RealType const& sqrt_precision, bool negative_matrix_is_error = true)
			: _internal(sqrt_precision, negative_matrix_is_error), _derived(derived_ptr) { };

		virtual ~GeneralResolvent() = default;

		bool dynamic_matrix_is_negative() {
			_derived->fill_M();
			if constexpr (is_complex<NumberType>()) {
				if (this->_internal.contains_negative(M.diagonal().real())) {
					return true;
				}
			}
			else {
				if (this->_internal.contains_negative(M.diagonal())) {
					return true;
				}
			}
			M = this->_internal.removeNoise(M);
			if (not matrix_wrapper<NumberType>::is_non_negative(M, this->_internal._sqrt_precision)) {
				return true;
			}

			return false;
		};

		template<int CheckHermitian = -1>
		std::vector<ResolventDataWrapper<RealType>> computeCollectiveModes(unsigned int LANCZOS_ITERATION_NUMBER)
		{
			std::chrono::time_point begin = std::chrono::steady_clock::now();
			std::chrono::time_point end = std::chrono::steady_clock::now();

			_derived->fillMatrices();
			_derived->createStartingStates();

			M = this->_internal.removeNoise(M);
			N = this->_internal.removeNoise(N);

			if constexpr (CheckHermitian > 0) {
				if ((M - M.adjoint()).norm() > constexprPower<-CheckHermitian, RealType, RealType>(10.)) {
					throw std::runtime_error("M is not Hermitian!");
				}
				if ((N - N.adjoint()).norm() > constexprPower<-CheckHermitian, RealType, RealType>(10.)) {
					throw std::runtime_error("N is not Hermitian!");
				}
			}

			end = std::chrono::steady_clock::now();
			std::cout << "Time for filling of M and N: "
				<< std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
			begin = std::chrono::steady_clock::now();

			Eigen::SelfAdjointEigenSolver<Matrix> M_solver(M);
			RealVector& evs_M = const_cast<RealVector&>(M_solver.eigenvalues());
			this->_internal.template applyMatrixOperation<IEOM_NONE>(evs_M);

			auto bufferMatrix = N * M_solver.eigenvectors();
			// = N * 1/M * N
			Matrix n_hacek = bufferMatrix
				* evs_M.unaryExpr([this](RealType x) { return abs(x) < this->_internal._precision ? 0 : 1. / x; }).asDiagonal()
				* bufferMatrix.adjoint();

			Eigen::SelfAdjointEigenSolver<Matrix> norm_solver(n_hacek);
			RealVector& evs_norm = const_cast<RealVector&>(norm_solver.eigenvalues());
			this->_internal.template applyMatrixOperation<IEOM_SQRT>(evs_norm);

			// n_hacek -> n_hacek^(-1/2)
			n_hacek = norm_solver.eigenvectors()
				* evs_norm.unaryExpr([this](RealType x) { return abs(x) < this->_internal._precision ? 0 : 1. / x; }).asDiagonal()
				* norm_solver.eigenvectors().adjoint();
			// Starting here M is the adjusted solver matrix (s s hackem)
			// n_hacek * M * n_hacek
			M = n_hacek * M_solver.eigenvectors() * evs_M.asDiagonal() * M_solver.eigenvectors().adjoint() * n_hacek;
			// Starting here N is the extra matrix that defines |a> (n s hackem N)
			N.applyOnTheLeft(n_hacek);

			// Starting here h_hacek is its own inverse (defining |b>)
			n_hacek = norm_solver.eigenvectors() * evs_norm.asDiagonal() * norm_solver.eigenvectors().adjoint();

			end = std::chrono::steady_clock::now();
			std::cout << "Time for adjusting of the matrices: "
				<< std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
			begin = std::chrono::steady_clock::now();

			const int N_RESOLVENT_TYPES = starting_states.size();
			std::vector<Resolvent<RealType, true>> resolvents(3 * N_RESOLVENT_TYPES);
			resolvents.resize(3 * N_RESOLVENT_TYPES);
			if(N_RESOLVENT_TYPES <= resolvent_names.size()) {
				for(size_t i = 0U; i < N_RESOLVENT_TYPES; ++i) {
					resolvents[3 * i].data.name = resolvent_names[i] + "_a";
					resolvents[3 * i + 1].data.name = resolvent_names[i] + "_a+b";
					resolvents[3 * i + 2].data.name = resolvent_names[i] + "_a+ib";
				}
			}

#pragma omp parallel for
			for (int i = 0; i < N_RESOLVENT_TYPES; i++)
			{
				ComplexVector a = N * starting_states[i];
				ComplexVector b = n_hacek * starting_states[i];

				resolvents[3 * i].setStartingState(a);
				resolvents[3 * i + 1].setStartingState(0.5 * (a + b));
				resolvents[3 * i + 2].setStartingState(0.5 * (a + std::complex<RealType>{ 0, 1 } *b));
			}
#pragma omp parallel for
			for (int i = 0; i < 3 * N_RESOLVENT_TYPES; i++)
			{
				resolvents[i].compute(M, LANCZOS_ITERATION_NUMBER);
			}

			end = std::chrono::steady_clock::now();
			std::cout << "Time for resolventes: "
				<< std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;

			std::vector<ResolventDataWrapper<RealType>> ret;
			ret.reserve(resolvents.size());
			for (const auto& re : resolvents)
			{
				ret.push_back(re.getData());
			}
			return ret;
		}

	private:
		ieom_internal<RealType> _internal;
		Derived* _derived;
	};
}