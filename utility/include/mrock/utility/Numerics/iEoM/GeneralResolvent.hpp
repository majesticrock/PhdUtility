#pragma once
#include "../Resolvent.hpp"
#include "_internal_functions.hpp"
#include "../../IsComplex.hpp"
#include "../../ConstexprPower.hpp"
#include "../PivotToBlockStructure.hpp"
#include <chrono>
#include <string>

namespace mrock::utility::Numerics::iEoM {
	template<class Derived, class NumberType>
	struct GeneralResolvent {
	public:
		using RealType = UnderlyingFloatingPoint_t<NumberType>;
		using Matrix = Eigen::Matrix<NumberType, Eigen::Dynamic, Eigen::Dynamic>;
		using Vector = Eigen::Vector<NumberType, Eigen::Dynamic>;
		using BlockedMatrix = BlockDiagonalMatrix<Matrix>;
		using ResolventType = Resolvent<BlockedMatrix, Eigen::Vector< std::complex<RealType>, Eigen::Dynamic >>;

	protected:
		Matrix M, N;
		std::vector<std::string> resolvent_names;
		std::vector<Vector> starting_states;

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
			M = this->_internal.remove_noise(M);
			if (not matrix_wrapper<Matrix>::is_non_negative(M, this->_internal._sqrt_precision)) {
				return true;
			}

			return false;
		};

		template<int CheckHermitian = -1>
		std::vector<ResolventDataWrapper<RealType>> compute_collective_modes(unsigned int LANCZOS_ITERATION_NUMBER)
		{
			using __matrix_wrapper__ = blocked_matrix_wrapper<NumberType>;
			std::chrono::time_point begin = std::chrono::steady_clock::now();
			std::chrono::time_point end = std::chrono::steady_clock::now();

			_derived->fillMatrices();
			_derived->createStartingStates();

			M = this->_internal.remove_noise(M);
			N = this->_internal.remove_noise(N);

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

			const auto pivot = pivot_to_block_structure(M);
			M = pivot.transpose() * M * pivot;
			const std::vector<HermitianBlock> blocks = identify_hermitian_blocks(M);
			N = pivot.transpose() * N * pivot;

			BlockedMatrix M_blocked(M, blocks);
			BlockedMatrix N_blocked(N, blocks);
			
			end = std::chrono::steady_clock::now();
			std::cout << "Time for pivoting: "
				<< std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
			begin = std::chrono::steady_clock::now();

			__matrix_wrapper__ M_solver = __matrix_wrapper__::solve_block_diagonal_matrix(M_blocked); // , blocks
			this->_internal.template apply_matrix_operation<IEOM_NONE>(M_solver.eigenvalues);

			end = std::chrono::steady_clock::now();
			std::cout << "Time for first solving: "
				<< std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
			begin = std::chrono::steady_clock::now();

			const auto bufferMatrix = N_blocked * M_solver.eigenvectors;
			// = N * 1/M * N
			BlockedMatrix n_hacek = bufferMatrix
				* M_solver.eigenvalues.unaryExpr([this](RealType x) { return abs(x) < this->_internal._precision ? 0 : 1. / x; }).asDiagonal()
				* bufferMatrix.adjoint();

			__matrix_wrapper__ norm_solver = __matrix_wrapper__::solve_block_diagonal_matrix(n_hacek); // , blocks
			this->_internal.template apply_matrix_operation<IEOM_SQRT>(norm_solver.eigenvalues);

			end = std::chrono::steady_clock::now();
			std::cout << "Time for norm solving: "
				<< std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
			begin = std::chrono::steady_clock::now();

			// n_hacek -> n_hacek^(-1/2)
			n_hacek = norm_solver.eigenvectors
				* norm_solver.eigenvalues.unaryExpr([this](RealType x) { return abs(x) < this->_internal._precision ? 0 : 1. / x; }).asDiagonal()
				* norm_solver.eigenvectors.adjoint();
			// Starting here M is the adjusted solver matrix (s s hackem)
			// n_hacek * M * n_hacek
			M_blocked = n_hacek * M_blocked * n_hacek; // M_solver.eigenvectors * M_solver.eigenvalues.asDiagonal() * M_solver.eigenvectors.adjoint()
			// Starting here N is the extra matrix that defines |a> (n s hackem N)
			N_blocked.applyOnTheLeft(n_hacek);

			// Starting here h_hacek is its own inverse (defining |b>)
			n_hacek = norm_solver.eigenvectors * norm_solver.eigenvalues.asDiagonal() * norm_solver.eigenvectors.adjoint();

			end = std::chrono::steady_clock::now();
			std::cout << "Time for adjusting of the matrices: "
				<< std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
			begin = std::chrono::steady_clock::now();

			const int N_RESOLVENT_TYPES = starting_states.size();
			std::vector< ResolventType > resolvents(3 * N_RESOLVENT_TYPES);
			resolvents.resize(3 * N_RESOLVENT_TYPES);
			if(N_RESOLVENT_TYPES <= resolvent_names.size()) {
				for (size_t i = 0U; i < N_RESOLVENT_TYPES; ++i) {
					resolvents[3 * i].data.name = resolvent_names[i] + "_a";
					resolvents[3 * i + 1].data.name = resolvent_names[i] + "_a+b";
					resolvents[3 * i + 2].data.name = resolvent_names[i] + "_a+ib";
				}
			}

#pragma omp parallel for
			for (int i = 0; i < N_RESOLVENT_TYPES; i++)
			{
				Vector a = N * (pivot.transpose() * starting_states[i]).eval();
				Vector b = n_hacek * (pivot.transpose() * starting_states[i]).eval();

				resolvents[3 * i].set_starting_state(a);
				resolvents[3 * i + 1].set_starting_state(0.5 * (a + b));
				resolvents[3 * i + 2].set_starting_state(0.5 * (a + std::complex<RealType>{ 0, 1 } *b));
			}
			//M = M_blocked.construct_matrix();
#pragma omp parallel for
			for (int i = 0; i < 3 * N_RESOLVENT_TYPES; i++)
			{
				resolvents[i].compute(M_blocked, LANCZOS_ITERATION_NUMBER);
			}

			end = std::chrono::steady_clock::now();
			std::cout << "Time for resolventes: "
				<< std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;

			std::vector<ResolventDataWrapper<RealType>> ret;
			ret.reserve(resolvents.size());
			for (const auto& re : resolvents)
			{
				ret.push_back(re.get_data());
			}
			return ret;
		}

	private:
		ieom_internal<RealType> _internal;
		Derived* _derived;
	};
}