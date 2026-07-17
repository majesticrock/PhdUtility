#ifndef MROCK_IEOM_INCLUDE_MROCK_IEOM_GENERALRESOLVENT_HPP
#define MROCK_IEOM_INCLUDE_MROCK_IEOM_GENERALRESOLVENT_HPP

#include <chrono>
#include <string>

#ifndef MROCK_IEOM_DO_NOT_PARALLELIZE
#ifndef _OPENMP
#define MROCK_IEOM_DO_NOT_PARALLELIZE
#else
#include <omp.h>
#endif  // ifndef _OPENMP
#endif  // ifndef MROCK_IEOM_DO_NOT_PARALLELIZE

#include "Resolvent.hpp"
#include "detail/PivotToBlockStructure.hpp"
#include "detail/UnderlyingRealType.hpp"
#include "detail/constexpr_power.hpp"
#include "detail/internal_functions.hpp"
#include "detail/is_complex.hpp"

namespace mrock::iEoM {
/**
 * @brief General resolvent manager for block-structured matrices.
 *
 * This CRTP-enabled struct prepares matrices M and N, applies block
 * diagonalization, and computes resolvent data for a set of starting states.
 * It supports both real and complex number types through the provided
 * NumberType template parameter.
 *
 * @tparam NumberType Numeric matrix scalar type, optionally complex.
 */
template <class NumberType>
struct GeneralResolvent {
protected:
    using RealType = detail::UnderlyingRealType_t<NumberType>;
    using Matrix = Eigen::Matrix<NumberType, Eigen::Dynamic, Eigen::Dynamic>;
    using Vector = Eigen::Vector<NumberType, Eigen::Dynamic>;
    using BlockedMatrix = detail::BlockDiagonalMatrix<Matrix>;
    using ResolventType = Resolvent<BlockedMatrix, Eigen::Vector<std::complex<RealType>, Eigen::Dynamic>>;

    Matrix M, N;
    std::vector<std::string> resolvent_names;
    std::vector<Vector> starting_states;

private:
    detail::iEoM_internal<RealType> _internal;

protected:
    // The following virtual functions must be implemented by the user
    // If one knows that certain functions are not needed, implementing an empty
    // function suffices.
    // An example would be: fill_M() is only needed by dynamic_matrix_is_negative()
    // Thus, if this function is not called, fill_M() is not required either.

    /**
     * @brief Must be implemented by the user. Must fill the starting states for the Laczos algorithm
     */
    virtual void create_starting_states() = 0;
    /**
     * @brief Must be implemented by the user. Must fill the matrices M and N
     */
    virtual void fill_matrices() = 0;
    /**
     * @brief Must be implemented by the user. Must fill the matrix M
     */
    virtual void fill_M() = 0;

public:
    /**
     * @brief Construct a GeneralResolvent instance.
     *
     * @param sqrt_precision Precision threshold used for numerical comparisons.
     * @param negative_matrix_is_error If true, negative matrix values are treated as errors.
     */
    GeneralResolvent(RealType const& sqrt_precision, bool negative_matrix_is_error = true)
        : _internal(sqrt_precision, negative_matrix_is_error){};

    /**
     * @brief Virtual default destructor.
     */
    virtual ~GeneralResolvent() = default;

    /**
     * @brief Checks whether the assembled dynamic matrix contains negative eigenvalues,
     * or its diagonal contains negative numbers.
     * These must be >= 0 in thermal equilibrium. Thus, if this function returns true,
     * one can deduce that the system under study is not in thermal equilibrium.
     * Note that the converse is not true.
     *
     * @return true when M contains negative eigenvalues.
     */
    bool dynamic_matrix_is_negative() {
        fill_M();
        if constexpr (detail::is_complex_v<NumberType>) {
            if (this->_internal.contains_negative(M.diagonal().real())) {
                return true;
            }
        } else {
            if (this->_internal.contains_negative(M.diagonal())) {
                return true;
            }
        }
        if (!detail::matrix_wrapper<Matrix>::is_non_negative(M, this->_internal._sqrt_precision)) {
            return true;
        }

        return false;
    };

    /**
     * @brief Compute resolvent data for all starting states.
     *
     * This function fills M and N, pivots them into block structure, performs
     * eigenvalue decompositions, and then computes Lanczos resolvents for each
     * prepared starting state.
     *
     * @tparam CheckHermitian If positive, enforces Hermiticity checks on M and N.
     * @param LANCZOS_ITERATION_NUMBER Number of Lanczos iterations to perform.
     * @return Vector of ResolventDataWrapper objects holding the resolvent results.
     */
    template <int CheckHermitian = -1>
    std::vector<ResolventDataWrapper<RealType>> compute_collective_modes(unsigned int LANCZOS_ITERATION_NUMBER) {
        using __matrix_wrapper__ = detail::blocked_matrix_wrapper<NumberType>;
        std::chrono::time_point begin = std::chrono::steady_clock::now();
        std::chrono::time_point end = std::chrono::steady_clock::now();

        fill_matrices();
        create_starting_states();

        if constexpr (CheckHermitian > 0) {
            if ((M - M.adjoint()).norm() > detail::constexpr_power<-CheckHermitian, RealType, RealType>(10.)) {
                throw std::runtime_error("M is not Hermitian!");
            }
            if ((N - N.adjoint()).norm() > detail::constexpr_power<-CheckHermitian, RealType, RealType>(10.)) {
                throw std::runtime_error("N is not Hermitian!");
            }
        }

        end = std::chrono::steady_clock::now();
        std::cout << "Time for filling of M and N: "
                  << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
        begin = std::chrono::steady_clock::now();

        const auto pivot = detail::pivot_to_block_structure(M);
        M = pivot.transpose() * M * pivot;
        const std::vector<detail::HermitianBlock> blocks = detail::identify_hermitian_blocks(M);
        N = pivot.transpose() * N * pivot;

        BlockedMatrix M_blocked(M, blocks);
        BlockedMatrix N_blocked(N, blocks);

        end = std::chrono::steady_clock::now();
        std::cout << "Time for pivoting: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()
                  << "[ms]" << std::endl;
        begin = std::chrono::steady_clock::now();

        __matrix_wrapper__ M_solver = __matrix_wrapper__::solve_block_diagonal_matrix(M_blocked);  // , blocks
        this->_internal.template apply_matrix_operation<detail::iEoM_operation::NONE>(M_solver.eigenvalues);

        end = std::chrono::steady_clock::now();
        std::cout << "Time for first solving: "
                  << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
        begin = std::chrono::steady_clock::now();

        const auto bufferMatrix = N_blocked * M_solver.eigenvectors;
        // = N * 1/M * N
        BlockedMatrix n_hacek =
            bufferMatrix *
            M_solver.eigenvalues
                .unaryExpr([this](RealType x) { return abs(x) < this->_internal._precision ? 0 : 1. / x; })
                .asDiagonal() *
            bufferMatrix.adjoint();

        __matrix_wrapper__ norm_solver = __matrix_wrapper__::solve_block_diagonal_matrix(n_hacek);  // , blocks
        this->_internal.template apply_matrix_operation<detail::iEoM_operation::SQRT>(norm_solver.eigenvalues);

        end = std::chrono::steady_clock::now();
        std::cout << "Time for norm solving: "
                  << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
        begin = std::chrono::steady_clock::now();

        // n_hacek -> n_hacek^(-1/2)
        n_hacek = norm_solver.eigenvectors *
                  norm_solver.eigenvalues
                      .unaryExpr([this](RealType x) { return abs(x) < this->_internal._precision ? 0 : 1. / x; })
                      .asDiagonal() *
                  norm_solver.eigenvectors.adjoint();
        // Starting here M is the adjusted solver matrix (M s hackem)
        // n_hacek * M * n_hacek
        M_blocked =
            n_hacek * M_blocked *
            n_hacek;  // M_solver.eigenvectors * M_solver.eigenvalues.asDiagonal() * M_solver.eigenvectors.adjoint()
        // Starting here N is the extra matrix that defines |a> ((N s hackem) * N)
        N_blocked.applyOnTheLeft(n_hacek);

        // Starting here h_hacek is its own inverse (defining |b>)
        n_hacek = norm_solver.eigenvectors * norm_solver.eigenvalues.asDiagonal() * norm_solver.eigenvectors.adjoint();

        end = std::chrono::steady_clock::now();
        std::cout << "Time for adjusting of the matrices: "
                  << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
        begin = std::chrono::steady_clock::now();

        const int N_RESOLVENT_TYPES = starting_states.size();
        std::vector<ResolventType> resolvents(3 * N_RESOLVENT_TYPES);
        resolvents.resize(3 * N_RESOLVENT_TYPES);
        if (N_RESOLVENT_TYPES <= resolvent_names.size()) {
            for (std::size_t i = 0U; i < N_RESOLVENT_TYPES; ++i) {
                resolvents[3 * i].data.name = resolvent_names[i] + "_a";
                resolvents[3 * i + 1].data.name = resolvent_names[i] + "_a+b";
                resolvents[3 * i + 2].data.name = resolvent_names[i] + "_a+ib";
            }
        }
#ifndef MROCK_IEOM_DO_NOT_PARALLELIZE
#pragma omp parallel for
#endif
        for (int i = 0; i < N_RESOLVENT_TYPES; i++) {
            Vector a = N * (pivot.transpose() * starting_states[i]).eval();
            Vector b = n_hacek * (pivot.transpose() * starting_states[i]).eval();

            resolvents[3 * i].set_starting_state(a);
            resolvents[3 * i + 1].set_starting_state(0.5 * (a + b));
            resolvents[3 * i + 2].set_starting_state(0.5 * (a + std::complex<RealType>{0, 1} * b));
        }
        // M = M_blocked.construct_matrix();
#ifndef MROCK_IEOM_DO_NOT_PARALLELIZE
#pragma omp parallel for
#endif
        for (int i = 0; i < 3 * N_RESOLVENT_TYPES; i++) {
            resolvents[i].compute(M_blocked, LANCZOS_ITERATION_NUMBER);
        }

        end = std::chrono::steady_clock::now();
        std::cout << "Time for resolventes: "
                  << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;

        std::vector<ResolventDataWrapper<RealType>> ret;
        ret.reserve(resolvents.size());
        for (const auto& re : resolvents) {
            ret.push_back(re.get_data());
        }
        return ret;
    }
};
}  // namespace mrock::iEoM
#endif  // MROCK_IEOM_INCLUDE_MROCK_IEOM_GENERALRESOLVENT_HPP
