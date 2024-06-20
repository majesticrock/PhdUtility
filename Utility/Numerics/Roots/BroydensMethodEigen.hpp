#pragma once
#include <Eigen/Dense>
#include <complex>
#include <type_traits>
#include <cstring>

namespace Utility::Numerics::Roots {
	template<typename RealType, int t_vector_size>
	class BroydensMethodEigen {
		static_assert(std::is_floating_point<RealType>::value, "You're data type must be a floating point number");

		using MatrixType = Eigen::Matrix<RealType, t_vector_size, t_vector_size>;
		using VectorType = Eigen::Vector<RealType, t_vector_size>;

		using error_type = std::conditional_t< sizeof(RealType) >= sizeof(double), double, RealType >;

		static void estimate_jacobian(MatrixType& J_new, const MatrixType& J_old,
			const VectorType& delta_x, const VectorType& delta_F)
		{
			J_new = (delta_x - J_old * delta_F) * (delta_x.transpose() * J_old);
			J_new /= (delta_x.dot(J_old * delta_F));
			J_new += J_old;
		};
	public:
		// the function must have the following signature void func(const VectorType& input, VectorType& output)
		template<class FunctionType>
		static bool compute(const FunctionType& func, VectorType& x0, const int MAX_ITER = 200)
		{
			size_t DIM = x0.rows();
			// You may play around with EPS_X and EPS_F to your desire
			// EPS_X is the minimum distance between x_i and x_i+1
			// EPS_F is the minimum f(x)

			constexpr error_type EPS_F = std::numeric_limits<error_type>::epsilon();
			constexpr error_type EPS_X = std::numeric_limits<error_type>::epsilon();
			RealType diff_x{ 100 };
			RealType diff_F{ 100 };
			int iter_num{};

			VectorType F_old{ VectorType::Zero(DIM) },
				F_new{ VectorType::Zero(DIM) },
				delta_x{ VectorType::Zero(DIM) },
				delta_F{ VectorType::Zero(DIM) };
			MatrixType J_old{ MatrixType::Zero(DIM, DIM) },
				J_new{ MatrixType::Identity(DIM, DIM) };
			func(x0, F_new);

			while (diff_x > EPS_X && diff_F > EPS_F && iter_num++ <= MAX_ITER && F_new.norm() > EPS_F) {
				delta_x = -J_new * F_new;
				x0 += delta_x;
				diff_x = delta_x.norm();
				F_old = F_new;
				func(x0, F_new);
				delta_F = F_new - F_old;
				diff_F = delta_F.norm();

				J_old = J_new;
				estimate_jacobian(J_new, J_old, delta_x, delta_F);
			}
			// This method returns true if convergence is achieved, in this case if |F(x_final)| < 1e-10
			return (F_new.norm() < 1e-10);
		}
	};

	template<typename RealType, int t_vector_size>
	class BroydensMethodEigen< std::complex<RealType>, t_vector_size > {
		static constexpr Eigen::Index double_size = t_vector_size == Eigen::Dynamic ? Eigen::Dynamic : 2 * t_vector_size;
		using VectorType = Eigen::Vector<std::complex<RealType>, t_vector_size>;
		using RealVector = Eigen::Vector<RealType, double_size>;

		using RealSolver = BroydensMethodEigen<RealType, t_vector_size>;
	public:
	template<class FunctionType>
		// the function must have the following signature void func(const VectorType& input, VectorType& output)
		static bool compute(const FunctionType& func, VectorType& x_complex, const int MAX_ITER = 200)
		{
			VectorType f_complex = x_complex;
			RealVector x0;
			if constexpr (t_vector_size == Eigen::Dynamic) {
				x0.setZero(2 * x_complex.size());
			}
			else {
				x0.setZero();
			}
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wclass-memaccess"
			auto call_f_from_real = [&](const RealVector& x_real, RealVector& f_real) {
				std::memcpy(x_complex.data(), x_real.data(), 2 * sizeof(RealType) * x_complex.size());
				func(x_complex, f_complex);
				std::memcpy(f_real.data(), f_complex.data(), 2 * sizeof(RealType) * x_complex.size());
				};

			std::memcpy(x0.data(), x_complex.data(), 2 * sizeof(RealType) * x_complex.size());
#pragma GCC diagnostic pop
			return RealSolver::compute(call_f_from_real, x0, MAX_ITER);
		}
	};
}