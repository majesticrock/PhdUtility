#pragma once
#include <Eigen/Dense>
#include <complex>
#include <type_traits>
#include <cstring>

namespace mrock::utility::Numerics::Roots {
	template<typename RealType, int t_vector_size>
	class BroydensMethodEigen {
		static_assert(std::is_floating_point<RealType>::value, "You're data type must be a floating point number");

		using MatrixType = Eigen::Matrix<RealType, t_vector_size, t_vector_size>;
		using VectorType = Eigen::Vector<RealType, t_vector_size>;

		using error_type = std::conditional_t< sizeof(RealType) >= sizeof(double), double, RealType >;

		VectorType F_old, F_new, delta_x, delta_F;
		MatrixType jacobian;

		VectorType tmp_vec;
		MatrixType tmp_mat;

	public:
		// the function must have the following signature void func(const VectorType& input, VectorType& output)
		template<class FunctionType>
		RealType compute(const FunctionType& func, VectorType& x0, const unsigned int MAX_ITER = 200)
		{
			const size_t DIM = x0.rows();
			// You may play around with EPS_X and EPS_F to your desire
			// EPS_X is the minimum distance between x_i and x_i+1
			// EPS_F is the minimum f(x)

			const error_type EPS_F = std::numeric_limits<error_type>::epsilon() * DIM;
			const error_type EPS_X = std::numeric_limits<error_type>::epsilon() * DIM;
			RealType diff_x{ 100 };
			RealType diff_F{ 100 };
			int iter_num{};

			F_old.setZero(DIM);
			F_new.setZero(DIM);
			delta_x.setZero(DIM);
			delta_F.setZero(DIM);
			tmp_vec.setZero(DIM);
			tmp_mat.setZero(DIM, DIM);
			jacobian.setIdentity(DIM, DIM);
			func(x0, F_new);

			while (diff_x > EPS_X && diff_F > EPS_F && iter_num++ <= MAX_ITER && F_new.norm() > EPS_F) {
				delta_x.noalias() = -jacobian * F_new;
				x0 += delta_x;
				diff_x = delta_x.norm();
				F_old = F_new;
				func(x0, F_new);
				delta_F.noalias() = F_new - F_old;
				diff_F = delta_F.norm();

				tmp_vec.noalias() = jacobian * delta_F;
				tmp_mat.noalias() = (delta_x - tmp_vec) * (delta_x.transpose() * jacobian) / (delta_x.dot(tmp_vec));
				// new estimate for the jacobian
				jacobian += tmp_mat;
			}
			// This method returns |F(x_final)|
			return F_new.norm();
		}

		void free_memory() {
			F_old.setZero(1);
			F_new.setZero(1);
			delta_x.setZero(1);
			delta_F.setZero(1);
			tmp_vec.setZero(1);
			tmp_mat.setZero(1,1);
			jacobian.setZero(1,1);
		}
	};

	template<typename RealType, int t_vector_size>
	class BroydensMethodEigen< std::complex<RealType>, t_vector_size > {
		static constexpr Eigen::Index double_size = t_vector_size == Eigen::Dynamic ? Eigen::Dynamic : 2 * t_vector_size;
		using VectorType = Eigen::Vector<std::complex<RealType>, t_vector_size>;
		using RealVector = Eigen::Vector<RealType, double_size>;

		using RealSolver = BroydensMethodEigen<RealType, t_vector_size>;

		RealSolver _solver;
	public:
	template<class FunctionType>
		// the function must have the following signature void func(const VectorType& input, VectorType& output)
		RealType compute(const FunctionType& func, VectorType& x_complex, const int MAX_ITER = 200)
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
			return _solver.compute(call_f_from_real, x0, MAX_ITER);
		}

		void free_memory() {
			_solver.free_memory();
		}
	};
}