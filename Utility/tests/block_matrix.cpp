#include "../include/mrock/utility/Numerics/PivotToBlockStructure.hpp"
#include <iostream>
#include <chrono>
//#include "../include/mrock/utility/print_type.hpp"

using namespace mrock::utility::Numerics;

Eigen::MatrixXd generateRandomHermitian(int blockSize) {
    Eigen::MatrixXd mat = Eigen::MatrixXd::Random(blockSize, blockSize);
    Eigen::MatrixXd hermitianMat = (mat + mat.adjoint()) / 2.0;
    return hermitianMat;
}

Eigen::MatrixXd constructBlockDiagonalMatrix(int blockSize, int numBlocks) {
    Eigen::MatrixXd blockDiagMatrix = Eigen::MatrixXd::Zero(blockSize * numBlocks, blockSize * numBlocks);

    for (int i = 0; i < numBlocks; ++i) {
        blockDiagMatrix.block(i * blockSize, i * blockSize, blockSize, blockSize) = generateRandomHermitian(blockSize);
    }

    return blockDiagMatrix;
}

inline double rel_error(const Eigen::MatrixXd& original, const Eigen::MatrixXd& other) {
    return (original - other).norm() / original.norm();
}

using BlockMatrix = BlockDiagonalMatrix<Eigen::MatrixXd>;

int main() {
    const int N = 250;
    const int n_blocks = 4;
    Eigen::MatrixXd toSolve = constructBlockDiagonalMatrix(N, n_blocks);

    std::vector<HermitianBlock> block_indizes = identify_hermitian_blocks(toSolve);
    BlockMatrix blocked_toSolve(toSolve, block_indizes);

    if((toSolve - blocked_toSolve.construct_matrix()).norm() > 1e-12) {
        return 1;
    }

    Eigen::MatrixXd second_matrix = constructBlockDiagonalMatrix(N, n_blocks);
    std::vector<HermitianBlock> second_block_indizes = identify_hermitian_blocks(second_matrix);
    BlockMatrix second_blocked(second_matrix, second_block_indizes);

    {
        Eigen::MatrixXd tester = toSolve + second_matrix;
        const auto blocked_tester = blocked_toSolve + second_blocked;
        const double error = rel_error(tester, blocked_tester.construct_matrix());
        if(error > 1e-12) {
            std::cerr << "Addition failed " << error << std::endl;
            return 2;
        }
    }
    {
        Eigen::MatrixXd tester = toSolve - second_matrix;
        const auto blocked_tester = blocked_toSolve - second_blocked;
        const double error = rel_error(tester, blocked_tester.construct_matrix());
        if(error > 1e-12) {
            std::cerr << "Substraction failed " << error << std::endl;
            return 3;
        }
    }
    {
        Eigen::MatrixXd tester = toSolve * second_matrix;
        const auto blocked_tester = blocked_toSolve * second_blocked;
        const double error = rel_error(tester, blocked_tester.construct_matrix());
        if(error > 1e-12) {
            std::cerr << "Multiplication failed " << error << std::endl;
            return 4;
        }
    }
    {
        Eigen::MatrixXd smaller = Eigen::MatrixXd::Random(N * n_blocks, 33);
        Eigen::MatrixXd result = toSolve * smaller;
        Eigen::MatrixXd result_blocked = blocked_toSolve * smaller;
        double error = rel_error(result, result_blocked);
        if(error > 1e-12) {
            std::cerr << "toSolve * smaller failed " << error << std::endl;
            return 5;
        }
        result = smaller.transpose() * toSolve;
        result_blocked = smaller.transpose() * blocked_toSolve;
        error = rel_error(result, result_blocked);
        if(error > 1e-12) {
            std::cerr << "smaller.transpose() * toSolve failed " << error << std::endl;
            return 5;
        }
    }

    auto start = std::chrono::high_resolution_clock::now();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen_solver(toSolve);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Diagonalization time (Eigen): " << duration.count() << " ms" << std::endl;

    start = std::chrono::high_resolution_clock::now();
    blocked_matrix_wrapper<double> blocked_solver = blocked_matrix_wrapper<double>::solve_block_diagonal_matrix(blocked_toSolve);
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Diagonalization time (Blocked): " << duration.count() << " ms" << std::endl;

    Eigen::MatrixXd recon_eigen = eigen_solver.eigenvectors() * eigen_solver.eigenvalues().asDiagonal() * eigen_solver.eigenvectors().adjoint();
    Eigen::MatrixXd recon_block = blocked_solver.reconstruct_matrix_as_eigen();
    const double error = rel_error(recon_eigen, recon_block);
    if(error > 1e-12) {
        std::cerr << "Diagonalization failed " << error << std::endl;
        return 6;
    }

    return 0;
}