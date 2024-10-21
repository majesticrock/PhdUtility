#include "../include/Utility/Numerics/PivotToBlockStructure.hpp"
#include <iostream>

using namespace Utility::Numerics;
using namespace Eigen;

MatrixXd generateRandomHermitian(int blockSize) {
    MatrixXd mat = MatrixXd::Random(blockSize, blockSize);
    MatrixXd hermitianMat = (mat + mat.adjoint()) / 2.0;
    return hermitianMat;
}

MatrixXd constructBlockDiagonalMatrix(int blockSize, int numBlocks) {
    MatrixXd blockDiagMatrix = MatrixXd::Zero(blockSize * numBlocks, blockSize * numBlocks);

    for (int i = 0; i < numBlocks; ++i) {
        blockDiagMatrix.block(i * blockSize, i * blockSize, blockSize, blockSize) = generateRandomHermitian(blockSize);
    }

    return blockDiagMatrix;
}

int main() {
    const int N = 100;
    const int n_blocks = 4;
    Eigen::MatrixXd toSolve = constructBlockDiagonalMatrix(N, n_blocks);

    std::vector<HermitianBlock> block_indizes = identify_hermitian_blocks(toSolve);
    BlockDiagonalMatrix<double> blocked_toSolve(toSolve, block_indizes);

    if((toSolve - blocked_toSolve.construct_matrix()).norm() > 1e-12) {
        return 1;
    }

    Eigen::MatrixXd second_matrix = constructBlockDiagonalMatrix(N, n_blocks);
    std::vector<HermitianBlock> second_block_indizes = identify_hermitian_blocks(second_matrix);
    BlockDiagonalMatrix<double> second_blocked(second_matrix, second_block_indizes);

    {
        Eigen::MatrixXd tester = toSolve + second_matrix;
        BlockDiagonalMatrix<double> blocked_tester = blocked_toSolve + second_blocked;
        if((tester - blocked_tester.construct_matrix()).norm() > 1e-12) {
            return 2;
        }
    }
    {
        Eigen::MatrixXd tester = toSolve - second_matrix;
        BlockDiagonalMatrix<double> blocked_tester = blocked_toSolve - second_blocked;
        if((tester - blocked_tester.construct_matrix()).norm() > 1e-12) {
            return 3;
        }
    }
    {
        Eigen::MatrixXd tester = toSolve * second_matrix;
        BlockDiagonalMatrix<double> blocked_tester = blocked_toSolve * second_blocked;
        if((tester - blocked_tester.construct_matrix()).norm() > 1e-12) {
            return 4;
        }
    }
    

    return 0;
}