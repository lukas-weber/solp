#include "solp.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>

namespace solp {

using SparseMatrix = Eigen::SparseMatrix<double>;

// Type definition to provide a QR Decomposition irrespective for both sparse and dense matrices.
// Eigen does not have something like this out of the box, as far as i know.
template<typename Matrix>
using QRDecomp =
    typename std::conditional_t<std::is_same_v<Matrix, SparseMatrix>,
                                Eigen::SparseQR<SparseMatrix, Eigen::COLAMDOrdering<int>>,
                                Eigen::ColPivHouseholderQR<Matrix>>;

const char *exception::what() const noexcept {
	switch(status) {
	case type::rank_deficient:
		return "constraint matrix does not have full row rank";
	case type::infeasible:
		return "constraints are incompatible/infeasible";
	case type::unbounded:
		return "problem is unbounded";
	case type::input_format:
		return "constraint arrays do not have matching dimensions";
	default:
		return "invalid exception";
	}
}

template<typename Matrix>
static Matrix drop_slacks(std::vector<int> &colperm, const Matrix &A) {
	Matrix Anew(A.rows(), A.cols() - A.rows());

	int j = 0;
	for(int i = 0; i < A.cols(); i++) {
		if(colperm[i] < A.rows()) {
			if(i < A.rows()) {
				// slack cannot be dropped
				throw exception{exception::type::infeasible};
			}
		} else {
			Anew.col(j) = A.col(i);
			colperm[j] = colperm[i] - A.rows();
			j++;
		}
	}
	colperm.resize(Anew.cols());
	assert(j == Anew.cols());
	return Anew;
}

// In the convention here, the current corner of the simplex is given by
// xb = B^-1 b
//
// where B is A.leftCols(nb). The algorithm assumes that we start in a corner in the positive cone
// (xb >= 0).
//
template<typename Matrix>
static result revised_simplex(Eigen::VectorXd &objective, Matrix &A, const Eigen::VectorXd &b,
                              std::vector<int> &colperm, const options &opts) {
	assert(A.cols() >= A.rows());

	int nb = A.rows();
	int nn = A.cols() - A.rows();

	Eigen::VectorXd xb(nb);
	Eigen::VectorXd d(nb);

	Eigen::VectorXd sn(nn);
	Eigen::VectorXd lambda(nn);

	while(true) {
		lambda = QRDecomp<Matrix>(A.leftCols(nb).transpose()).solve(objective.head(nb));
		sn = objective.tail(nn) - A.rightCols(nn).transpose() * lambda;

		int q{-1};
		double smin = INFINITY;
		for(int i = 0; i < nn; i++) {
			if(sn(i) < -opts.tolerance && sn[i] < smin) {
				q = nb + i;
				smin = sn[i];
			}
		}

		auto Bfac = QRDecomp<Matrix>(A.leftCols(nb));
		xb = Bfac.solve(b);

		if(q < 0) {
			break;
		}

		d = Bfac.solve(A.col(q));

		int p = -1;
		double rmin = INFINITY;
		for(int i = 0; i < nb; i++) {
			if(d(i) > opts.tolerance && xb(i) / d(i) < rmin) {
				p = i;
				rmin = xb(i) / d(i);
			}
		}

		if(p < 0) {
			throw exception{exception::type::unbounded};
		}

		std::swap(colperm[p], colperm[q]);
		// XXX: Eigen does not provide .swap() for sparse matrices. This is a dirty workaround to
		// still get an in-place swap.
		A.col(p) += A.col(q);
		A.col(q) = A.col(p) - A.col(q);
		A.col(p) -= A.col(q);
		std::swap(objective(p), objective(q));
	}

	result res;
	res.x.resize(A.cols());
	for(int i = 0; i < nb; i++) {
		res.x[colperm[i]] = xb[i];
		if(xb[i] < -opts.tolerance) {
			throw exception{exception::type::infeasible};
		}
	}
	return res;
}

// The method used here is a two-phase revised simplex method. In the presolve step, slack variables
// are introduced so we are sure we start from a feasible solution. Then, we try to find a solution
// where they are zero using LP.
//
// If it exists, we can drop them and start from that solution to solve our original objective.
//
template<typename Matrix>
static result solve_common(const std::vector<double> &objective, Matrix &&A,
                           const Eigen::VectorXd &b, const options &opts) {
	Eigen::VectorXd presolve_obj = Eigen::VectorXd::Zero(A.cols());
	presolve_obj.head(A.rows()) = Eigen::VectorXd::Ones(A.rows());

	if(QRDecomp<Matrix>(A.rightCols(objective.size())).rank() != A.rows()) {
		throw exception{exception::type::rank_deficient};
	}

	std::vector<int> colperm(A.cols());
	for(int i = 0; i < A.cols(); i++) {
		colperm[i] = i;
	}

	revised_simplex(presolve_obj, A, b, colperm, opts);

	A = drop_slacks(colperm, A);

	Eigen::VectorXd obj(objective.size());
	for(int i = 0; i < A.cols(); i++) {
		obj[i] = objective[colperm[i]];
	}

	return revised_simplex(obj, A, b, colperm, opts);
}

result solve(const std::vector<double> &objective, const std::vector<constraint> &constraints,
             const options &opts) {
	Eigen::MatrixXd A(constraints.size(), constraints.size() + objective.size());
	Eigen::VectorXd b(constraints.size());
	A.leftCols(A.rows()) = Eigen::MatrixXd::Identity(A.rows(), A.rows());

	for(int i = 0; i < A.rows(); i++) {
		if(constraints[i].coeff.size() != objective.size()) {
			throw exception{exception::type::input_format};
		}

		A.rightCols(objective.size()).row(i) = Eigen::Map<const Eigen::RowVectorXd>(
		    constraints[i].coeff.data(), constraints[i].coeff.size());
		b(i) = constraints[i].rhs;

		if(b(i) < 0) {
			A(i, i) *= -1;
		}
	}

	return solve_common(objective, std::move(A), b, opts);
}

result solve(const std::vector<double> &objective,
             const std::vector<constraint_sparse> &constraints, const options &opts) {
	SparseMatrix A(constraints.size(), constraints.size() + objective.size());
	Eigen::VectorXd b(constraints.size());
	Eigen::VectorXd presolve_obj = Eigen::VectorXd::Zero(A.cols());

	std::vector<Eigen::Triplet<double>> triplets;
	triplets.reserve(A.rows() * 5);

	for(int i = 0; i < A.rows(); i++) {
		b(i) = constraints[i].rhs;

		triplets.push_back({i, i, b(i) < 0 ? -1.0 : 1.0});
		if(constraints[i].coeff.size() != constraints[i].idx.size()) {
			throw exception{exception::type::input_format};
		}

		for(int j = 0; j < static_cast<int>(constraints[i].coeff.size()); j++) {
			double coeff = constraints[i].coeff[j];
			int idx = constraints.size() + constraints[i].idx[j];
			if(idx > A.cols()) {
				throw exception{exception::type::input_format};
			}

			triplets.push_back({i, idx, coeff});
		}
	}

	A.setFromTriplets(triplets.begin(), triplets.end());

	return solve_common(objective, std::move(A), std::move(b), opts);
}

}
