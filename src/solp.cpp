#include "solp.h"
#include <Eigen/Dense>
#include <iostream>

namespace solp {

const char *exception::what() const noexcept {
	switch(status) {
	case type::rank_deficient:
		return "constraint matrix does not have full row rank";
	case type::infeasible:
		return "constraints are incompatible/infeasible";
	case type::unbounded:
		return "problem is unbounded";
	default:
		return "invalid exception";
	}
}

static Eigen::MatrixXd drop_slacks(std::vector<int> &colperm, Eigen::MatrixXd &A) {
	Eigen::MatrixXd Anew(A.rows(), A.cols() - A.rows());

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
static result revised_simplex(Eigen::VectorXd &objective, Eigen::MatrixXd &A, Eigen::VectorXd &b,
                              std::vector<int> &colperm, const options &opts) {
	assert(A.cols() >= A.rows());

	int nb = A.rows();
	int nn = A.cols() - A.rows();

	Eigen::VectorXd xb(nb);
	Eigen::VectorXd d(nb);

	Eigen::VectorXd sn(nn);
	Eigen::VectorXd lambda(nn);

	while(true) {
		lambda = A.leftCols(nb).transpose().colPivHouseholderQr().solve(objective.head(nb));
		sn = objective.tail(nn) - A.rightCols(nn).transpose() * lambda;

		int q{-1};
		int minidx = A.cols() + 1;
		for(int i = 0; i < nn; i++) {
			if(sn(i) < -opts.tolerance && sn(i) < minidx) {
				q = nb + i;
				minidx = sn(i); //colperm[q];
			}
		}

		auto Bfac = A.leftCols(nb).colPivHouseholderQr();
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
		A.col(p).swap(A.col(q));
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
result solve(const std::vector<double> &objective, const std::vector<constraint> &constraints,
             const options &opts) {
	Eigen::VectorXd obj(objective.size());

	Eigen::MatrixXd A(constraints.size(), constraints.size() + obj.size());
	Eigen::VectorXd b(constraints.size());
	A.leftCols(A.rows()) = Eigen::MatrixXd::Identity(A.rows(), A.rows());
	Eigen::VectorXd presolve_obj = Eigen::VectorXd::Zero(A.cols());

	for(int i = 0; i < A.rows(); i++) {
		A.rightCols(obj.size()).row(i) = Eigen::Map<const Eigen::RowVectorXd>(
		    constraints[i].coeff.data(), constraints[i].coeff.size());
		b(i) = constraints[i].rhs;
		presolve_obj(i) = 1;

		if(b(i) < 0) {
			A(i, i) *= -1;
		}
	}

	if(A.rightCols(obj.size()).colPivHouseholderQr().rank() != A.rows()) {
		throw exception{exception::type::rank_deficient};
	}

	std::vector<int> colperm(A.cols());
	for(int i = 0; i < A.cols(); i++) {
		colperm[i] = i;
	}

	revised_simplex(presolve_obj, A, b, colperm, opts);

	A = drop_slacks(colperm, A);
	for(int i = 0; i < A.cols(); i++) {
		obj[i] = objective[colperm[i]];
	}

	return revised_simplex(obj, A, b, colperm, opts);
}

}
