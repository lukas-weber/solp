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

// find a basis for the column space of A
// returns indices of linearly independent columns
static std::vector<int> find_column_basis(Eigen::MatrixXd &A) {
	auto QR = A.fullPivHouseholderQr();
	auto &R = QR.matrixQR().template triangularView<Eigen::Upper>();
	auto &P = QR.colsPermutation();

	std::vector<int> basis(QR.rank());
	for(size_t i = 0; i < basis.size(); i++) {
		for(int j = i; j < R.cols(); j++) {
			if(R(i, j) != 0) {
				basis[i] = P.indices()(j);
				break;
			}
		}
	}

	return basis;
}

// reorders the problem so that A starts with linearly independent columns.
static std::vector<int> simplex_initialize(Eigen::VectorXd &objective, Eigen::MatrixXd &A) {
	auto basis = find_column_basis(A);

	if(static_cast<Eigen::Index>(basis.size()) != A.rows()) {
		throw exception{exception::type::rank_deficient};
	}

	std::vector<int> idxs = basis;
	for(int i = 0; i < A.cols(); i++) {
		if(std::find(basis.begin(), basis.end(), i) == basis.end()) {
			idxs.push_back(i);
		}
	}

	auto Atmp = A;
	auto objectivetmp = objective;
	for(int i = 0; i < A.cols(); i++) {
		A.col(i) = Atmp.col(idxs[i]);
		objective(i) = objectivetmp(idxs[i]);
	}
	return idxs;
}

static result revised_simplex(Eigen::VectorXd &objective, Eigen::MatrixXd &A, Eigen::VectorXd &b) {
	std::vector<int> idxs = simplex_initialize(objective, A);
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
			if(sn(i) < 0 && idxs[nb + i] < minidx) {
				q = nb + i;
				minidx = idxs[q];
			}
		}

		auto Bfac = A.leftCols(nb).colPivHouseholderQr();
		xb = Bfac.solve(b.head(nb));

		if(q < 0) {
			break;
		}

		d = Bfac.solve(A.col(q));

		int p = -1;
		double rmin = INFINITY;
		for(int i = 0; i < nb; i++) {
			if(d(i) > 0 && xb(i) / d(i) < rmin) {
				p = i;
				rmin = xb(i) / d(i);
			}
		}
		if(p < 0) {
			throw exception{exception::type::unbounded};
		}

		std::swap(idxs[p], idxs[q]);
		A.col(p).swap(A.col(q));
		std::swap(objective(p), objective(q));
	}

	result res;
	res.x.resize(A.cols());
	for(int i = 0; i < nb; i++) {
		res.x[idxs[i]] = xb[i];
		if(xb[i] < 0) {
			throw exception{exception::type::infeasible};
		}
	}
	return res;
}

result solve(const std::vector<double> &objective, const std::vector<constraint> &constraints) {
	Eigen::VectorXd obj = Eigen::Map<const Eigen::VectorXd>(objective.data(), objective.size());

	Eigen::MatrixXd A(constraints.size(), obj.size());
	Eigen::VectorXd b(constraints.size());

	for(int i = 0; i < A.rows(); i++) {
		A.row(i) = Eigen::Map<const Eigen::RowVectorXd>(constraints[i].coeff.data(),
		                                                constraints[i].coeff.size());
		b(i) = constraints[i].rhs;
	}

	return revised_simplex(obj, A, b);
}

}
