#include "solp.h"
#include <iostream>
#include <Eigen/Dense>

namespace solp {

result solve(const std::vector<double> &objective, const std::vector<constraint> &constraints) {
	Eigen::VectorXd obj = Eigen::Map<const Eigen::VectorXd>(objective.data(), objective.size());

	Eigen::MatrixXd A(constraints.size(), obj.size());
	Eigen::VectorXd b;

	for(int i = 0; i < A.rows(); i++) {
		A.row(i) = Eigen::Map<const Eigen::RowVectorXd>(constraints[i].coeff.data(), constraints[i].coeff.size());
		b(i) = constraints[i].rhs;
	}

	result r;
	return r;

}

static std::vector<int> find_column_basis(Eigen::MatrixXd &A) {
	auto &QR = A.colPivHouseholderQr();
	auto &R = QR.matrixR().template triangularView<Eigen::Upper>();
	auto &P = QR.colsPermutation();

	assert(QR.rank() == A.rows());

	std::vector<int> basis(R.rows());
	for(int i = 0; i < R.rows(); i++) {
		for(int j = i; j < R.cols(); j++) {
			if(R(i,j) != 0) {
				basis[i] = P.indices()(j);
				break;
			}
		}
	}

	return basis;
}

static result canonical_simplex(Eigen::VectorXd &objective, Eigen::MatrixXd &A, Eigen::VectorXd &b) {
	assert(A.cols() > A.rows());

	int nb = A.rows();
	int nn = A.cols()-A.rows();

	Eigen::VectorXd xb(nb);
	Eigen::VectorXd d(nb);

	Eigen::VectorXd sn(nn);
	Eigen::VectorXd lambda(nn);

	//assert(A.leftCols(nb).isApprox(Eigen::MatrixXd::Identity(nb,nb))); // problem canonical?

	std::vector<int> idxs(A.cols());
	for(size_t i = 0; i < idxs.size(); i++) {
		idxs[i] = i;
	}


	while(true) {
		lambda = A.leftCols(nb).transpose().colPivHouseholderQr().solve(objective.head(nb));
		sn = objective.tail(nn) - A.rightCols(nn).transpose()*lambda;

		int q{-1};
		int minidx = A.cols()+1;
		for(int i = 0; i < nn; i++) {
			if(sn(i) < 0 && idxs[nb+i] < minidx) {
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
			if(d(i) > 0 && xb(i)/d(i) < rmin) {
				p = i;
				rmin = xb(i)/d(i);
			}
		}
		if(p < 0) {
			throw std::runtime_error("problem unbounded");
		}
		std::cout << p << "<-" << q << "\n";

		std::swap(idxs[p], idxs[q]);
		A.col(p).swap(A.col(q));
		std::swap(objective(p), objective(q));
	}

	result res;
	res.x.resize(A.cols());
	for(int i = 0; i < nb; i++) {
		res.x[idxs[i]] = xb[i];
		//res.basis[i] = idxs[i];	
	}
	return res;
}

}
