#define CATCH_CONFIG_MAIN
#include <catch.hpp>

#include "solp.cpp"

TEST_CASE("canonical problems") {
	// example from wikipedia
	Eigen::VectorXd objective(5);
	objective << 0,0,-2,-3,-4;
	Eigen::MatrixXd A(2,5);
	A << 1, 0, 3, 2, 1,
	     0, 1, 2, 5, 3;
	Eigen::VectorXd b(2);
	b << 10, 15;

	auto res = solp::revised_simplex(objective, A, b);

	Eigen::VectorXd x(5);
	x << 5,0,0,0,5;
	Eigen::Map<Eigen::VectorXd> xres(res.x.data(),5);
	REQUIRE(xres.isApprox(x));
}

TEST_CASE("basic infeasible") {
	Eigen::VectorXd objective(2);
	objective << -1, -1;

	SECTION("negative region") {
		Eigen::MatrixXd A(1,2);
		A << 1, 1;
		Eigen::VectorXd b(1);
		b << -1;

		CHECK_THROWS(solp::revised_simplex(objective, A, b));
		try {
			solp::revised_simplex(objective, A, b);
		} catch(const solp::exception &e) {
			CHECK(e.status == solp::exception::type::infeasible);
		}
	}

	SECTION("empty region") {
		Eigen::MatrixXd A(2,2);
		A << 1, 1,
		     1, 1;
		Eigen::VectorXd b(2);
		b << 2, 3;

		CHECK_THROWS(solp::revised_simplex(objective, A, b));
		try {
			solp::revised_simplex(objective, A, b);
		} catch(const solp::exception &e) {
			e.what();
			CHECK(e.status == solp::exception::type::rank_deficient);
		}

	}
}

TEST_CASE("invalid exception") {
	std::cout << solp::exception{static_cast<solp::exception::type>(9000)}.what();
}

TEST_CASE("column basis") {
	static const int cols = 5;
	static const int rows = 2;

	int c1 = GENERATE(range(0,cols));
	int c2 = (c1 + GENERATE(range(1,cols)))%cols;

	
		
	Eigen::MatrixXd A(rows,cols);
	A.col(c1) << 1,1.2;
	A.col(c2) << -1,3.3;
	for(int i = 0; i < A.cols(); i++) {
		if(i != c1 && i != c2) {
			A.col(i) = A.col(c1)*GENERATE(-1.2,0.2) + A.col(c2)*GENERATE(-0.3,0.,2.);
		}
	}


	auto b = solp::find_column_basis(A);
	Eigen::MatrixXd B(A.rows(),A.rows());
	for(int i = 0; i < A.rows(); i++) {
		B.col(i) = A.col(b[i]);
	}

	int rank = B.colPivHouseholderQr().rank();
	REQUIRE(rank == A.rows());
}
