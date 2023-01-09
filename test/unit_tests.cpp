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

	auto res = solp::canonical_simplex(objective, A, b);

	Eigen::VectorXd x(5);
	x << 5,0,0,0,5;
	Eigen::Map<Eigen::VectorXd> xres(res.x.data(),5);
	//REQUIRE(res.x == std::vector<double>{5,0,0,0,5});
	REQUIRE(xres.isApprox(x));
}

TEST_CASE("noncanonical") {
	Eigen::VectorXd objective(5);
	objective << -2,-3,-4;
	Eigen::MatrixXd A(2,5);
	A << 1, 0, 3, 2, 1,
	     0, 1, 2, 5, 3;
	Eigen::VectorXd b(2);
	b << 10, 15;

}

TEST_CASE("column basis") {
	Eigen::MatrixXd A(3,5);

	A << 2, 1, 3, 5, 2,
	     1, 3, 4, 2, 1,
	     2, 0, 2, 3, 4;

	auto b = find_column_basis(A);
	auto B = A(Eigen::All, b);
	auto &QR = 
}
