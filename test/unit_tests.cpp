#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "solp.cpp"

TEST_CASE("canonical problem") {
	// example from wikipedia
	std::vector<double> objective = {0, 0, -2, -3, -4};
	std::vector<solp::constraint> A = {
	    {{1, 0, 3, 2, 1}, 10},
	    {{0, 1, 2, 5, 3}, 15},
	};

	auto res = solp::solve(objective, A);

	Eigen::VectorXd x(5);
	x << 5, 0, 0, 0, 5;
	Eigen::Map<Eigen::VectorXd> xres(res.x.data(), 5);
	REQUIRE(xres.isApprox(x));
}

TEST_CASE("basic infeasible") {
	std::vector<double> objective = {-1, -1};

	SECTION("negative region") {
		std::vector<solp::constraint> A{
		    {{1, 1}, -1},
		};

		CHECK_THROWS(solp::solve(objective, A));
		try {
			solp::solve(objective, A);
		} catch(const solp::exception &e) {
			CHECK(e.status == solp::exception::type::infeasible);
		}
	}

	SECTION("empty region") {
		std::vector<solp::constraint> A{
		    {{1, 1}, 2},
		    {{1, 1}, 3},
		};

		CHECK_THROWS(solp::solve(objective, A));
		try {
			solp::solve(objective, A);
		} catch(const solp::exception &e) {
			e.what();
			CHECK(e.status == solp::exception::type::rank_deficient);
		}
	}
}

TEST_CASE("example problem") {
	const double pCoff = 0.2, pCook = 0.5;
	const double vCook = 0.5;
	const double B = 20;
	const double V = 15;

	const double aCoff = 2, aCook = 5.5;
	std::vector<double> objective = {-aCoff, -aCook, 0};
	std::vector<solp::constraint> constraints = {
	    {{pCoff, pCook, 0}, B},
	    {{0, vCook, 1}, V},
	};

	CHECK_NOTHROW(solp::solve(objective, constraints));
	solp::result res = solp::solve(objective, constraints);
	CHECK(res.x[0] == Approx(25));
	CHECK(res.x[1] == Approx(30));
	CHECK(res.x[2] == Approx(0));
}

TEST_CASE("sse example") {
	std::vector<solp::constraint> constraints = {
	    {{1, 1, 1, 0, 0, 0}, 0.75}, {{0, 1, 0, 1, 1, 0}, 0}, {{0, 0, 1, 0, 1, 1}, 0}};

	std::vector<double> objective = {1, 0, 0, 1, 0, 1};

	solp::result res = solp::solve(objective, constraints);
	Eigen::VectorXd sol(6);
	sol << 0.75, 0, 0, 0, 0, 0;
	CHECK(sol.isApprox(Eigen::Map<Eigen::VectorXd>(res.x.data(), res.x.size())));
}

TEST_CASE("cycling") {
	std::vector<double> obj = {1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0};

	
	std::vector<solp::constraint> constraints = {
	{{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 1.0},
	{{0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 1.0},
	{{0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 1.0},
	{{0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 1.0},
	{{0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 0.0},
	{{0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0}, 0.0},
	{{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0}, 0.0},
	{{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0}, 0.0}};

	CHECK_NOTHROW(solp::solve(obj, constraints));
}

TEST_CASE("slack dropping") {
	Eigen::MatrixXd A(2, 5);
	A << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10;

	std::vector<int> idxs1 = {1, 3, 2, 0, 4};
	CHECK_THROWS(solp::drop_slacks(idxs1, A));
	std::vector<int> idxs2 = {2, 3, 0, 4, 1};

	auto Anew = solp::drop_slacks(idxs2, A);
	Eigen::MatrixXd Acheck(2, 3);
	Acheck << 1, 2, 4, 6, 7, 9;
	CHECK(idxs2 == std::vector<int>{0, 1, 2});
	CHECK(Acheck == Anew);
}
