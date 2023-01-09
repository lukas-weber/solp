#ifndef solp_h_INCLUDED
#define solp_h_INCLUDED

#include <stdexcept>
#include <vector>
namespace solp {

class exception : public std::exception {
public:
	enum class type {
		rank_deficient,
		infeasible,
		unbounded,
	};

	explicit exception(type t) : status{t} {}

	type status;

	virtual const char *what() const noexcept override;

private:
	constexpr const char *msg(type t);
};

struct result {
	std::vector<double> x;
};

// constraint represents an equality constraint on the variables x
//
//    coeff * x = rhs
//
struct constraint {
	std::vector<double> coeff;
	double rhs;
};

// solve the linear programming problem
//
//     Minimize objective * x
//     with
//         x >= 0
//         and the given constraints.
//
// Note that constraints only encode equality constraints. However,
// every linear programming problem can be brought into this standard form
// by shifting x and introducing additional slack variables.
//
// If there is no solution or the problem is unbounded, a solp::exception is thrown.
//
result solve(const std::vector<double> &objective, const std::vector<constraint> &constraints);

}
#endif // solp_h_INCLUDED
