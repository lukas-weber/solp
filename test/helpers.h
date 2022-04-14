#ifndef helpers_h_INCLUDED
#define helpers_h_INCLUDED

#include <algorithm>
#include <catch2/catch.hpp>
#include <solp.h>

// To test the sparse implementation, we are cheeky and just use the same test cases
// as in the dense case and convert it into sparse form (even if it is actually a dense problem).

inline std::vector<solp::constraint_sparse>
    convert_to_sparse(const std::vector<solp::constraint> &constraints) {
	std::vector<solp::constraint_sparse> result(constraints.size());
	std::transform(constraints.begin(), constraints.end(), result.begin(), [](auto c) {
		solp::constraint_sparse cs;
		for(int i = 0; i < static_cast<int>(c.coeff.size()); i++) {
			double coeff = c.coeff[i];
			if(coeff != 0) {
				cs.coeff.push_back(coeff);
				cs.idx.push_back(i);
			}
			cs.rhs = c.rhs;
		}
		return cs;
	});

	return result;
}

inline Catch::Generic::PredicateMatcher<solp::exception>
    MatchSolpException(solp::exception::type type) {
	return Catch::Matchers::Predicate<solp::exception>(
	    [=](const auto &e) { return e.status == type; }, solp::exception{type}.what());
}

#endif // helpers_h_INCLUDED
