#!/usr/bin/env python

import numpy as np
import scipy.optimize as spo
import argparse
import textwrap

testfile_header='''// generated by lpgen.py
#include <catch2/catch.hpp>
#include <solp.h>
#include <Eigen/Dense>
#include "helpers.h"
TEST_CASE("random matrices") {
'''
types = {
    'valid' : 0,
    'infeasible' : 2,
    'unbounded' : 3,
}
    

def lpgen(m, n, prob_type):
    while True:
        A = np.random.normal(size=(m,n))
        b = np.random.normal(size=(m))
        c = -np.random.exponential(size=n)

        if np.linalg.matrix_rank(A) != A.shape[0]:
            continue

        res = spo.linprog(c, A_eq=A, b_eq=b, options={'tol':1e-11})

        if res.status != types[prob_type]:
            continue

        return A, b, c, res.x

def eigen_comma(a):
    return textwrap.fill(', '.join([str(x) for x in a.flatten()]),120)

def print_cpp_tests(outfile, cases):
    with open(outfile, 'w') as f:
        f.write(testfile_header)
        for case in cases:
            f.write('''SECTION("{type} {m}x{n}: {N}") {{
            solp::result res;
            std::vector<solp::constraint> A;
            std::vector<double> c;

            '''.format(**case))
            
            for n in range(case['N']):
                A, b, c, x = lpgen(case['m'], case['n'], case['type'])
                Astr = '{'+',\n'.join(['{{'+eigen_comma(A[i,:])+'},'+str(b[i])+'}' for i in range(len(b))])+'}'
                f.write('''A = {a};
                c = {{{c}}};

                '''.format(a = Astr,
                                    c = eigen_comma(c)))

                for A in ['A', 'convert_to_sparse(A)']:

                    if case['type'] == 'valid':
                        f.write('''res = solp::solve(c, {A});
                        CHECK_THAT(res.x, Catch::Matchers::Approx(std::vector<double>{{{solution}}}).margin(1e-10));
                        '''.format(A = A, solution = eigen_comma(x)))
                    elif case['type'] == 'unbounded':
                        f.write('''CHECK_THROWS_MATCHES(solp::solve(c, {A}), solp::exception, MatchSolpException(solp::exception::type::unbounded));
                        '''.format(A=A))
                    elif case['type'] == 'infeasible':
                        f.write('''CHECK_THROWS_MATCHES(solp::solve(c, {A}), solp::exception, MatchSolpException(solp::exception::type::infeasible));
                        '''.format(A=A))

                
            f.write('}\n')
        f.write('}\n')

parser = argparse.ArgumentParser(description='generator for linear programming tests on random matrices')
parser.add_argument('outfile')
args = parser.parse_args()

np.random.seed(0)

print_cpp_tests(args.outfile, [
    {'m': 2, 'n': 4, 'type': 'valid', 'N': 10},
    {'m': 3, 'n': 5, 'type': 'valid', 'N': 20},
    {'m': 10, 'n': 20, 'type': 'valid', 'N': 1},
    {'m': 3, 'n': 5, 'type': 'unbounded', 'N': 5},
    {'m': 3, 'n': 5, 'type': 'infeasible', 'N': 5},
])

