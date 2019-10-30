#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# txmn: A tool for inferring endosymbiont biology from metagenome data.
# Copyright (C) 4-11-2018 Mark Juers

#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.

import numpy as np
import pytest

from txmn.base import polymorphic_to_dict
from txmn.fixtures import single_replicate_meta
from txmn.fixtures import double_replicate_meta
from txmn.fixtures import double_replicate_meta_exclude_pop_1
from txmn.fixtures import double_replicate_meta_3_popl_exclude_pop_2

# These must match the args in txmn.fixtures.<>_replicate
num_samples = 4
num_pop = 2

# single replicate
@pytest.mark.parametrize(
    "bias, expected",
    [
        (
            False,
            np.array(
                [[0.375, 0.375, 0.375, 0.5, 0.375], [0.0, 0.0, 0.0, 0.0, 0.5]]
            ),
        ),
        (
            True,
            4.0
            / (4.0 - 1.0)
            * np.array(
                [[0.375, 0.375, 0.375, 0.5, 0.375], [0.0, 0.0, 0.0, 0.0, 0.5]]
            ),
        ),
    ],
)
def test_h_no_average_one_rep(bias, expected, single_replicate_meta):
    assert np.isclose(
        single_replicate_meta.h(bias=bias, by_population=True)[0], expected
    ).all()


@pytest.mark.parametrize(
    "bias, expected",
    [
        (False, np.array([[0.400], [0.1]])),
        (True, 4.0 / (4.0 - 1.0) * np.array([[0.400], [0.100]])),
    ],
)
def test_h_average_one_rep(bias, expected, single_replicate_meta):
    assert np.isclose(
        single_replicate_meta.h(bias=bias, by_population=True, average=True),
        expected,
    ).all()


@pytest.mark.parametrize(
    "polymorphic_opts, expected",
    [
        ({"threshold": 0}, np.array([[2.727272727], [0.54545454]])),
        ({"threshold": 1}, np.array([[0.54545454], [0.54545454]])),
    ],
)
def test_theta_w_one_rep(polymorphic_opts, expected, single_replicate_meta):
    assert np.isclose(
        expected, single_replicate_meta.theta_w(True, polymorphic_opts)
    ).all()


# double replicate
@pytest.mark.parametrize(
    "bias, expected",
    [
        (
            False,
            [
                np.array(
                    [
                        [0.375, 0.375, 0.375, 0.5, 0.375],
                        [0.0, 0.0, 0.0, 0.0, 0.5],
                    ]
                ),
                np.array(
                    [
                        [0.375, 0.375, 0.375, 0.0, 0.375, 0.0, 0.0],
                        [0.375, 0.0, 0.0, 0.5, 0.375, 0.375, 0.375],
                    ]
                ),
            ],
        ),
        (
            True,
            [
                4
                / (4 - 1)
                * np.array(
                    [
                        [0.375, 0.375, 0.375, 0.5, 0.375],
                        [0.0, 0.0, 0.0, 0.0, 0.5],
                    ]
                ),
                4
                / (4 - 1)
                * np.array(
                    [
                        [0.375, 0.375, 0.375, 0.0, 0.375, 0.0, 0.0],
                        [0.375, 0.0, 0.0, 0.5, 0.375, 0.375, 0.375],
                    ]
                ),
            ],
        ),
    ],
)
def test_h_no_average_two_reps(bias, expected, double_replicate_meta):
    heterozygosity = double_replicate_meta.h(by_population=True, bias=bias)
    out = [
        np.isclose(expected_local, heterozygosity[i]).all()
        for i, expected_local in enumerate(expected)
    ]
    assert all(out)


@pytest.mark.parametrize(
    "bias, expected",
    [
        (False, np.array([[0.4, 0.214285714], [0.1, 0.285714286]])),
        (
            True,
            4.0
            / (4.0 - 1.0)
            * np.array([[0.4, 0.214285714], [0.1, 0.285714286]]),
        ),
    ],
)
def test_h_average_two_reps(bias, expected, double_replicate_meta):
    heterozygosity = double_replicate_meta.h(
        by_population=True, bias=bias, average=True
    )
    out = [
        np.isclose(heterozygosity[i], expected_local).all()
        for i, expected_local in enumerate(expected)
    ]
    assert all(out)


@pytest.mark.parametrize(
    "by_population, polymorphic_opts, expected",
    [
        (
            True,
            {"threshold": 0},
            np.array([[2.72727273, 2.18181818], [0.54545455, 2.72727273]]),
        ),
        (
            True,
            {"threshold": 1},
            np.array([[0.54545455, 1.09090909], [0.54545455, 0.54545455]]),
        ),
        (False, {"threshold": 0}, np.array([[1.928374656, 2.699724518]])),
        (False, {"threshold": 1}, np.array([[0.771349862, 1.157024793]])),
    ],
)
def test_theta_w_two_reps(
    by_population, polymorphic_opts, expected, double_replicate_meta
):
    out = np.isclose(
        double_replicate_meta.theta_w(by_population, polymorphic_opts),
        expected,
    )
    assert out.all()


@pytest.mark.parametrize(
    "method, expected",
    [
        ("tajima", np.array([1.714285714, 2.571428571])),
        ("nei", np.array([1.5, 2.25])),
    ],
)
def test_pi_two_reps(method, expected, double_replicate_meta):
    out = np.isclose(double_replicate_meta.pi(pi_method=method), expected)
    assert out.all()


@pytest.mark.parametrize(
    "expected",
    [(np.array([[1, 1, 1, 2, 3]]), np.array([[4, 1, 1, 2, 4, 1, 1]]))],
)
def test_num_mutants(expected, double_replicate_meta):
    test = double_replicate_meta.num_mutants()
    assert all(np.array_equal(test[i], expected[i]) for i in range(len(test)))


@pytest.mark.parametrize(
    "by_population, threshold, expected",
    [
        (
            False,
            0,
            (
                {"which": (np.zeros(5), np.arange(5)), "num": np.array([5])},
                {"which": (np.zeros(7), np.arange(7)), "num": np.array([7])},
            ),
        ),
        (
            False,
            1,
            (
                {
                    "which": (np.array([0, 0]), np.array([3, 4])),
                    "num": np.array([2]),
                },
                {
                    "which": (np.array([0, 0, 0]), np.array([0, 3, 4])),
                    "num": np.array([3]),
                },
            ),
        ),
        (
            True,
            0,
            (
                {
                    "which": (
                        np.array([0, 0, 0, 0, 0, 1]),
                        np.array([0, 1, 2, 3, 4, 4]),
                    ),
                    "num": np.array([5, 1]),
                },
                {
                    "which": (
                        np.array([0, 0, 0, 0, 1, 1, 1, 1, 1]),
                        np.array([0, 1, 2, 4, 0, 3, 4, 5, 6]),
                    ),
                    "num": np.array([4, 5]),
                },
            ),
        ),
        (
            True,
            1,
            (
                {
                    "which": (np.array([0, 1]), np.array([3, 4])),
                    "num": np.array([1, 1]),
                },
                {
                    "which": (np.array([0, 0, 1]), np.array([0, 4, 3])),
                    "num": (np.array([2, 1])),
                },
            ),
        ),
    ],
)
def test_polymorphic_two_reps(
    by_population, threshold, expected, double_replicate_meta
):
    test = tuple(double_replicate_meta.polymorphic(by_population, threshold))
    assert hash(frozenset(polymorphic_to_dict(expected))) == hash(
        frozenset(polymorphic_to_dict(test))
    )


@pytest.mark.parametrize("expected", [np.array([5, 7])])
def test_segsites_two_reps(expected, double_replicate_meta):
    test = double_replicate_meta.segsites()
    assert all(test == expected)


@pytest.mark.parametrize(
    "by_population, expected",
    [
        (
            True,
            [
                np.array([[0.375, 0.375, 0.375, 0.5, 0.375]]),
                np.array([[0.375, 0.375, 0.375, 0, 0.375, 0, 0]]),
            ],
        ),
        (
            False,
            [
                np.array([[0.375, 0.375, 0.375, 0.5, 0.375]]),
                np.array([[0.375, 0.375, 0.375, 0, 0.375, 0, 0]]),
            ],
        ),
    ],
)
def test_h_no_average_two_reps_exclude_pop_1(
    by_population, expected, double_replicate_meta_exclude_pop_1
):
    out = []
    h = double_replicate_meta_exclude_pop_1.h(
        bias=False, by_population=by_population
    )
    for i, x in enumerate(h):
        out.append(np.isclose(x, expected[i]).all())
    assert all(out)


@pytest.mark.parametrize("expected", [np.array([[2.72727273, 2.18181818]])])
def test_theta_w_exclude_pop_1(expected, double_replicate_meta_exclude_pop_1):
    theta_w = double_replicate_meta_exclude_pop_1.theta_w(by_population=True)
    assert np.isclose(expected, theta_w).all()


@pytest.mark.parametrize(
    "expected",
    [
        (
            np.array(
                [
                    0.14285714,
                    0.0,
                    0.14285714,
                    0.14285714,
                    0.14285714,
                    0.14285714,
                    0.06666667,
                    0.14285714,
                    0.14285714,
                    0.14285714,
                ]
            ),
            np.array([0.14285714, 0.0, 0.14285714, 0.14285714]),
        )
    ],
)
def test_fst_exclude_pop_2(
    expected, double_replicate_meta_3_popl_exclude_pop_2
):
    out = double_replicate_meta_3_popl_exclude_pop_2.fst(
        average_sites=False, h_opts={"bias": False}
    )
    test_results = tuple(
        True if np.isclose(x, y).all() else False
        for x, y in zip(out, expected)
    )
    assert all(test_results)


@pytest.mark.parametrize(
    "expected",
    [
        (
            (
                {
                    "which": (np.array([1]), np.array([4])),
                    "num": np.array([1]),
                },
                {
                    "which": (
                        np.array([1, 1, 1, 1, 1]),
                        np.array([0, 3, 4, 5, 6]),
                    ),
                    "num": np.array([5]),
                },
            )
        )
    ],
)
def test_polymorphic_exclude_pop_1(
    expected, double_replicate_meta_exclude_pop_1
):
    out = double_replicate_meta_exclude_pop_1.polymorphic(output="num")
    assert hash(frozenset(polymorphic_to_dict(expected))) == hash(
        frozenset(polymorphic_to_dict(out))
    )


@pytest.mark.parametrize(
    "expected, method",
    [(np.array([2.0, 1.5]), "nei"), (np.array([2.66666667, 2.0]), "tajima")],
)
def test_pi_exclude_pop_1(
    expected, method, double_replicate_meta_exclude_pop_1
):
    out = np.isclose(
        expected, double_replicate_meta_exclude_pop_1.pi(pi_method=method)
    )
    assert out.all()


@pytest.mark.parametrize(
    "expected",
    [
        (
            np.array(
                [
                    [1, 0, 1, 1, 1, 0, 1, 0, 0, 0],
                    [1, 1, 1, 1, 1, 0, 0, 0, 0, 0],
                    [1, 0, 1, 1, 1, 0, 1, 0, 0, 0],
                    [1, 0, 1, 1, 1, 0, 1, 0, 0, 0],
                    [1, 0, 1, 1, 1, 0, 1, 0, 0, 0],
                    [1, 0, 1, 1, 1, 0, 1, 0, 0, 0],
                    [1, 1, 1, 1, 1, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 1, 0, 1, 1, 1],
                ],
                dtype=np.uint8,
            ),
            np.array(
                [
                    [0, 1, 0, 0],
                    [0, 1, 0, 0],
                    [0, 0, 0, 0],
                    [0, 0, 0, 0],
                    [1, 1, 1, 0],
                    [0, 1, 0, 1],
                    [0, 0, 0, 0],
                    [0, 0, 0, 0],
                ],
                dtype=np.uint8,
            ),
        )
    ],
)
def test_filter_variants(expected, double_replicate_meta_3_popl_exclude_pop_2):
    test = []
    for rep in double_replicate_meta_3_popl_exclude_pop_2.filter_variants():
        out = []
        for site in rep:
            out.append(site.genotypes)
        test.append(np.array(out).T)
    assert all((x == y).all() for x, y in zip(test, expected))
