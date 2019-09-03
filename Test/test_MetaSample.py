#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# transmission: A tool for inferring endosymbiont biology from metagenome data.
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

from transmission.fixtures import single_replicate_meta
from transmission.fixtures import double_replicate_meta
from transmission.fixtures import double_replicate_meta_exclude_pop_1
from transmission.fixtures import double_replicate_meta_3_popl

# These must match the args in transmission.fixtures.<>_replicate
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
    "threshold, expected",
    [
        (0, np.array([[2.727272727], [0.54545454]])),
        (1, np.array([[0.54545454], [0.54545454]])),
    ],
)
def test_theta_w_one_rep(threshold, expected, single_replicate_meta):
    assert np.isclose(
        expected,
        single_replicate_meta.theta_w(by_population=True, threshold=threshold),
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
    "by_population, threshold, expected",
    [
        (
            True,
            0,
            np.array([[2.72727273, 2.18181818], [0.54545455, 2.72727273]]),
        ),
        (
            True,
            1,
            np.array([[0.54545455, 1.09090909], [0.54545455, 0.54545455]]),
        ),
        (False, 0, np.array([[1.928374656, 2.699724518]])),
        (False, 1, np.array([[0.771349862, 1.157024793]])),
    ],
)
def test_theta_w_two_reps(
    by_population, threshold, expected, double_replicate_meta
):
    out = np.isclose(
        double_replicate_meta.theta_w(
            by_population=by_population, threshold=threshold
        ),
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
    "expected",
    [
        {
            "which": (
                (np.zeros(5, int), np.arange(5)),
                (np.zeros(7, int), np.arange(7)),
            ),
            "num": (5, 7),
        },
        {
            "which": (
                (np.zeros(2, int), np.array([3, 4])),
                (np.zeros(3, int), np.array([0, 3, 4])),
            ),
            "num": (2, 3),
        },
    ],
)
def test_polymorphic_two_reps(expected, double_replicate_meta):
    test = double_replicate_meta.polymorphic()
    assert hash(frozenset(expected)) == hash(frozenset(test))


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
def test_fst_exclude_pop_1(expected, double_replicate_meta_3_popl):
    out = double_replicate_meta_3_popl.fst(
        average_sites=False, h_opts={"bias": False}
    )
    test_results = tuple(
        True if np.isclose(x, y).all() else False
        for x, y in zip(out, expected)
    )
    assert all(test_results)


@pytest.mark.parametrize("expected", [np.array([[5, 4]])])
def test_polymorphic_exclude_pop_1(
    expected, double_replicate_meta_exclude_pop_1
):
    out = double_replicate_meta_exclude_pop_1.polymorphic(output="num")
    assert (expected == out).all()
