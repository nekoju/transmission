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

import msprime as ms
import numpy as np
import pytest
from transmission.transmission import Sample


@pytest.fixture
def single_replicate():
    """
    Generate a single replicate of two populations with msprime.
    """
    population_config = [ms.PopulationConfiguration(sample_size=4)
                         for _ in (0, 1)]
    migration = np.full((2, 2), 10.)
    np.fill_diagonal(migration, 0.)
    out = ms.simulate(
        population_configurations=population_config,
        migration_matrix=migration,
        mutation_rate=0.5,
        random_seed=3
        )
    return Sample(out)

# single_replicate genotypes
# array([[1, 1, 0, 1, 0],
#        [0, 0, 0, 0, 0],
#        [0, 0, 0, 0, 1],
#        [0, 0, 1, 1, 0],
#        [0, 0, 0, 0, 0],
#        [0, 0, 0, 0, 0],
#        [0, 0, 0, 0, 1],
#        [0, 0, 0, 0, 1]])


@pytest.fixture
def double_replicate():
    """
    Generate a double replicate of two populations with msprime.
    """
    population_config = [ms.PopulationConfiguration(sample_size=4)
                         for _ in (0, 1)]
    migration = np.full((2, 2), 10.)
    np.fill_diagonal(migration, 0.)
    out = ms.simulate(
        population_configurations=population_config,
        migration_matrix=migration,
        mutation_rate=0.5,
        num_replicates=2,
        random_seed=3
        )
    return Sample(out)

# double_replicate genotypes
# [array([[1, 1, 0, 1, 0],
#         [0, 0, 0, 0, 0],
#         [0, 0, 0, 0, 1],
#         [0, 0, 1, 1, 0],
#         [0, 0, 0, 0, 0],
#         [0, 0, 0, 0, 0],
#         [0, 0, 0, 0, 1],
#         [0, 0, 0, 0, 1]]),
#  array([[1, 0, 0, 0, 1, 0, 0],
#         [1, 0, 0, 0, 1, 0, 0],
#         [0, 0, 0, 1, 0, 0, 0],
#         [0, 1, 1, 0, 0, 0, 0],
#         [0, 0, 0, 1, 0, 0, 0],
#         [0, 1, 1, 0, 0, 0, 0],
#         [1, 0, 0, 0, 1, 0, 0],
#         [0, 0, 0, 0, 0, 1, 1]])]


# single replicate
@pytest.mark.parametrize("bias, expected", [
    (False, np.array([0.21875, 0.21875, 0.21875, 0.375, 0.46875])),
    (True,
     8. / (8. - 1.) * np.array([0.21875, 0.21875, 0.21875, 0.375, 0.46875])
     )
    ])
def test_h_no_average_one_rep(bias, expected, single_replicate):
    assert np.isclose(single_replicate.h(bias=bias), expected).all()


@pytest.mark.parametrize("bias, expected", [
    (False, 0.300),
    (True, 8. / (8. - 1.) * 0.300)
    ])
def test_h_average_one_rep(bias, expected, single_replicate):
    assert np.isclose(
        single_replicate.h(bias=bias, average=True),
        expected
        )


# double replicate
@pytest.mark.parametrize("bias, expected", [
    (False, [
        np.array([0.21875, 0.21875, 0.21875, 0.375, 0.46875]),
        np.array([0.46875, 0.375, 0.375, 0.375, 0.46875, 0.21875, 0.21875])
        ]),
    (True, [
        (8. / (8. - 1.) *
         np.array([0.21875, 0.21875, 0.21875, 0.375, 0.46875])),
        (8. / (8. - 1) *
         np.array([0.46875, 0.375, 0.375, 0.375, 0.46875, 0.21875, 0.21875]))
        ]),
    ])
def test_h_no_average_two_reps(bias, expected, double_replicate):
    heterozygosity = double_replicate.h(bias=bias)
    out = [np.isclose(expected, heterozygosity[i]).all()
           for i, expected in enumerate(expected)]
    assert all(out)


@pytest.mark.parametrize("bias, expected", [
    (False, [0.300, 0.357142857]),
    (True, [(8. / (8. - 1.) * 0.300), (8. / (8. - 1) * 0.357142857)]),
    ])
def test_h_average_two_reps(bias, expected, double_replicate):
    heterozygosity = double_replicate.h(bias=bias, average=True)
    out = [np.isclose(expected, heterozygosity[i]).all()
           for i, expected in enumerate(expected)]
    assert all(out)
