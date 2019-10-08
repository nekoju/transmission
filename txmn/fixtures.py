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


import msprime as ms
import numpy as np
import pytest

from txmn.classes import Sample
from txmn.classes import MetaSample


@pytest.fixture
def single_replicate(num_samples=4):
    """
    Generate a single replicate of two populations with msprime.
    """
    num_pop = 2
    population_config = [
        ms.PopulationConfiguration(sample_size=num_samples)
        for _ in range(num_pop)
    ]
    migration = np.full((num_pop, num_pop), 10.0)
    np.fill_diagonal(migration, 0.0)
    out = ms.simulate(
        population_configurations=population_config,
        migration_matrix=migration,
        mutation_rate=0.5,
        random_seed=3,
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
def double_replicate(num_samples=4):
    """
    Generate a double replicate of two populations with msprime.
    """
    num_pop = 2
    population_config = [
        ms.PopulationConfiguration(sample_size=num_samples)
        for _ in range(num_pop)
    ]
    migration = np.full((num_pop, num_pop), 10.0)
    np.fill_diagonal(migration, 0.0)
    out = ms.simulate(
        population_configurations=population_config,
        migration_matrix=migration,
        mutation_rate=0.5,
        num_replicates=2,
        random_seed=3,
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
#         [0, 1, 1, 0, 0, 0, 0],
#         [1, 0, 0, 0, 1, 0, 0],
#         [1, 0, 0, 0, 1, 0, 0],
#         [0, 0, 0, 0, 0, 1, 1],
#         [1, 0, 0, 0, 1, 0, 0],
#         [0, 0, 0, 1, 0, 0, 0],
#         [0, 0, 0, 1, 0, 0, 0]])


@pytest.fixture
def double_replicate_meta_exclude_pop_1(
    num_samples=4,
    populations=np.repeat([0, 1], 4),
    keep_populations=np.array([0]),
):
    """
    Generate a double replicate of two populations with msprime.
    """
    num_pop = len(set(populations))
    population_config = [
        ms.PopulationConfiguration(sample_size=num_samples)
        for _ in range(num_pop)
    ]
    migration = np.full((num_pop, num_pop), 10.0)
    np.fill_diagonal(migration, 0.0)
    out = ms.simulate(
        population_configurations=population_config,
        migration_matrix=migration,
        mutation_rate=0.5,
        num_replicates=2,
        random_seed=3,
    )
    return MetaSample(out, populations, keep_populations)


@pytest.fixture
def single_replicate_meta(num_samples=4, populations=np.repeat([0, 1], 4)):
    """
    Generate a single replicate of two populations with msprime.
    """
    num_pop = len(set(populations))
    population_config = [
        ms.PopulationConfiguration(sample_size=num_samples)
        for _ in range(num_pop)
    ]
    migration = np.full((num_pop, num_pop), 10.0)
    np.fill_diagonal(migration, 0.0)
    out = ms.simulate(
        population_configurations=population_config,
        migration_matrix=migration,
        mutation_rate=0.5,
        random_seed=3,
    )
    return MetaSample(out, populations)


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
def double_replicate_meta(num_samples=4, populations=np.repeat([0, 1], 4)):
    """
    Generate a double replicate of two populations with msprime.
    """
    num_pop = len(set(populations))
    population_config = [
        ms.PopulationConfiguration(sample_size=num_samples)
        for _ in range(num_pop)
    ]
    migration = np.full((num_pop, num_pop), 10.0)
    np.fill_diagonal(migration, 0.0)
    out = ms.simulate(
        population_configurations=population_config,
        migration_matrix=migration,
        mutation_rate=0.5,
        num_replicates=2,
        random_seed=3,
    )
    return MetaSample(out, populations)


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
#         [0, 1, 1, 0, 0, 0, 0],
#         [1, 0, 0, 0, 1, 0, 0],
#         [1, 0, 0, 0, 1, 0, 0],
#         [0, 0, 0, 0, 0, 1, 1],
#         [1, 0, 0, 0, 1, 0, 0],
#         [0, 0, 0, 1, 0, 0, 0],
#         [0, 0, 0, 1, 0, 0, 0]])


@pytest.fixture
def double_replicate_meta_3_popl(
    num_samples=4,
    populations=np.repeat([0, 1, 2], 4),
    keep_populations=np.array([0, 1]),
):
    """
    Generate a triple replicate of two populations with msprime.
    """
    num_pop = len(set(populations))
    population_config = [
        ms.PopulationConfiguration(sample_size=num_samples)
        for _ in range(num_pop)
    ]
    migration = np.full((num_pop, num_pop), 10.0)
    np.fill_diagonal(migration, 0.0)
    out = ms.simulate(
        population_configurations=population_config,
        migration_matrix=migration,
        mutation_rate=0.1,
        num_replicates=2,
        random_seed=3,
    )
    out = MetaSample(out, populations, keep_populations)
    return out


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
#         [0, 1, 1, 0, 0, 0, 0],
#         [1, 0, 0, 0, 1, 0, 0],
#         [1, 0, 0, 0, 1, 0, 0],
#         [0, 0, 0, 0, 0, 1, 1],
#         [1, 0, 0, 0, 1, 0, 0],
#         [0, 0, 0, 1, 0, 0, 0],
#         [0, 0, 0, 1, 0, 0, 0]])
