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


@pytest.mark.parametrize("bias, replace, expected", [
    (False, True, np.array([0.21875, 0.21875, 0.21875, 0.375, 0.46875])),
    (True, True,
     8. / (8. - 1.) * np.array([0.21875, 0.21875, 0.21875, 0.375, 0.46875])
     ),
    (False, False, np.array([0.25, 0.25, 0.25, 0.428571429, 0.535714286])),
    (True, False,
     8. / (8. - 1.) * np.array([0.25, 0.25, 0.25, 0.428571429, 0.535714286])
     ),
    ])
def test_h_no_average_one_rep(bias, replace, expected, single_replicate):
    assert np.isclose(single_replicate.h(bias=bias, replace=replace),
                      expected).all()
