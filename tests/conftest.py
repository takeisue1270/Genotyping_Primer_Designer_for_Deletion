"""Shared fixtures.

Every test runs against a synthetic sequence built from a fixed seed, so the
expected coordinates are exact and nothing here touches the network.
"""

from __future__ import annotations

import random

import pytest

from genotyping_primers.config import DesignConfig
from genotyping_primers.models import Target


def synthetic_sequence(length: int, seed: int, gc: float = 0.48) -> str:
    """A reproducible sequence with a realistic base composition."""
    rng = random.Random(seed)
    bases = []
    for _ in range(length):
        roll = rng.random()
        if roll < gc / 2:
            bases.append("G")
        elif roll < gc:
            bases.append("C")
        elif roll < gc + (1 - gc) / 2:
            bases.append("A")
        else:
            bases.append("T")
    return "".join(bases)


@pytest.fixture
def make_sequence():
    """The generator itself, for tests that want their own sequence."""
    return synthetic_sequence


@pytest.fixture
def sequence() -> str:
    return synthetic_sequence(8_000, seed=7)


@pytest.fixture
def target(sequence: str) -> Target:
    return Target(name="synthetic", sequence=sequence, source="fixture")


@pytest.fixture
def plasmid() -> Target:
    return Target(
        name="synthetic_plasmid",
        sequence=synthetic_sequence(4_000, seed=11),
        circular=True,
        source="fixture",
    )


@pytest.fixture
def genomic_target(sequence: str) -> Target:
    """The same sequence, positioned on a chromosome like an Ensembl fetch."""
    return Target(
        name="region",
        sequence=sequence,
        source="Ensembl test",
        species="mus_musculus",
        assembly="GRCm39",
        region="11",
        region_start=69_000_001,
    )


@pytest.fixture
def config() -> DesignConfig:
    """Defaults, with the specificity check kept inside the fixture sequence."""
    return DesignConfig(specificity="local", alternatives=2)
