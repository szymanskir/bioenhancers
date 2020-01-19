#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `bioenhancers` package."""
import pytest
from typing import Dict, List

from bioenhancers.features import create_kmer_features, extract_features


def test_create_kmer_features():
    kmer_features: List[str] = create_kmer_features(4)
    assert len(kmer_features) == 136
    assert "AAAA" in kmer_features
    assert "TTTT" not in kmer_features


def test_extract_features():
    test_seq: str = "AAAATTTT"

    result: Dict[str, int] = extract_features(test_seq, k=4, step=4)
    assert result["AAAA"] == 2 / 8
