import numpy as np

from Bio.Seq import reverse_complement, Seq
from itertools import product
from typing import Dict, List, Tuple


def get_nucleotides() -> List[str]:
    return ["A", "T", "C", "G"]


def should_kmer_be_reversed(kmer: str, rev_kmer: str):
    return kmer <= rev_kmer


def create_kmer_combinations(k: int) -> List[str]:
    nucleotides: List[str] = get_nucleotides()
    nucleotides_permutations: List[Tuple[str]] = product(nucleotides, repeat=k)
    kmers = ["".join(np) for np in nucleotides_permutations]
    return kmers


def create_kmer_features(k: int) -> List[str]:
    kmer_combinations: List[str] = create_kmer_combinations(k=4)
    kmer_features: List[str] = [
        kc for kc in kmer_combinations
        if should_kmer_be_reversed(kc, reverse_complement(kc))]
    return kmer_features


def extract_features(seq: Seq, k: int, step: int) -> Dict[str, int]:
    kmer_features: List[str] = create_kmer_features(k=k)
    init_vals: List[int] = [0] * len(kmer_features)
    feature_dict: Dict[str, int] = dict(zip(kmer_features, init_vals))

    seq_len: int = len(seq)

    for ind in range(0, (k - 1) * (seq_len // k), step):
        cur_kmer = seq[ind:(ind+k)]
        rev_cur_kmer = reverse_complement(cur_kmer)
        kmer_to_count = cur_kmer if should_kmer_be_reversed(cur_kmer, rev_cur_kmer) \
            else rev_cur_kmer
        feature_dict[kmer_to_count] += 1

    for kmer in feature_dict:
        feature_dict[kmer] /= seq_len

    return feature_dict


class KmerVectorizer():
    def __init__(self, k: int, step: int):
        self._k = k
        self._step = step

    def _extract_features_single(self, seq: str):
        features: Dict[str, int] = extract_features(
            seq, k=self._k, step=self._k)

        return list(features.values())

    def transform(self, seq_list: List[str]):
        features: List[Dict[str, int]] = [
            self._extract_features_single(seq) for seq in seq_list]

        return np.array(features)
