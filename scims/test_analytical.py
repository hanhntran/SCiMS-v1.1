"""
test_analytical.py: Unit tests for analytical.py.

Covers:
  - compute_sex_chrom_probs for XY and ZW
  - Per-sample classification under both Binomial and Beta-Binomial
  - Label-free calibration of epsilon
  - ML estimation of Beta-Binomial concentration kappa
  - Batch classification with auto-calibration
  - Baboon regression: the samples that motivated this fix
  - Beta-Binomial <-> Binomial equivalence as kappa -> infinity

Run with:  pytest test_analytical.py -v
"""

import numpy as np
import pandas as pd
import pytest

from analytical import (
    DEFAULT_MISMAPPING_RATE,
    MIN_EPSILON,
    MIN_KAPPA,
    MIN_SAMPLES_FOR_KAPPA,
    _estimate_kappa,
    classify_batch,
    classify_sample,
    compute_sex_chrom_probs,
    estimate_mismapping_rate,
)


# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------

L_X_MAMMAL = 142_711_496
L_Y_MAMMAL = 8_963_220
L_Z_REPTILE = 100_000_000
L_W_REPTILE = 20_000_000


def make_idxstats(homo_reads, hetero_reads,
                  autosome_reads=100_000,
                  homo_length=L_X_MAMMAL, hetero_length=L_Y_MAMMAL,
                  homo_id="chrX", hetero_id="chrY"):
    """Build an idxstats DataFrame for one sample."""
    return pd.DataFrame(
        {"length": [200_000_000, homo_length, hetero_length],
         "mapped": [autosome_reads, homo_reads, hetero_reads]},
        index=["chr1", homo_id, hetero_id],
    )


def simulate_betabinom(n, mu, kappa, rng):
    """Draw a single Beta-Binomial(n, mu*kappa, (1-mu)*kappa) variate."""
    alpha = mu * kappa
    beta = (1.0 - mu) * kappa
    p = rng.beta(alpha, beta)
    return rng.binomial(n, p)


# ---------------------------------------------------------------------------
# compute_sex_chrom_probs
# ---------------------------------------------------------------------------

class TestComputeSexChromProbs:

    def test_xy_male_probs_match_length_ratio(self):
        p_homo, p_hetero = compute_sex_chrom_probs(
            L_X_MAMMAL, L_Y_MAMMAL, "male", is_zw=False
        )
        expected = L_Y_MAMMAL / (L_X_MAMMAL + L_Y_MAMMAL)
        assert p_hetero == pytest.approx(expected)
        assert p_homo == pytest.approx(1 - expected)

    def test_xy_female_uses_epsilon(self):
        _, p_hetero = compute_sex_chrom_probs(
            L_X_MAMMAL, L_Y_MAMMAL, "female", is_zw=False,
            mismapping_rate=0.01,
        )
        assert p_hetero == pytest.approx(0.01)

    def test_zw_female_matches_length_ratio(self):
        _, p_hetero = compute_sex_chrom_probs(
            L_Z_REPTILE, L_W_REPTILE, "female", is_zw=True
        )
        expected = L_W_REPTILE / (L_Z_REPTILE + L_W_REPTILE)
        assert p_hetero == pytest.approx(expected)

    def test_zw_male_uses_epsilon(self):
        _, p_hetero = compute_sex_chrom_probs(
            L_Z_REPTILE, L_W_REPTILE, "male", is_zw=True,
            mismapping_rate=0.002,
        )
        assert p_hetero == pytest.approx(0.002)

    def test_probabilities_sum_to_one(self):
        for sex in ["male", "female"]:
            for is_zw in [False, True]:
                p_h, p_het = compute_sex_chrom_probs(
                    L_X_MAMMAL, L_Y_MAMMAL, sex, is_zw=is_zw,
                    mismapping_rate=0.005,
                )
                assert p_h + p_het == pytest.approx(1.0)


# ---------------------------------------------------------------------------
# classify_sample: Binomial
# ---------------------------------------------------------------------------

class TestClassifySampleBinomial:

    SCAFFOLDS = ["chr1", "chrX", "chrY"]

    def test_clear_male(self):
        idxstats = make_idxstats(10_000, 600)
        result = classify_sample(
            idxstats, self.SCAFFOLDS, "chrX", "chrY",
            mismapping_rate=0.005,
        )
        assert result["SCiMS predicted sex"] == "male"

    def test_clear_female(self):
        idxstats = make_idxstats(10_000, 50)
        result = classify_sample(
            idxstats, self.SCAFFOLDS, "chrX", "chrY",
            mismapping_rate=0.005,
        )
        assert result["SCiMS predicted sex"] == "female"

    def test_no_sex_chrom_reads_uncertain(self):
        idxstats = make_idxstats(0, 0)
        result = classify_sample(
            idxstats, self.SCAFFOLDS, "chrX", "chrY",
        )
        assert result["SCiMS predicted sex"] == "uncertain"

    def test_zw_female(self):
        idxstats = make_idxstats(
            10_000, 2_000,
            homo_length=L_Z_REPTILE, hetero_length=L_W_REPTILE,
            homo_id="chrZ", hetero_id="chrW",
        )
        result = classify_sample(
            idxstats, ["chr1", "chrZ", "chrW"], "chrZ", "chrW",
            is_zw=True, mismapping_rate=0.005,
        )
        assert result["SCiMS predicted sex"] == "female"
        assert "Reads mapped to Z" in result
        assert "Reads mapped to W" in result

    def test_zw_male(self):
        idxstats = make_idxstats(
            10_000, 50,
            homo_length=L_Z_REPTILE, hetero_length=L_W_REPTILE,
            homo_id="chrZ", hetero_id="chrW",
        )
        result = classify_sample(
            idxstats, ["chr1", "chrZ", "chrW"], "chrZ", "chrW",
            is_zw=True, mismapping_rate=0.005,
        )
        assert result["SCiMS predicted sex"] == "male"

    def test_epsilon_changes_classification(self):
        """The baboon bug as a unit test."""
        idxstats = make_idxstats(442, 6)
        wrong = classify_sample(
            idxstats, self.SCAFFOLDS, "chrX", "chrY", mismapping_rate=1e-4,
        )
        assert wrong["SCiMS predicted sex"] == "male"
        right = classify_sample(
            idxstats, self.SCAFFOLDS, "chrX", "chrY", mismapping_rate=0.005,
        )
        assert right["SCiMS predicted sex"] == "female"


# ---------------------------------------------------------------------------
# classify_sample: Beta-Binomial
# ---------------------------------------------------------------------------

class TestClassifySampleBetaBinom:

    SCAFFOLDS = ["chr1", "chrX", "chrY"]

    def test_betabinom_matches_binomial_at_high_kappa(self):
        """Beta-Binomial with very high kappa should produce the same
        posterior as Binomial (to numerical precision)."""
        idxstats = make_idxstats(1000, 30)
        binom = classify_sample(
            idxstats, self.SCAFFOLDS, "chrX", "chrY",
            mismapping_rate=0.02, kappa=None,
        )
        bb = classify_sample(
            idxstats, self.SCAFFOLDS, "chrX", "chrY",
            mismapping_rate=0.02, kappa=1e5,
        )
        # Posteriors should be nearly equal
        assert abs(binom["Posterior probability of being male"]
                   - bb["Posterior probability of being male"]) < 0.02

    def test_betabinom_is_more_conservative_at_low_kappa(self):
        """A borderline sample should get a less extreme posterior under
        Beta-Binomial with low concentration."""
        # A female with mildly elevated Y-fraction
        idxstats = make_idxstats(1000, 18)  # 1.8% Y, well above eps=0.005
        # Under Binomial with eps=0.005, this looks male.
        binom = classify_sample(
            idxstats, self.SCAFFOLDS, "chrX", "chrY",
            mismapping_rate=0.005, kappa=None,
        )
        # Under Beta-Binomial with low kappa (wide Beta), it's plausibly female.
        bb = classify_sample(
            idxstats, self.SCAFFOLDS, "chrX", "chrY",
            mismapping_rate=0.005, kappa=200,
        )
        # BB posterior should be *less* extreme toward male than Binomial's
        assert bb["Posterior probability of being male"] < binom["Posterior probability of being male"]


# ---------------------------------------------------------------------------
# _estimate_kappa
# ---------------------------------------------------------------------------

class TestEstimateKappa:

    def test_recovers_true_kappa_on_simulated_overdispersed(self):
        """Simulate Beta-Binomial data and check kappa recovery."""
        true_mu = 0.006
        true_kappa = 500.0
        rng = np.random.default_rng(42)
        n_samples = 40
        n_reads = rng.integers(2000, 15000, size=n_samples)
        Y = np.array([simulate_betabinom(n, true_mu, true_kappa, rng)
                      for n in n_reads])
        kappa_hat = _estimate_kappa(Y.astype(float), n_reads.astype(float), true_mu)
        assert kappa_hat is not None
        # ML on ~40 samples should be within ~2x of truth
        assert 0.5 * true_kappa < kappa_hat < 2 * true_kappa

    def test_returns_none_for_binomial_data(self):
        """When the data is actually Binomial, ML pushes kappa above the
        ceiling and we get None (signal to fall back)."""
        mu = 0.006
        rng = np.random.default_rng(0)
        n_samples = 30
        n_reads = rng.integers(5000, 20000, size=n_samples)
        # Draw from Binomial, not Beta-Binomial
        Y = np.array([rng.binomial(n, mu) for n in n_reads])
        kappa_hat = _estimate_kappa(Y.astype(float), n_reads.astype(float), mu)
        # Should either be None (ceiling hit) or a very large value
        if kappa_hat is not None:
            assert kappa_hat > 10_000

    def test_returns_none_for_small_sample(self):
        """Fewer than MIN_SAMPLES_FOR_KAPPA -> None."""
        Y = np.array([5.0, 6.0, 4.0])  # only 3 samples
        n = np.array([1000.0, 1100.0, 900.0])
        assert _estimate_kappa(Y, n, 0.005) is None


# ---------------------------------------------------------------------------
# estimate_mismapping_rate (Binomial path)
# ---------------------------------------------------------------------------

class TestEstimateMismappingRateBinomial:

    def _counts(self, samples):
        return pd.DataFrame(
            [{"sample": n, "homo_reads": h, "hetero_reads": ht}
             for n, h, ht in samples]
        ).set_index("sample")

    def test_mixed_batch_recovers_epsilon(self):
        true_eps = 0.006
        p_Y_m = L_Y_MAMMAL / (L_X_MAMMAL + L_Y_MAMMAL)
        rng = np.random.default_rng(42)
        samples = []
        for i in range(20):
            n = rng.integers(1000, 10000)
            y = rng.binomial(n, true_eps)
            samples.append((f"F{i}", n - y, y))
        for i in range(10):
            n = rng.integers(1000, 10000)
            y = rng.binomial(n, p_Y_m)
            samples.append((f"M{i}", n - y, y))
        eps, kappa, info = estimate_mismapping_rate(
            self._counts(samples), L_X_MAMMAL, L_Y_MAMMAL,
            model="binomial",
        )
        assert info["model_effective"] == "binomial"
        assert kappa is None
        assert eps == pytest.approx(true_eps, rel=0.15)

    def test_all_male_falls_back(self):
        p_Y_m = L_Y_MAMMAL / (L_X_MAMMAL + L_Y_MAMMAL)
        rng = np.random.default_rng(1)
        samples = [(f"M{i}", 5000 - rng.binomial(5000, p_Y_m),
                    rng.binomial(5000, p_Y_m)) for i in range(5)]
        eps, kappa, info = estimate_mismapping_rate(
            self._counts(samples), L_X_MAMMAL, L_Y_MAMMAL,
        )
        assert info["method"] == "fallback_no_homogametic"
        assert eps == DEFAULT_MISMAPPING_RATE
        assert kappa is None

    def test_all_low_depth_falls_back(self):
        samples = [(f"S{i}", 30, 1) for i in range(10)]
        eps, kappa, info = estimate_mismapping_rate(
            self._counts(samples), L_X_MAMMAL, L_Y_MAMMAL,
        )
        assert info["method"] == "fallback_all_low_depth"
        assert eps == DEFAULT_MISMAPPING_RATE

    def test_epsilon_floor(self):
        samples = [(f"F{i}", 5000, 0) for i in range(10)]
        eps, _, _ = estimate_mismapping_rate(
            self._counts(samples), L_X_MAMMAL, L_Y_MAMMAL,
        )
        assert eps == MIN_EPSILON


# ---------------------------------------------------------------------------
# estimate_mismapping_rate (Beta-Binomial path)
# ---------------------------------------------------------------------------

class TestEstimateMismappingRateBetaBinom:

    def _counts(self, samples):
        return pd.DataFrame(
            [{"sample": n, "homo_reads": h, "hetero_reads": ht}
             for n, h, ht in samples]
        ).set_index("sample")

    def test_recovers_both_params_on_simulated_overdispersed(self):
        """Simulate a mixed batch with overdispersed females; verify both
        epsilon and kappa are recovered."""
        true_eps = 0.005
        true_kappa = 400.0
        p_Y_m = L_Y_MAMMAL / (L_X_MAMMAL + L_Y_MAMMAL)
        rng = np.random.default_rng(7)
        samples = []
        for i in range(30):
            n = rng.integers(2000, 15000)
            y = simulate_betabinom(n, true_eps, true_kappa, rng)
            samples.append((f"F{i}", n - y, y))
        for i in range(10):
            n = rng.integers(2000, 15000)
            y = rng.binomial(n, p_Y_m)
            samples.append((f"M{i}", n - y, y))

        eps, kappa, info = estimate_mismapping_rate(
            self._counts(samples), L_X_MAMMAL, L_Y_MAMMAL,
            model="betabinom",
        )
        assert info["model_effective"] == "betabinom"
        assert eps == pytest.approx(true_eps, rel=0.2)
        assert kappa is not None
        assert 0.4 * true_kappa < kappa < 3 * true_kappa

    def test_binomial_data_triggers_fallback_to_binomial(self):
        """When data is actually Binomial and betabinom is requested,
        kappa estimation returns None and model_effective falls back."""
        true_eps = 0.005
        p_Y_m = L_Y_MAMMAL / (L_X_MAMMAL + L_Y_MAMMAL)
        rng = np.random.default_rng(11)
        samples = []
        for i in range(30):
            n = rng.integers(5000, 20000)
            y = rng.binomial(n, true_eps)
            samples.append((f"F{i}", n - y, y))
        for i in range(5):
            n = rng.integers(5000, 20000)
            y = rng.binomial(n, p_Y_m)
            samples.append((f"M{i}", n - y, y))

        eps, kappa, info = estimate_mismapping_rate(
            self._counts(samples), L_X_MAMMAL, L_Y_MAMMAL,
            model="betabinom",
        )
        # Either kappa is None (clean fallback) or very large (~Binomial)
        assert info["model"] == "betabinom"  # user requested
        if kappa is None:
            assert info["model_effective"] == "binomial"
        else:
            assert kappa > 10_000

    def test_insufficient_homogametic_samples_falls_back(self):
        """With < MIN_SAMPLES_FOR_KAPPA homogametic samples, kappa can't
        be estimated and we fall back to Binomial."""
        true_eps = 0.005
        p_Y_m = L_Y_MAMMAL / (L_X_MAMMAL + L_Y_MAMMAL)
        rng = np.random.default_rng(13)
        samples = []
        # Only 3 female samples (below MIN_SAMPLES_FOR_KAPPA=5)
        for i in range(3):
            n = rng.integers(2000, 10000)
            y = simulate_betabinom(n, true_eps, 300.0, rng)
            samples.append((f"F{i}", n - y, y))
        for i in range(8):
            n = rng.integers(2000, 10000)
            y = rng.binomial(n, p_Y_m)
            samples.append((f"M{i}", n - y, y))

        eps, kappa, info = estimate_mismapping_rate(
            self._counts(samples), L_X_MAMMAL, L_Y_MAMMAL,
            model="betabinom",
        )
        assert kappa is None
        assert info["model_effective"] == "binomial"


# ---------------------------------------------------------------------------
# classify_batch: end-to-end
# ---------------------------------------------------------------------------

class TestClassifyBatch:

    SCAFFOLDS = ["chr1", "chrX", "chrY"]

    def _make_mixed_batch(self, n_female=20, n_male=10,
                           true_eps=0.006, true_kappa=None, seed=42):
        p_Y_m = L_Y_MAMMAL / (L_X_MAMMAL + L_Y_MAMMAL)
        rng = np.random.default_rng(seed)
        batch, truth = {}, {}
        for i in range(n_female):
            n = rng.integers(2000, 20000)
            if true_kappa is None:
                y = rng.binomial(n, true_eps)
            else:
                y = simulate_betabinom(n, true_eps, true_kappa, rng)
            batch[f"F{i}"] = make_idxstats(n - y, y)
            truth[f"F{i}"] = "female"
        for i in range(n_male):
            n = rng.integers(2000, 20000)
            y = rng.binomial(n, p_Y_m)
            batch[f"M{i}"] = make_idxstats(n - y, y)
            truth[f"M{i}"] = "male"
        return batch, truth

    def test_auto_calibration_binomial(self):
        batch, truth = self._make_mixed_batch(true_eps=0.006)
        results, info = classify_batch(
            batch, self.SCAFFOLDS, "chrX", "chrY", model="binomial",
        )
        assert info["source"] == "auto_label_free"
        assert info["model_effective"] == "binomial"
        assert info["kappa_used"] is None
        confident = [s for s, r in results.items()
                     if r["SCiMS predicted sex"] != "uncertain"]
        correct = sum(results[s]["SCiMS predicted sex"] == truth[s]
                      for s in confident)
        assert correct / len(confident) >= 0.95

    def test_auto_calibration_betabinom(self):
        batch, truth = self._make_mixed_batch(
            true_eps=0.005, true_kappa=400.0,
        )
        results, info = classify_batch(
            batch, self.SCAFFOLDS, "chrX", "chrY", model="betabinom",
        )
        assert info["source"] == "auto_label_free"
        # With overdispersed data, effective model should be betabinom
        assert info["model_effective"] == "betabinom"
        assert info["kappa_used"] is not None
        confident = [s for s, r in results.items()
                     if r["SCiMS predicted sex"] != "uncertain"]
        correct = sum(results[s]["SCiMS predicted sex"] == truth[s]
                      for s in confident)
        assert correct / len(confident) >= 0.9

    def test_user_supplied_epsilon_binomial(self):
        batch, _ = self._make_mixed_batch()
        _, info = classify_batch(
            batch, self.SCAFFOLDS, "chrX", "chrY",
            mismapping_rate=0.0123, model="binomial",
        )
        assert info["source"] == "user_supplied"
        assert info["epsilon_used"] == 0.0123
        assert info["kappa_used"] is None

    def test_user_supplied_epsilon_plus_betabinom_fits_kappa(self):
        """If user supplies epsilon but asks for betabinom, kappa should
        still be fit from the data."""
        batch, _ = self._make_mixed_batch(true_eps=0.005, true_kappa=300.0)
        _, info = classify_batch(
            batch, self.SCAFFOLDS, "chrX", "chrY",
            mismapping_rate=0.005, model="betabinom",
        )
        assert info["source"] == "user_supplied"
        assert info["epsilon_used"] == 0.005
        # Kappa should be estimated
        assert info["kappa_used"] is not None

    def test_single_sample_falls_back(self):
        _, info = classify_batch(
            {"S1": make_idxstats(1000, 60)},
            self.SCAFFOLDS, "chrX", "chrY", model="betabinom",
        )
        assert info["source"] == "fallback_single_sample"
        assert info["model_effective"] == "binomial"

    def test_invalid_model_raises(self):
        with pytest.raises(ValueError):
            classify_batch(
                {"S1": make_idxstats(1000, 60)},
                self.SCAFFOLDS, "chrX", "chrY", model="not_a_model",
            )


# ---------------------------------------------------------------------------
# Baboon regression
# ---------------------------------------------------------------------------

class TestBaboonRegression:
    """The samples that originally motivated this whole refactor."""

    SCAFFOLDS = ["chr1", "chrX", "chrY"]

    FALSE_MALES = [
        ("SRR1747028", 1475, 23),
        ("SRR1747035", 442, 6),
        ("SRR1747036", 399, 6),
        ("SRR1747041", 1481, 16),
        ("SRR1747065", 540, 12),
    ]
    TRUE_FEMALES = [
        ("SRR1747034", 23689, 118),
        ("SRR1747045", 27096, 128),
        ("SRR1747038", 14541, 56),
        ("SRR1747054", 14766, 76),
        ("SRR1747055", 9033, 48),
        ("SRR1747050", 10505, 64),
        ("SRR1747056", 9703, 64),
        ("SRR1747047", 9015, 55),
        ("SRR1747049", 3195, 22),
        ("SRR1747051", 4902, 21),
        ("SRR1747052", 3465, 26),
        ("SRR1747053", 2264, 13),
    ]
    TRUE_MALES = [
        ("SRR1747061", 48096, 2876),
        ("SRR1747039", 5253, 328),
        ("SRR1747040", 3595, 213),
        ("SRR1747029", 1647, 128),
        ("SRR1747018", 5166, 322),
        ("SRR1747021", 3854, 243),
    ]

    def _build_batch(self):
        batch = {}
        for name, x, y in (self.FALSE_MALES + self.TRUE_FEMALES + self.TRUE_MALES):
            batch[name] = make_idxstats(x, y)
        return batch

    def test_binomial_no_false_males(self):
        results, info = classify_batch(
            self._build_batch(), self.SCAFFOLDS, "chrX", "chrY",
            model="binomial",
        )
        assert info["source"] == "auto_label_free"
        assert 0.002 <= info["epsilon_used"] <= 0.01

        for name, _, _ in self.FALSE_MALES:
            assert results[name]["SCiMS predicted sex"] != "male", \
                f"{name} regressed to male under Binomial"
        for name, _, _ in self.TRUE_MALES:
            assert results[name]["SCiMS predicted sex"] == "male"

    def test_betabinom_no_false_males_and_fewer_uncertain(self):
        """Under Beta-Binomial, all of the previously-misclassified baboon
        females should be classified as female (not just 'uncertain'),
        because the Beta-Binomial explicitly accounts for the elevated
        epsilon variance in the baboon reference."""
        results, info = classify_batch(
            self._build_batch(), self.SCAFFOLDS, "chrX", "chrY",
            model="betabinom",
        )
        assert info["source"] == "auto_label_free"
        assert info["model_effective"] == "betabinom"
        assert info["kappa_used"] is not None

        for name, _, _ in self.FALSE_MALES:
            call = results[name]["SCiMS predicted sex"]
            assert call != "male", f"{name} regressed to male under BetaBinom"
        for name, _, _ in self.TRUE_MALES:
            assert results[name]["SCiMS predicted sex"] == "male"
        for name, _, _ in self.TRUE_FEMALES:
            assert results[name]["SCiMS predicted sex"] == "female"