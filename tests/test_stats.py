import dreamcoat as dc
import numpy as np


def test_std_bias_correction():
    # Test values taken from table in section Bias correction
    # https://en.wikipedia.org/wiki/Unbiased_estimation_of_standard_deviation
    # on 4th May 2024
    sample_size = np.array([2, 3, 4, 5, 6, 7, 8, 9, 10, 100])
    c4 = dc.stats.std_bias_correction(sample_size)
    c4_wikipedia = np.array(
        [
            0.7978845608,  # 2
            0.8862269255,  # 3
            0.9213177319,  # 4
            0.9399856030,  # 5
            0.9515328619,  # 6
            0.9593687891,  # 7
            0.9650304561,  # 8
            0.9693106998,  # 9
            0.9726592741,  # 10
            0.9974779761,  # 100
            # Can't use sample_sizes below here because factorial(k - 1) returns inf
            # 0.9997497811,  # 1000
            # 0.9999749978,  # 10000
        ]
    )
    diff = np.round(c4, decimals=10) - c4_wikipedia
    # We get the final (10th) decimal place slightly wrong for sample_sizes 7 and 9
    assert np.all(np.abs(diff) < 1e-9)


# test_std_bias_correction()
