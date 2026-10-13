import numpy as np

from scipy.stats import lognorm

  
  

 ============================================================

# TCR koff distribution

 ============================================================

  

KOFF_VALUES = np.array([

    0.057,

    0.090,

    0.150,

    0.063,

    0.360,

    0.440,

    0.020,

    0.068,

    0.039,

    0.146,

    0.200,

    0.600,

    0.057,

    0.360,

    5.000,

    0.064,

    0.140,

], dtype=float)

  
  

def fit_koff_distribution(koff_values=KOFF_VALUES, shift_factor=0.3, sigma_factor=6.0):

    """Fit a log-normal distribution to empirical koff values.

    shift_factor > 1

    → moves the entire distribution right

    shift_factor < 1

    → moves it left

    sigma_factor > 1

    → increases heterogeneity and right tail

    sigma_factor < 1

    → reduces heterogeneity and right tail

    """

  

    shape, loc, scale = lognorm.fit(koff_values, floc=0)

  

    mu_log = np.log(scale)

    sigma_log = shape

  

    mu_log += np.log(shift_factor)

    sigma_log *= sigma_factor

  

    print("\nkoff log-normal parameters")

    print("-" * 40)

    print(f"mu_log     = {mu_log:.6f}")

    print(f"sigma_log  = {sigma_log:.6f}")

  

    return mu_log, sigma_log

  
  

def generate_koff_values(n_cells, rng,

                         koff_values=KOFF_VALUES,

                         shift_factor: float = 0.3,

                         sigma_factor: float = 6.0

                         ):

    """Generate one baseline koff value for each thymocyte."""

    mu_log, sigma_log = fit_koff_distribution(

        koff_values,

        shift_factor=shift_factor,

        sigma_factor=sigma_factor

        )

  

    return rng.lognormal(

        mean=mu_log,

        sigma=sigma_log,

        size=n_cells

    )

  
  

# ============================================================

# ZAP70 / proximal signalling-capacity distribution

# ============================================================

  

def generate_zap70_values(

    n_cells,

    rng,

    mean=1.0,

    cv=0.30

):

    """

    Generate relative ZAP70-dependent proximal signalling capacity.

  

    Mean is normalized to 1.0.

    CV controls cell-to-cell heterogeneity.

    """

  

    sigma_log = np.sqrt(np.log(1.0 + cv**2))

    mu_log = np.log(mean) - 0.5 * sigma_log**2

  

    return rng.lognormal(

        mean=mu_log,

        sigma=sigma_log,

        size=n_cells

    )

  

def summarize_distribution(values, name="Distribution"):

    """Print summary statistics for a simulated population."""

  

    values = np.asarray(values)

  

    print(f"\n{name}")

    print("-" * 40)

    print(f"N          = {len(values)}")

    print(f"Mean       = {np.mean(values):.6f}")

    print(f"Median     = {np.median(values):.6f}")

    print(f"Std        = {np.std(values, ddof=1):.6f}")

    print(f"Min        = {np.min(values):.6f}")

    print(f"Max        = {np.max(values):.6f}")

    print(f"5th pct    = {np.percentile(values, 5):.6f}")

    print(f"25th pct   = {np.percentile(values, 25):.6f}")

    print(f"75th pct   = {np.percentile(values, 75):.6f}")

    print(f"95th pct   = {np.percentile(values, 95):.6f}")

  

if __name__ == "__main__":

    rng = np.random.default_rng(42)

  

    N_CELLS = 1000

  

    koff_population = generate_koff_values(

            n_cells=N_CELLS,

            rng=rng

    )

  

    zap70_population = generate_zap70_values(

            n_cells=N_CELLS,

            rng=rng,

            mean=1.0,

            cv=0.30

    )

  

    summarize_distribution(

            koff_population,

            name="Simulated koff_base (s^-1)"

    )

  

    summarize_distribution(

            zap70_population,

            name="Simulated ZAP70 relative capacity"

    )