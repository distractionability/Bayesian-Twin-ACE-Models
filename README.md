## README: Stan Code for Bayesian-Twin-ACE-Models
This repository contains the Stan code used for the twin ACE (Additive Genetic, Common Environment, Unique Environment) model analyses detailed in **Social change, modernization, and the heritability of educational attainment**. All models are used to estimate the heritability of educational attainment, but they vary along two dimensions:

1.  **Inference strategy** - The baseline model treats educational attainment as an **ordinal categorical outcome"** that reflects a latent educational liability, and estimates the MZ and DZ correlations in this latent liability. These estimated correlations are then used to infer the implied biometric parameters. The "cardinal" model treats the normed education years of each attainment level as a continuous cardinal variable (analogous to height or income). The "constrained" model uses the liability approach, but expresses the model using the biometric parameters directly - ensuring that estimates satisfy theoretical constraints (e.g., A, C and E all non-negative and summing to 1).\
2.  **Temporal resolution** - Each model variant comes in three subvariants that differ in how biometric parameters are allowed to differ across birth cohorts. These models are nested versions of each other: One model has parameters fixed across all birth cohorts. The next model groups birth cohorts into sets of 20 cohorts and lets parameters for these blocks vary around the global mean. The last model groups birth cohorts into sets of 5 cohorts and lets parameters for these blocks differ around the mean of their 20-cohort "parent block".

The logic of the base model in brief:

For each birth cohort, the observed distribution of educational attainment is modelled by inferring a set of cutpoints that split a standard normal distribution of latent educational liability into attainment bins. The probability that a randomly drawn individual would get some educational attainment is the probability that the latent liability falls within that category's attainment bin.

Next, liability for any twin is the sum of a twinpair-specific shared component and a unique individual component, with the relative size of the shared and unique components differing by twin type and potentially (depending on model) birth cohort.

To calculate the likelihood we take each twin in a pair and use an inverse probit to identify the probability that the twin's liability will fall within the binned range reflecting their observed educational attainment - conditional on their shared component.

### 1. General Model Structure

All Stan files follow a standard structure for estimating variance components of an ordinal outcome in a twin-pair design.

| Stan Block | Purpose | Key Components |
|:-----------------------|:-----------------------|:-----------------------|
| `data` | Defines observed data and meta-data. | Zygosity-specific sample sizes (`mz_N`, `dz_N`), number of distinct educational categories (`outcomes_N`), and the number of distinct **cohorts** (`cohort_N`). |
| `transformed data` | Calculates necessary constants. | `cutpoints_N` (always $outcomes\_N - 1$). |
| `parameters` | Variables to be estimated. | **Cutpoints** (via `outcome_shares` simplex) , and the **zygosity-specific shared factor loadings** (`mz_shared_mu`, `dz_shared_mu`) , and the **shared latent twin-pair effects** (`twinpair_std`). |
| `transformed parameters` | Calculates derived parameters used in the likelihood function. | **Cohort-specific cutpoints** (transformed from the simplex using `inv_Phi` and `cumulative_sum`) , the **zygosity-specific shared variance** (MZ_shared, DZ_shared, transformed using `inv_logit`) , and the corresponding **factor coefficients (**$\beta_Z$) and residual standard deviations ($\sigma_\epsilon$) (`MZ_coef`, `DZ_coef`, `MZ_sigma`, `DZ_sigma`)]. |

### 2. Core Model Mechanics (The `model` Block)

The `model` block calculates the log-likelihood by integrating the latent normal variable between two cutpoints, conditioned on the twin-pair's shared latent factor. This method is mathematically equivalent to the ordered probit model.

#### Latent Variable and Correlation

1.  **Shared Latent Factor (**$Z$): Each twin pair has a shared, unscaled latent factor $Z \sim N(0, 1)$ (`twinpair_std`) that captures the combined effect of the $\mathbf{A}$ and $\mathbf{C}$ factors.
2.  **Individual Liability (**$L$): The latent liability for a single twin is $L = \beta_Z Z + \epsilon$, where $\epsilon \sim N(0, \sigma_\epsilon^2)$ and $\sigma_\epsilon^2 = 1 - \beta_Z^2$. The coefficient $\beta_Z$ is either `MZ_coef` or `DZ_coef`.
3.  **Within-Pair Correlation:** The within-pair correlation is simply $\rho = \beta_Z^2$, which corresponds to `MZ_shared` or `DZ_shared`.

#### Log-Likelihood Calculation

The observed ordinal outcome $Y=k$ corresponds to the latent liability $L$ falling between two cutpoints, $\tau_{k-1} < L \le \tau_k$.

-   The probability for a single twin $T_1$ having outcome $Y_1=k$ (conditioned on the shared factor $Z$) is $\Phi(\tau_k | \beta_Z Z, \sigma_\epsilon) - \Phi(\tau_{k-1} | \beta_Z Z, \sigma_\epsilon)$.
-   This is calculated efficiently in Stan using the **Log-Cumulative Distribution Function (LCDF)** for the normal distribution (`normal_lcdf`). The difference of CDFs is calculated on the log scale as $\log(\Phi(\tau_k) - \Phi(\tau_{k-1}))$.

#### Handling Twin Pairs

The model accumulates the log-likelihood (`temp_sum`) for both twins in a pair:

-   **Identical Outcomes (**$Y_1 = Y_2 = k$): If both twins share the same outcome (`twin_same == 1` ), the log-likelihood contribution for the outcome $k$ is added **twice**. This is because, conditioned on the shared factor $Z$, the unique environmental effects ($\epsilon_1$ and $\epsilon_2$) are assumed to be independent.
-   **Different Outcomes (**$Y_1 \ne Y_2$): The log-likelihood contribution for $Y_1$ is added, and then the log-likelihood for $Y_2$ is added separately.

### 3. Derived Variance Components (The `generated quantities` Block)

The `generated quantities` block calculates the variance component estimates (A, C, E) from the posterior samples of the zygosity-specific correlations (`MZ_shared` and `DZ_shared`).

#### Standard ACE Decomposition (Falconer's Formulas)

-   **A (Additive Genetic):** $A = 2 \times (r_{MZ} - r_{DZ})$
-   **C (Common Environment):** $C = r_{MZ} - A$
-   **E (Unique Environment):** $E = 1 - r_{MZ}$

#### Assortative Mating Adjustment

Models that include "assortative" in the filename (e.g., `twin_ace_corr_assortative.stan`) also compute $\mathbf{A}_{AM}$ and $\mathbf{C}_{AM}$ using a different assumed DZ genetic correlation ($\frac{r_A}{2} = 0.68$) to account for non-random mating:

$$\mathbf{A}_{AM} = \frac{1}{1 - 0.68} \times (r_{MZ} - r_{DZ})$$ $$\mathbf{C}_{AM} = r_{MZ} - \mathbf{A}_{AM}$$
