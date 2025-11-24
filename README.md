# Bidding Under Uncertain Demand — Project Guide

> The code implements the equilibrium model in **Luo & Takahashi (RAND Journal of Economics, 2025)**, “Bidding for Contracts under Uncertain Demand: Skewed Bidding and Risk Sharing.” Key objects used below—inner portfolio choice, pseudo‑costs, score ODE, and entry fixed point—correspond directly to their Equations (4), (5), (7), (8) for UP, and (10)–(11) for FP, plus the observables scaling in (18). 

------

## 1) What this repository does

- **Simulates UP and FP auctions** with risk‑averse (CARA) bidders and ex‑post quantity risk.
- **Builds equilibrium bidding strategies**:
  - **UP**: two‑loop structure—inner portfolio of itemized bids and an outer score rule that solves an ODE. 
  - **FP**: risk is borne by the contractor; the outer problem is a standard **first‑price** auction in pseudo‑costs. 
- **Endogenizes entry** via a fixed point in entry probability using the win‑probability–weighted expected utility. 
- **Runs counterfactuals**:
  - UP→FP switches (with or without re‑optimizing entry).
  - **Fixed‑entry** comparisons (hold *n* fixed across formats).
- **Outputs  summaries** (auction‑level and bidder‑level; pooled and by risk state) with “homogenized” measures consistent with the paper’s observable rescaling . 

------

## 2) Repository layout

```
.
├── config/                         # Hyper-parameters & switches (sim sizes, solver tolerances, scaling)
│   └── config_up_1d.m
├── utils/                          # General utilities
│   ├── build_Ff_splines.m          # KDE/spline law of pseudo-costs (UP)
│   ├── build_Ff_fp_normal.m        # Normal law of pseudo-costs (FP)
│   ├── expected_utility_integral_weighted.m
│   ├── sample_c_via_sobol.m        # Low-variance draws for pseudo-costs (UP)
│   └── solve_s_via_ode.m           # ODE solver wrapper for UP score rule
├── strategy/                       # Equilibrium strategies
│   ├── build_s_handle.m            # UP: s(c;n) via Eq. (7) (outer loop) 
│   └── build_bfp_handle.m          # FP: b(c;n) 
├── sim/                            # Auction simulators
│   ├── simulate_one_auction_full.m     # UP
│   ├── simulate_one_auction_full_fp.m  # FP
│   ├── compute_entry_probability_bisect.m     # UP entry (Eq. (8))
│   └── compute_entry_probability_bisect_fp.m  # FP entry (Eq. (8) analogue)
├── stats/                          # Summaries
│   └── compute_and_save_summaries.m
├── scripts/                        # Reproducible runners
│   ├── run_up_1d.m
│   ├── run_fp_1d.m
│   ├── run_cf_main.m               # UP vs FP (endogenous entry)
│   └── run_cf_fixed_entry.m        # UP vs FP (fixed n)
└── output/                         # CSVs, MAT files, and logs
```



------

## 3) Typical workflow

1. **Choose a runner** in `scripts/`:
   - `run_up_1d.m` (UP only), `run_fp_1d.m` (FP only), or counterfactual `run_cf_main.m` / `run_cf_fixed_entry.m`.
2. **Load config** (`config/*.m`), which sets:
   - `cfg.sim.*` (number of auctions, NN, risk schedule or probabilities),
   - `cfg.model.*` (baseline θ0,θ1,α,σL,σH\theta_0,\theta_1,\alpha,\sigma_L,\sigma_H),
   - `cfg.scale.*` (activate Eq. (18) rescaling), and solver tolerances (`ode`, `integral`).
3. **Build the pseudo‑cost law** for each risk state:
   - **UP**: draw pseudo‑costs via Eq. (5) and smooth (F,fF,f) using KDE splines.
   - **FP**: form cfc_f analytically under CARA+Normal; use a **Normal** F,fF,f for stability.
4. **Precompute entry** per risk state by solving the fixed point δ∗\delta^* via `compute_entry_probability_bisect*.m` (optional; otherwise computed per‑auction on the fly).
6. **Simulate each auction**:
   - Draw entry outcomes with δ∗\delta^*, construct the entrant pool nn.
   - Draw bidder types, compute pseudo‑costs, evaluate s(c;n)s(c;n) or b(c;n)b(c;n), determine winner and payment rule:
     - **UP**: payment = (lump‑sum price) + (non‑lump‑sum price) × realized quantity; **overrun** may be nonzero.
     - **FP**: payment = bid; **overrun = 0** by construction.
7. **Write outputs**:
   - `auctions.csv` (auction‑level), `bids.csv` (bidder‑level panel), `dgp_results.mat`.
   - `compute_and_save_summaries.m` produces pooled and by‑risk tables, plus homogenized versions (dividing by the project scale ss per Eq. (18)). 

------

## 4) Key modules

### A) Strategy builders

- `build_s_handle.m` (UP): wraps **Eq. (7)**. Inputs: (F,f,supp,n,α)(F,f,\text{supp},n,\alpha). Returns a PCHIP‑interpolated handle c↦s(c;n)c\mapsto s(c;n), enforcing s(cˉ)=cˉs(\bar c)=\bar c. 
- `build_bfp_handle.m` (FP): computes W(c)=[1−F(c)]n−1W(c)=[1-F(c)]^{n-1}, builds its cumulative ∫W\int W via `cumtrapz`, and returns b(c)=c+∫ccˉWW(c)b(c)=c+\frac{\int_c^{\bar c}W}{W(c)} with underflow guards.

### B) Entry solvers

- `compute_entry_probability_bisect.m` (UP) and `compute_entry_probability_bisect_fp.m` (FP) implement the δ\delta fixed point of **Eq. (8)** with **log‑domain** weights over nn to avoid overflow. Per‑nn expectations call:

  - `expected_utility_integral_weighted.m`, which integrates

    Un  =  ∫suppu(payoff(c;n))  [1−F(c)]n−1⏟win prob  f(c) dc.  U_n \;=\; \int_{\text{supp}} u\big(\text{payoff}(c;n)\big)\;\underbrace{[1-F(c)]^{n-1}}_{\text{win prob}}\;f(c)\,dc.

    **Important:** the integrand must use the **raw CARA utility** u(x)=−exp⁡(−αx)u(x)=-\exp(-\alpha x) so that U<0U<0 and CE=−1αlog⁡(−U)CE=-\frac{1}{\alpha}\log(-U) is valid. See **Pitfall P1**.

### C) Simulators

- `simulate_one_auction_full*.m` unify:
  - Entry draw →\to get nn and strategy handles,
  - Draw bidder types →\to pseudo‑costs,
  - Evaluate ss or bb, rank, pick winner,
  - Compute payments:
    - **UP:** payment incorporates realized adjustments; aligns with the paper’s Eq. (2)–(3). 
    - **FP:** payment equals bid (Eq. (9)); no ex‑post adjustments on contracted items. 

------

## 5) Outputs you should expect (and why)

- **Under FP**, pseudo‑costs (and bids) **increase with risk** through the common risk premium, so mean payments rise in σ\sigma. 
- **Under UP**, more risk **dampens skewing**, homogenizing pseudo‑costs via the mean‑variance tradeoff in the inner portfolio (Eq. (4) & (5)). In the numerical illustration of the paper, expected UP payments fall as risk rises while FP payments rise. The summaries mirror these channels; aggregated patterns can be compared to Figure 1 and Figure 2 qualitatively. 

------

## 6) Configuration 

- **Project scale (Eq. (18))**: turn on `cfg.scale.useProjectScale=true` to obtain effective parameters (θ0s,θ1s,α/s)(\theta_0 s,\theta_1 s,\alpha/s). Homogenized outputs then divide by ss, aligning cross‑project comparability as in the paper’s econometric specification. 
- **Risk states** {L,H}\{L,H\}: the code follows the paper’s two‑state mixture intuition (Section 7), letting you schedule states or draw them with a probability PP. 
- **Solver tolerances**: `cfg.solver.integral.RelTol/AbsTol`, `cfg.solver.ode` control numerical accuracy of (7) and the utility integrals.

------

##  References 

- Luo, Y., & Takahashi, H. (2025). *Bidding for Contracts under Uncertain Demand: Skewed Bidding and Risk Sharing*. **RAND Journal of Economics**, 56(3), 325–343. (Core model, Eqs. (2)–(11), scoring ODE (7), entry condition (8), rescaling (18), and comparative statics, including **Figures 1–2**.) 

------

Portions of this project's annotations and documentation were prepared with assistance from ChatGPT 5.0.
