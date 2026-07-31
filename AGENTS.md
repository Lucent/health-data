# Working instructions for this repository

Read these files before extending the analysis:

1. `new_findings/README.md` — current audit status and report index.
2. `new_findings/01_claim_registry.md` — what is supported, exploratory, contradicted, or quarantined.
3. `new_findings/10_current_action_plan.md` — current practical priorities and safety boundaries.
4. `new_findings/11_P4b_forward_model_audit.md` — valid uses and limits of the opt-in forward model.
5. `new_findings/15_long_run_trajectory_audit.md` and `new_findings/16_calendar_and_extreme_clustering_audit.md` — direct regime shifts, missingness, and short-run continuation versus failed long-horizon warning claims.
6. `new_findings/17_energy_density_and_nutrient_leverage_audit.md` — why retrospective whole-diet energy density is unavailable and where current meal leverage remains.
7. `CLAUDE.md` — data layout, commands, and known pitfalls.
8. `BACKGROUND.md`, `FINDINGS.md`, and `RECOMMENDATIONS.md` — legacy hypotheses and analyses. Treat their causal and medical claims as hypotheses unless the claim registry says otherwise.

## Evidence boundaries

- Preserve the distinction between raw observations, deterministic transformations, imputations, latent-state estimates, and external evidence.
- The uninterrupted daily intake ledger is a major strength. It does not make model-derived fat mass, TDEE, set point, drug sensitivity, or tachyphylaxis direct observations.
- Do not turn an in-sample association, parameter sweep, or retrospectively smoothed state into a prediction claim.
- Medication schedules and supplement doses require current authoritative evidence and clinician/pharmacist review. Never operationalize click-count dosing from this repository.
- Treat `intake/diet_epochs.csv` as overlapping annotations, not a complete partition.
- Do not treat unweighed foods as negative for the Fazzino sodium clusters. The current exact-mass coverage is too sparse for a complete hyper-palatable-food index; `analysis/EJ_hyperpalatable_food_audit.py` uses the fat+sugar cluster as its primary exposure.
- Do not call exact-mass-subset kcal/g whole-diet energy density. Exact mass covers only 17.0% of solid calories overall and 10.8% in the final 90 days, with nonrandom missingness.
- Treat long weigh-in gaps as informative missingness, not evidence that weighing itself causes weight control. Do not interpolate a causal path across sparse intervals.

## Analysis standard

- Use time-ordered, locked validation for predictive claims. Fit preprocessing, feature selection, and hyperparameters on the training period only.
- Export observed/imputed/source/device flags and uncertainty alongside derived estimates.
- Do not set random states unless the user explicitly asks. Prefer deterministic methods when possible.
- New Python scripts should take no command-line arguments and keep configuration variables near the top after imports.
- Use `uv` for Python environments and dependencies. A basic audit runs with:

  `uv run --with pandas python analysis/EA_repository_audit.py`

- Rebuild the observation-first cardiometabolic/nutrition snapshot with:

  `uv run --with pandas python analysis/EH_cardiometabolic_snapshot.py`

- Rebuild the exact-calendar sleep/activity re-audit with:

  `uv run --with pandas --with numpy python analysis/EI_sleep_activity_reaudit.py`

- Rebuild the direct long-run, calendar-continuation, and nutrient-leverage audits with:

  `uv run --with numpy --with pandas python analysis/EL_long_run_trajectory_audit.py`

  `uv run --with numpy --with pandas python analysis/EM_calendar_extreme_clustering.py`

  `uv run --with numpy --with pandas python analysis/EN_energy_density_nutrient_leverage.py`

- P4b is opt-in and does not replace P4 or downstream consumers. Rebuild and test it with:

  `uv run --with numpy --with pandas python analysis/P4b_forward_state_space.py`

  `uv run --with numpy --with pandas python analysis/test_P4b_forward_state_space.py`

- The 2025–2026 P4b interval has already been opened. If P4b is tuned in response, label that interval calibration and reserve genuinely new data for validation.
- Put material new conclusions in `new_findings/` and update `new_findings/01_claim_registry.md` when a claim changes status.
- Keep legacy results reproducible, but add a warning where a result has been invalidated rather than silently rewriting history.
