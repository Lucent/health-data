# Current menu leverage: protein without more calories

Date: 2026-07-30
Status: direct 90-day food-ledger analysis; intended for practical menu design

## Bottom line

The final 90 days do contain low-hanging fruit—not because the subject is secretly eating huge unlogged meals, but because a large share of logged calories is nutritionally inefficient for the current treatment phase.

- Mean intake: 2,004 kcal/day
- Mean protein: 54.6 g/day
- Recurring foods below 3 g protein per 100 kcal: 914 kcal/day, or 45.6% of all calories
- The top 10 of those sources alone: 428 kcal/day
- Breakfast: 364 kcal and 6.4 g protein on average
- Snacks: 316 kcal and 4.1 g protein on logged snack days

This creates a calorie-neutral path from roughly 55 g to roughly 80 g protein: substitute, do not simply add.

## Reproducer

Command:

```bash
uv run --with numpy --with pandas python analysis/EF_current_protein_menu.py
```

Artifacts: `analysis/EF_current_food_sources.csv`, `analysis/EF_current_protein_by_meal.csv`, `analysis/EF_current_swap_candidates.csv`, `analysis/EF_current_menu_summary.csv`

The window is 2025-12-01 through 2026-02-28. Food strings are conservatively normalized to combine portion variants.

## Where the calories are

The largest recurring sources below 3 g protein per 100 kcal were:

| Familiar food cluster | Days | 90-day kcal | Protein density |
|---|---:|---:|---:|
| Homemade brown-sugar cinnamon pop-tarts | 16 | 5,397 | 0.9 g/100 kcal |
| Fritos chili-cheese package | 23 | 5,280 | 1.3 |
| A Taste of Thai peanut noodles | 7 | 4,410 | 1.9 |
| Baklava | 21 | 4,320 | 1.3 |
| Praline | 13 | 3,403 | 0.7 |
| Liège waffle | 12 | 3,370 | 1.7 |
| Thumbprint cookies | 12 | 3,365 | 1.0 |
| Second Fritos chili-cheese entry | 16 | 3,200 | 1.2 |
| Eggless waffles | 10 | 2,970 | 2.2 |
| Momofuku sweet-and-spicy noodles | 8 | 2,800 | 2.9 |

Other recurring entries include peppermint bark, heavy cream, chocolate bites, cinnamon rolls, corn nuts, butter cookies, chips, and popcorn.

This is not a verdict that any one item is forbidden. It shows that “no fast food” and “everything weighed” do not guarantee a protein- or nutrient-dense diet. Many precisely logged, delicious foods can still occupy calories needed for protein, fiber, and micronutrients when appetite is pharmacologically reduced.

## Familiar protein sources already in the ledger

The best recent protein densities among foods eaten on at least four days were:

| Food | Days | Protein density |
|---|---:|---:|
| Nasoya plant-based chick'n bites | 4 | 12.9 g/100 kcal |
| Basic seitan | 4 | 11.1 |
| Gardein plant-based chick'n filets | 4 | 7.1 |

Extra-firm tofu appeared on three days at 10.0 g/100 kcal, and vegetarian chili appeared on 21 days at 6.6 g/100 kcal. The chili was the largest single protein source in the window: 397 g total, about 19 g per logged serving day.

These are particularly useful because they are already familiar. That preserves the menu-stability signal rather than solving protein by creating a constantly changing “health food” project.

## The calorie-neutral substitution

The math is unusually favorable:

- 300 kcal from a 1.5 g/100-kcal source supplies about 4.5 g protein.
- 300 kcal from a 10 g/100-kcal source supplies about 30 g.
- The substitution adds roughly 25.5 g protein without changing calories.
- Applied to the last-90 mean, that moves 54.6 g to about 80 g/day.

The exact swap need not be 300 kcal or happen all at once. A gradual 150-kcal swap at breakfast plus 150 kcal at snacks is easier to tolerate and directly targets the two lowest-protein eating occasions.

Possible familiar templates:

- seitan or tofu plus a preferred breakfast carbohydrate;
- a Nasoya or Gardein portion replacing part of chips/sweets/noodles;
- vegetarian chili as a default lunch or supper;
- dairy protein such as Greek yogurt/cottage cheese if desired and tolerated;
- a tolerable protein drink when early satiety makes solid food difficult.

The choice should respect vegetarian preferences, GI effects, kidney/medical context, enjoyment, and budget. Protein powders are optional tools, not a requirement.

## A two-stage plan

### Stage 1: improve composition at the current calorie level

For two to four weeks:

1. Keep the approximate 2,000-kcal level rather than cutting harder.
2. Replace 250–350 kcal/day of recurring low-protein foods with familiar protein-dense foods.
3. Aim for a clinician/dietitian-reviewed protein range, plausibly 80–100 g/day.
4. Raise fiber gradually from the current 14 g/day using tolerable legumes, vegetables, fruit, whole grains, nuts, and seeds.
5. Log GI symptoms, hunger, training performance, and whether the menu remains enjoyable.

This stage tests whether protein and nutrient adequacy can improve without destabilizing the successful weight-maintenance pattern.

### Stage 2: if further weight loss is desired

Only after Stage 1 is easy and the prescriber agrees:

- remove a modest additional 150–250 kcal/day from the remaining low-protein calorie sources;
- keep protein and strength training stable;
- use weekly weight trends, waist, function, and clinical markers rather than daily noise;
- pause escalation if fatigue, strength, nutrition, or GI tolerance worsens.

A modest change is more defensible than returning to 1,200-kcal restriction, trying to out-exercise intake, or improvising the medication schedule.

## Dessert and chips can stay—deliberately

The target is not purity. A stable plan can include portion-controlled treats. The leverage comes from deciding which treats are worth the calories, reducing duplicate default snack sources, and protecting protein first.

For example, keep one preferred dessert event and one preferred savory snack rather than cycling among many large packages. That also reduces the food-novelty pattern associated with extreme-intake days.

## Travel version

Travel historically adds roughly 600 kcal/day. The highest-value preparation is not packing every supplement; it is preserving the new defaults:

- identify a protein-forward breakfast before arrival;
- know where the first chili/tofu/seitan/vegetarian-protein meal will come from;
- bring one familiar protein option if practical;
- choose planned novelty for one meal rather than letting every meal become unstructured.

## Caveats

- MFP nutrient entries can be wrong, and protein quality/digestibility is not captured.
- The 3 g and 7 g per 100-kcal cutoffs are menu-design rules, not biological thresholds.
- Food clusters with fewer than four days were omitted from the “familiar” list, so additional good options exist.
- The analysis identifies available substitution calories. It does not prove that a specific food caused obesity.
