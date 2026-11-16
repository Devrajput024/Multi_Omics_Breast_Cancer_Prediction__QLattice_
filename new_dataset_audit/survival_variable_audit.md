# Survival variable audit

## Observed variables and counts

Dataset 2 includes four relevant fields:

| Field | Observed values | Count/summary |
|---|---|---|
| Overall Survival (Months) | numeric | 780 nonmissing; range 1.233333333–337.0333333; 739 distinct values |
| Overall Survival Status | Living / Deceased | 373 / 407 |
| Relapse Free Status (Months) | numeric | 780 nonmissing; range 1.22–298.88; 743 distinct values |
| Relapse Free Status | Not Recurred / Recurred | 461 / 319 |
| Patient's Vital Status | Living / Died of Disease / Died of Other Causes | 373 / 256 / 151 |

The two time fields contain plausible fractional values and are labeled in months. This is the only evidence in the archive for their units; no codebook or source documentation states the time origin, censoring rule, event definition, last follow-up definition, rounding, or whether values have been transformed. Do not claim these details are verified.

## Coding and cross-checks

`Overall Survival Status` codes vital outcome coarsely as `Living` or `Deceased`. It agrees exactly with the more detailed `Patient's Vital Status`: all 373 `Living` cases are Living; all 407 `Deceased` cases are split into 256 `Died of Disease` and 151 `Died of Other Causes`.

For recurrence, `Relapse Free Status` uses `Recurred` and `Not Recurred`, but it is not an overall survival status. Among the 319 Recurred cases, 256 are `Died of Disease`, 50 are Living, and 13 died of other causes. Among the 461 Not Recurred cases, 323 are Living and 138 died of other causes. The existence of a `Not Recurred` label among deceased cases is compatible with other-cause mortality but cannot be interpreted further without endpoint definitions.

## Documentation conclusion

The names assert months and provide categorical coding, and internal consistency checks support the OS status mapping above. The archive has no documentation verifying units, endpoint definitions, time origin, censoring, or relapse event semantics. Treat survival data as undocumented time/status fields pending confirmation from the original data source; do not silently equate recurrence-free status with censoring/event indicators or infer censoring from `Living`/`Deceased` alone.
