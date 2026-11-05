# Changelog

## 1.4.0 (2026-11-05)

The models are renamed. Their coefficients, standardization, features, thresholds and scores are identical to 1.3.1.

| Before 1.4.0 | From 1.4.0 | Model file | Use when the cell type has |
|---|---|---|---|
| CisLink-E2G (`e2g`) | **CisLink-Extended** (`extended`) | `src/cislink/models/extended.json` | accessibility, histone marks, RNA-seq and Hi-C |
| CisLink-V2G (`v2g`) | **CisLink-Base** (`base`) | `src/cislink/models/base.json` | accessibility and RNA-seq only |

The new names follow the data each model needs, as with ENCODE-rE2G and ENCODE-rE2G Extended. Both models score
element–gene pairs and every variant–gene pair. Score columns, `load_model` keys, `cislink score --model` choices and
the score columns of the shipped data tables use the new names; the old names are no longer accepted.

## 1.3.1 (2026-10-27)

Zenodo software archive; no change to the models.

## 1.3.0 (2026-10-27)

Models with cross-epigenome activity–expression correlations from the Roadmap activity reference and distance ×
accessibility; `cislink reference download` and `cislink add-f1`.
