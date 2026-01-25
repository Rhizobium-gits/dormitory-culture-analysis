# Data Directory

## Data Availability

Due to privacy considerations for survey participants, the raw data used in this study cannot be made publicly available.

## Data Access

Researchers interested in accessing the data for research purposes should contact the corresponding author:

- **Tsubasa Sato**: satotsubasa@keio.jp
- **Misaki Iio**: misaki.iio@keio.jp

## Data Description

The survey data includes responses from 153 domestic students residing in H village dormitory at Keio University. The data was collected through five longitudinal surveys conducted in May, June, July, October, and November 2025.

### Variables

| Variable | Description | Type |
|----------|-------------|------|
| `id` | Participant identifier (anonymized) | Character |
| `building` | Building name (Basil, Turmeric, Rosemary, Paprika) | Factor |
| `building_type` | Building type (Male-only, Female-only, Co-ed) | Factor |
| `sex` | Participant sex (Male, Female) | Factor |
| `resident_type` | Resident type (New, Incumbent) | Factor |
| `survey_month` | Survey month (May, June, July, October, November) | Factor |
| `Q1` | Easy to communicate with people in the dormitory (1-5) | Integer |
| `Q2` | Feel barriers when talking with opposite sex (1-5) | Integer |
| `Q3` | Good friendships within my unit (1-5) | Integer |
| `Q4` | Keep my room clean (1-5) | Integer |
| `Q5` | Want to actively participate in dormitory events (1-5) | Integer |
| `Q6` | Interested in presentations/public speaking (1-5) | Integer |
| `Q7` | Can confide personal matters to others (1-5) | Integer |
| `Q8` | Feel attachment to my building (1-5) | Integer |

### Sample Size by Building

| Building | Gender Composition | n | Percentage |
|----------|-------------------|---|------------|
| Basil | Male-only | 53 | 34.6% |
| Turmeric | Female-only | 29 | 19.0% |
| Rosemary | Co-ed (same-floor) | 29 | 19.0% |
| Paprika | Co-ed (floor-separated) | 42 | 27.5% |
| **Total** | | **153** | **100%** |

## Ethical Considerations

This survey was approved and conducted with informed consent from all participants. The consent statement assured participants that:

1. Responses would be statistically processed
2. Data would be used for improving event management and living conditions at H-village
3. Information identifying specific individuals would not be published
4. Participation was voluntary

## Simulated Data

For testing purposes, the analysis script (`R/analysis_main.R`) includes a function `create_demo_data()` that generates simulated data with similar structure and statistical properties to the actual data. This allows researchers to run and understand the analysis pipeline without access to the original data.
