# Variant-pair LD by population

<!-- workflow-badges:start -->
[![Docker Image CI](https://github.com/AoU-Multiomics-Analysis/AnalysisTools/actions/workflows/docker-image.yml/badge.svg)](https://github.com/AoU-Multiomics-Analysis/AnalysisTools/actions/workflows/docker-image.yml)
[![R lint](https://github.com/AoU-Multiomics-Analysis/AnalysisTools/actions/workflows/r-lint.yml/badge.svg)](https://github.com/AoU-Multiomics-Analysis/AnalysisTools/actions/workflows/r-lint.yml)
[![Update README workflow badges](https://github.com/AoU-Multiomics-Analysis/AnalysisTools/actions/workflows/update-readme-badges.yml/badge.svg)](https://github.com/AoU-Multiomics-Analysis/AnalysisTools/actions/workflows/update-readme-badges.yml)
[![WDL validation](https://github.com/AoU-Multiomics-Analysis/AnalysisTools/actions/workflows/wdl-validation.yml/badge.svg)](https://github.com/AoU-Multiomics-Analysis/AnalysisTools/actions/workflows/wdl-validation.yml)
<!-- workflow-badges:end -->

This WDL 1.0 workflow calculates unphased r² with PLINK 2.
It uses genotype dosages when they are present. Otherwise, it uses hard calls.
It produces one row for each requested pair in the full cohort and each population.

## Inputs

Use `workflows/pair_ld.wdl`. Supply these inputs:

| Input | Required content |
| --- | --- |
| `pgen`, `pvar`, `psam` | The three files from one PLINK 2 dataset |
| `pgi` | Optional external PGEN index; supply it when the PGEN needs one |
| `variant_pairs` | TSV with `variant_a` and `variant_b` columns |
| `sample_populations` | TSV with `FID`, `IID`, and `population` columns |
| `docker_image` | Published image built from `envs/Dockerfile` |

Use exact variant IDs from the PVAR file. Supply two different IDs per pair.
The workflow retains the pair order and duplicate or reversed pairs.
The PVAR can be plain text or Zstandard compressed (`.pvar.zst`).
The PSAM must be plain text and must have a header.

Each PSAM sample must have exactly one row in the population table.
The table must contain no extra samples. Match samples with both FID and IID.
Use FID `0` when the PSAM has no FID column. Duplicate FID/IID combinations are not supported, including combinations that differ only by SID.
The population label `ALL` is reserved for the full cohort.
The three-row population example shows the format only. Replace it with the complete cohort table.

The workflow supports human autosomal, biallelic variants.
It reports a status for other chromosomes and allele sets.
It does not combine alternate alleles at a multiallelic site.
Variants on different autosomes are allowed; their r² measures cross-chromosome correlation.

## Calculation

One task selects all eligible requested variants into a small PGEN dataset.
This step preserves dosages. For each enabled group, PLINK runs one LD calculation
with an anchor list from the requested `variant_a` IDs.
The script reads the report in sequence and retains only the exact requested pairs.
Input order, duplicates, and reversed pairs remain in the output.

For pairs on the same chromosome, the distance limit covers the largest requested
distance, rounded up to a positive whole number of kilobases.
PLINK can report additional nearby pairs; the script discards them.
If any eligible pair spans chromosomes, the batch uses `inter-chr` and no distance
limit. This also retains distant pairs on the same chromosome in a mixed request.
In both cases, the negative r² threshold retains zero, low, and undefined LD values.
The task removes each temporary LD report after reading it. The logs remain in the
exported archive. Cross-chromosome batches can produce many extra comparisons;
increase `extra_disk_gb` if one report needs more task storage.

Every sample in the PSAM contributes to the full-cohort group.
The task creates a temporary PSAM without parent columns to include samples with known parents.
It preserves the sample order and does not remove related samples.
If you need unrelated samples only, filter the input dataset first.

By default, a group needs at least 50 samples for LD calculation.
Set `allow_small_groups` to `true` to permit groups with 3–49 samples through PLINK's `--bad-ld` option.
These rows have `small_group=true`. Groups with fewer than 3 samples remain disabled.

## Outputs

`ld_results` is a TSV with these columns:

| Column | Meaning |
| --- | --- |
| `pair_index` | Input row number, starting at 1 after the header |
| `variant_a`, `variant_b` | Requested variant IDs |
| `population` | Population label, or `ALL` for the full cohort |
| `n_samples` | Number of selected samples in the group |
| `r2` | Unphased dosage r², or `NA` |
| `status` | Calculation status |
| `small_group` | Whether the group has fewer than 50 samples |

Status values are `ok`, `missing_variant`, `unsupported_chromosome`,
`unsupported_alleles`, `too_few_samples`, and `undefined_ld`.
PLINK reports undefined LD for cases such as a monomorphic variant or insufficient nonmissing genotype data.
An unexpected missing pair report or PLINK error fails the task.
Missing or unsupported variants take status precedence over sample count.

`n_samples` is the group size. It is not the count of samples with nonmissing data at both variants.
The workflow also returns sample counts, the PLINK version, PLINK logs, and task stdout and stderr.

## Terra use

1. Use the published image from `examples/inputs.json`. That example pins the exact tested image by digest.
2. Upload the genotype files and the two TSV files to your workspace bucket.
3. Add the WDL as a Terra method.
4. Set the inputs with `examples/inputs.json` as a guide. Replace every example path and image name.
5. Supply `pgi` as an additional File input if the dataset has an external index.

The JSON example is a Terra workflow input map. It is not a script argument wrapper.
WDL passes explicit File inputs through safely quoted named CLI arguments.
The task creates all sample lists and variant lists after localization.
The script checks input readability and fails on an unresolved cloud URI.
The workflow creates no files at workflow scope.

The default task uses 2 CPUs and 8 GB of memory.
Its disk request is twice the input file size plus 10 GB.
Increase `extra_disk_gb` if the dataset needs more temporary space.
Use positive CPU and memory values.

The complete workflow has **not been tested on Terra**.
No cloud job has been submitted.

## Small checks and GitHub Actions

Source repository: [AoU-Multiomics-Analysis/AnalysisTools](https://github.com/AoU-Multiomics-Analysis/AnalysisTools).
The `.github` workflows use paths relative to the repository root.

`wdl-validation.yml` checks WDL syntax, workflow-scope file creation, and command-time localization and quoting.
It also checks two input error cases. It builds no image.

`docker-image.yml` builds the Linux image and runs the synthetic LD regressions and input checks.
The tests use small synthetic genotype files. They test dosages, missing dosages,
population selection, missing variants, monomorphic variants, indels, reversed and
duplicate pairs, distant pairs, cross-chromosome pairs, and small sample groups.
A regression test compares the batch results with real per-anchor PLINK calculations
from the prior method. It also checks the PLINK logs for one LD run per enabled group.
Automatic image builds run only when `scripts/pair_ld.py` changes.
Use the manual workflow trigger after a Dockerfile or dependency change.
After the checks pass on `main`, the action publishes the tested image to GitHub Container Registry.
It publishes a `latest` tag and a `sha-<full commit SHA>` tag. Pull requests run the checks without publishing.
Use the digest from the action summary as the WDL `docker_image` input to select the exact tested image.
The action uses the repository's `GITHUB_TOKEN`; no separate registry secret is required.
The [published LD image](https://github.com/orgs/AoU-Multiomics-Analysis/packages/container/package/plink-pair-ld) is public.
Anonymous manifest access was verified for the digest in `examples/inputs.json`.
Terra can pull this image without registry credentials.
No local container build is needed for these checks.

Run the checks directly when PLINK 2 and miniwdl are available:

```bash
miniwdl check workflows/pair_ld.wdl
python tests/check_wdl.py
python -m unittest discover -s tests -v
```

Local checks passed with PLINK `2.0.0-a.6.9` and miniwdl `1.14.2`.
The [GitHub Actions image build](https://github.com/AoU-Multiomics-Analysis/AnalysisTools/actions/runs/38347449373)
passed all five tests and published the batched container pinned in `examples/inputs.json`.

## References

- [PLINK 2 LD measures and filters](https://www.cog-genomics.org/plink/2.0/ld)
- [PLINK 2 file formats](https://www.cog-genomics.org/plink/2.0/formats#vcor)
- [PLINK 2 input files and indexes](https://www.cog-genomics.org/plink/2.0/input#pgen)
