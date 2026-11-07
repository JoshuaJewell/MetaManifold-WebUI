<!-- SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk> -->
<!-- SPDX-License-Identifier: CC-BY-SA-4.0 -->

# MiSeq data: published terms and verified source

Reviewed 7 November 2026. **Published licence basis: CC BY 4.0.** This records
the basis already adopted upstream; it supplies no new permission, scientific
sign-off or release approval.

## Original-publisher and upstream evidence

The mothur website's footer states: **“All content on this site is available
under the CC-BY 4.0 license.”** See the publisher's pinned
[footer](https://github.com/mothur/mothur.github.io/blob/4ea040bd2b5c3f5564ab5eb8dc22fd81c218ce73/_includes/footer.html#L25),
[licence](https://github.com/mothur/mothur.github.io/blob/4ea040bd2b5c3f5564ab5eb8dc22fd81c218ce73/LICENSE.md)
and [author configuration](https://github.com/mothur/mothur.github.io/blob/4ea040bd2b5c3f5564ab5eb8dc22fd81c218ce73/_config.yml#L19),
which names Patrick D. Schloss. The [original SOP](https://mothur.org/wiki/miseq_sop/)
links the mothur-hosted archive inspected here.

Upstream already assigns `CC-BY-4.0` specifically to
`data/MiSeq_SOP/**/*.fastq.gz` in
[REUSE.toml](https://github.com/JoshuaJewell/MetaManifold-WebUI/blob/89be534ab5cf5e286b34a8e0bf9310527c492ec5/REUSE.toml#L31)
and its [README](https://github.com/JoshuaJewell/MetaManifold-WebUI/blob/89be534ab5cf5e286b34a8e0bf9310527c492ec5/README.md#L133).
JoshuaJewell introduced this in [commit 4dce274](https://github.com/JoshuaJewell/MetaManifold-WebUI/commit/4dce274580ad88230b7d97818919ceb6997aae43),
before the statistics port. Preserve the existing `2019 Patrick D. Schloss`
notice; its year is an upstream attribution, not a newly verified date of
the experiment.

**Assessment:** the publisher's site-wide notice, exact file matches and
upstream's data-specific annotation establish a documented public-licence
basis. Applying the publisher's notice to its own SOP download is a scope
inference. The archive contains no separately named licence/readme/terms/notice
entry and has no archive comment. No dataset-specific correspondence was
obtained. A second notice is not inherently required for a public licence;
publisher clarification remains available if the maintainer questions that scope.

The [Galaxy tutorial deposit](https://doi.org/10.5281/zenodo.800651), created by
Saskia Hiltemann, also records CC BY 4.0, and its six corresponding declared
MD5 values match these inputs. This is corroborative distribution history,
not proof that a third-party depositor can grant the original publisher's rights.

## Input identity

The SOP's [original ZIP](https://mothur.s3.us-east-2.amazonaws.com/wiki/miseqsopdata.zip),
fetched 7 November 2026, has 45 entries and SHA-256
`728f7544824f8ca0223618b18bb60fce6a5dc20ae3d476267dbd5cf5594a35d4`.
Each repository gzip decompresses to exactly the corresponding original
`MiSeq_SOP/<filename>.fastq` entry. Gzip packaging and assignment to two run
directories changed; raw FASTQ contents did not.

| Input under `data/MiSeq_SOP` | Records | Decompressed SHA-256 |
| --- | ---: | --- |
| `run_A/F3D0_S188_L001_R1_001.fastq.gz` | 7,793 | `f908bf940a2413f9d026615afca8cacfbcb36bf3c91ffa5ad809738b0f7064d3` |
| `run_A/F3D0_S188_L001_R2_001.fastq.gz` | 7,793 | `1b3feb32e2bef0a573382cdd2b4eb833af5a3fb093e3544f6882da371db64937` |
| `run_A/F3D1_S189_L001_R1_001.fastq.gz` | 5,869 | `27f8791c021cbdbc6e772a608f3491ac62dfe6b4d4d148871b09b6716621286d` |
| `run_A/F3D1_S189_L001_R2_001.fastq.gz` | 5,869 | `0ea14696f4ededb4c693bb5d552da13ff0953ff3e8fb5b5ba456dc8ec3e929a3` |
| `run_B/F3D2_S190_L001_R1_001.fastq.gz` | 19,620 | `216b03461b02b333d693c52fbfc6aba595b59c0da57939ff3dafc41fc8d9b655` |
| `run_B/F3D2_S190_L001_R2_001.fastq.gz` | 19,620 | `2c6126d632a1a42c254d8a6455e2096f7f8107d5e25a22d00266865d622d1467` |

The review attachment retains the source manifest, archive inventory and
input-identity CSV. Record counts per direction are not biological replicate counts.

## Reuse and attribution

[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/), subject to its
[legal text](https://creativecommons.org/licenses/by/4.0/legalcode.en), allows
copying and adaptation, including commercial use. Preserve supplied credits
and notices, link the source and licence, identify changes, and avoid implying
endorsement or adding restrictions to the licensed material. It gives no
warranty and licenses only rights the licensor can grant. The application's
code licences remain separate.

Suggested attribution for this subset and its derived evidence:

> MiSeq SOP example reads, Patrick D. Schloss / mothur project. Kozich JJ,
> Westcott SL, Baxter NT, Highlander SK, Schloss PD (2013), *Development of a
> dual-index sequencing strategy and curation pipeline for analyzing amplicon
> sequence data on the MiSeq Illumina sequencing platform*, Applied and
> Environmental Microbiology 79(17):5112–5120,
> https://doi.org/10.1128/AEM.01043-13. Source:
> https://mothur.org/wiki/miseq_sop/ (accessed 7 November 2026).
> Published licence basis: CC BY 4.0,
> https://creativecommons.org/licenses/by/4.0/.
> Changes: F3D0/F3D1/F3D2 paired reads were selected, gzip packaged and divided
> into two runs, with original FASTQ contents unchanged. Derived counts and
> descriptive results were produced by the recorded MetaManifold pipeline and
> analysis; they are not original-author findings or endorsements.

Preserve the existing raw-data copyright notice alongside distribution.
This attribution does not relabel the entire application or every evidence
file. No change is made here to README licensing/acknowledgements, `LICENSES`
or `CITATION.cff`.

## Remaining decisions

The [publisher's provenance clarification](https://forum.mothur.org/t/where-is-the-miseq-sop-example-data-from/21335)
ties the sequencing to Kozich et al.; the SOP identifies these samples as
time points from Female 3. Their use does not establish independent biological
replication, PR2's suitability for bacterial 16S, or broad biological performance.
The data-use basis is documented. Use the [review record](issue-65-review.md)
for actual scientific and maintainer decisions.
