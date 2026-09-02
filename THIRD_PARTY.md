# Third-party software

## FFPErase

The adapter in `bin/run_ffperase_single_picard_pileup_nf.sh` invokes the
third-party `papaemmelab/nf-ffperase` workflow at runtime. The adapter does not
vendor or redistribute that workflow's source code, container image, or model
files.

- Upstream project: <https://github.com/papaemmelab/nf-ffperase>
- Source revision used by default:
  [`b0dd56cbd0a939896a966b9ce30c4d719b158170`](https://github.com/papaemmelab/nf-ffperase/commit/b0dd56cbd0a939896a966b9ce30c4d719b158170)
- Terms at that revision:
  [`LICENSE`](https://github.com/papaemmelab/nf-ffperase/blob/b0dd56cbd0a939896a966b9ce30c4d719b158170/LICENSE)
- Optional runtime container:
  `docker://papaemmelab/nf-ffperase:v1.0.0`
- Optional runtime models:
  <https://huggingface.co/papaemmelab/ffperase>

Nextflow fetches the pinned upstream source revision when the adapter runs. A
configured container runtime may fetch the container, and the adapter may fetch
missing models when model downloading is enabled. Those runtime downloads are
stored in user-selected caches or output directories and are not part of this
repository.

### Upstream terms

The upstream file is titled **FFPErase Terms of Use**. It is a restrictive
custom agreement, not a standard open-source license. At the pinned revision,
its material conditions include:

- FFPErase, its underlying content, and its output may be used only for
  personal, academic-research, and noncommercial purposes, including teaching
  and research at educational institutions.
- Publishing the upstream content or research results requires express written
  permission from Memorial Sloan Kettering Cancer Center (MSK).
- Redistributing or sharing the upstream content, in whole or in part, requires
  MSK's express written permission.
- Commercial use, commercial research, commercial products or services, and
  commercial distribution require express written permission and an
  appropriate license from MSK.
- Diagnosis, treatment, patient care, medical services, and generation of
  reports in medical, laboratory, hospital, or other patient-care settings are
  prohibited.
- Personally identifiable information (PII) and protected health information
  (PHI) must not be submitted to MSK in connection with the software.
- Information or data submitted to MSK is subject to the broad use rights
  described in the upstream agreement, and the submitter must have authority
  to provide it and grant those rights.
- Use of the upstream content must retain the appropriate original copyright
  notices required by MSK.
- The software and output are provided without warranties and are not a
  substitute for professional medical judgment.

This is only a practical summary. Users must read and comply with the complete
upstream terms before accessing or running FFPErase. The upstream terms direct
permission and licensing inquiries to `papaemme@mskcc.org`,
`domenicd@mskcc.org`, or `arangooj@mskcc.org`.

The source revision is pinned, but the container tag and the Hugging Face
`main` model URLs are not immutable digest or commit pins. Operators who require
fully reproducible third-party artifacts should pin or mirror those artifacts
separately, subject to the upstream terms.

This notice does not grant a license to FFPErase and does not add or imply an
overall license for this repository.
