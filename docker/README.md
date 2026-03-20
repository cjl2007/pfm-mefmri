# Container Notes

This repository now supports containerized launches without changing the native pipeline entrypoint.

## 1) Build an image

Two Dockerfiles are provided:

- `docker/Dockerfile`: minimal scaffold (copy repo only)
- `docker/Dockerfile.runtime`: production-candidate runtime image

Build the runtime image:

```bash
bash bin/build_container_image.sh
```

Optional custom tag/base:

```bash
bash bin/build_container_image.sh --tag mefmri-pipeline:runtime --base nipreps/fmriprep:24.1.1
```

Run a smoke test after build:

```bash
bash bin/smoke_test_container_image.sh mefmri-pipeline:runtime
```

## 2) Run with wrapper (recommended)

```bash
bash bin/mefmri_pipeline_container.sh \
  --engine docker \
  --image mefmri-pipeline:runtime \
  --fs-license /path/to/license.txt \
  /path/to/study/ME06 \
  /path/to/config/mefmri_wrapper_config_container.sh
```

Use `--engine apptainer` and a `.sif` image on HPC systems where Docker is not available.

## 3) FreeSurfer license behavior

The wrapper bind-mounts the license into `/licenses/` and sets both:

- `FS_LICENSE`
- `FS_LICENSE_FILE`

The pipeline resolves license in this order:

1. `FS_LICENSE`
2. `FS_LICENSE_FILE`
3. `${FREESURFER_HOME}/license.txt`

`anat_hcp` stage now fails fast if no readable license is found.

## 4) Notes on optional tools

- SimNIBS/CHARM is not installed by default in `Dockerfile.runtime`.
- If you run `anat_charm`, either:
  - build a derived image that installs SimNIBS and places `charm` on PATH, or
  - set `CHARM_BIN` in your config to a valid in-container path.
