# Changelog

## tEDM 1.3

## tEDM 1.2

CRAN release: 2026-01-15

#### enhancements

- Support specifying library units via `lib` parameter in
  `multispatialccm` generic
  ([\#161](https://github.com/stscl/tEDM/issues/161)).

- Permit `simplex` and `ic` generics to accept varying E, k, and tau
  inputs ([\#148](https://github.com/stscl/tEDM/issues/148)).

- Enforce strict floating-point comparisons across cpp sources
  ([\#143](https://github.com/stscl/tEDM/issues/143)).

- Unify font specification in S3 plotting method for cross-mapping
  results ([\#135](https://github.com/stscl/tEDM/issues/135)).

#### breaking changes

- Use consistent masking for library and prediction indices in cross
  mapping parameter selection
  ([\#157](https://github.com/stscl/tEDM/issues/157)).

- Correct library and prediction indices handling in cross mapping
  ([\#145](https://github.com/stscl/tEDM/issues/145)).

## tEDM 1.1

CRAN release: 2025-08-25

#### new

- Introduce configurable distance metrics for cross mapping
  ([\#120](https://github.com/stscl/tEDM/issues/120)).

#### enhancements

- Safeguard transient removal logic in logistic map to prevent index
  errors ([\#96](https://github.com/stscl/tEDM/issues/96)).

#### breaking changes

- Adjust time-delay embedding logic to support flexible tau
  configurations ([\#111](https://github.com/stscl/tEDM/issues/111)).

- Refine `lib` and `pred` index filtering logic for cross mapping
  ([\#102](https://github.com/stscl/tEDM/issues/102)).

## tEDM 1.0

CRAN release: 2025-07-15

#### new

- First stable release
  ([\#93](https://github.com/stscl/tEDM/issues/93)).
