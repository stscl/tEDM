# tEDM 1.3

### new

* Correct maintainer surname spelling from `Lv` to `Lyu` for pinyin compliance (#175).

# tEDM 1.2

### enhancements

* Support specifying library units via `lib` parameter in `multispatialccm` generic (#161).

* Permit `simplex` and `ic` generics to accept varying E, k, and tau inputs (#148).

* Enforce strict floating-point comparisons across cpp sources (#143).

* Unify font specification in S3 plotting method for cross-mapping results (#135).

### breaking changes

* Use consistent masking for library and prediction indices in cross mapping parameter selection (#157).

* Correct library and prediction indices handling in cross mapping (#145).

# tEDM 1.1

### new

* Introduce configurable distance metrics for cross mapping (#120).

### enhancements

* Safeguard transient removal logic in logistic map to prevent index errors (#96).

### breaking changes

* Adjust time-delay embedding logic to support flexible tau configurations (#111).

* Refine `lib` and `pred` index filtering logic for cross mapping (#102).

# tEDM 1.0

### new

* First stable release (#93).
