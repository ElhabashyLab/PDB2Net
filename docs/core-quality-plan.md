# Verification status of the scientific baseline

This baseline completed 617 tests under Python 3.12/Linux without skipped or
expected-failure cases. The frozen dataset and polymer-length corrections are
included. The results describe local verification, not a published release.
See [validation](core-quality-audit.md) and [development](development.md).

## Remaining work

- Preserve fully modified supported chains through initial selection and resolve
  otherwise unknown types from consistent polymer metadata.
- Correct BLAST alignment-span coverage when qcovs is unavailable and align
  readiness checks with the actual query-quality filter.
- Add negative search-boundary tests and correct the independent identity-edge
  oracle to accept a spanning tree instead of all cross-structure pairs.
- Exercise scientific calculations with the installed wheel outside the checkout.
- Check native Windows file labels and Windows/macOS process behavior.
- Confirm remote CI for the eventual submitted revision before release work.

This status records remaining work; it does not authorize publication,
reference downloads, store rebuilds or deployment.
