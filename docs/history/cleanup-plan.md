# Clean release architecture history

The clean-release work established the standalone Core and its external CLI
integration boundary. This summary preserves architectural rationale, not an
implementation checklist or permission to resume previous work.

## Decisions at that baseline

- Raw-file analysis remains the scientific reference path.
- Offline `precompute` is the only store writer; `assemble` consumes complete
  fixed-profile entries without resolving archives or downloading structures.
- Public capability, CLI, web-output and artifact contracts are explicit.
  External applications install the package and use those contracts.
- Headless CX2 export is independent of the optional live Cytoscape extra.
- Web access, scheduling, databases and deployment are owned by external
  applications and are not Core responsibilities.
- Removed runtime lazy population and earlier compatibility/release-management
  mechanisms are not part of the supported architecture.

The former Core-owned store has since been replaced by the prepared-input API.
Current behavior is specified in [Architecture](../architecture.md),
[Prepared structures](../prepared_structures.md) and the
[backend contract](../server_backend_usage.md). Remaining verification work is
tracked in [current status](../core-quality-plan.md).

Earlier published revisions may still contain the original implementation plan.
This document does not redefine the behavior of those historical revisions.
