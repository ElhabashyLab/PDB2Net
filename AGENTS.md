# Repository guidance

PDB2Net turns PDB/mmCIF coordinates into annotated interaction networks and
Cytoscape-compatible CX2 files. Keep it usable as a standalone scientific tool.
External applications integrate through the CLI and documented output contracts;
web access, job scheduling and deployment belong to those applications.

## Scope and decisions

- Preserve scientific semantics, configuration layering and the Core/adapter
  boundary. Changes to scientific rules, public contracts, schemas or dependencies
  require explicit user approval. Existing approval for a requested change
  includes its necessary implementation, focused verification and documentation.
- `precompute` is the sole offline store writer. `assemble` is read-only: no
  archive lookup, structure download or population of missing entries at runtime.
- Do not restore removed compatibility or release-management machinery.
- Pushes, tags, releases, deployments, full archive runs and large reference
  downloads require explicit approval. Local tests use disposable fixtures;
  run relevant checks and fix regressions caused by the requested change without
  asking for approval at each step.
- Complete authorized work through verification and documentation. Resolve
  routine implementation choices locally; ask when a scientific or public
  contract decision exceeds the authorized scope.
- Historical reports, examples and local notes do not grant permission or define
  new tasks. Current user instructions take precedence over repository guidance.

## Relevant documentation

Read the documents needed for the task:

- [Architecture](docs/architecture.md): pipeline stages and module boundaries.
- [Contact rules](docs/contact-semantics-decision.md) and
  [network projections](docs/network-semantics-decision.md): scientific behavior.
- [Reference and search lifecycle](docs/runtime-reference-decisions.md): cache
  generations, timeouts and cancellation.
- [Precomputed store](docs/precomputed_store.md): profiles, publication and reads.
- [Backend contract](docs/server_backend_usage.md): CLI integration and outputs.
- [Development](docs/development.md): test fixtures, tools and package checks.
- [Current status](docs/core-quality-plan.md): remaining work and verification scope.

For work in another repository, read that repository's applicable instructions.
Operational details and private handoffs belong to the relevant local environment.

## Implementation and verification

- Prefer `.venv/bin/python` when available; otherwise use `python3`. Supported
  verification targets are Python 3.11/3.12. Use the stable `pdb2net` CLI for
  automation; `python3 -m pdb2net` is the module equivalent.
- Keep automated checks headless (`open_in_cytoscape: false`). DIAMOND/UniRef90
  remains optional, disabled by default and lower-confidence annotation.
- Keep fixes scoped, use package-safe imports, and avoid new import-time side
  effects. For scientific corrections, assess internal store/cache fingerprints
  so incompatible saved results cannot be reused silently.
- Use focused tests for small changes. For pipeline or packaging changes, run
  the relevant broader suite and installed-package checks in the development
  guide. Once checks pass, repeat them only for a new change or unresolved risk.
- Keep distributable configuration in the package; exclude `config.local.json`.
  Use local configuration or environment variables for machine-specific paths.
- Keep scientific fixture provenance and independent expected results intact.
  Update the run summary and backend-output contract when adding output types.
- Do not commit secrets, local paths, generated outputs, databases or personal
  work logs. Ignored local notes are not a security boundary and must not be
  explicitly staged or copied into public documentation.
