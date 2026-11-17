# RS-utils 0.11.0rc1 release checklist

This checklist separates local REL-01 evidence from release actions that must
be performed by a human with deployment authority. A Git tag, package
publication, and cluster deployment are intentionally outside this task.

## REL-01 evidence

- [x] version updated consistently
- [x] wheel metadata reports `0.11.0rc1`
- [x] changelog/release notes updated
- [x] HT quick-start added
- [x] HT documentation linked from Sphinx
- [x] README points to the HT quick-start and readiness matrix
- [x] legacy-vs-HT UppASD distinction documented
- [x] stale documentation targets removed/fixed
- [x] cluster deployment guidance added
- [x] release checklist added
- [x] full source tests pass
- [x] installed-wheel tests pass
- [x] documentation builds
- [x] HT support matrix reviewed
- [x] wheel artifact retained
- [x] CLI/campaign smoke test passes
- [x] scientific source code unchanged (only non-scientific Python 3.9 typing compatibility and version metadata changed)

## Recorded verification

- `import rsutil` and wheel metadata both report `0.11.0rc1`.
- Wheel: `rsutil-0.11.0rc1-cp311-cp311-macosx_26_0_arm64.whl`.
- Source suite: **179 passed, 16 warnings**.
- Installed-wheel suite in an external environment: **179 passed, 16 warnings**;
  `pip check` passed and both package/native-extension paths resolved under
  `site-packages`.
- CLI smoke: `rsutil --help`, `rsutil campaign --help`, and fixture
  `validate -> create -> status` passed.
- Documentation: `make clean html` with warnings-as-errors under Python 3.11
  passed; required HT pages rendered and stale targets were removed.

## Human release and deployment steps

- [x] release commit identified
- [ ] Git tag created
- [ ] wheel published/copied to deployment location
- [ ] cluster smoke test passes

## Publishing to PyPI (trusted publishing)

Releases are built and published by `.github/workflows/release.yaml` without
API tokens. One-time setup (done by a maintainer on the websites):

1. On <https://pypi.org/manage/account/publishing/> add a *pending publisher*
   for the project name `rslmto-utils`: owner `aibergman`, repository
   `RS-utils`, workflow `release.yaml`, environment `pypi`. Repeat on
   <https://test.pypi.org/manage/account/publishing/> with environment
   `testpypi`.
2. In the GitHub repository settings create the environments `pypi` and
   `testpypi` (optionally require a reviewer for `pypi`).

Per release:

- [ ] `rsutil/_version.py` and `CHANGELOG.md` updated; tests green on `main`
- [ ] `NOTICE` reviewed (provenance of derived sources is current)
- [ ] Dry run: *Actions → Release → Run workflow* with target `testpypi`, then
      `pip install --pre -i https://test.pypi.org/simple/ rslmto-utils` in a clean
      environment and run `rsutil --help`
- [ ] Tag `vX.Y.Z` matching the package version and push it; the workflow
      refuses a tag that does not match `rsutil/_version.py`
- [ ] A final (non-`rc`) version is needed for plain `pip install rslmto-utils`;
      pre-releases require `--pre`
