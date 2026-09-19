# gTranslate

Predicts the genetic translation table (GTT) used by prokaryotic genomes from Prodigal
gene calls and an ensemble classifier.

## Releasing a new version

The version is defined in one place, `gtranslate/__init__.py`. `pyproject.toml` reads it
dynamically (`version = {attr = "gtranslate.__version__"}`), as does `docs/src/conf.py`,
so it must not be duplicated anywhere else.

Whenever the version is bumped, **both** of these documentation files must be updated in
the same change:

* `docs/src/changelog.rst` — a section for the new version listing every user-visible
  change, grouped under a heading such as `Major Changes:` or `Bug Fixes:`.
* `docs/src/announcements.rst` — a section for the new version with the month and year,
  summarising what users need to know or act on.

Newest version first in both files. Neither file should be left behind when the version
changes, even for a patch release.

## Tests

```
python -m unittest discover -s tests
```

Running gTranslate itself needs `prodigal` on `PATH` and the classifier models, located
via `--custom_model_path` or the `GTRANSLATE_MODEL_PATH` environment variable.
