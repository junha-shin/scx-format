# Contributing

## Writing

Use [SimpleEnglish and the shared text rules](https://github.com/junha-shin/agent-rules/blob/main/AGENTS.md).
The [SCX project rules](https://github.com/junha-shin/agent-rules/blob/main/projects/scx-format.md)
are stored in the same public repository.

Use short sentences and common words. Keep API names, data types, numbers,
requirements, limits, and data-safety steps exact.
A prose edit must not change runnable code or recorded output unless requested.
Write README files in ASCII English and save text as UTF-8 with LF line endings.

A link does not load rules into an agent automatically. Ask the agent to read
those files for SCX work. For another local project, review and copy the shared
AGENTS.md to a place where your tool reads instructions. Do not overwrite an
existing rule file without merging it first.

## Checks

Run these checks before committing:

```sh
python scripts/check_readmes.py
python scripts/build_standalone.py --check
```

When changing tutorial prose, keep code cells and results intact.
When changing executable examples, run the tutorial checks:

```sh
python scripts/check_tutorials.py --output .compat-results/my-tutorial-check
```

Use a new output path. CI checks Python/R examples on Linux, Windows, and macOS.
Keep published release tags and assets unchanged. A release-file fix needs a new version.
