# Repository rules

- Write all README.md files in plain ASCII English, saved as UTF-8 with LF line endings.
- Use `-`, `x`, `->`, `<->`, and straight quotes instead of typographic symbols.
- Always specify encoding="utf-8" when reading or writing text in maintenance scripts.
- Never rewrite text using the operating system default encoding.
- Run `python scripts/check_readmes.py` before committing documentation changes.
- Preserve historical release tags and assets.

## SimpleEnglish

Apply these rules to the repository description, README files, quickstarts,
and new or revised user guides:

- Use short sentences and common words.
- State what the user can do, then show the steps.
- Remove repeated text and details that do not help the current task.
- Link to detailed guides instead of repeating their options.
- Keep exact API names, data types, numbers, and technical meaning.
- Keep requirements, limits, data sources, and steps needed to protect data.
- Use the same terms in Python and R guides when their meaning is the same.
- A prose edit must not change runnable code or recorded output unless requested.

The approved repository description is in docs/repository-description.txt.
Keep that file and the GitHub About description in sync when changing either.
