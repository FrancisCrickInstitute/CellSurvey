# Working agreement: what Claude may do without asking

This file records how Claude Code and Dave work together in this repo. It sets out what Claude may do on its own, what needs a "yes" first, and what Claude should never do.

The hard limits are also enforced by `.claude/settings.json` (permission `allow` / `ask` / `deny` rules), which Claude Code applies whatever the model decides. This file covers the judgment calls that can't be written as command patterns. **Update both as we go**: when Dave says "you can just do X from now on" or "never do Y", add it here (and to `settings.json` if it's a command).

Sections 1–4 are generic and meant to be copied to other repos (or to `~/.claude/CLAUDE.md` for all repos). Section 5 is specific to this project.

## 1. How we work

- **Planning is not permission.** When we're discussing or planning a feature, or Dave has answered design questions, summarise the plan and wait for an explicit "go ahead", "make the changes" or "implement it" before editing files. (Learned 2026-11-06: InstanSeg support was implemented during what was meant to be planning.)
- **Claude prepares changes; Dave commits and pushes.** Claude edits files, checks them as far as it can, and reports what changed. Dave reviews, commits and pushes, unless he explicitly asks Claude to.
- **After an interruption or a rejected tool call, stop.** Don't retry or carry on with the next step; say where things stand and wait. A rejection may be deliberate even if the terminal looked frozen.
- **Report results faithfully.** Say what was actually run and what wasn't (e.g. "compiles, but not run end to end"). Don't call something tested unless it was run.
- **Decisions that change results** (defaults, units, which features clustering uses) are Dave's: lay out the options with a recommendation and ask.
- **Keep the docs in step.** When behaviour, flags or plans change, update `AGENTS.md` (and `docs/` / `README.md` if users see it) in the same change.

## 2. OK without asking

- **Reading**: any file in the repo; searching (Glob / Grep); `git status`, `log`, `diff`, `show`, `branch` (listing), `fetch`, `check-ignore`, `ls-files`, `rev-parse`.
- **Editing code and docs in the repo** once Dave has asked for the change (see "Planning is not permission"). No new edits after the requested change is done without checking.
- **Python files (`.py`) in the repo: OK to create or edit without asking** (added 2026-11-21), while working on something Dave has asked for. This covers the code itself; "Planning is not permission" still applies to *starting* new work.
- **Markdown files (`.md`) in the repo: OK to create or edit without asking** (added 2026-11-21; extends the earlier `AGENTS.md`/`CLAUDE.md` rule), e.g. `AGENTS.md`, `CLAUDE.md`, `README.md`, `docs/*.md`: recording findings, test results, plans, notes and doc updates. Still say what changed. **Exception: this file (`WORKING_AGREEMENT.md`)**: changing its rules needs a "yes" first (see §3), except recording a rule Dave has just stated.
- **Scratch work**: writing and running throwaway scripts, test data and outputs in the session's scratchpad directory (never in the repo).
- **Checks**: compiling (`python -m py_compile`) and running small local test scripts in the existing environments without installing anything (`pixi run --no-install …`).
- **Web lookups**: reading documentation and package pages (PyPI, GitHub, project docs) to check APIs.
- **Memory**: saving or updating Claude's memory notes about how we work.

## 3. Ask first

- **Git that changes state**: `add`, `commit`, `push`, `pull`, `switch` / `checkout`, `merge`, `rebase`, `stash`, creating a branch, tagging.
  - Before switching branches, check `git status` and stop if there are uncommitted changes or untracked files that would be overwritten or carried along.
- **Environments and dependencies**: `pixi install` / `lock` / `add` / `remove` / `update` / `clean`, editing `pixi.toml`, or anything else that downloads and installs packages.
- **Moving, renaming or deleting files in the repo**, including untracked ones (e.g. a stray `pixi.lock`). Back up to the scratchpad rather than deleting, and say so.
- **Long or heavy runs**: the full pipeline (`run.py`) or anything expected to take more than a few minutes, or to use a lot of memory or disk.
- **Anything outward-facing**: publishing (Artifacts, docs, gists), posting to GitHub (issues, PRs, comments), sending anything to an external service, or contacting servers on Dave's behalf (e.g. OMERO).
- **Changing Claude Code settings**: `.claude/settings*.json`, hooks, or this file's rules (except recording a rule Dave has just stated).

## 4. Never (unless Dave explicitly asks for that specific action)

- `git push --force` / `-f`, `git reset --hard`, `git clean`, `git branch -d` / `-D`, rewriting published history (`rebase` / `commit --amend` on pushed commits).
- Recursive or bulk deletes (`rm -rf`, `Remove-Item -Recurse`), or deleting anything outside the scratchpad.
- Hand-editing generated files (`pixi.lock`).
- Disabling safety checks: `--no-verify`, skipping hooks, turning off TLS verification to get round certificate errors.
- Putting credentials on a command line or in files (e.g. OMERO passwords); use the existing login session.
- Changing another user's state, e.g. an OMERO session's group (`setGroupForSession`).

## 5. This project

- **Local machine**: Windows, CPU only. The `default` pixi env currently has a broken shapely build; use the `instanseg` env for local checks (`pixi run --no-install -e instanseg …`). pixi isn't on PATH in Claude's shell (`~\.pixi\bin`). pixi can't reach conda-forge from Claude's shell (certificate error), so Dave runs installs.
- **HPC** (linux-64, GPU): all real runs happen there, run by Dave. Claude writes the commands and the expected results, then records Dave's results in `AGENTS.md`.
- **Output locations**: test outputs on the HPC go under `out/`; temporary files on the local machine go in the scratchpad.
