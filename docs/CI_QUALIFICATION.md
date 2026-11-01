# CI qualification and wheel-install evidence

This document distinguishes an implemented GitHub Actions workflow, local
tests, and a **GitHub-hosted run with an actual conclusion**. Only the last of
these establishes cloud CI at an exact commit. A green local test suite does
not substitute for a remote workflow run.

## Workflow contract

`.github/workflows/test.yml` supports push, pull-request, and manual
`workflow_dispatch` events. The token is limited to `contents: read`; the
checkout action does not persist Git credentials. It does not use the
privileged `pull_request_target` event. The test job uses Ubuntu and CPython
3.11, 3.12, and 3.13, runs the repository's standard-library unit suite, and
compiles the source. The package job runs only after that matrix passes.

The package job provisions setuptools/wheel first (which **may download**
build tools from PyPI). It then builds the project wheel without build
isolation or runtime-dependency resolution, creates a separate virtualenv,
installs the produced wheel with `--no-index --no-deps`, changes outside the
source checkout, imports packaged code and `topic_manifest.json`, and invokes
the installed `calc-atlas` entry point. This is an installed-artifact check,
not a PyPI publication or a proof of portability beyond the runner.

## How to get a real GitHub-hosted result

After the new workflow is on the default branch:

1. Open [Python mathematics checks](https://github.com/thebrazenbeard/calc/actions/workflows/test.yml)
   in GitHub Actions and choose **Run workflow**, branch `main`. GitHub's
   documented `workflow_dispatch` event requires the workflow file on the
   default branch.
2. Or, with a suitably authorized GitHub CLI session:

   ```bash
   gh workflow run test.yml --repo thebrazenbeard/calc --ref main
   gh run list --repo thebrazenbeard/calc --workflow test.yml --branch main --limit 10
   gh run view RUN_ID --repo thebrazenbeard/calc --log
   ```

3. Record the run ID, the **exact run head SHA**, the event, and each Python
   matrix and package-job conclusion. Confirm that the head SHA is the
   intended release commit, not merely an earlier passing ancestor.
4. If a run fails, inspect the failing job and fix its mechanism; do not
   treat a manual rerun as proof the original failure was flaky.

Useful read-only API paths for independent confirmation are the repository's
Actions **workflow runs** and commit **check runs**. Do not infer GitHub
Actions status from the legacy combined commit-status API alone.

## Historical observation and evidence limitation

At `c2137b844b5d6a7268ef98a4502c6f78218dc0ff` (the Wave 9 receipt
commit), the available GitHub connector reported an empty legacy `statuses`
array and an empty result from a workflow-runs endpoint explicitly restricted
to **pull-request-triggered** runs. Neither result can establish that no
push-triggered Actions checks existed. GitHub distinguishes Actions **checks**
from older **commit statuses**, and Actions normally produces the former.
There is no independently confirmed hosted run verdict in this document.

The repository's local Python 3.12 tests passed 757/757 before this CI
change, and the CI contract checks cover the workflow configuration itself.
Those are local-test evidence, not cloud-execution evidence. The interactive
mathematical explorer and future calculus modules must be verified separately
as they are integrated.

## Sources

- [GitHub: events that trigger workflows](https://docs.github.com/en/actions/reference/workflows-and-actions/events-that-trigger-workflows)
- [GitHub: manually running a workflow](https://docs.github.com/en/actions/how-tos/manage-workflow-runs/manually-run-a-workflow?tool=webui)
- [GitHub: workflow permissions](https://docs.github.com/en/actions/reference/workflows-and-actions/workflow-syntax#permissions)
- [GitHub: status-check types](https://docs.github.com/en/enterprise-cloud@latest/pull-requests/reference/status-checks)

No external source code was copied into this workflow.
