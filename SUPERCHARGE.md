# Supercharge Plan

Action items distilled from Mateusz Pusz's ["Your Library is Great and Nobody
Uses It, Here's Why"](https://mpusz.github.io/mp-units/latest/blog/2026/09/03/your-library-is-great-nobody-uses-it-heres-why/)
and his "Supercharge Your C++ Project: 10 Tips to Elevate from Repo to
Professional Product" talk (using std::cpp, 2026-03-17). Flat list, most
foundational first; each item is an action plus its rationale, followed by
sub-items that are concretely actionable for **this** project
(`FrancoisCarouge/TypedLinearAlgebra`, CMake project
`fcarouge-typed-linear-algebra`, target `fcarouge-typed-linear-algebra::tlinalg`).
Several tips are already substantially implemented here (extensive CI badges,
a hand-written README reference, a "Lessons Learned" design-rationale section,
a CppNow 2026 talk) — those sub-items focus on closing the remaining gap
rather than starting from zero.

1. **Treat the repository as a product, not just code.** Technical excellence
   and "project excellence" (discovery, docs, community, packaging) are
   different skills — a brilliant library with none of the latter still fails
   to get adopted. Deliberately engineer every stage of the user journey
   (discover → evaluate → onboard → contribute → advocate).
   - a. Add a short "Roadmap"/"Vision" note (a `ROADMAP.md`, or a top section
     in `README.md`) framing the path toward `1.0` and the `std::linalg`
     integration ambitions already implied by the "Use Cases" section, so the
     product narrative is explicit rather than only implied by code quality.
   - b. Once per release, review `README.md`, `CONTRIBUTING.md`,
     `INSTALL.md`, and `AGENTS.md` together to confirm they still describe
     the same onboarding funnel instead of drifting independently (they
     currently already overlap on build commands — see item 14).
   - c. When triaging the backlog, explicitly weight adoption/onboarding work
     (docs, CI, packaging) against feature work in the same list instead of
     treating it as a lower-priority bucket.

2. **Give the project a unique, branded name** — consistent across GitHub
   repo, CMake target, package-manager name, namespace, and include path.
   Generic names (`units`, `linalg`, `http`) create discovery and branding
   collisions and make the project hard to search for or discuss.
   - a. The repo (`TypedLinearAlgebra`), CMake project
     (`fcarouge-typed-linear-algebra`), and namespace (`fcarouge`) are already
     distinct and consistent — document explicitly in `INSTALL.md` that the
     linked alias `tlinalg` (`fcarouge-typed-linear-algebra::tlinalg`) is a
     deliberate short form for ergonomics, so a first-time reader doesn't
     mistake it for an inconsistency.
   - b. When preparing the Conan/vcpkg recipes (item 26), use the exact
     `fcarouge-typed-linear-algebra` name so the package name matches the
     CMake project name on every distribution channel.
   - c. Standardize on "TypedLinearAlgebra" (no spaces) in new conference-talk
     titles, blog posts, and social bios — the README title currently reads
     "Typed Linear Algebra" (with spaces) while the repo/package name has
     none, which is a minor but avoidable search-fragmentation risk.

3. **Add status badges to the README** (build, coverage, license, version).
   They signal "the lights are on" within the first few seconds a stranger
   spends evaluating the project.
   - a. `README.md`'s "Continuous Integration & Deployment Actions" section
     already has Pipeline, Sanitizer, Format, ClangTidy, CppCheck, Doxygen,
     Valgrind, License, License Scan, OpenSSF, and coverage badges — add a
     release/version badge (e.g. `shields.io/github/v/release/...`) next to
     them, since the only version indicator today is the hardcoded `0.3.0` in
     `CMakeLists.txt`.
   - b. Add a GitHub Discussions badge/link once Discussions is enabled (item
     12), matching the existing badge-row convention.
   - c. Add a "Try it on Compiler Explorer" badge/link once the library is
     added to Godbolt (item 22).

4. **Write a specific, differentiated pitch**, not a generic one-liner. A
   vague pitch ("a units library for C++") fails to communicate unique value;
   specific, verifiable claims are what convert an evaluator into a user.
   - a. Rewrite `README.md`'s opening sentence to lead with the concrete claim
     already buried in "Use Cases" — that a mismatched unit, out-of-order
     axis, or mixed-up reference frame *fails to compile* instead of
     producing a silently wrong number — instead of the current generic
     "brings type safety to matrix operations."
   - b. Pull the concrete failure framing already referenced by
     `documentation/.../mars_climate_orgiter_99.png` into the opening
     paragraph; a real historical failure story is more convincing than
     abstract phrasing.
   - c. Name the supported backends (Eigen, `std::linalg`/`mdspan`, Au,
     mp-units, nholthaus/units, Kokkos) directly in the pitch — "works with
     the numerical library you already use" is a real differentiator versus
     a from-scratch linear-algebra implementation.

5. **Put a runnable code example above the fold**, visible without scrolling.
   Users decide whether to keep reading in about five seconds.
   - a. The `state x{...}; x * transposed(x)` block is already first in
     `README.md` — once item 22 lands, link it directly to a Compiler
     Explorer shortlink instead of only to the `sample/` directory.
   - b. Add a one-line comment above that block naming which backend it uses
     (it currently reads as backend-agnostic pseudocode), so a first-time
     reader isn't left guessing which `#include`s are implied.

6. **Add zero-install "try it live" links** wherever code is shown. Zero-
   install prototyping removes the single biggest barrier to a stranger's
   first contact with the code.
   - a. File the Compiler Explorer library-addition request described in item
     22 first — this tip depends on it.
   - b. Once added, replace the `sample/` link in `README.md`'s first code
     block with a `godbolt.org/z/...` shortlink.
   - c. Add the same kind of link to every row of the backend-integration
     table in `README.md` (Au, Eigen, Kokkos, mp-units, …) so each backend is
     independently explorable, not just the default example.

7. **Track GitHub star history against comparable projects.** Stars aren't
   just vanity metrics — they proxy for reach and community trust, and
   inflection points reveal which marketing actions worked.
   - a. Add a `star-history.com` chart to a "Project Health" subsection of
     `README.md`, tracking `FrancoisCarouge/TypedLinearAlgebra` against
     related projects such as `mpusz/mp-units` and `aurora-opensource/au`
     (both already dependencies/acknowledged in "Third Party
     Acknowledgement").
   - b. Re-check that chart after each release announcement (item 30) or
     conference talk to see which channel actually moved the needle, and
     jot the finding somewhere durable (e.g. a short note in `documentation/`).

8. **Keep `master` green at all times.** A green build is a promise made to
   a stranger; one broken build destroys the trust that gets someone to
   invest time evaluating the library.
   - a. `pipeline.yml` already runs on a `cron: '0 0 * * */5'` schedule to
     catch upstream drift — extend the same nightly schedule to
     `sanitizer.yml` and `valgrind_memory.yml` if they don't already have one,
     so drift in those checks doesn't wait for the next PR to surface.
   - b. Treat a red scheduled run (e.g. a FetchContent'd dependency changed
     upstream) as a same-day fix, not a backlog item — the badge row in
     `README.md` is a live promise to every visitor who lands on it.

9. **Build a CI matrix spanning every supported compiler, architecture, and
   configuration axis.** Feature support varies wildly across the C++
   compiler landscape; users are spread across every cell of that matrix.
   - a. `pipeline.yml` already covers Clang 20/21 + libc++, GCC 14/15, and
     MSVC (Debug/Release) on Ubuntu/Windows — add Apple Clang on a macOS
     runner, since macOS is the one toolchain named in `CLAUDE.md`
     ("Clang with libc++, GCC 14+, or MSVC") that the current six-job matrix
     doesn't exercise.
   - b. Verify every job actually builds the *full* backend list from
     `support/support.cmake` (`au_eigen`, `au_std`, `chrono_eigen`,
     `chrono_std`, `eigen`, `eigexed`, `mp_units_eigen`, `mp_units_std`,
     `nested_typed_eigen`, `nholthaus_eigen`, `nholthaus_std`), not a subset —
     for this project the backend list *is* the real "feature matrix," more
     than compiler flags.
   - c. Label the GCC 14 / Clang 20 legs explicitly as "minimum supported" in
     the job `name` field, so a failure there reads immediately as "we broke
     the floor," distinct from a failure on the newest compiler.

10. **Fuzz the build matrix with a reproducible random seed** instead of
    testing every combination, when the combinatorics get too large for
    exhaustive CI.
    - a. The compiler/OS matrix here is small (6 legs) and doesn't need this
      yet — but the backend × test-category combinatorics
      (`support/support.cmake`'s ~11 backends × the 32 `test/` categories)
      are the actual combinatorial risk; if `ctest` ever starts timing out in
      CI, add a `workflow_dispatch` input that samples a random subset of
      backend/category pairs, with `seed=0` meaning "randomize" and any other
      value meaning "reproduce exactly."
    - b. If adopted, print the chosen seed to the job summary
      (`$GITHUB_STEP_SUMMARY`) so a failing nightly run can be reproduced
      locally with the same `ctest --test-dir build -R <regex>` pattern
      already documented in `AGENTS.md`.

11. **Chain a multi-stage quality gate before merge** (pre-commit → build
    matrix → static analysis/tests → integration → merge). Independent
    validation layers mean no single point of failure, and "all checks
    passed" becomes a claim users can actually trust.
    - a. The chain already exists across `pipeline.yml`, `format.yml`,
      `clang_tidy.yml`, `cppcheck.yml`, `codeql.yml`, `sanitizer.yml`,
      `valgrind_memory.yml`, and `coverage.yml` — document it explicitly as a
      single ordered list in `CONTRIBUTING.md` ("what must go green before
      merge"), since today a contributor discovers each check one PR comment
      at a time.
    - b. Confirm every one of those workflows is configured as a *required*
      status check in the `master` branch-protection rule (GitHub repo
      Settings → Branches), not just running informationally — a check that
      can be merged around isn't actually a gate.

12. **Separate Issues (actionable bugs, tasks, specific feature requests,
    security reports) from Discussions** (Q&A, show-and-tell, ideas, RFCs).
    Enabling Discussions is free and keeps the issue tracker actionable.
    - a. Enable GitHub Discussions (Settings → Features) and seed it with
      categories fitting this project: "Q&A" for backend-integration
      questions, "Show and tell" (complementing `README.md`'s existing
      "Projects" list, currently just `Kalman`), and "Ideas" for new
      index-type/backend proposals.
    - b. Add a fourth bullet to `CONTRIBUTING.md`, alongside the existing bug
      template / feature template / security policy links: "Have a question
      or an idea? Start a Discussion instead."
    - c. Migrate any currently-open issues that are really questions or
      design debates into Discussions, so the tracker stays a clean,
      actionable backlog.

13. **Adopt issue templates that require a description, repro steps, a
    Compiler Explorer link, and build config.** This eliminates the "it
    doesn't work" / "what doesn't work?" loop.
    - a. `.github/ISSUE_TEMPLATE/bug_report.md` already asks for OS,
      compiler, version, and commit — add a required **Compiler Explorer
      link** field to its "To Reproduce" section once item 6/22 makes that
      possible.
    - b. Add a **Backend/Integration** field (values matching the identifiers
      in `support/support.cmake`, e.g. `eigen`, `mp_units_eigen`, `au_std`) —
      "which of the ~11 backend integrations" is exactly the environment
      detail this project needs that a generic template doesn't ask for.
    - c. Convert `bug_report.md`/`feature_request.md` from the legacy
      front-matter format to GitHub's YAML `.github/ISSUE_TEMPLATE/*.yml`
      forms, which support a required dropdown for the backend field instead
      of an easy-to-skip free-text prompt.

14. **Write a `CONTRIBUTING.md` with an explicit contribution path.** "Just
    look at the code" and no build instructions are classic symptoms of a
    "fortress" library nobody can get into.
    - a. Add a "Building & Testing Locally" section to `CONTRIBUTING.md`
      inlining the same commands already in `AGENTS.md`
      (`cmake -S . -B build -G Ninja`, `cmake --build build --parallel`,
      `ctest --test-dir build --parallel --verbose`) — today that
      information exists only in the AI-agent guidance file, not the human
      contributor guide.
    - b. Add a "Code Style" section pointing at the exact
      `clang-format-22 -i -style=file` / `cmake-format -i` commands and the
      `.pre-commit-config.yaml` hooks, so a contributor doesn't discover the
      formatting requirement only after a failed `format.yml` run.
    - c. Cross-link `CONTRIBUTING.md` and `AGENTS.md` explicitly ("Human
      contributors: this file. AI coding agents: see `AGENTS.md`.") so the
      two guides stay one source of truth instead of silently diverging.

15. **Label a curated set of small issues "good first issue."** It turns
    "figure it out yourself" into a concrete, low-risk entry point.
    - a. Create the `good first issue` label and apply it to 3-5 small,
      well-scoped tasks — e.g. "add a `sample/` for the Kokkos backend" or
      "add a compile-fail test for mismatched Au units."
    - b. Mention the label explicitly in `CONTRIBUTING.md`'s bug/feature
      bullets: "New to the project? Look for issues labeled
      `good first issue`."

16. **Treat communication as a feature of the library, and respond with
    empathy over ego.** How an issue is handled defines the community around
    it.
    - a. Add a short "How we respond" note to `CONTRIBUTING.md` (e.g. "always
      thank the reporter first, then ask for missing repro details"),
      complementing `CODE_OF_CONDUCT.md`'s conduct rules with concrete triage
      etiquette.
    - b. When closing an issue as "won't fix" or "works as intended," link to
      the relevant entry in `README.md`'s "Lessons Learned" section instead
      of a bare rejection — that design-rationale material already exists,
      it just isn't being cross-linked from triage yet.

17. **Welcome users as future contributors** and lower the friction between
    using and contributing. Adoption and contribution are connected stages of
    the same funnel.
    - a. Add a closing line to `CONTRIBUTING.md` inviting bug reporters to
      submit the fix themselves, pointing at the `good first issue` label
      (item 15) for anyone who wants a warm-up task first.
    - b. In `README.md`'s "Projects" section, add a line inviting maintainers
      of dependent projects (currently just `Kalman`) to also become
      contributors, not only downstream users.

18. **Keep the auto-generated reference where it already does its job, and
    stop duplicating it in prose.** Auto-generated references are "inventory
    documentation" — fine for the reference quadrant, useless for teaching.
    - a. `documentation/Doxyfile` + `deploy_doxygen.yml` already publish a
      full API reference to GitHub Pages — keep it exactly as the Diátaxis
      "Reference" quadrant, and trim `README.md`'s hand-written "Reference"
      section to what Doxygen genuinely can't express (design intent), rather
      than re-describing every member.
    - b. Explicitly label `README.md`'s existing "Lessons Learned" section as
      Diátaxis "Explanation" content (a one-line framing sentence is enough)
      — it already *is* that category, it's just not recognized as one yet.

19. **Structure documentation with the Diátaxis framework** (tutorials,
    how-to guides, explanation, reference). Doxygen-only docs address only
    the reference quadrant.
    - a. Write a dedicated **Tutorial**: a `documentation/tutorial.md` (or a
      new "Tutorial" section in `README.md`) walking a first-time reader from
      `#include` through building a small state vector with one backend
      (Eigen) end to end — distinct from the "Declaration / Template
      Parameters / Member Functions" style walkthrough currently under
      "Reference."
    - b. Write **How-To guides** for the two most common undocumented tasks:
      "How to add a new backend integration plug-in" (design knowledge exists
      in `support/*/fcarouge/*.hpp` but no step-by-step guide) and "How to
      add a new strongly-typed index."
    - c. Add a short navigation line at the top of `README.md` pointing to
      all four categories (Tutorial, How-To, Explanation = "Lessons
      Learned", Reference = Doxygen + README "Reference"), mirroring the
      mp-units navigation-bar pattern from the reference talk.

20. **Use an LLM as a documentation-drafting aid**, with mandatory human
    review for technical accuracy. This can cut prose-writing time roughly
    10x without giving up quality control.
    - a. For the next new backend plug-in or index type, draft its first
      tutorial/how-to prose by feeding an LLM the relevant
      `support/<backend>/fcarouge/*.hpp` header plus the matching design note
      from `AGENTS.md`, then have the maintainer review before merging.
    - b. Keep the LLM out of `README.md`'s "Lessons Learned" section
      specifically — that's first-person design rationale from the
      maintainer's own experience and should stay human-authored (see item
      34's primary-author principle).

21. **Treat writing documentation as a design exercise.** If a feature is
    hard to explain in prose, the API is probably wrong; writing docs
    surfaces corner cases before users do.
    - a. Before merging any new public API surface (a new operator, alias, or
      backend), require a one-paragraph "how would I explain this in
      README's Reference section" draft in the PR description — if it's hard
      to write, treat that as a signal to revisit the design, the way the
      existing "Lessons Learned" entries (e.g. the lvalue-reference-
      assignment tradeoff) clearly document happened before.
    - b. Write the next "Lessons Learned" entry *before* declaring the
      related feature done, not after, so the friction is captured while
      it's fresh.

22. **Get the library added to Compiler Explorer (Godbolt).** It gives
    zero-install prototyping, exact bug-report reproductions, and interactive
    doc embeds, and header-only libraries are trivial to add.
    - a. Prepare a PR/issue against `compiler-explorer/infra`'s library list:
      header-only, single include root (`include/`), noting it needs to pair
      with an already-supported backend (Eigen or mp-units are both already
      on Godbolt) to be useful standalone.
    - b. Pick Eigen as the default paired backend for the initial Godbolt
      entry (it's the most widely available on Godbolt already), and note in
      the request that Au/mp-units/`std::linalg` pairings can follow once the
      base entry works.

23. **Prove "zero overhead" claims with generated assembly** instead of
    asserting them. A concrete asm diff is far more convincing than a claim.
    - a. Add a "Zero-overhead" subsection to `README.md` near "Use Cases"
      showing a Compiler Explorer link comparing a raw Eigen/`std::linalg`
      matrix operation against the equivalent `typed_matrix` call, once item
      22 makes that link possible — mirroring the mp-units `ttg_s`/`ttg`
      example from the reference talk.
    - b. Reuse an existing isolated operation from `benchmark/` for the
      comparison if one already fits, rather than authoring a new example
      from scratch.

24. **Provide a `devcontainer.json` / GitHub Codespaces setup** with the
    toolchain, dependencies, and editor extensions pre-installed. This
    removes the "first install Clang 18, CMake 3.28…" barrier that kills
    drive-by contributions before they start.
    - a. Add `.devcontainer/devcontainer.json` pinning CMake ≥ 4.3,
      `clang++-20`/`clang++-21` + `libc++`, `g++-14`/`g++-15`,
      `clang-format-22`, `clang-tidy-21`, and `cmake-format` — the exact
      toolchain versions already pinned in `AGENTS.md` and
      `.github/workflows/*.yml`, so the container matches CI instead of
      drifting from it.
    - b. Add the C/C++ and CMake Tools VS Code extensions to
      `customizations.vscode.extensions`, and set a `postCreateCommand` that
      runs the configure command from `AGENTS.md`
      (`cmake -S . -B build -G Ninja`).
    - c. Add an "Open in GitHub Codespaces" badge to `README.md`'s
      "Installation & Usage" section and to `CONTRIBUTING.md`, offering it as
      the recommended path for first-time contributors alongside `INSTALL.md`'s
      native-toolchain instructions.

25. **Reuse the same container image locally and in CI.** This guarantees
    local == CI, kills "works on my machine," and still allows fully offline
    development.
    - a. Once `.devcontainer/devcontainer.json` exists (item 24), keep its
      Clang/GCC versions numerically identical to `pipeline.yml`'s matrix
      entries, and note the pairing explicitly in a comment in both files so
      future version bumps happen together.
    - b. Document in `CONTRIBUTING.md`/`AGENTS.md` that the devcontainer is
      optional — contributors on macOS/Windows can still build natively per
      `INSTALL.md` — so it's offered as a shortcut, not a hard requirement.

26. **Publish official packages on Conan Center and/or vcpkg, and support
    CPM/FetchContent, `add_subdirectory()`, and install + `find_package()`.**
    The more consumption paths supported, the fewer reasons a user has to say
    no.
    - a. Author a `conanfile.py` for `fcarouge-typed-linear-algebra`
      (header-only `package_type`, options mirroring the backend list in
      `support/support.cmake`) and submit it to Conan Center, keeping the
      package name identical to the CMake project name (item 2).
    - b. Investigate a `vcpkg.json` port for Microsoft's vcpkg registry —
      MSVC is already a first-class CI target (`pipeline.yml`'s
      `windows-2025`/`cl` jobs), making vcpkg the natural complement for
      Windows users who don't want to `FetchContent` from source.
    - c. Add a "Package Managers" subsection to `INSTALL.md` documenting the
      Conan/vcpkg snippets once available, alongside the existing
      FetchContent and manual-clone instructions.

27. **Expose the library only through modern, namespaced CMake imported
    targets**, never `include_directories()` or global `CMAKE_CXX_FLAGS`.
    - a. Both `README.md` and `INSTALL.md` already correctly document
      `target_link_libraries(your_target PRIVATE fcarouge-typed-linear-algebra::tlinalg)`
      — audit `cmake/` and `include/CMakeLists.txt` to confirm nothing
      outside the `tlinalg` INTERFACE target's own
      `target_include_directories`/`target_compile_features` leaks into a
      consumer's build.
    - b. Keep the MSVC `/std:c++latest` workaround scoped to
      `support/CMakeLists.txt` only, exactly as `support/support.cmake`'s
      comment and `AGENTS.md` already specify, so the root `tlinalg` target
      stays compiler-agnostic for consumers.

28. **Expose configuration switches for the axes that matter, and test every
    combination in CI.** Users are spread across every point in a library's
    supported configuration space.
    - a. This library targets C++26 only (no std-version axis); the real
      configuration axis is backend choice. Confirm every backend in
      `support/support.cmake`'s full list (`au_eigen`, `au_std`,
      `chrono_eigen`, `chrono_std`, `eigen`, `eigexed`, `mp_units_eigen`,
      `mp_units_std`, `nested_typed_eigen`, `nholthaus_eigen`,
      `nholthaus_std`) is exercised on every CI matrix leg, not a default
      subset, so "all green" means "every integration works everywhere."
    - b. Extend `pipeline.yml`'s nightly cron to also build against each
      FetchContent'd dependency's *latest* release (Eigen, mp-units, Au,
      Kokkos, mdspan), catching upstream breakage before a user reports it.

29. **Publish real prose release notes for every version.** A git tag is not
    user communication; releases are marketing events.
    - a. Only `0.1.0` and `0.2.0` are tagged while `CMakeLists.txt` is
      already at `0.3.0` — write a GitHub Release with prose notes (why, 2-3
      highlights, breaking changes) for `0.3.0`, using `CITATION.cff`'s
      version metadata as the anchor for "what's new."
    - b. Add a `CHANGELOG.md` (none currently exists) accumulating the same
      why/highlights/breaking-changes structure per release, browsable
      without digging through the GitHub Releases UI.
    - c. For future releases, write the notes as part of the same PR that
      bumps `VERSION` in `CMakeLists.txt`, so it's a required step, not an
      afterthought.

30. **Actively distribute release announcements** to the channels where
    users actually are. Don't be a hidden gem.
    - a. Post `0.3.0`'s (and future) release notes to r/cpp and relevant
      C++ community Discord/Slack channels — the project already has strong
      announcement material (the CppNow 2026 talk, the `Kalman` dependent
      project) not yet explicitly tied to a release announcement.
    - b. Whenever a release incorporates feedback from the CppNow 2026 talk
      already linked in `README.md`, close the loop publicly ("the community
      asked for X at CppNow, it's in 0.4.0").

31. **Maintain a contributors list and credit contributors by name.** Many
    contributors are motivated by resume-worthy proof of impact.
    - a. `README.md`'s "Sponsors" section recognizes financial supporters but
      has no equivalent for code contributors — add a `CONTRIBUTORS.md` (or a
      parallel "Contributors" README subsection) with GitHub-linked names and
      what they contributed.
    - b. Seed the initial list with `git shortlog -sne`, then keep it updated
      as part of the PR-merge checklist.
    - c. Once release notes exist (item 29), always credit the specific
      contributor by `@handle` for the specific change, not just in the
      standing contributors list.

32. **Offer recommendations, blog-post help, or co-presenting opportunities**
    to significant contributors. This deepens the win-win beyond a name in a
    list.
    - a. For any contributor whose PR meaningfully shapes the library (a new
      backend integration, a new index-safety mechanism), proactively offer a
      LinkedIn recommendation or a co-author credit on a future conference
      submission — the project already has an active speaking track (the
      2026 CppNow talk) a contributor could realistically join.

33. **Publish an explicit AI-contribution policy.** LLM-authored PRs are
    rising fast and risk silently degrading a bleeding-edge codebase.
    - a. Add an "AI-Generated Contributions" section to `CONTRIBUTING.md`
      (the project already has a detailed `AGENTS.md` for AI *agents working
      in* the repo, but nothing yet governs AI-assisted *external
      contributions*): require disclosure of AI assistance, require a
      rationale for any AI-suggested design change, keep diffs scoped to the
      actual fix.
    - b. Note explicitly that the strict `.clang-tidy` (`Checks: '*'`, zero
      `NOLINT` comments per `AGENTS.md`) and zero-warning CI already filter
      out a lot of low-effort AI-generated noise — but a stated policy still
      sets expectations up front rather than relying only on CI rejection.

34. **Apply a "primary author" principle**: don't let an AI- or contributor-
    suggested redesign override an existing decision without asking the
    original author for their rationale first.
    - a. Formalize the existing implicit practice — `README.md`'s "Lessons
      Learned" section already documents *why* designs like strongly-typed
      memory storage and lvalue-reference assignment were rejected — as an
      explicit reviewer rule: check "Lessons Learned" before accepting a
      suggested redesign of `typed_matrix` or anything under
      `include/fcarouge/typed_linear_algebra_internal/algorithm/`.
    - b. When rejecting a suggested redesign for this reason, add the
      rationale as a new "Lessons Learned" entry if it isn't already
      covered, growing that section into the project's living design-
      decision log.

35. **Adopt welcoming-culture communication standards** ("we" language,
    empathy over ego, GitHub Community Standards). This closes the feedback
    loop that sustains a project long-term.
    - a. Audit compile-fail diagnostics and any custom `static_assert`
      messages under `include/` and `test/*_fail.cpp` for tone — a compile
      error is often a user's very first interaction with the library's
      voice, so messages should explain ("this requires an integer index
      because…") rather than blame.
    - b. Confirm the GitHub repo settings page (not just the files) reports
      full Community Standards: Code of Conduct ✅, Contributing ✅, Issue
      templates ✅, plus License, Security Policy, and a filled-in
      description/topics.

36. **Proactively evangelize the project** (talks, meetups, community chat).
    Passive maintenance leads to project death.
    - a. Beyond the existing CppNow 2026 talk, submit a proposal to at least
      one more conference or local meetup (CppCon, C++ on Sea, a regional
      user group) specifically on the index-safety angle — the project's
      most novel claim.
    - b. Participate regularly in the cpplang Slack and in the mp-units/Au
      maintainer communities — the project already depends on and integrates
      with both, making them the highest-yield audiences for adoption and
      contribution.

37. **Write educational content about the underlying domain**, not just the
    API, for niche/domain-specific projects. This broadens the addressable
    audience beyond people who already know they need the library.
    - a. Turn one bullet from `README.md`'s existing "Use Cases" list (e.g.
      "Guidance, navigation, and control," or the Mars-orbiter framing
      already referenced by the `mars_climate_orgiter_99.png` asset) into a
      standalone write-up on dimensional-analysis pitfalls in engineering,
      independent of any API discussion.
    - b. Cross-link that write-up from the "Resources" section's existing
      academic bibliography, giving domain newcomers and API newcomers two
      clearly separate on-ramps.

38. **Proactively solicit production-usage feedback**, and treat "boring"
    complaints (compile time, header bloat) as valuable signal. Real
    deployment problems reveal barriers enthusiastic early adopters won't
    mention unprompted.
    - a. Once Discussions is enabled (item 12), open a pinned "Looking for
      production users" thread asking who has integrated the library beyond
      the already-listed `Kalman` project, and what friction they hit.
    - b. `deploy_time_trace.yml` already tracks compile-time; surface/link
      its output prominently in `README.md` so compile-time trends are
      visible to both the maintainer and prospective production users
      without having to ask.

39. **Automate the repetitive maintainer burden.** This frees time for the
    human-judgment work that can't be automated.
    - a. `.pre-commit-config.yaml` already runs gitleaks, shellcheck,
      cpplint, end-of-file-fixer, and trailing-whitespace — add a stale-
      issue/stale-PR workflow (`actions/stale`), since none of the current
      `.github/workflows/*.yml` handle staleness yet.
    - b. Have an LLM draft the `CHANGELOG.md`/release-notes prose (item 29)
      from the commit log between two tags as a starting point for the
      human-reviewed final version — the same LLM-as-drafting-aid pattern
      already proposed for documentation in item 20.

40. **Delegate triage permissions to trusted contributors and document the
    decision-making process.** It's okay to say "not for this release."
    - a. Once a second regular contributor emerges (tracked via the
      `CONTRIBUTORS.md` from item 31), grant them `triage` role on the
      GitHub repo so they can label/close issues without needing full write
      access.
    - b. Treat `README.md`'s "Lessons Learned" section as documented
      decision-making precedent (item 34), so triage responses can point to
      it as prior art instead of the maintainer re-explaining the same
      rationale each time.

41. **Publicly and transparently acknowledge the project's current
    shortcomings.** This builds trust with evaluators and surfaces real
    problems instead of hiding them.
    - a. Add a short "Known Limitations" subsection to `README.md` (or fold
      it into "Lessons Learned") candidly listing current gaps — no
      `CONTRIBUTORS.md` yet, no Conan/vcpkg package yet, pre-1.0 API per
      `SECURITY.md`'s supported-versions table — so evaluators see honesty,
      not only polish.
    - b. Surface `SECURITY.md`'s existing "< 1.0: security updates provided
      while developing the first release" framing directly in `README.md`,
      so the pre-1.0 stability caveat is visible to a first-time reader
      instead of only someone who goes looking for the security policy.
