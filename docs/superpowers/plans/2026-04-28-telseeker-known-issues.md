# TelSeeker Known Issues Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix the known TelSeeker contract and consistency issues found during review.

**Architecture:** Keep changes local to TelSeeker modules. Add focused unittest coverage for metadata parsing, direct connection metadata, Step2 completeness detection, and TRC consistency.

**Tech Stack:** Python 3, unittest, Biopython-dependent project modules.

---

### Task 1: Metadata Contract

**Files:**
- Modify: `bin/TelSeeker.py`
- Modify: `bin/TelSeekerPart2.py`
- Test: `tests/test_telseeker_contracts.py`

- [ ] Write failing tests for `linker.info` parsing and direct metadata writing.
- [ ] Implement robust key-value parsing in `TelSeeker._parse_linker_info`.
- [ ] Preserve `method=direct` in `finalize_chr_end_result`.
- [ ] Run `python3 -m unittest tests.test_telseeker_contracts`.

### Task 2: Step2 Completeness

**Files:**
- Modify: `bin/TelSeeker.py`
- Test: `tests/test_telseeker_contracts.py`

- [ ] Write failing test where `part2.chr.end.job` has only one of two required ends.
- [ ] Make `_check_step2_complete` read `need_extension_chr_end.txt` and require a terminal artifact per expected end.
- [ ] Treat no required ends as complete without running Part2.
- [ ] Run `python3 -m unittest tests.test_telseeker_contracts`.

### Task 3: TRC Consistency

**Files:**
- Modify: `bin/TelSeekerPart1.py`
- Test: `tests/test_telseeker_contracts.py`

- [ ] Write failing test that Part1 and TelSeekerCheck calculate identical TRC for the same input.
- [ ] Remove the divergent local TRC implementation in Part1 by delegating to `TelSeekerCheck`.
- [ ] Update the generated worker script to use non-overlapping matching semantics.
- [ ] Run `python3 -m unittest tests.test_telseeker_contracts`.

### Task 4: Final Verification

**Files:**
- Verify: `bin/TelSeeker.py`
- Verify: `bin/TelSeekerCheck.py`
- Verify: `bin/TelSeekerPart1.py`
- Verify: `bin/TelSeekerPart2.py`
- Verify: `bin/TelSeekerVisualizer.py`
- Verify: `bin/DEGAP.py`

- [ ] Run `python3 -m py_compile bin/TelSeeker.py bin/TelSeekerCheck.py bin/TelSeekerPart1.py bin/TelSeekerPart2.py bin/TelSeekerVisualizer.py bin/DEGAP.py`.
- [ ] Run `python3 -m unittest tests.test_telseeker_contracts`.
