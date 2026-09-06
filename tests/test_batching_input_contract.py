import multiprocessing
import os
from pathlib import Path
import time

import pytest

from pdb2net import batching
from pdb2net.batching import batch_run
from pdb2net.input_contract import InputValidationError


def _wait_for_grandchild(_files) -> None:
    grandchild = multiprocessing.Process(target=time.sleep, args=(10,))
    grandchild.start()
    grandchild.join()


class FakeProcess:
    def __init__(
        self,
        *,
        exitcode: int | None,
        alive: bool = False,
        terminate_stops: bool = True,
        pid: int | None = None,
    ) -> None:
        self.exitcode = exitcode
        self.pid = pid
        self._alive = alive
        self._terminate_stops = terminate_stops
        self.kill_called = False

    def start(self) -> None:
        pass

    def join(self, timeout: float | None = None) -> None:
        del timeout

    def is_alive(self) -> bool:
        return self._alive

    def terminate(self) -> None:
        if self._terminate_stops:
            self._alive = False

    def kill(self) -> None:
        self.kill_called = True
        self._alive = False


def _write_batch_inputs(tmp_path: Path, count: int = 2, size: int = 1024) -> Path:
    input_dir = tmp_path / "inputs"
    input_dir.mkdir()
    for index in range(count):
        (input_dir / f"{index + 1}abc.cif").write_bytes(b"x" * size)
    return input_dir


def _process_factory(outcomes: list[tuple[int | None, bool, bool]]):
    remaining = iter(outcomes)
    created: list[FakeProcess] = []

    def factory(*, target, args):
        del target, args
        exitcode, alive, terminate_stops = next(remaining)
        process = FakeProcess(
            exitcode=exitcode,
            alive=alive,
            terminate_stops=terminate_stops,
        )
        created.append(process)
        return process

    factory.created = created
    return factory


def test_batch_run_rejects_missing_input_folder_before_processing(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setitem(batching.config, "output_path", str(tmp_path / "outputs"))
    missing = tmp_path / "missing"

    with pytest.raises(InputValidationError) as exc_info:
        batch_run(str(missing), lambda files: None)

    assert exc_info.value.code == "INPUT_PATH_NOT_FOUND"
    assert list((tmp_path / "outputs").glob("20*/manifest.json"))


def test_batch_run_rejects_file_input_before_processing(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setitem(batching.config, "output_path", str(tmp_path / "outputs"))
    file_path = tmp_path / "single.cif"
    file_path.write_text("data_single\n", encoding="utf-8")

    with pytest.raises(InputValidationError) as exc_info:
        batch_run(str(file_path), lambda files: None)

    assert exc_info.value.code == "INPUT_PATH_NOT_DIRECTORY"


def test_batch_run_rejects_folder_without_supported_structure_files(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setitem(batching.config, "output_path", str(tmp_path / "outputs"))
    (tmp_path / "notes.txt").write_text("not a structure\n", encoding="utf-8")

    with pytest.raises(InputValidationError) as exc_info:
        batch_run(str(tmp_path), lambda files: None)

    assert exc_info.value.code == "NO_VALID_INPUT_FILES"


@pytest.mark.parametrize("open_in_cytoscape", [False, True])
def test_batch_run_continues_after_any_child_failure(
    tmp_path: Path, monkeypatch, open_in_cytoscape: bool
) -> None:
    input_dir = _write_batch_inputs(tmp_path)
    output_dir = tmp_path / "outputs"
    monkeypatch.setitem(batching.config, "output_path", str(output_dir))
    monkeypatch.setitem(batching.config, "open_in_cytoscape", open_in_cytoscape)
    process_factory = _process_factory([(7, False, True), (0, False, True)])

    complete = batch_run(
        str(input_dir),
        lambda _files: None,
        size_limit_kb=1,
        process_factory=process_factory,
    )

    assert complete is False
    assert len(process_factory.created) == 2
    log = (output_dir / "error_in_batch_log" / "skipped_batches.txt").read_text(
        encoding="utf-8"
    )
    assert "exitcode=7" in log


def test_batch_run_returns_false_after_timeout(tmp_path: Path, monkeypatch) -> None:
    input_dir = _write_batch_inputs(tmp_path)
    output_dir = tmp_path / "outputs"
    monkeypatch.setitem(batching.config, "output_path", str(output_dir))
    monkeypatch.setitem(batching.config, "open_in_cytoscape", False)

    process_factory = _process_factory([(None, True, False), (0, False, True)])
    complete = batch_run(
        str(input_dir),
        lambda _files: None,
        timeout_minutes=1,
        size_limit_kb=1,
        process_factory=process_factory,
    )

    assert complete is False
    assert len(process_factory.created) == 2
    assert process_factory.created[0].kill_called is True
    assert "timeout>1min" in (
        output_dir / "error_in_batch_log" / "skipped_batches.txt"
    ).read_text(encoding="utf-8")


def test_batch_run_returns_true_when_every_batch_succeeds(
    tmp_path: Path, monkeypatch
) -> None:
    input_dir = _write_batch_inputs(tmp_path)
    monkeypatch.setitem(batching.config, "output_path", str(tmp_path / "outputs"))
    monkeypatch.setitem(batching.config, "open_in_cytoscape", False)

    assert batch_run(
        str(input_dir),
        lambda _files: None,
        size_limit_kb=1,
        process_factory=_process_factory([(0, False, True), (0, False, True)]),
    ) is True


def test_batch_run_returns_false_when_input_is_too_large(
    tmp_path: Path, monkeypatch
) -> None:
    input_dir = _write_batch_inputs(tmp_path, count=1, size=2048)
    output_dir = tmp_path / "outputs"
    monkeypatch.setitem(batching.config, "output_path", str(output_dir))
    monkeypatch.setitem(batching.config, "open_in_cytoscape", False)

    complete = batch_run(
        str(input_dir),
        lambda _files: None,
        size_limit_kb=1,
        process_factory=lambda **_kwargs: (_ for _ in ()).throw(
            AssertionError("no process should start")
        ),
    )

    assert complete is False
    assert "size_kb=2>1" in (
        output_dir / "error_in_batch_log" / "skipped_batches.txt"
    ).read_text(encoding="utf-8")


def test_batch_run_returns_false_when_no_batch_remains(
    tmp_path: Path, monkeypatch
) -> None:
    input_dir = _write_batch_inputs(tmp_path, count=1)
    monkeypatch.setitem(batching.config, "output_path", str(tmp_path / "outputs"))
    monkeypatch.setitem(batching.config, "open_in_cytoscape", False)
    monkeypatch.setattr(
        batching,
        "create_batches_streaming",
        lambda *_args, **_kwargs: iter(()),
    )

    assert batch_run(str(input_dir), lambda _files: None) is False


def test_batch_run_returns_false_after_success_plus_oversized_input(
    tmp_path: Path, monkeypatch
) -> None:
    input_dir = _write_batch_inputs(tmp_path, count=1, size=1024)
    (input_dir / "2abc.cif").write_bytes(b"x" * 2048)
    output_dir = tmp_path / "outputs"
    monkeypatch.setitem(batching.config, "output_path", str(output_dir))
    monkeypatch.setitem(batching.config, "open_in_cytoscape", False)

    complete = batch_run(
        str(input_dir),
        lambda _files: None,
        size_limit_kb=1,
        process_factory=_process_factory([(0, False, True)]),
    )

    assert complete is False
    assert "size_kb=2>1" in (
        output_dir / "error_in_batch_log" / "skipped_batches.txt"
    ).read_text(encoding="utf-8")


def test_batch_run_returns_false_after_success_plus_stat_error(
    tmp_path: Path, monkeypatch
) -> None:
    input_dir = _write_batch_inputs(tmp_path)
    output_dir = tmp_path / "outputs"
    monkeypatch.setitem(batching.config, "output_path", str(output_dir))
    monkeypatch.setitem(batching.config, "open_in_cytoscape", False)
    getsize = batching.os.path.getsize

    def fail_second_file(path):
        if str(path).endswith("2abc.cif"):
            raise OSError("file disappeared")
        return getsize(path)

    monkeypatch.setattr(batching.os.path, "getsize", fail_second_file)

    complete = batch_run(
        str(input_dir),
        lambda _files: None,
        size_limit_kb=1,
        process_factory=_process_factory([(0, False, True)]),
    )

    assert complete is False
    assert "size_error=OSError" in (
        output_dir / "error_in_batch_log" / "skipped_batches.txt"
    ).read_text(encoding="utf-8")


def test_child_rejects_same_path_input_replacement(
    tmp_path: Path,
) -> None:
    source = tmp_path / "1abc.cif"
    source.write_text("original\n", encoding="utf-8")
    expected = batching.input_file_signature(source)
    source.unlink()
    source.write_text("replacement\n", encoding="utf-8")
    with pytest.raises(InputValidationError) as error:
        batching._validate_batch_signatures(
            [str(source)],
            {str(source): expected},
        )

    assert error.value.code == "INPUT_CHANGED_DURING_PROCESSING"


@pytest.mark.skipif(os.name != "posix", reason="POSIX process-group behavior")
def test_real_timeout_stops_descendant_process_tree(
    tmp_path: Path, monkeypatch
) -> None:
    input_dir = _write_batch_inputs(tmp_path, count=1)
    monkeypatch.setitem(batching.config, "output_path", str(tmp_path / "outputs"))
    monkeypatch.setitem(batching.config, "open_in_cytoscape", False)
    monkeypatch.setattr(
        batching.multiprocessing,
        "Process",
        multiprocessing.get_context("spawn").Process,
    )
    started = time.monotonic()

    complete = batch_run(
        str(input_dir),
        _wait_for_grandchild,
        timeout_minutes=0.01,
        size_limit_kb=1,
    )

    assert complete is False
    assert time.monotonic() - started < 3


@pytest.mark.skipif(os.name != "posix", reason="POSIX process-group behavior")
def test_process_group_readiness_race_falls_back_to_direct_termination(
    monkeypatch,
) -> None:
    process = FakeProcess(exitcode=None, alive=True, pid=12345)
    monkeypatch.setattr(
        batching.os,
        "killpg",
        lambda _pid, _signal: (_ for _ in ()).throw(ProcessLookupError()),
    )

    assert batching._stop_batch_process(process, terminate_tree=True) is True
    assert process.is_alive() is False


def test_windows_tree_kill_failure_is_not_reported_as_confirmed(
    monkeypatch,
) -> None:
    process = FakeProcess(exitcode=None, alive=True, pid=12345)
    monkeypatch.setattr(batching.os, "name", "nt")
    monkeypatch.setattr(
        batching.subprocess,
        "run",
        lambda *_args, **_kwargs: type("Result", (), {"returncode": 1})(),
    )

    assert batching._stop_batch_process(process, terminate_tree=True) is False
    assert process.is_alive() is False
