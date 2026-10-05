"""Real bounded processes for search timeout/cancellation, without databases."""

import os
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import signal
import subprocess
import sys
import threading
import time

import pytest

from pdb2net import uniprot_matcher as matcher
from pdb2net import search_process
from pdb2net.search_process import SearchAborted, SearchProcesses


@pytest.mark.parametrize("code", [0, 3])
def test_disabled_limit_preserves_completed_process_output(code):
    result = SearchProcesses(None).run(
        [sys.executable, "-c", f"import sys; print('out'); print('err',file=sys.stderr); sys.exit({code})"], tool="probe"
    )
    assert (result.returncode, result.stdout.strip(), result.stderr.strip()) == (code, "out", "err")


@pytest.mark.skipif(os.name != "posix", reason="Controlled SIGINT after child acquisition")
def test_interrupt_immediately_after_popen_still_reaps_child(monkeypatch):
    real_popen, real_clock = subprocess.Popen, time.monotonic
    children = []
    interrupted = False
    def record(*args, **kwargs):
        child = real_popen(*args, **kwargs)
        children.append(child)
        return child
    def clock():
        nonlocal interrupted
        if children and not interrupted:
            interrupted = True
            os.kill(os.getpid(), signal.SIGINT)
        return real_clock()
    monkeypatch.setattr(search_process.subprocess, "Popen", record)
    monkeypatch.setattr(search_process.time, "monotonic", clock)
    try:
        with pytest.raises(KeyboardInterrupt):
            SearchProcesses(None).run([sys.executable, "-c", "import time; time.sleep(3)"], tool="probe")
        assert children and children[0].poll() is not None
        assert children[0].stdout.closed and children[0].stderr.closed
    finally:
        for child in children:
            if child.poll() is None:
                child.kill()
            child.communicate()


def test_cancel_during_process_launch_stops_child_and_blocks_later_launches(monkeypatch):
    started, release, cancelling = threading.Event(), threading.Event(), threading.Event()
    processes = []
    real_popen = subprocess.Popen

    def delayed_start(*args, **kwargs):
        child = real_popen(*args, **kwargs)
        processes.append(child)
        started.set()
        assert release.wait(3)
        return child

    monkeypatch.setattr(search_process.subprocess, "Popen", delayed_start)
    execution = SearchProcesses(None)
    command = [sys.executable, "-c", "import time; time.sleep(3)"]
    with ThreadPoolExecutor(max_workers=2) as pool:
        result = pool.submit(execution.run, command, tool="probe")
        try:
            assert started.wait(3)
            def cancel():
                cancelling.set()
                execution.cancel()
            cancelled = pool.submit(cancel)
            assert cancelling.wait(1)
        finally:
            release.set()
        cancelled.result(timeout=2)
        with pytest.raises(SearchAborted):
            result.result(timeout=2)
    assert len(processes) == 1 and processes[0].poll() is not None
    with pytest.raises(SearchAborted):
        execution.run(command, tool="must-not-start")
    assert len(processes) == 1


@pytest.mark.skipif(os.name != "posix", reason="SIGTERM-ignoring test process")
def test_timeout_escalates_and_reaps_term_ignoring_child(tmp_path, monkeypatch):
    marker = tmp_path / "pid"
    monkeypatch.setattr(search_process, "TERMINATION_GRACE_SECONDS", .15)
    command = [sys.executable, "-c", "import os,signal,time; from pathlib import Path; "
               "signal.signal(signal.SIGTERM,signal.SIG_IGN); "
               f"Path({str(marker)!r}).write_text(str(os.getpid())); time.sleep(3)"]
    with pytest.raises(SearchAborted, match="time limit"):
        SearchProcesses(.3).run(command, tool="probe")
    assert marker.exists()
    with pytest.raises(ProcessLookupError):
        os.kill(int(marker.read_text()), 0)


@pytest.mark.skipif(os.name != "posix", reason="Fake backend executable uses a POSIX shebang")
@pytest.mark.parametrize("backend", ["blast", "diamond"])
def test_search_timeout_is_an_error_and_removes_query_files(tmp_path, monkeypatch, backend):
    marker = tmp_path / "started"
    tool = tmp_path / "fake-search"
    tool.write_text(f"#!{sys.executable}\n" + f'''
import sys,time
from pathlib import Path
flag = '--out' if '--out' in sys.argv else '-out'
output = Path(sys.argv[sys.argv.index(flag)+1])
Path({str(marker)!r}).write_text(str(output.parent))
time.sleep(1.8)
output.write_text('')
''')
    tool.chmod(0o755)
    monkeypatch.setitem(matcher.config, "sequence_search_timeout_seconds", 1)
    monkeypatch.setitem(matcher.config, "blastp_executable", str(tool))
    monkeypatch.setitem(matcher.config, "blast_db_path", str(tmp_path))
    monkeypatch.setitem(matcher.config, "diamond", {
        "enabled": True, "executable": str(tool), "uniref90_db_path": str(tmp_path / "unused.dmnd"),
        "temp_dir": str(tmp_path),
    })
    search = matcher._run_blastp_swissprot_batch if backend == "blast" else matcher._run_diamond_uniref90_batch
    with pytest.raises(matcher.SequenceSearchError, match="time limit of 1 seconds"):
        search([("probe", "A" * 40, "probe")])
    assert marker.exists()
    assert not Path(marker.read_text()).exists()


@pytest.mark.skipif(os.name != "posix", reason="Fake backend executable uses a POSIX shebang")
def test_later_timeout_cancels_earlier_retry_without_waiting_for_map_order(tmp_path, monkeypatch):
    tool = tmp_path / "fake-blast"
    retry = tmp_path / "retry"
    pending = tmp_path / "pending"
    tool.write_text(f"#!{sys.executable}\n" + f'''
import sys,time
from pathlib import Path
seq = Path(sys.argv[sys.argv.index('-query')+1]).read_text().splitlines()[1]
fmt = sys.argv[sys.argv.index('-outfmt')+1]
if seq.startswith('A') and 'qcovs' in fmt:
    time.sleep(.7)
    print('qcovs unsupported',file=sys.stderr)
    sys.exit(1)
if seq.startswith('A'):
    Path({str(retry)!r}).write_text(str(time.monotonic()))
if seq.startswith('T'):
    Path({str(pending)!r}).touch()
time.sleep(4)
Path(sys.argv[sys.argv.index('-out')+1]).write_text('')
''')
    tool.chmod(0o755)
    monkeypatch.setitem(matcher.config, "blastp_executable", str(tool))
    monkeypatch.setitem(matcher.config, "blast_db_path", str(tmp_path))
    monkeypatch.setitem(matcher.config, "blast", {"batch_max_sequences": 1})
    monkeypatch.setitem(matcher.config, "sequence_search_timeout_seconds", 1)
    cancellations = []
    original = SearchProcesses.cancel
    def observed(self, reason="Sequence search interrupted."):
        cancellations.append(time.monotonic())
        return original(self, reason)
    monkeypatch.setattr(SearchProcesses, "cancel", observed)
    with pytest.raises(matcher.SequenceSearchError, match="time limit of 1 seconds"):
        matcher._run_blastp_swissprot_batch([(s, s * 40, s) for s in "AGT"], max_workers=2)
    assert retry.exists() and cancellations
    # The second input times out while the first input's compatibility retry
    # still has most of its own second budget left. map waits for input one.
    assert 0 <= min(cancellations) - float(retry.read_text()) < .8
    assert not pending.exists()


@pytest.mark.skipif(os.name != "posix", reason="Controlled SIGINT to the probe parent only")
def test_parent_interrupt_during_threaded_search_cleans_up_promptly(tmp_path):
    marker = tmp_path / "marker.json"
    child_code = f'''
import json,os,signal,sys,tempfile,threading,time
from pathlib import Path
from pdb2net import config_loader, uniprot_matcher as search
marker=Path({str(marker)!r})
with tempfile.TemporaryDirectory() as temp:
    tool=Path(temp)/'fake-blast'
    tool.write_text('#!'+sys.executable+'\\nimport json,os,sys,time\\nfrom pathlib import Path\\n'
        + 'Path('+repr(str(marker))+').write_text(json.dumps([os.getpid(),sys.argv[sys.argv.index("-out")+1]]))\\n'
        + 'time.sleep(3)\\n')
    tool.chmod(0o755)
    cfg=dict(config_loader.get_config())
    cfg.update(blastp_executable=str(tool),blast_db_path=temp,sequence_search_timeout_seconds=None)
    config_loader.activate_config(cfg)
    def interrupt():
        deadline=time.monotonic()+3
        while not marker.exists() and time.monotonic()<deadline: time.sleep(.01)
        if marker.exists(): os.kill(os.getpid(),signal.SIGINT)
    watcher=threading.Thread(target=interrupt)
    watcher.start()
    started=time.monotonic()
    try:
        search._run_blastp_swissprot_batch([('probe','A'*40,'probe')])
    except KeyboardInterrupt:
        print(json.dumps({{"interrupted":True,"seconds":time.monotonic()-started}}))
    finally:
        watcher.join(timeout=3)
'''
    result = subprocess.run([sys.executable, "-c", child_code], capture_output=True, text=True, timeout=10)
    assert result.returncode == 0, result.stderr
    observation = json.loads(result.stdout)
    assert observation["interrupted"] and observation["seconds"] < 2
    pid, output = json.loads(marker.read_text())
    with pytest.raises(ProcessLookupError):
        os.kill(pid, 0)
    assert not Path(output).parent.exists()


def test_timeout_aborts_annotation_without_cache_or_diamond(monkeypatch):
    def fail(*args, **kwargs):
        raise matcher.SequenceSearchError("BLASTP exceeded the configured time limit of 1 seconds")
    def forbidden(*args, **kwargs):
        raise AssertionError("A timeout must not be cached or forwarded to DIAMOND")
    monkeypatch.setattr(matcher, "_cache_get", lambda key: (False, None))
    monkeypatch.setattr(matcher, "_cache_put", forbidden)
    monkeypatch.setattr(matcher, "_run_blastp_swissprot_batch", fail)
    monkeypatch.setattr(matcher, "_run_diamond_uniref90_batch", forbidden)
    monkeypatch.setitem(matcher.config, "diamond", {"enabled": True})
    data = [{"pdb_id": "LOCAL", "file_path": "local.pdb", "atom_data": [
        {"chain_id": "A", "unique_chain_id": "LOCAL:A", "uniprot_id": None,
         "molecule_type": "Protein", "residues": [{"residue_name": "ALA"}] * 40}]}]
    with pytest.raises(matcher.SequenceSearchError, match="time limit"):
        matcher.parallel_blast_search(data, max_workers=1)


def test_runtime_limit_does_not_change_store_or_sequence_cache_signature(monkeypatch, tmp_path):
    from pdb2net.precomputed import profile_id

    monkeypatch.setitem(matcher.config, "reference_manifest_id", "unchanged-generation")
    monkeypatch.setitem(matcher.config, "blast_db_path", str(tmp_path))
    monkeypatch.setitem(matcher.config, "uniprot_fasta_path", str(tmp_path / "missing.fasta"))
    monkeypatch.setitem(matcher.config, "structure_model_policy", "first")
    monkeypatch.setitem(matcher.config, "sequence_search_timeout_seconds", None)
    before = profile_id(), matcher._db_signature()
    monkeypatch.setitem(matcher.config, "sequence_search_timeout_seconds", 30)
    assert (profile_id(), matcher._db_signature()) == before
