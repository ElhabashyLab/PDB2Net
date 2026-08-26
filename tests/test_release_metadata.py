import re
import tomllib
from pathlib import Path

from pdb2net import __version__


ROOT = Path(__file__).resolve().parents[1]


def _citation_value(name: str) -> str:
    content = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
    match = re.search(rf'^{re.escape(name)}:\s*["\']?([^"\'\n]+)', content, re.MULTILINE)
    assert match is not None
    return match.group(1).strip()


def test_release_version_metadata_is_synchronized() -> None:
    project = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    assert project["project"]["version"] == __version__
    assert _citation_value("version") == __version__

    changelog_version = __version__.replace("rc", "-rc")
    changelog = (ROOT / "CHANGELOG.md").read_text(encoding="utf-8")
    assert re.search(rf"^## v{re.escape(changelog_version)}\b", changelog, re.MULTILINE)


def test_unreleased_software_citation_does_not_claim_a_manuscript_release() -> None:
    citation = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
    assert "date-released:" not in citation
    assert "preferred-citation:" not in citation
    assert "manuscript" not in citation.lower()


def test_cytoscape_is_optional_and_quality_tools_are_development_only() -> None:
    project = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))[
        "project"
    ]
    assert all("py4cytoscape" not in dependency for dependency in project["dependencies"])
    assert project["optional-dependencies"]["cytoscape"] == [
        "py4cytoscape==1.11.0"
    ]
    assert any(
        dependency.startswith("ruff")
        for dependency in project["optional-dependencies"]["dev"]
    )
