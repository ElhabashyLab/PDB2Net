import tomllib
from pathlib import Path


def test_config_json_files_are_included_as_package_data() -> None:
    pyproject = tomllib.loads(Path("pyproject.toml").read_text(encoding="utf-8"))

    package_data = pyproject["tool"]["setuptools"]["package-data"]

    assert pyproject["tool"]["setuptools"]["include-package-data"] is False
    assert package_data["pdb2net"] == [
        "configs/config.base.json",
        "configs/config.darwin.json",
        "configs/config.linux.json",
        "configs/config.windows.json",
        "configs/config.local.example.json",
    ]
    assert pyproject["tool"]["setuptools"]["exclude-package-data"]["pdb2net"] == [
        "configs/config.local.json",
    ]
