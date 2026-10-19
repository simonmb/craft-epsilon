"""Locating xtb: release asset selection (mocked GitHub release), cache folders, explicit path."""
import sys

import pytest

from craft_epsilon import conformers

RELEASE = {"tag_name": "v6.7.1", "assets": [{"name": n, "browser_download_url": "https://example.invalid/" + n}
                                             for n in ("xtb-6.7.1-linux-x86_64.tar.xz",
                                                       "xtb-6.7.1-linux-x86_64.tar.xz.sha256",
                                                       "xtb-6.7.1-source.tar.xz",
                                                       "xtb-6.7.1pre-windows-x86_64.zip",
                                                       "xtb-6.7.1pre-windows-x86_64.zip.sig")]}


@pytest.mark.parametrize("os_name, arch, expected", [
    ("linux", "x86_64", "xtb-6.7.1-linux-x86_64.tar.xz"),
    ("windows", "amd64", "xtb-6.7.1pre-windows-x86_64.zip"),
    ("darwin", "arm64", None),                 # no macOS binary in the release
    ("linux", "aarch64", None),
])
def test_select_xtb_asset(os_name, arch, expected):
    asset = conformers.select_xtb_asset(RELEASE, os_name, arch)
    assert (asset["name"] if asset else None) == expected


def fake_release(root, tag, name="xtb.exe" if sys.platform == "win32" else "xtb"):
    binary = root / tag / f"xtb-{tag.lstrip('v')}" / "bin" / name
    binary.parent.mkdir(parents=True)
    binary.write_text("", encoding="utf-8")
    return binary


def test_cached_binary_prefers_newest_release(tmp_path):
    fake_release(tmp_path, "v6.6.1")
    newest = fake_release(tmp_path, "v6.7.1")
    assert conformers.cached_xtb_binary([tmp_path]) == str(newest)


def test_cached_binary_ignores_files_outside_bin(tmp_path):
    stray = tmp_path / "v6.7.1" / ("xtb.exe" if sys.platform == "win32" else "xtb")
    stray.parent.mkdir(parents=True)
    stray.write_text("", encoding="utf-8")
    assert conformers.cached_xtb_binary([tmp_path]) is None


def test_environment_variable_overrides_cache_folders(tmp_path, monkeypatch):
    monkeypatch.setenv("CRAFT_EPSILON_XTB_DIR", str(tmp_path))
    assert conformers.xtb_cache_dirs() == [tmp_path]
    assert conformers.writable_xtb_cache_dir() == tmp_path


def test_default_cache_is_inside_the_package(monkeypatch):
    monkeypatch.delenv("CRAFT_EPSILON_XTB_DIR", raising=False)
    package_dir, user_dir = conformers.xtb_cache_dirs()
    assert package_dir.parent.name == "craft_epsilon" and package_dir.name == "xtb"
    assert user_dir.parts[-2:] == (".craft_epsilon", "xtb")


def test_find_uses_cache_when_xtb_is_not_on_path(tmp_path, monkeypatch):
    binary = fake_release(tmp_path, "v6.7.1")
    monkeypatch.setenv("CRAFT_EPSILON_XTB_DIR", str(tmp_path))
    monkeypatch.setattr(conformers.shutil, "which", lambda name: None)
    assert conformers.find_xtb_binary(download=False) == str(binary)


def test_find_without_xtb_and_without_download_raises(tmp_path, monkeypatch):
    monkeypatch.setenv("CRAFT_EPSILON_XTB_DIR", str(tmp_path))
    monkeypatch.setattr(conformers.shutil, "which", lambda name: None)
    with pytest.raises(FileNotFoundError, match="not found"):
        conformers.find_xtb_binary(download=False)


def test_explicit_binary(tmp_path):
    binary = fake_release(tmp_path, "v6.7.1")
    assert conformers.find_xtb_binary(str(binary)) == str(binary.resolve())
    with pytest.raises(FileNotFoundError):
        conformers.find_xtb_binary(str(tmp_path / "missing" / "xtb"))
