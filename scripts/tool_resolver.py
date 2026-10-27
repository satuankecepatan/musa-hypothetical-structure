"""Automated binary detection and resolver for IQ-TREE 2 and MAFFT across Windows and Google Colab."""

from pathlib import Path
import os
import shutil
import subprocess
import sys
import urllib.request
import zipfile

# Ensure workspace root is in sys.path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from scripts.config import BIN_DIR, IS_COLAB


def find_binary(binary_name: str) -> str | None:
    """Finds binary executable either in BIN_DIR or system PATH."""
    # Check project bin directory first
    for candidate in [
        BIN_DIR / binary_name,
        BIN_DIR / f"{binary_name}.exe",
        BIN_DIR / f"{binary_name}.bat",
    ]:
        if candidate.is_file():
            return str(candidate)

    # Check system PATH
    found = shutil.which(binary_name)
    if found:
        return found
    if sys.platform == "win32":
        return shutil.which(f"{binary_name}.exe")
    return None


def resolve_iqtree() -> str | None:
    """Resolves or installs IQ-TREE 2."""
    binary = find_binary("iqtree2") or find_binary("iqtree")
    if binary:
        return binary

    if IS_COLAB or sys.platform.startswith("linux"):
        print("[Tool Resolver] Attempting apt-get install for iqtree...")
        try:
            subprocess.run(["apt-get", "update", "-qq"], check=True)
            subprocess.run(["apt-get", "install", "-y", "-qq", "iqtree"], check=True)
            return find_binary("iqtree") or find_binary("iqtree2")
        except Exception as e:
            print(f"[Tool Resolver] Failed to install iqtree via apt: {e}")

    elif sys.platform == "win32":
        # Attempt downloading precompiled standalone Windows binary
        windows_url = "https://github.com/iqtree/iqtree2/releases/download/v2.3.6/iqtree-2.3.6-Windows.zip"
        zip_path = BIN_DIR / "iqtree_windows.zip"
        print(f"[Tool Resolver] Downloading standalone IQ-TREE 2 binary from {windows_url}...")
        try:
            urllib.request.urlretrieve(windows_url, zip_path)
            with zipfile.ZipFile(zip_path, "r") as zip_ref:
                for member in zip_ref.namelist():
                    if member.endswith("iqtree2.exe"):
                        with zip_ref.open(member) as source, open(BIN_DIR / "iqtree2.exe", "wb") as target:
                            target.write(source.read())
                        break
        except Exception as e:
            print(f"[Tool Resolver] Could not auto-download IQ-TREE 2 on Windows: {e}")
        finally:
            if zip_path.exists():
                try:
                    zip_path.unlink()
                except OSError:
                    pass
        return find_binary("iqtree2")

    return None


def resolve_mafft() -> str | None:
    """Resolves or installs MAFFT."""
    binary = find_binary("mafft")
    if binary:
        return binary

    # Check for WSL mafft on Windows (recommended by MAFFT authors for Windows)
    if sys.platform == "win32":
        try:
            res = subprocess.run(["wsl", "which", "mafft"], capture_output=True, text=True)
            if res.returncode == 0 and res.stdout.strip():
                return "wsl mafft"
        except Exception:
            pass

    # Check for mafft.bat inside bin/mafft-win
    mafft_win_candidate = BIN_DIR / "mafft-win" / "mafft.bat"
    if mafft_win_candidate.is_file():
        return str(mafft_win_candidate)

    if IS_COLAB or sys.platform.startswith("linux"):
        print("[Tool Resolver] Attempting apt-get install for mafft...")
        try:
            subprocess.run(["apt-get", "update", "-qq"], check=True)
            subprocess.run(["apt-get", "install", "-y", "-qq", "mafft"], check=True)
            return find_binary("mafft")
        except Exception as e:
            print(f"[Tool Resolver] Failed to install mafft via apt: {e}")

    elif sys.platform == "win32":
        mafft_url = "https://mafft.cbrc.jp/alignment/software/mafft-7.526-win64-signed.zip"
        zip_path = BIN_DIR / "mafft_win.zip"
        print(f"[Tool Resolver] Downloading MAFFT Windows package from {mafft_url}...")
        try:
            # Set User-Agent to avoid HTTP 403 on some mirrors
            req = urllib.request.Request(mafft_url, headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(req) as resp, open(zip_path, "wb") as out_f:
                out_f.write(resp.read())
            with zipfile.ZipFile(zip_path, "r") as zip_ref:
                zip_ref.extractall(BIN_DIR)
            if (BIN_DIR / "mafft-win" / "mafft.bat").is_file():
                return str(BIN_DIR / "mafft-win" / "mafft.bat")
        except Exception as e:
            print(f"[Tool Resolver] Could not auto-download MAFFT on Windows: {e}")
        finally:
            if zip_path.exists():
                try:
                    zip_path.unlink()
                except OSError:
                    pass

    return None


def get_tool_status() -> dict[str, str | None]:
    """Returns mapping of available binary paths."""
    return {
        "iqtree": resolve_iqtree(),
        "mafft": resolve_mafft(),
    }


if __name__ == "__main__":
    status = get_tool_status()
    print("Tool Status:")
    for tool, path in status.items():
        print(f"  {tool}: {path or 'Not Available'}")
