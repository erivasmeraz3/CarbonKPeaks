"""
Build script for CarbonKPeaks application.
Creates a standalone executable with embedded icon.
"""

import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

def install_pyinstaller():
    """Ensure PyInstaller is installed."""
    try:
        import PyInstaller
        print(f"PyInstaller version: {PyInstaller.__version__}")
    except ImportError:
        print("Installing PyInstaller...")
        subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'pyinstaller', '-q'])

def build():
    """Build the executable."""
    install_pyinstaller()

    base_dir = Path(__file__).parent
    main_script = base_dir / "c1s_peak_viewer_gui_final.py"
    icon_file = base_dir / "carbonpeaks.ico"

    # Splash image for PyInstaller bootloader (shown during extraction before
    # Python starts). Prefer the 512px render: the bootloader window is not
    # DPI-aware, so scaled displays upscale the image — a larger source
    # keeps it sharp.
    splash_image = None
    for size in [512, 256, 128, 64]:
        candidate = base_dir / f"carbonpeaks_{size}.png"
        if candidate.exists():
            splash_image = candidate
            break

    # Data files to include (icon PNGs for runtime)
    data_files = []
    for png in base_dir.glob("carbonpeaks_*.png"):
        data_files.append(f"--add-data={png};.")
    if icon_file.exists():
        data_files.append(f"--add-data={icon_file};.")

    # Hidden imports that PyInstaller might miss
    hidden_imports = [
        '--hidden-import=PIL._tkinter_finder',
        '--hidden-import=scipy.special._cdflib',
        '--hidden-import=scipy._lib.array_api_compat.numpy.fft',
        '--hidden-import=lmfit',
        '--hidden-import=lmfit.models',
    ]

    # Build on the local disk: cloud-synced folders (Google Drive/OneDrive)
    # lock freshly written files mid-build and make PyInstaller's copyfile
    # step fail. The finished exe is copied back into dist/ afterwards.
    work_root = Path(tempfile.gettempdir()) / "carbonkpeaks_build"
    work_path = work_root / "build"
    temp_dist = work_root / "dist"

    # PyInstaller command
    cmd = [
        sys.executable, '-m', 'PyInstaller',
        '--name=CarbonKPeaks',
        '--onefile',  # Single executable
        '--windowed',  # No console window
        f'--icon={icon_file}',
        '--clean',  # Clean cache before building
        '--noconfirm',  # Overwrite without asking
        f'--workpath={work_path}',
        f'--distpath={temp_dist}',
        f'--specpath={work_root}',
    ]

    # Add splash screen (shown by bootloader during extraction, before Python starts)
    if splash_image is not None:
        cmd.append(f'--splash={splash_image}')

    cmd += data_files + hidden_imports + [str(main_script)]

    print("\nBuilding CarbonKPeaks...")
    print(f"Command: {' '.join(cmd)}\n")

    result = subprocess.run(cmd, cwd=str(base_dir))

    if result.returncode == 0:
        exe_path = base_dir / "dist" / "CarbonKPeaks.exe"
        exe_path.parent.mkdir(exist_ok=True)
        # Copy under a temp name, then rename: cloud sync (Google Drive)
        # blocks direct writes to a just-synced filename, but renames pass.
        tmp_target = exe_path.with_name("CarbonKPeaks.exe.part")
        shutil.copy2(temp_dist / "CarbonKPeaks.exe", tmp_target)
        if exe_path.exists():
            exe_path.unlink()
        os.replace(tmp_target, exe_path)
        print(f"\n{'='*60}")
        print(f"Build successful!")
        print(f"Executable: {exe_path}")
        print(f"{'='*60}")
    else:
        print(f"\nBuild failed with return code {result.returncode}")

if __name__ == '__main__':
    build()
