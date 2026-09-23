"""Per-user background supervisor. It never starts SecureCRT or submits jobs."""
import json
import msvcrt
import os
from pathlib import Path
import subprocess
import sys
import time
import urllib.request
import urllib.error
import webbrowser

ROOT = Path(__file__).resolve().parent
URL = 'http://127.0.0.1:8765/'


def health():
    try:
        with urllib.request.urlopen(URL + 'api/health', timeout=3) as response:
            return json.load(response).get('app') == 'gwas-local-console'
    except (OSError, ValueError):
        return False


def main():
    if '--open' in sys.argv:
        subprocess.Popen([str(Path(sys.executable).with_name('pythonw.exe')), str(ROOT / 'launcher.py')],
                         creationflags=subprocess.CREATE_NO_WINDOW, cwd=ROOT)
        for _ in range(40):
            if health():
                webbrowser.open(URL)
                return
            time.sleep(.5)
        raise RuntimeError('GWAS background did not start; inspect service.log')
    os.chdir(ROOT)
    guard = open(ROOT / 'launcher.lock', 'a+b')
    if guard.tell() == 0:
        guard.write(b'0'); guard.flush()
    guard.seek(0)
    try:
        msvcrt.locking(guard.fileno(), msvcrt.LK_NBLCK, 1)
    except OSError:
        return  # One supervisor per local project.
    while True:
        if health():
            time.sleep(5)
            continue
        env = dict(os.environ, GWAS_WEB_PORT='8765', PYTHONIOENCODING='utf-8')
        with open(ROOT / 'service.log', 'ab') as output:
            process = subprocess.Popen([sys.executable, str(ROOT / 'server.py')], cwd=ROOT,
                                       stdin=subprocess.DEVNULL, stdout=output, stderr=output,
                                       creationflags=subprocess.CREATE_NO_WINDOW, env=env)
            process.wait()
        time.sleep(5)


if __name__ == '__main__':
    main()
