process DOWNLOAD_OMAMER_DB {
    label 'exogap_python'
    output:
    path("LUCA.h5")

    script:
    """
    python - <<'PYDOWNLOAD'
from urllib.request import urlopen
import shutil, hashlib
with urlopen('${params.omamer_download_url}', timeout=120) as source, open('LUCA.h5', 'wb') as destination:
    shutil.copyfileobj(source, destination)
with open('LUCA.h5', 'rb') as source:
    digest = hashlib.md5()
    for block in iter(lambda: source.read(1024 * 1024), b''):
        digest.update(block)
if digest.hexdigest() != '${params.omamer_download_md5}':
    raise ValueError('OMA database checksum does not match the pinned download')
PYDOWNLOAD
    """

    stub:
    """
    touch LUCA.h5
    """
}
