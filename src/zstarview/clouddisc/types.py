from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional
import datetime as dt

@dataclass
class CloudMeta:
    satellite: str         # "G16"|"G18"|"HIMAWARI"
    product: str           # "CMIPF-C13"|"HSD-B13"|"ISatSS-B13"
    time_utc: dt.datetime  # 10分丸め
    src_paths: List[Path]

class VisibilityError(RuntimeError): ...
class DataNotFoundError(RuntimeError): ...
class DownloadError(RuntimeError): ...
class RenderError(RuntimeError): ...
