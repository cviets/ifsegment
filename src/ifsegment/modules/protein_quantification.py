from aicsimageio import AICSImage
import numpy as np
from typing import Tuple

def quantify_czi(
        czi_image: AICSImage, 
        nuclear_channel: int, 
        cyto_channel: int
        ) -> Tuple[np.ndarray[np.bool_], ]: