from .cyto_segment import cyto_segment_czi
from aicsimageio import AICSImage
import numpy as np
from .normalizations import minmax_percentile
from .io_utils import get_czi_in_folder, read_czi, get_well_from_file, write_tiff
from tqdm import tqdm

def nucleus_segment_czi(czi_image: AICSImage, nuclear_channel: int, cyto_channel: int) -> np.ndarray[np.bool_]:
    """
    Generates nuclear mask from input AICS Image (from read_czi)

    Parameters 
    -----------
    czi_image : AICSImage
        AICSImage object containing CZI image information
    nuclear_channel : int
        Nuclear marker channel number (0-indexed) (eg, Hoechst)
    cyto_channel : int
        Cytoplasmic marker channel number (0-indexed), used to mask nuclei with associated cytoplasm

    Returns
    ---------
    masks : np.ndarray[np.bool_]
        Segmentation array containing cytoplasmic and nuclear masks
    """
    cyto_mask = cyto_segment_czi(czi_image, cyto_channel)

    img_dask = czi_image.get_image_dask_data("ZYX", T=0, C=nuclear_channel)
    img_numpy = img_dask.compute()

    # take max z stack
    img_numpy = np.max(img_numpy, axis=0)
    img_numpy = minmax_percentile(img_numpy, 2, 98)

    # remove nuclei not in cell cytoplasm (usually means no neurites)
    img_numpy[np.logical_not(cyto_mask)] = -1
    # TODO: finish implementation here :-)

def segment_folder(path_to_folder: str, output_folder: str, nuclear_channel: int, cyto_channel: int) -> None:
    czi_paths = get_czi_in_folder(path_to_folder)
    for czi_path in tqdm(czi_paths, desc="Masking cells"):
        img = read_czi(czi_path)
        mask = nucleus_segment_czi(img, nuclear_channel, cyto_channel)

        well_name = get_well_from_file(czi_path)
        write_tiff(mask, output_folder, well_name)