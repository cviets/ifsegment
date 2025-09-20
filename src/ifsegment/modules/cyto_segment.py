import numpy as np
from .io_utils import read_czi, get_czi_in_folder, write_tiff, get_well_from_file
from .normalizations import minmax_percentile
from aicsimageio import AICSImage
from tqdm import tqdm
from skimage.morphology import remove_small_objects, remove_small_holes
from skimage.measure import label
from scipy.ndimage import binary_closing

def cyto_segment_array(image: np.ndarray[np.float64]) -> np.ndarray[np.bool_]:
    """
    Create cytoplasmic segmentation of input 2D image

    Parameters 
    -----------
    image: np.ndarray[np.float64]
        Input image as numpy array
    
    Returns 
    -----------
    mask : np.ndarray[np.bool_]
        Mask of image
    """
    assert image.ndim == 2
    mu = np.mean(image[image < 1])
    std = np.std(image[image < 1])
    mask = image > mu + 5*std
    mask = remove_small_objects(mask, min_size=600, connectivity=2)
    mask = remove_small_holes(mask, 100, connectivity=2)
    return mask

def cyto_segment_czi(czi_image: AICSImage, cyto_channel: int) -> np.ndarray[np.bool_]:
    """
    Generates cytoplasmic mask from input AICS image (from read_czi)

    Parameters 
    ----------
    czi_image : AICSImage
        Input AICSImage object (TCZYX) (will access dask arrays from this)
    cyto_channel : int
        Channel containing cytoplasmic information
    """
    img_dask = czi_image.get_image_dask_data("ZYX", T=0, C=cyto_channel)
    img_numpy = img_dask.compute()

    # take z-stack
    img_numpy = np.max(img_numpy, axis=0)
    img_numpy = minmax_percentile(img_numpy, 2, 98)
    return cyto_segment_array(img_numpy)

def cyto_segment_folder(path_to_folder: str, output_folder:str, cyto_channel: int) -> None:
    czi_paths = get_czi_in_folder(path_to_folder)
    for czi_path in tqdm(czi_paths, desc="Cytoplasm masks"):
        img = read_czi(czi_path)
        mask = cyto_segment_czi(img, cyto_channel)

        well_name = get_well_from_file(czi_path)
        write_tiff(mask, output_folder, well_name)