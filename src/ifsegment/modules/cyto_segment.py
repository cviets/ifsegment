import numpy as np
from .io_utils import read_czi, get_czi_in_folder, write_tiff, get_well_from_file
from .normalizations import minmax_percentile, clip_512
from aicsimageio import AICSImage
from tqdm import tqdm
from skimage.morphology import remove_small_objects, remove_small_holes, diamond
from skimage.measure import label
from scipy.ndimage import binary_closing, binary_dilation, binary_erosion

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

    # we are interested in the non-one values 
    # (want to ignore super outlier-y bright spots)
    mu = np.mean(image[image < 1])
    std = np.std(image[image < 1])

    # hard coded: mask according to mean +/- 5 std for non-bright spots
    thresh = min(0.75, mu + 5*std)
    mask = image > thresh

    diamond_strel = diamond(1)
    mask = binary_dilation(mask, diamond_strel, iterations=3)
    mask = remove_small_holes(mask, 1200, connectivity=2)
    mask = binary_erosion(mask, diamond_strel, iterations=3)

    mask = remove_small_objects(mask, min_size=600, connectivity=2)
    
    return mask

def czi_preprocess(czi_image: AICSImage, cyto_channel: int, mode:str="max") -> np.ndarray[np.float64]: 
    mode = mode.lower()
    assert mode in {"max", "average", "avg", "mean"}

    img_dask = czi_image.get_image_dask_data("ZYX", T=0, C=cyto_channel)
    img_numpy = img_dask.compute()

    # take z-stack
    if mode == "max":
        img_numpy = np.max(img_numpy, axis=0)
    else:
        img_numpy = np.mean(img_numpy, axis=0)
    
    # HARD CODED: normalize 0 percentile to -1 and 95 percentile to 1
    img_numpy = clip_512(img_numpy)
    img_numpy = minmax_percentile(img_numpy, 0, 95)
    return img_numpy

def cyto_segment_czi(czi_image: AICSImage, cyto_channel: int, mode:str="max") -> np.ndarray[np.bool_]:
    preprocessed = czi_preprocess(czi_image, cyto_channel, mode)
    return cyto_segment_array(preprocessed)

def cyto_segment_folder(path_to_folder: str, output_folder:str, cyto_channel: int, mode:str="max") -> None:
    czi_paths = get_czi_in_folder(path_to_folder)
    for czi_path in tqdm(czi_paths, desc="Cytoplasm masks"):
        img = read_czi(czi_path)
        well_name = get_well_from_file(czi_path)
        mask = cyto_segment_czi(img, cyto_channel, mode)
        
        write_tiff(mask, output_folder, well_name)